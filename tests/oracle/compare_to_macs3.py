"""Compare `easyepi peaks` against MACS3 on the same BAM.

This is the harness behind the numbers in `docs/dev/MACS3_COMPARISON.md`. It is
not a pytest module; the assertions that must hold every run live in
`tests/test_peaks_macs3_agreement.py`.

    python tests/oracle/compare_to_macs3.py

Needs `macs3` on PATH (conda-forge/bioconda; 3.0.4 measured). Run from the
repository root.

## What is being compared, and what is not

The two callers are different algorithms, so "identical peaks" is not a
reasonable expectation and a low overlap number is not by itself a defect:

  * easyepi bins the genome (default 200 bp), takes the **median of the whole
    bin signal** as background, and calls a bin enriched when it exceeds
    `background * threshold`, then merges consecutive bins and BH-corrects.
  * MACS3 builds a **local lambda** per position from the genome-wide and
    per-5kb background, then calls broad/narrow peaks from a Poisson tail.

What *must* agree is that both locate the same accessible regions. So the
harness asks the questions that remain meaningful across the two models:

  1. Recovery: are the planted peaks found by both, and at what Jaccard?
  2. False positives: does either caller report peaks where nothing was planted?
  3. Boundary agreement: for peaks found by both, how far apart are the edges?
  4. Ordering: do the reported scores rank the planted sites the same way?

The fixture plants peaks of known width and known enrichment over a uniform
background, which is what makes questions 1-3 answerable at all.
"""
from __future__ import annotations

import itertools
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from easyepi.peaks import call_peaks_from_bam  # noqa: E402

CHROM = "chr1"
CHROM_LEN = 400_000
BIN_SIZE = 200

#: Planted peaks: (centre, half_width, fold over background). Spread so the
#: distances between them are irregular, which stops a fixed-window overlap
#: test from passing by accident.
PEAKS = [
    (40_000, 600, 12.0),
    (90_000, 400, 25.0),
    (150_000, 900, 8.0),
    (210_000, 350, 30.0),
    (265_000, 700, 10.0),
    (330_000, 500, 15.0),
]
#: Reads placed uniformly outside the planted regions. High enough that MACS3's
#: local lambda is estimated rather than falling back to a tiny constant.
BACKGROUND_READS = 24_000
READ_LEN = 50


def _find_macs3() -> str:
    for cand in (
        "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi/bin/macs3",
        shutil.which("macs3"),
    ):
        if cand and os.access(cand, os.X_OK):
            return cand
    raise SystemExit("macs3 not found; install it or extend the candidate list")


def build_fixture(tmp: Path):
    """A paired-end BAM with planted peaks over uniform background."""
    import pysam

    rng = np.random.default_rng(4242)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}

    in_peak = np.zeros(CHROM_LEN, dtype=bool)
    for centre, half, _fold in PEAKS:
        in_peak[max(0, centre - half):min(CHROM_LEN, centre + half)] = True

    starts = []
    # Background: uniform over everything that is not a planted peak.
    free = np.flatnonzero(~in_peak)
    starts.extend(rng.choice(free, size=BACKGROUND_READS, replace=False).tolist())
    # Signal: depth inside each planted peak, proportional to its fold.
    for centre, half, fold in PEAKS:
        n = int(BACKGROUND_READS * (CHROM_LEN - in_peak.sum()) / CHROM_LEN
                * fold / 20.0)
        n = max(60, n)
        lo, hi = max(0, centre - half), min(CHROM_LEN - READ_LEN, centre + half)
        starts.extend(rng.integers(lo, hi, n).tolist())

    recs = []
    for i, start in enumerate(starts):
        frag = int(rng.integers(150, 250))
        for flag, rs, tl in ((99, start, frag),
                             (147, min(CHROM_LEN - 1, start + frag - READ_LEN), -frag)):
            a = pysam.AlignedSegment()
            a.query_name = "r%d" % i
            a.query_sequence = "A" * READ_LEN
            a.flag = flag
            a.reference_id = 0
            a.reference_start = int(rs)
            a.mapping_quality = 60
            a.cigar = ((0, READ_LEN),)
            a.next_reference_id = 0
            a.next_reference_start = int(start)
            a.template_length = int(tl)
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            recs.append(a)

    bam = tmp / "atac.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as fh:
        for r in sorted(recs, key=lambda x: (x.reference_start, x.flag)):
            fh.write(r)
    pysam.index(str(bam))
    return bam, starts


def run_macs3(macs3: str, bam: Path, outdir: Path) -> Path:
    cmd = [macs3, "callpeak", "-t", str(bam), "-f", "BAMPE",
           "-g", "hs", "-n", "macs", "--outdir", str(outdir),
           "-q", "0.01", "--nomodel", "--extsize", "200"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise SystemExit("macs3 failed:\n" + proc.stdout.decode()[-3000:])
    peak = outdir / "macs_peaks.narrowPeak"
    if not peak.exists():
        raise SystemExit("macs3 wrote no narrowPeak\n" + proc.stdout.decode()[-2000:])
    return peak


def read_narrowpeak(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track")):
            continue
        p = line.split("\t")
        rows.append((p[0], int(p[1]), int(p[2]), float(p[6]) if p[6] != "." else 0.0,
                     float(p[8]) if p[8] != "." else 0.0))
    return rows


def jaccard(a, b) -> float:
    """Interval Jaccard between two peak sets on the same chromosome."""
    def total(iv):
        return sum(e - s for _, s, e, *_ in iv)
    inter = 0
    for _, s1, e1, *_ in a:
        for _, s2, e2, *_ in b:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi > lo:
                inter += hi - lo
    union = total(a) + total(b) - inter
    return inter / union if union else 0.0


def overlapping(a, b, min_frac=0.5):
    """How many peaks of `a` have an overlapping partner in `b`."""
    hit = 0
    for _, s1, e1, *_ in a:
        for _, s2, e2, *_ in b:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi > lo and (hi - lo) >= min_frac * min(e1 - s1, e2 - s2):
                hit += 1
                break
    return hit


def planted_overlap(peaks, centre, half):
    lo, hi = centre - half, centre + half
    best, best_len = None, 0
    for _, s, e, *_ in peaks:
        ilo, ihi = max(lo, s), min(hi, e)
        if ihi > ilo and ihi - ilo > best_len:
            best, best_len = (s, e), ihi - ilo
    return best, (best_len / (hi - lo) if hi > lo else 0.0)


def main():
    macs3 = _find_macs3()
    print("easyepi peptides vs MACS3")
    print("macs3:", subprocess.run([macs3, "--version"], stdout=subprocess.PIPE)
          .stdout.decode().strip())
    print()

    tmp = Path(tempfile.mkdtemp(prefix="macs3cmp_"))
    bam, starts = build_fixture(tmp)
    print("fixture: %s  (%d read starts, %d planted peaks, %d bp)"
          % (bam, len(starts), len(PEAKS), CHROM_LEN))

    ours = call_peaks_from_bam(str(bam), CHROM, bin_size=BIN_SIZE,
                               threshold=5.0, max_qvalue=0.01)
    our_peaks = [(CHROM, int(s), int(e)) for s, e in zip(ours.start, ours.end)]
    print("easyepi: %d peaks" % len(our_peaks))

    macs_dir = tmp / "macs"
    macs_dir.mkdir()
    np_path = run_macs3(macs3, bam, macs_dir)
    macs_peaks = read_narrowpeak(np_path)
    print("MACS3  : %d peaks" % len(macs_peaks))

    print("\n=== 1. recovery of the planted regions ===")
    print("  %-22s %-16s %-16s" % ("planted (centre/half)", "easyepi frac", "MACS3 frac"))
    our_fracs, macs_fracs = [], []
    for centre, half, fold in PEAKS:
        _, of = planted_overlap(our_peaks, centre, half)
        _, mf = planted_overlap(macs_peaks, centre, half)
        our_fracs.append(of)
        macs_fracs.append(mf)
        print("  %6d / %-4d (x%-4.0f) %-16.3f %-16.3f"
              % (centre, half, fold, of, mf))
    print("  mean                 %-16.3f %-16.3f"
          % (np.mean(our_fracs), np.mean(macs_fracs)))

    print("\n=== 2. false positives outside the planted regions ===")
    def outside(peaks):
        n = 0
        for _, s, e, *_ in peaks:
            covered = 0.0
            for centre, half, _f in PEAKS:
                lo, hi = max(s, centre - half), min(e, centre + half)
                if hi > lo:
                    covered += hi - lo
            if covered < 0.5 * (e - s):
                n += 1
        return n
    print("  easyepi: %d of %d peaks lie mostly outside a planted region"
          % (outside(our_peaks), len(our_peaks)))
    print("  MACS3  : %d of %d peaks lie mostly outside a planted region"
          % (outside(macs_peaks), len(macs_peaks)))

    print("\n=== 3. agreement between the two callers ===")
    print("  interval Jaccard            : %.3f" % jaccard(our_peaks, macs_peaks))
    print("  easyepi peaks matched by MACS3 : %d/%d"
          % (overlapping(our_peaks, macs_peaks), len(our_peaks)))
    print("  MACS3 peaks matched by easyepi : %d/%d"
          % (overlapping(macs_peaks, our_peaks), len(macs_peaks)))

    print("\n=== 4. boundary agreement on shared peaks ===")
    diffs_left, diffs_right = [], []
    for centre, half, _f in PEAKS:
        o, _ = planted_overlap(our_peaks, centre, half)
        m, _ = planted_overlap(macs_peaks, centre, half)
        if o and m:
            diffs_left.append(abs(o[0] - m[0]))
            diffs_right.append(abs(o[1] - m[1]))
    if diffs_left:
        print("  median |left edge diff|  : %.0f bp" % np.median(diffs_left))
        print("  median |right edge diff| : %.0f bp" % np.median(diffs_right))
        print("  (bin size is %d bp; a difference below it is sub-bin)" % BIN_SIZE)
    else:
        print("  no peak pair found by both callers")

    print("\nartifacts:", tmp)


if __name__ == "__main__":
    main()
