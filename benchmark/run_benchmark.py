#!/usr/bin/env python3
"""Benchmark easyepi against the tools it wraps.

Purpose: show that easyepi is orchestration over published tools, not a
reimplementation. Three comparisons, one script:

  1. peaks   ``easyepi.specific.chip.narrow_peaks`` vs ``macs3 callpeak``
  2. motifs  ``easyepi.cli motifs --backend homer``  vs ``findMotifs.pl``
  3. qc      ``easyepi.cli chip-qc`` (FRiP)          vs a direct pysam scan

Run it inside the easyepi environment:

    micromamba run -n easyepi python benchmark/run_benchmark.py

Design notes
------------
* One file on purpose: a benchmark that needs five scripts is a benchmark
  nobody runs.
* Inputs are generated when absent, so the script is reproducible from scratch.
* Comparisons use set overlap and top-N rank agreement rather than bitwise
  identity, because the wrapped tools are independently implemented.
* Every comparison is isolated: one failure does not hide the others.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

CHROM = "chr1"
SLICE_LEN = 1_000_000
N_PEAKS = 300
PEAK_WIDTH = 500
FRAGS_PER_PEAK = 40
N_BACKGROUND = 8000
MOTIF = "AGGTCAAGGTCA"  # RAR:RXR DR0, planted at every peak centre
FRAG_MIN, FRAG_MAX = 150, 400
READ_LEN = 50


# --------------------------------------------------------------------- utilities
def run(cmd, *, check=True):
    """Run a command and return (returncode, combined output)."""
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    out = (proc.stdout or "") + (proc.stderr or "")
    if check and proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(map(str, cmd))}\n{out[-1200:]}")
    return proc.returncode, out


def tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"required tool not found on PATH: {name}")
    return path


class IntervalIndex:
    """Sorted-per-chromosome interval lookup with binary search."""

    def __init__(self, intervals):
        by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in intervals:
            by_chrom[chrom].append((start, end))
        for chrom in by_chrom:
            by_chrom[chrom].sort()
        self._by_chrom = dict(by_chrom)

    def contains(self, chrom: str, pos: int) -> bool:
        arr = self._by_chrom.get(chrom)
        if not arr:
            return False
        lo, hi = 0, len(arr) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            start, end = arr[mid]
            if pos < start:
                hi = mid - 1
            elif pos >= end:
                lo = mid + 1
            else:
                return True
        return False

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        arr = self._by_chrom.get(chrom)
        if not arr:
            return False
        lo, hi = 0, len(arr) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            s, e = arr[mid]
            if end <= s:
                hi = mid - 1
            elif start >= e:
                lo = mid + 1
            else:
                return True
        return False


def bed_intervals(path: Path) -> list[tuple[str, int, int]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue
        try:
            rows.append((f[0], int(f[1]), int(f[2])))
        except ValueError:
            continue
    return rows


# --------------------------------------------------------------------- fixtures
def build_inputs(work: Path, reference: Path) -> dict:
    """Create reference slice, truth peaks, peak FASTA, background, BAM and GTF."""
    import numpy as np
    import pysam

    work.mkdir(parents=True, exist_ok=True)
    fasta = work / "chr1_1Mb.fa"
    peaks_bed = work / "truth_peaks.bed"
    peaks_fa = work / "peaks.fa"
    background = work / "background.fa"
    bam = work / "sample.bam"
    genes = work / "genes.gtf"

    if not fasta.exists():
        with pysam.FastaFile(str(reference)) as fa:
            seq = fa.fetch(CHROM, 0, SLICE_LEN).upper()
        with fasta.open("w") as fh:
            fh.write(f">{CHROM}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
        pysam.faidx(str(fasta))
    with pysam.FastaFile(str(fasta)) as fa:
        seq = fa.fetch(CHROM, 0, SLICE_LEN)
    length = len(seq)

    rng = np.random.default_rng(23)

    if not peaks_bed.exists():
        starts = np.sort(rng.choice(np.arange(1000, length - PEAK_WIDTH - 2000), size=N_PEAKS, replace=False))
        with peaks_bed.open("w") as fh:
            for i, s in enumerate(starts):
                fh.write(f"{CHROM}\t{int(s)}\t{int(s) + PEAK_WIDTH}\ttp{i}\n")
    peaks = [tuple(map(int, line.split("\t")[1:3])) for line in peaks_bed.read_text().splitlines()]

    if not peaks_fa.exists():
        half = len(MOTIF) // 2
        with peaks_fa.open("w") as fh:
            for i, (s, e) in enumerate(peaks):
                centre = (s + e) // 2
                frag = list(seq[centre - 150:centre + 150])
                frag[150 - half:150 - half + len(MOTIF)] = list(MOTIF)
                fh.write(f">peak_{i}\n{''.join(frag)}\n")

    if not background.exists():
        centres = [(s + e) // 2 for s, e in peaks]
        written = 0
        with background.open("w") as fh:
            while written < N_PEAKS:
                c = int(rng.integers(500, length - 500))
                if any(abs(c - m) < 1000 for m in centres):
                    continue
                fh.write(f">bg_{written}\n{seq[c - 150:c + 150]}\n")
                written += 1

    if not bam.exists():
        # Distinct fragment coordinates: MACS3 deduplicates identical fragments,
        # so synthetic data must not reuse them or signal collapses to nothing.
        frags = []
        for s, e in peaks:
            picks = rng.choice(np.arange(s, max(s + 2, e - FRAG_MIN)), size=FRAGS_PER_PEAK, replace=False)
            for st in picks:
                frags.append((int(st), int(st) + int(rng.integers(FRAG_MIN, FRAG_MAX))))
        for st in rng.choice(np.arange(0, length - 600), size=N_BACKGROUND, replace=False):
            frags.append((int(st), int(st) + int(rng.integers(FRAG_MIN, FRAG_MAX))))
        frags.sort()

        # Collect every record, then emit in coordinate order: read2 starts at
        # end - READ_LEN and can precede the next fragment's read1, so writing
        # them pairwise leaves the BAM unsorted and pysam.index refuses it.
        records = []
        for i, (fs, fe) in enumerate(frags):
            for flag, start, tlen, mate in (
                (99, fs, fe - fs, fe - READ_LEN),
                (147, fe - READ_LEN, -(fe - fs), fs),
            ):
                records.append((start, flag, i, tlen, mate))
        records.sort()

        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": length}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
            for start, flag, i, tlen, mate in records:
                rec = pysam.AlignedSegment(out.header)
                rec.query_name = f"f{i}"
                rec.query_sequence = "".join(rng.choice(list("ACGT"), READ_LEN))
                rec.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
                rec.flag = flag
                rec.reference_id = 0
                rec.reference_start = start
                rec.mapping_quality = 60
                rec.cigarstring = f"{READ_LEN}M"
                rec.next_reference_id = 0
                rec.next_reference_start = mate
                rec.template_length = tlen
                out.write(rec)
        pysam.index(str(bam))

    if not genes.exists():
        with genes.open("w") as fh:
            for i in range(200):
                s = int(rng.integers(1000, length - 20000))
                e = s + int(rng.integers(5000, 15000))
                fh.write(f'{CHROM}\tsim\tgene\t{s + 1}\t{e}\t.\t{"+" if i % 2 == 0 else "-"}\t.\tgene_id "G{i:04d}";\n')

    return {"fasta": fasta, "peaks_bed": peaks_bed, "peaks_fa": peaks_fa,
            "background": background, "bam": bam, "genes": genes, "length": length}


# --------------------------------------------------------------------- parsing
def homer_known_motifs(result_dir: Path) -> list[tuple[str, str, float]]:
    """Return [(name, consensus, log_p)] from HOMER knownResults, best first."""
    rows = []
    for f in sorted(result_dir.glob("known*.motif")):
        first = next((ln for ln in f.read_text().splitlines() if ln.startswith(">")), None)
        if not first:
            continue
        parts = first[1:].split("\t")
        if len(parts) < 4:
            continue
        try:
            log_p = float(parts[3])
        except ValueError:
            continue
        rows.append((parts[1], parts[0], log_p))
    rows.sort(key=lambda r: r[2])
    return rows


# --------------------------------------------------------------------- comparisons
def bench_peaks(work: Path, inp: dict, results: dict) -> None:
    easyepi_dir = work / "bench" / "peaks_easyepi"
    direct_dir = work / "bench" / "peaks_macs3"
    manifest = work / "manifest.tsv"
    manifest.write_text(f"sample_id\tassay\tbam\ns1\tChIP\t{inp['bam']}\n")

    run([sys.executable, "-m", "easyepi.specific.chip.narrow_peaks",
         "--project", manifest, "--out", easyepi_dir,
         "--genome-size", inp["length"], "--threads", "4"])
    easyepi_peaks = sorted(easyepi_dir.rglob("*_peaks.narrowPeak"))[0]

    direct_dir.mkdir(parents=True, exist_ok=True)
    run([tool("macs3"), "callpeak", "-t", inp["bam"], "-f", "BAMPE",
         "-g", inp["length"], "-n", "direct", "--outdir", direct_dir, "-q", "0.01"])
    direct_peaks = direct_dir / "direct_peaks.narrowPeak"

    a = bed_intervals(easyepi_peaks)
    b = bed_intervals(direct_peaks)
    b_idx = IntervalIndex(b)
    a_idx = IntervalIndex(a)

    a_in_b = sum(1 for c, s, e in a if b_idx.overlaps(c, s, e))
    b_in_a = sum(1 for c, s, e in b if a_idx.overlaps(c, s, e))
    identical = set(a) == set(b)

    results["peaks"] = {
        "easyepi_n": len(a),
        "macs3_n": len(b),
        "easyepi_overlapping_macs3": a_in_b,
        "macs3_overlapping_easyepi": b_in_a,
        "identical_interval_sets": identical,
        "pass": identical,
    }


def bench_motifs(work: Path, inp: dict, results: dict) -> None:
    easyepi_dir = work / "bench" / "motifs_easyepi"
    direct_dir = work / "bench" / "motifs_homer"

    run([sys.executable, "-m", "easyepi.cli", "motifs",
         "--peaks", inp["peaks_bed"], "--fasta", inp["peaks_fa"],
         "--background", inp["background"], "--backend", "homer",
         "--out", easyepi_dir])

    direct_dir.mkdir(parents=True, exist_ok=True)
    run([tool("findMotifs.pl"), inp["peaks_fa"], "fasta", direct_dir,
         "-fasta", inp["background"], "-p", "4"])

    ep_rows = homer_known_motifs(easyepi_dir / "knownResults")
    dir_rows = homer_known_motifs(direct_dir / "knownResults")

    ep_top = [r[1] for r in ep_rows[:10]]
    dir_top = [r[1] for r in dir_rows[:10]]
    shared = len(set(ep_top) & set(dir_top))
    planted_top3 = any("RAR" in r[1] or "RAR" in r[0] for r in ep_rows[:3])

    results["motifs"] = {
        "easyepi_n_known": len(ep_rows),
        "homer_cli_n_known": len(dir_rows),
        "easyepi_top10": ep_top,
        "homer_cli_top10": dir_top,
        "shared_in_top10": shared,
        "planted_RAR_in_top3": planted_top3,
        "pass": shared >= 8 and planted_top3,
    }


def bench_qc(work: Path, inp: dict, results: dict) -> None:
    import pysam

    out = work / "bench" / "qc_easyepi"
    peak_file = sorted((work / "bench" / "peaks_easyepi").rglob("*_peaks.narrowPeak"))[0]
    manifest = work / "manifest.tsv"

    run([sys.executable, "-m", "easyepi.cli", "chip-qc", "--project", manifest,
         "--peak", peak_file, "--out", out, "--max-reads", "200000"])

    payload = json.loads((out / "chip_seq_qc.json").read_text())
    samples = payload["samples"]
    row = samples[0] if isinstance(samples, list) else next(iter(samples.values()))
    easyepi_frip = (row.get("assay_qc") or {}).get("frip")

    index = IntervalIndex(bed_intervals(peak_file))
    total = in_peaks = 0
    with pysam.AlignmentFile(str(inp["bam"]), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            total += 1
            if index.contains(rec.reference_name, rec.reference_start):
                in_peaks += 1
    direct_frip = in_peaks / total if total else 0.0
    delta = abs((easyepi_frip or 0.0) - direct_frip)

    results["qc"] = {
        "easyepi_frip": easyepi_frip,
        "direct_frip": direct_frip,
        "abs_delta": delta,
        "reads_total": total,
        "reads_in_peaks": in_peaks,
        "pass": delta < 0.02,
    }


# --------------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description="Benchmark easyepi against macs3 and HOMER")
    ap.add_argument("--workdir", default="/home/houoweibo/e2e")
    ap.add_argument("--reference", default="/data1/GaoLab/houoweibo/reference/human/genome.fa")
    ap.add_argument("--skip", nargs="*", default=[], choices=["peaks", "motifs", "qc"])
    args = ap.parse_args()

    work = Path(args.workdir)
    print(f"[bench] workdir   : {work}")
    print(f"[bench] reference : {args.reference}")

    inp = build_inputs(work, Path(args.reference))
    print(f"[bench] inputs ready (bam {inp['bam'].stat().st_size // 1024} KB)")

    results: dict = {}
    for name, fn in (("peaks", bench_peaks), ("motifs", bench_motifs), ("qc", bench_qc)):
        if name in args.skip:
            print(f"[bench] {name:7s} SKIPPED")
            continue
        print(f"[bench] {name:7s} running ...")
        try:
            fn(work, inp, results)
        except Exception as exc:  # keep the other comparisons running
            results[name] = {"error": f"{type(exc).__name__}: {exc}", "pass": False}
        print(f"[bench] {name:7s} {'PASS' if results[name].get('pass') else 'FAIL'}")

    out_dir = work / "benchmark"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results.json").write_text(json.dumps(results, indent=2))
    print(f"[bench] results   : {out_dir / 'results.json'}")

    print("\n=== SUMMARY ===")
    for name in ("peaks", "motifs", "qc"):
        if name not in results:
            continue
        r = results[name]
        print(f"  {name:7s} {'PASS' if r.get('pass') else 'FAIL'}")
        for key, value in r.items():
            if key != "pass":
                print(f"          {key}: {value}")

    return 0 if all(results.get(n, {}).get("pass") for n in ("peaks", "motifs", "qc")) else 1


if __name__ == "__main__":
    raise SystemExit(main())
