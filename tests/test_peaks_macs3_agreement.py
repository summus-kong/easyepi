"""Agreement with MACS3 on the peaks both callers are meant to find.

The harness that produces the full table is `tests/oracle/compare_to_macs3.py`.
This module asserts the part that must hold every run, so a regression in the
peak caller shows up as a failing test rather than as a number nobody re-reads.

The two callers use different models -- easyepi takes the median of the whole
binned signal as background, MACS3 builds a local lambda -- so identical peak
sets is not the expectation and is not asserted. What is asserted is that they
locate the same accessible regions and that neither invents regions the fixture
does not contain.

Skips cleanly when macs3 is absent, matching how `tests/oracle/` handles R.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("pysam")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from easyepi.peaks import call_peaks_from_bam  # noqa: E402

MACS3_CANDIDATES = (
    "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi/bin/macs3",
    "macs3",
)


def _macs3() -> str:
    for cand in MACS3_CANDIDATES:
        found = cand if Path(cand).is_file() else shutil.which(cand)
        if found and Path(found).exists():
            return found
    pytest.skip("macs3 not installed; the MACS3 oracle needs it")


CHROM = "chr1"
CHROM_LEN = 200_000
BIN_SIZE = 200
PEAKS = [(30_000, 500, 20.0), (90_000, 400, 25.0), (150_000, 600, 12.0)]
BACKGROUND_READS = 12_000
READ_LEN = 50


def _fixture(tmp: Path) -> Path:
    import pysam

    rng = np.random.default_rng(99)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    in_peak = np.zeros(CHROM_LEN, dtype=bool)
    for centre, half, _f in PEAKS:
        in_peak[max(0, centre - half):min(CHROM_LEN, centre + half)] = True

    starts = list(rng.choice(np.flatnonzero(~in_peak), size=BACKGROUND_READS,
                             replace=False))
    for centre, half, fold in PEAKS:
        n = max(50, int(BACKGROUND_READS * fold / 25.0))
        starts.extend(rng.integers(max(0, centre - half),
                                   min(CHROM_LEN - READ_LEN, centre + half), n).tolist())

    # Collect first, then write in coordinate order: a BAM must be sorted before
    # it can be indexed, and emitting both mates per iteration does not do that.
    recs = []
    for i, start in enumerate(starts):
        frag = 200
        for flag, rs in ((99, start), (147, min(CHROM_LEN - 1, start + frag - READ_LEN))):
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
            a.template_length = frag if flag == 99 else -frag
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            recs.append(a)

    bam = tmp / "atac.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as fh:
        for r in sorted(recs, key=lambda x: (x.reference_start, x.flag)):
            fh.write(r)
    pysam.index(str(bam))
    return bam


def _macs3_peaks(bam: Path, tmp: Path):
    outdir = tmp / "macs"
    outdir.mkdir()
    cmd = [_macs3(), "callpeak", "-t", str(bam), "-f", "BAMPE", "-g", "hs",
           "-n", "m", "--outdir", str(outdir), "-q", "0.01",
           "--nomodel", "--extsize", "200"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert proc.returncode == 0, proc.stdout.decode()[-2000:]
    path = outdir / "m_peaks.narrowPeak"
    assert path.exists(), "MACS3 wrote no narrowPeak"
    peaks = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track")):
            continue
        p = line.split("\t")
        peaks.append((int(p[1]), int(p[2])))
    return peaks


def _covers(peaks, centre, half, frac=0.5):
    lo, hi = centre - half, centre + half
    for s, e in peaks:
        ilo, ihi = max(lo, s), min(hi, e)
        if ihi > ilo and (ihi - ilo) >= frac * (hi - lo):
            return True
    return False


@pytest.fixture(scope="module")
def both_callers():
    tmp = Path(tempfile.mkdtemp(prefix="macs3t_"))
    bam = _fixture(tmp)
    ours = call_peaks_from_bam(str(bam), CHROM, bin_size=BIN_SIZE,
                               threshold=5.0, max_qvalue=0.01)
    our_peaks = [(int(s), int(e)) for s, e in zip(ours.start, ours.end)]
    return our_peaks, _macs3_peaks(bam, tmp)


def test_every_planted_peak_is_found_by_both(both_callers):
    our_peaks, macs_peaks = both_callers
    for centre, half, _fold in PEAKS:
        assert _covers(our_peaks, centre, half), (
            "easyepi missed the planted peak at %d" % centre)
        assert _covers(macs_peaks, centre, half), (
            "MACS3 missed the planted peak at %d" % centre)


def test_neither_caller_invents_regions(both_callers):
    """A peak must lie mostly on a planted region, not in the background."""
    our_peaks, macs_peaks = both_callers

    def mostly_outside(peaks):
        bad = []
        for s, e in peaks:
            covered = 0.0
            for centre, half, _f in PEAKS:
                lo, hi = max(s, centre - half), min(e, centre + half)
                if hi > lo:
                    covered += hi - lo
            if covered < 0.5 * (e - s):
                bad.append((s, e))
        return bad

    assert not mostly_outside(our_peaks), mostly_outside(our_peaks)
    assert len(mostly_outside(macs_peaks)) <= 1, mostly_outside(macs_peaks)


def test_the_two_callers_overlap_positions(both_callers):
    """Different models, same accessible regions: check the positional overlap."""
    our_peaks, macs_peaks = both_callers
    matched = 0
    for s1, e1 in our_peaks:
        for s2, e2 in macs_peaks:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi > lo and (hi - lo) >= 0.5 * min(e1 - s1, e2 - s2):
                matched += 1
                break
    frac = matched / max(1, len(our_peaks))
    assert frac >= 0.8, "only %d/%d easyepi peaks have a MACS3 partner" % (
        matched, len(our_peaks))


def test_boundaries_agree_within_a_bin(both_callers):
    """Edges cannot agree better than the bin width the caller works in."""
    our_peaks, macs_peaks = both_callers
    diffs = []
    for centre, half, _f in PEAKS:
        o = next(((s, e) for s, e in our_peaks
                  if max(s, centre - half) < min(e, centre + half)), None)
        m = next(((s, e) for s, e in macs_peaks
                  if max(s, centre - half) < min(e, centre + half)), None)
        if o and m:
            diffs.append(abs(o[0] - m[0]))
    assert diffs, "no peak pair was found by both callers"
    assert np.median(diffs) <= 2 * BIN_SIZE, (
        "median left-edge difference %.0f bp exceeds 2 bins" % np.median(diffs))
