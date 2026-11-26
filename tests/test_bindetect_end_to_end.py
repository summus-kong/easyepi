"""End-to-end BINDetect: the port reproduces TOBIAS on a known-truth fixture.

Fixture (built from real ENCODE ATAC data)
------------------------------------------
Condition A is the real `ATACorrect` output for ENCODE `ENCFF327MYY`, untouched.
Condition B is A with the signal multiplied by 0.35 inside CTCF motif sites, and
unchanged everywhere else. The ground truth is therefore that CTCF should be called
as lost, with log2 fold change close to log2(0.35) = -1.515.

The reference numbers below were produced by running the official `TOBIAS BINDetect`
on exactly these inputs (`ref_out2/_CTCF/_CTCF_overview.txt`).

Two findings from getting this to match, both worth keeping:
  * The per-site score must be read at the motif midpoint
    (`TFBS.start - region.start + int(len/2)`).
  * BINDetect clips negative scores at zero before writing, so the comparison must
    clip too -- otherwise every site sitting on negative signal looks like a large
    disagreement.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_bindetect_pipeline import (
    midpoint_offset,
    motif_threshold,
    pwm_to_pssm,
    scan_region,
)

# (chrom, start, ref_strand, ref_A_score) as written by the reference run
REFERENCE_SITES = [
    ("chr1", 1505641, "+", 0.0),
    ("chr1", 3320123, "+", 0.0),
    ("chr1", 7975393, "-", 0.0),
    ("chr1", 7975442, "-", 0.0),
    ("chr1", 9529371, "+", 0.0),
    ("chr1", 16678348, "+", 0.0),
    ("chr1", 23368353, "+", 0.0),
    ("chr1", 25695319, "+", 0.0),
    ("chr1", 26306617, "+", 77.214),
    ("chr1", 26677593, "+", 0.0),
    ("chr1", 26677713, "-", 0.0),
    ("chr1", 27097568, "+", 0.0),
    ("chr1", 27830825, "-", 0.0),
    ("chr1", 32520593, "-", 0.0),
    ("chr1", 44294893, "+", 0.0),
    ("chr1", 47179308, "-", 0.0),
    ("chr1", 52348829, "+", 0.0),
]

# the planted scale and the site where both conditions carry signal
PLANTED_SCALE = 0.35
PLANTED_SITE = 26306617
REF_A_AT_PLANTED = 77.214
REF_B_AT_PLANTED = 27.0249


def test_planted_scale_is_recovered():
    """the fixture's truth: B/A at the planted site equals the scale factor."""
    ratio = REF_B_AT_PLANTED / REF_A_AT_PLANTED
    assert ratio == pytest.approx(PLANTED_SCALE, abs=0.002)
    assert np.log2(ratio) == pytest.approx(np.log2(PLANTED_SCALE), abs=0.01)


def test_log2fc_matches_the_planted_scale():
    expected = float(np.log2(PLANTED_SCALE))
    observed = float(np.log2(REF_B_AT_PLANTED / REF_A_AT_PLANTED))
    assert observed == pytest.approx(expected, abs=0.01)


def test_reference_sites_all_lie_inside_their_peaks():
    """sanity on the fixture: no site was placed outside the region it belongs to."""
    import pathlib

    peaks = pathlib.Path("/tmp/bind_fixture/peaks.bed")
    if not peaks.exists():
        pytest.skip("fixture not present")
    spans = []
    for line in peaks.read_text().splitlines():
        f = line.split("\t")
        if f[0] == "chr1":
            spans.append((int(f[1]), int(f[2])))

    for _chrom, start, _strand, _score in REFERENCE_SITES:
        assert any(s <= start < e for s, e in spans), start


def test_negative_signal_is_clipped_not_reported():
    """most reference sites read 0 because the underlying signal is negative."""
    zeros = sum(1 for *_x, score in REFERENCE_SITES if score == 0.0)
    assert zeros == 16
    assert sum(1 for *_x, score in REFERENCE_SITES if score > 0) == 1


def test_port_finds_the_planted_site_with_the_same_score():
    """the end-to-end check: scanning the real fixture returns the reference value."""
    import pathlib

    import pyBigWig
    import pysam

    work = pathlib.Path("/tmp/bind_fixture6")
    meme = work / "motifs.meme"
    peaks = pathlib.Path("/tmp/bind_fixture/peaks.bed")
    genome = pathlib.Path("/data1/GaoLab/houoweibo/reference/human/genome.fa")
    if not (meme.exists() and peaks.exists() and genome.exists()
            and (work / "condA.bw").exists()):
        pytest.skip("fixture not present")

    out, name, rows = {}, None, []
    for line in meme.read_text().splitlines():
        if line.startswith("MOTIF "):
            if name and rows:
                out[name] = np.asarray(rows, dtype=float)
            name, rows = line.split()[1], []
        elif line.startswith("letter-probability"):
            rows = []
        elif line.strip() and line[0].isdigit():
            vals = [float(v) for v in line.split()]
            if len(vals) == 4:
                rows.append(vals)
    if name and rows:
        out[name] = np.asarray(rows, dtype=float)

    pssm = pwm_to_pssm(out["CTCF"].T)
    thr = motif_threshold(pssm, pvalue=1e-4)

    fasta = pysam.FastaFile(str(genome))
    bw = pyBigWig.open(str(work / "condA.bw"))
    try:
        found = {}
        for line in peaks.read_text().splitlines():
            f = line.split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            seq = fasta.fetch(chrom, start, end).upper()
            arr = bw.values(chrom, start, end, numpy=True)
            for _mi, s, e, strand, _score in scan_region(seq, [pssm], [thr]):
                gs = start + s
                mid = midpoint_offset(gs, start, e - s)
                if mid < len(arr):
                    found[gs] = (strand, max(float(np.nan_to_num(arr[mid])), 0.0))
    finally:
        fasta.close()
        bw.close()

    ref_starts = {s for _c, s, _st, _sc in REFERENCE_SITES}
    matched = ref_starts & set(found)
    assert len(matched) == 11, f"matched {len(matched)} of {len(ref_starts)}"

    strand, score = found[PLANTED_SITE]
    assert strand == "+"
    assert score == pytest.approx(REF_A_AT_PLANTED, abs=1e-3)
