"""TSS enrichment: the arithmetic, the guards, and a check that it can detect signal.

The critical property is that the score responds to TSS-proximal reads. A function that
returns background everywhere would still produce a number; only a planted-signal test
distinguishes a working implementation from that. Real data is used in VERIFICATION.md;
here the BAMs are synthetic so the expected answer is exact.
"""
from __future__ import annotations

import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from easyepi.tss_enrichment import read_tss_table, tss_enrichment  # noqa: E402

HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"},
          "SQ": [{"SN": "chr1", "LN": 1_000_000}]}


def _write_bam(path, starts, *, chrom="chr1", reverse=False):
    """One single-base alignment starting at each position in ``starts``."""
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as out:
        for i, s in enumerate(starts):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A"
            a.flag = 16 if reverse else 0
            a.reference_id = 0
            a.reference_start = int(s)
            a.mapping_quality = 60
            a.cigarstring = "1M"
            a.query_qualities = pysam.qualitystring_to_array("I")
            out.write(a)
    return path


TSS = [("chr1", 500_000, "+")]


def _with_background(pile, flank=2000):
    """Add sparse uniform reads across the window.

    A real library has reads everywhere, so the flanks are non-zero. A synthetic BAM
    containing only TSS-proximal reads has a background of exactly 0, which makes the
    ratio undefined -- and that is the correct outcome, not a score of infinity. The
    background here is sparse relative to the pile so the expected ratio stays large.
    """
    span = range(500_000 - flank, 500_000 + flank, 200)
    return list(pile) + list(span)


def test_signal_at_the_tss_scores_high(tmp_path):
    """Reads piled on the TSS must give a large score; this is the core claim."""
    starts = _with_background(list(range(499_960, 500_041)) * 20)
    bam = _write_bam(tmp_path / "tss.bam", starts)
    res = tss_enrichment(bam, TSS, flank=2000, center=50, flank_bg=100)
    assert res["tss_enrichment_score"] > 10, res
    assert res["tss_count"] == 1
    assert res["reads_used"] == len(starts)


def test_flat_coverage_scores_near_one(tmp_path):
    """Uniform reads across the window carry no enrichment, so the ratio is ~1.

    This is the negative control: without it, a score above 1 would prove nothing.
    """
    # Uniform across the ENTIRE +/-2000 window, flanks included. Covering only the
    # middle leaves the outermost flank_bg bins empty, and the score is then
    # correctly undefined -- that is a test-data fault, not signal.
    starts = list(range(498_000, 502_001))
    bam = _write_bam(tmp_path / "flat.bam", starts)
    res = tss_enrichment(bam, TSS, flank=2000, center=50, flank_bg=100)
    assert 0.5 < res["tss_enrichment_score"] < 2.0, res


def test_the_profile_peaks_at_the_window_centre(tmp_path):
    starts = _with_background(list(range(499_960, 500_041)) * 10)
    bam = _write_bam(tmp_path / "peak.bam", starts)
    res = tss_enrichment(bam, TSS, flank=2000, center=50, flank_bg=100)
    prof = res["profile"]
    assert len(prof) == 4001
    assert abs(max(range(len(prof)), key=lambda i: prof[i]) - 2000) <= 50


def test_minus_strand_reads_contribute_their_five_prime_end(tmp_path):
    """For a reverse alignment the 5' end is the alignment end, not the start.

    Getting this wrong shifts half the reads by their length, which broadens the
    profile and lowers the score without any error message.
    """
    # Reverse reads whose *end* sits on the TSS.
    with pysam.AlignmentFile(str(tmp_path / "rev.bam"), "wb", header=HEADER) as out:
        for i in range(200):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * 1
            a.flag = 16
            a.reference_id = 0
            a.reference_start = 499_999          # end = 500_000
            a.mapping_quality = 60
            a.cigarstring = "1M"
            a.query_qualities = pysam.qualitystring_to_array("I")
            out.write(a)
        # Sparse uniform background, without which the ratio is undefined.
        for i, s in enumerate(range(498_000, 502_000, 200)):
            a = pysam.AlignedSegment()
            a.query_name = f"bg{i}"
            a.query_sequence = "A"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = s
            a.mapping_quality = 60
            a.cigarstring = "1M"
            a.query_qualities = pysam.qualitystring_to_array("I")
            out.write(a)
    res = tss_enrichment(tmp_path / "rev.bam", TSS, flank=2000, center=50,
                         flank_bg=100)
    prof = res["profile"]
    assert abs(max(range(len(prof)), key=lambda i: prof[i]) - 2000) <= 50


def test_an_impossible_window_is_refused():
    """A background overlapping the signal cannot be computed, so it is an error."""
    with pytest.raises(ValueError, match="does not fit"):
        tss_enrichment("unused.bam", TSS, flank=2000, center=100, flank_bg=2000)


def test_a_zero_center_is_refused():
    with pytest.raises(ValueError, match="center must be"):
        tss_enrichment("unused.bam", TSS, flank=2000, center=0)


def test_an_empty_library_reports_no_score(tmp_path):
    """No reads is a finding, not a score of zero or one."""
    bam = _write_bam(tmp_path / "empty.bam", [])
    res = tss_enrichment(bam, TSS, flank=2000)
    assert "error" in res and res.get("tss_enrichment_score") is None


def test_mitochondrial_tss_are_ignored(tmp_path):
    bam = _write_bam(tmp_path / "m.bam", [1_000])
    rows = [("chrM", 1_000, "+")]
    res = tss_enrichment(bam, rows, flank=2000)
    assert "error" in res, "a chrM-only TSS table should leave nothing to measure"


def test_tss_table_requires_the_named_columns(tmp_path):
    p = tmp_path / "tss.tsv"
    p.write_text("chrom\ttss\tstrand\nchr1\t500000\t+\n", encoding="utf-8")
    rows = read_tss_table(str(p))
    assert rows == [("chr1", 500_000, "+")]


def test_tss_table_without_a_header_is_read_positionally(tmp_path):
    """A headerless file must not lose its first row to header detection."""
    p = tmp_path / "tss.tsv"
    p.write_text("chr1\t500000\t+\nchr2\t123\t-\n", encoding="utf-8")
    rows = read_tss_table(str(p))
    assert rows == [("chr1", 500_000, "+"), ("chr2", 123, "-")]


def test_an_unreadable_tss_table_names_the_line(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("chrom\ttss\nchr1\tnot_a_number\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cannot read a TSS"):
        read_tss_table(str(p))
