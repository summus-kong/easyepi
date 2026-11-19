"""Methylation read directly from BAM: tag detection, parsing, beta aggregation."""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.methylation_bam import (XM_CODES, MethCall, calls_to_dataframe,
                                     cpg_level_methylation,
                                     has_methylation_tags, read_methylation_calls)

pysam = pytest.importorskip("pysam", reason="pysam is required to build test BAMs")


def _bam(path, records, chrom="chr1", length=1000):
    """Write a tiny BAM.

    ``seq`` must use A/C/G/T/N: pysam normalises other characters to ``N`` on
    write, so a fixture written with "." would be read back as "N" and the XM
    offsets would no longer line up with what the test assumed.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": chrom, "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, spec in enumerate(sorted(records, key=lambda r: r["start"])):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = spec["seq"]
            a.flag = spec.get("flag", 0)
            a.reference_id = 0
            a.reference_start = spec["start"]
            a.mapping_quality = spec.get("mapq", 60)
            a.cigartuples = [(0, len(spec["seq"]))]
            a.query_qualities = pysam.qualitystring_to_array("I" * len(spec["seq"]))
            for tag, value in spec.get("tags", {}).items():
                a.set_tag(tag, value)
            out.write(a)
    pysam.index(str(path))
    return path


# --------------------------------------------------------------------------
# Tag detection
# --------------------------------------------------------------------------
def test_detects_xm_tag(tmp_path):
    bam = _bam(tmp_path / "x.bam", [
        {"start": 100, "seq": "CNNA", "tags": {"XM": "Z..h"}}])
    probe = has_methylation_tags(bam)
    assert probe["has_xm"] is True and probe["has_mm"] is False


def test_reports_absence_of_any_methylation_tag(tmp_path):
    bam = _bam(tmp_path / "plain.bam", [{"start": 100, "seq": "CGAA"}])
    probe = has_methylation_tags(bam)
    assert probe["has_xm"] is False and probe["has_mm"] is False


def test_missing_tags_raise_rather_than_implying_zero_methylation(tmp_path):
    bam = _bam(tmp_path / "plain.bam", [{"start": 100, "seq": "CGAA"}])
    with pytest.raises(ValueError, match="no methylation calls"):
        read_methylation_calls(bam)


# --------------------------------------------------------------------------
# XM parsing
# --------------------------------------------------------------------------
def test_xm_codes_cover_the_documented_alphabet():
    # every documented character must map to a (modified, context) pair
    assert XM_CODES["Z"] == (True, "CpG")
    assert XM_CODES["z"] == (False, "CpG")
    assert XM_CODES["X"] == (True, "CHG") and XM_CODES["x"] == (False, "CHG")
    assert XM_CODES["H"] == (True, "CHH") and XM_CODES["h"] == (False, "CHH")


def test_xm_forward_read_yields_calls_at_reference_positions(tmp_path):
    """Bismark writes a code only at cytosines; every other base is '.'.

    Sequence CGAC therefore carries codes at offsets 0 and 3 only, so the read
    yields exactly two CpG calls: position 100 methylated, position 103 unmethylated.
    """
    bam = _bam(tmp_path / "f.bam", [
        {"start": 100, "seq": "CNAC", "tags": {"XM": "Z..z"}}])
    calls, stats = read_methylation_calls(bam, contexts=("CpG", "CHG", "CHH"))
    assert stats["encoding"] == "XM"
    assert stats["reads_used"] == 1
    by_pos = {(c.pos, c.context, c.methylated) for c in calls}
    assert by_pos == {(100, "CpG", True), (103, "CpG", False)}


def test_xm_reverse_read_is_labelled_minus(tmp_path):
    bam = _bam(tmp_path / "r.bam", [
        {"start": 100, "seq": "CNAC", "flag": 16, "tags": {"XM": "z..h"}}])
    calls, _ = read_methylation_calls(bam, contexts=("CpG", "CHH"))
    assert calls, "reverse-strand calls were dropped"
    assert all(c.strand == "-" for c in calls)


def test_context_filter_is_applied(tmp_path):
    # one CpG cytosine (offset 0) and one CHH cytosine (offset 3)
    bam = _bam(tmp_path / "c.bam", [
        {"start": 100, "seq": "CNNC", "tags": {"XM": "Z..h"}}])
    only_cpg, _ = read_methylation_calls(bam, contexts=("CpG",))
    all_ctx, _ = read_methylation_calls(bam, contexts=("CpG", "CHH"))
    assert all(c.context == "CpG" for c in only_cpg)
    assert len(all_ctx) >= len(only_cpg)


def test_unmapped_and_low_mapq_records_are_skipped(tmp_path):
    bam = _bam(tmp_path / "s.bam", [
        {"start": 100, "seq": "CNNA", "tags": {"XM": "Z..h"}},
        {"start": 300, "seq": "CNNA", "mapq": 0, "tags": {"XM": "Z..h"}}])
    calls, stats = read_methylation_calls(bam, min_mapq=30)
    assert stats["reads_seen"] == 2
    assert stats["reads_used"] == 1
    assert all(c.pos < 300 for c in calls)


def test_duplicates_are_skipped_by_default(tmp_path):
    bam = _bam(tmp_path / "d.bam", [
        {"start": 100, "seq": "CNNA", "flag": 1024, "tags": {"XM": "Z..h"}}])
    _, stats = read_methylation_calls(bam)
    assert stats["reads_used"] == 0
    _, stats2 = read_methylation_calls(bam, skip_duplicates=False)
    assert stats2["reads_used"] == 1


def test_unknown_chromosome_is_reported(tmp_path):
    bam = _bam(tmp_path / "x.bam", [{"start": 10, "seq": "CG", "tags": {"XM": "ZZ"}}])
    with pytest.raises(ValueError, match="not in the BAM header"):
        read_methylation_calls(bam, chrom="chrZ")


# --------------------------------------------------------------------------
# Tabulation and beta
# --------------------------------------------------------------------------
def test_calls_to_dataframe_has_bed_like_columns():
    calls = [MethCall("chr1", 10, "+", "CpG", True),
             MethCall("chr1", 10, "-", "CpG", False)]
    frame = calls_to_dataframe(calls, {"chr1": 1000})
    assert list(frame.columns) == ["chrom", "start", "end", "strand", "context",
                                   "methylated"]
    assert frame["end"].iloc[0] == frame["start"].iloc[0] + 1


def test_empty_calls_produce_an_empty_frame_with_columns():
    frame = calls_to_dataframe([], {"chr1": 100})
    assert frame.empty
    assert "methylated" in frame.columns


def test_beta_is_methylated_over_coverage():
    calls = [MethCall("chr1", 10, "+", "CpG", True),
             MethCall("chr1", 10, "-", "CpG", True),
             MethCall("chr1", 20, "+", "CpG", False),
             MethCall("chr1", 20, "-", "CpG", False)]
    frame = calls_to_dataframe(calls, {"chr1": 1000})
    sites, stats = cpg_level_methylation(frame, min_coverage=1)
    by_pos = dict(zip(sites["start"], sites["beta"]))
    assert by_pos[10] == pytest.approx(1.0)
    assert by_pos[20] == pytest.approx(0.0)


def test_min_coverage_drops_thin_sites_and_reports_the_count():
    calls = [MethCall("chr1", 10, "+", "CpG", True),
             MethCall("chr1", 20, "+", "CpG", False),
             MethCall("chr1", 20, "-", "CpG", False),
             MethCall("chr1", 20, "+", "CpG", False)]
    frame = calls_to_dataframe(calls, {"chr1": 1000})
    sites, stats = cpg_level_methylation(frame, min_coverage=2)
    assert stats["dropped_below_coverage"] == 1
    assert 10 not in set(sites["start"])
    assert 20 in set(sites["start"])


def test_single_strand_sites_are_counted_not_hidden():
    calls = [MethCall("chr1", 10, "+", "CpG", True),
             MethCall("chr1", 10, "+", "CpG", False),
             MethCall("chr1", 20, "+", "CpG", True),
             MethCall("chr1", 20, "-", "CpG", True)]
    frame = calls_to_dataframe(calls, {"chr1": 1000})
    _sites, stats = cpg_level_methylation(frame, min_coverage=1)
    assert stats["sites_single_strand"] == 1


def test_empty_frame_returns_empty_beta():
    import pandas as pd
    sites, stats = cpg_level_methylation(pd.DataFrame(
        columns=["chrom", "start", "end", "strand", "context", "methylated"]))
    assert sites.empty and stats["sites"] == 0


# --------------------------------------------------------------------------
# Large-call behaviour
# --------------------------------------------------------------------------
def test_read_count_cap_is_honoured(tmp_path):
    recs = [{"start": 100 + 10 * i, "seq": "CNNA", "tags": {"XM": "Z..h"}}
            for i in range(20)]
    bam = _bam(tmp_path / "many.bam", recs, length=1000)
    _, stats = read_methylation_calls(bam, max_reads=5)
    assert stats["reads_used"] == 5
