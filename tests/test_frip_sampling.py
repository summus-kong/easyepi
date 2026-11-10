"""Regression tests for FRiP sampling correctness.

A BAM is coordinate-sorted, so capping the scan at N records samples only the
first chromosome and biases FRiP badly (47x high on a real ENCODE H3K4me3 file).
These tests pin the whole-file striding fix.
"""

from __future__ import annotations

import pytest

pytest.importorskip("pysam")

import pysam  # noqa: E402

from easyepi.assay_metrics import frip  # noqa: E402

CHROMS = [("chr1", 300_000), ("chr2", 300_000)]
# Peaks sit on the SECOND chromosome only: a head-truncated scan sees chr1 first.
PEAKS_CHR2 = [(150_000, 152_000), (200_000, 202_000)]
CHR2_INDEX = 1
READ_LEN = 50


def _records(reads_per_chrom=2000, in_peak_chr2=400):
    rows = []
    for i in range(reads_per_chrom):
        rows.append((0, 1000 + i * 100))
    for i in range(reads_per_chrom):
        rows.append((1, 1000 + i * 100))
    for i in range(in_peak_chr2):
        s, e = PEAKS_CHR2[i % len(PEAKS_CHR2)]
        rows.append((1, s + (i % max(1, e - s - READ_LEN))))
    rows.sort()
    return rows


def _build_bam(path, rows):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": n} for c, n in CHROMS]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for idx, (ref, start) in enumerate(rows):
            rec = pysam.AlignedSegment(out.header)
            rec.query_name = f"r{idx}"
            rec.query_sequence = "A" * READ_LEN
            rec.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            rec.flag = 0
            rec.reference_id = ref
            rec.reference_start = start
            rec.mapping_quality = 60
            rec.cigarstring = f"{READ_LEN}M"
            out.write(rec)
    pysam.index(str(path))


def _build_peaks(path):
    with open(path, "w", encoding="utf-8") as fh:
        for i, (s, e) in enumerate(PEAKS_CHR2):
            fh.write(f"chr2\t{s}\t{e}\tpk{i}\n")


def _expected_hits(rows):
    return sum(1 for ref, start in rows
               if ref == CHR2_INDEX and any(ps <= start < pe for ps, pe in PEAKS_CHR2))


def _make_inputs(tmp_path):
    rows = _records()
    bam, peaks = tmp_path / "s.bam", tmp_path / "p.bed"
    _build_bam(bam, rows)
    _build_peaks(peaks)
    return bam, peaks, rows


def test_frip_exhaustive_matches_expected(tmp_path):
    bam, peaks, rows = _make_inputs(tmp_path)
    res = frip(str(bam), str(peaks), max_reads=10**9)
    expected = _expected_hits(rows) / len(rows)
    assert res["stride"] == 1
    assert res["sampling"] == "exhaustive"
    assert res["reads_inspected"] == len(rows)
    assert res["reads_in_peaks"] == _expected_hits(rows)
    assert res["frip"] == pytest.approx(expected, rel=0.02)


def test_truncated_scan_does_not_collapse_to_zero(tmp_path):
    bam, peaks, _rows = _make_inputs(tmp_path)
    full = frip(str(bam), str(peaks), max_reads=10**9)["frip"]
    trunc = frip(str(bam), str(peaks), max_reads=500)
    assert trunc["sampling"] == "whole_file_stride"
    assert trunc["stride"] > 1
    assert trunc["reads_scanned"] > trunc["reads_inspected"]
    assert trunc["frip"] == pytest.approx(full, rel=0.25)
    assert trunc["frip"] > 0.01


def test_frip_reports_scan_provenance(tmp_path):
    bam, peaks, _rows = _make_inputs(tmp_path)
    res = frip(str(bam), str(peaks), max_reads=1000)
    for key in ("reads_scanned", "reads_inspected", "stride", "sampling", "frip"):
        assert key in res


def test_frip_no_peaks_returns_zero(tmp_path):
    bam, _peaks, _rows = _make_inputs(tmp_path)
    other = tmp_path / "none.bed"
    other.write_text("chr9\t100\t200\tx\n", encoding="utf-8")
    assert frip(str(bam), str(other), max_reads=10**9)["frip"] == 0.0
