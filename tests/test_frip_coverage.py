"""FRiP must not present a region-restricted value as a genome-wide QC metric."""
from __future__ import annotations

import pytest

from easyepi.assay_metrics import _frip_coverage_check


def _intervals(rows):
    out = {}
    for chrom, s, e in rows:
        out.setdefault(chrom, []).append((s, e))
    for chrom in out:
        out[chrom].sort()
    return out


def test_uniformly_gc_and_window_functions_exist():
    # guard against a refactor silently removing the check
    assert callable(_frip_coverage_check)


def test_region_limited_peaks_are_flagged(tmp_path):
    """Peaks over 5 Mb of a 3 Gb genome cannot produce a genome-wide FRiP."""
    bam = tmp_path / "tiny.bam"
    bam.write_bytes(b"")          # unreadable: read span falls back to 0
    intervals = _intervals([("chr1", 0, 5_000_000)])
    detail = _frip_coverage_check(str(bam), intervals, 1000, 10)
    # with no readable BAM the function must not claim interpretability
    assert detail["interpretable"] is None
    assert detail["peak_span_bp"] == 5_000_000
    assert detail["peak_bp"] == 5_000_000


def test_coverage_check_reports_the_peak_extent():
    intervals = _intervals([("chr1", 100, 200), ("chr1", 1000, 1200),
                            ("chr2", 5, 105)])
    detail = _frip_coverage_check("/nonexistent.bam", intervals, 100, 5)
    assert detail["peak_bp"] == 100 + 200 + 100
    # global min/max across chromosomes: 5 .. 1200
    assert detail["peak_span_bp"] == 1200 - 5


def test_coverage_check_without_intervals_is_not_interpretable():
    detail = _frip_coverage_check("/nonexistent.bam", {}, 0, 0)
    assert detail["interpretable"] is None
    assert detail["peak_bp"] == 0


def test_panel_note_is_always_present():
    detail = _frip_coverage_check("/nonexistent.bam",
                                  _intervals([("chr1", 0, 10)]), 1, 0)
    assert "same genome" in detail["panel"]
