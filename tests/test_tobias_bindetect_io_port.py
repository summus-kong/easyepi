"""BINDetect's sampling and output schema, checked against the reference rules.

Covers what the previous test file left out: reading footprint values at the TFBS
midpoint, the deterministic background sample, the per-site BED column order, and
the per-TF summary numbers.
"""
from __future__ import annotations

import random

import numpy as np
import pytest

from easyepi.specific.atac.tobias_bindetect_pipeline import (
    RANDOM_WINDOW,
    background_positions,
    midpoint_offset,
    observed_statistics,
    sample_signal,
    score_tfbs,
    tfbs_bed_row,
)


class FakeBigWig:
    """Minimal stand-in exposing the same ``values`` contract as pyBigWig."""

    def __init__(self, fill=np.nan, every=7, value=2.0):
        self.fill, self.every, self.value = fill, every, value

    def values(self, chrom, start, end, numpy=True):
        out = np.full(end - start, self.fill, dtype=float)
        out[::self.every] = self.value
        return out


# --- midpoint ---------------------------------------------------------------

def test_midpoint_matches_reference_formula():
    for tfbs_start, region_start, mlen in [(1000, 1000, 10), (1500, 1000, 20),
                                           (1234, 1000, 8)]:
        expected = tfbs_start - region_start + int(mlen / 2)
        assert midpoint_offset(tfbs_start, region_start, mlen) == expected


def test_score_tfbs_indexes_region_relative_array():
    """Hits are region-relative, so the array is indexed without subtracting again.

    Before the fix the region start was subtracted a second time, reading every
    score from the wrong position.
    """
    values = np.arange(500, dtype=float) * 10
    hits = [(0, 100, 110, "+", 0.0)]
    row = score_tfbs(hits, [10], {"c1": values})[0]

    assert row["mid"] == 105
    assert row["c1"] == pytest.approx(1050.0)


def test_score_tfbs_clamps_to_array_bounds():
    values = np.arange(50, dtype=float)
    hits = [(0, 48, 58, "+", 0.0)]      # midpoint would run past the end
    row = score_tfbs(hits, [10], {"c1": values})[0]
    assert row["c1"] == pytest.approx(49.0)


# --- background sampling ----------------------------------------------------

def test_background_positions_match_reference_sampling():
    mine = background_positions(4000)
    random.seed(4000)
    ref = random.sample(range(4000), max(1, int(4000 / RANDOM_WINDOW)))
    assert mine == ref


def test_background_positions_are_order_independent():
    """Seeded by region length, so interleaving other regions changes nothing."""
    first = background_positions(4000)
    _ = background_positions(1234)
    second = background_positions(4000)
    assert first == second


def test_background_positions_always_sample_at_least_one():
    assert len(background_positions(10)) == 1


# --- signal sampling --------------------------------------------------------

def test_sample_signal_keeps_nan_for_missing_values():
    """Missing bases must stay NaN; treating them as zero would deflate means."""
    arr = sample_signal(FakeBigWig(), "chr1", 0, 700)
    assert len(arr) == 700
    assert np.isnan(arr).any()
    assert np.nanmean(arr) == pytest.approx(2.0)


def test_sample_signal_rejects_empty_interval():
    assert sample_signal(FakeBigWig(), "chr1", 100, 100).size == 0


# --- output schema ----------------------------------------------------------

def test_bed_row_column_order():
    row = {"chrom": "chr1", "start": 100, "end": 110, "strand": "+",
           "c1_score": 1.23456, "c2_score": 0.5,
           "c1_bound": 1, "c2_bound": 0, "c1_c2_log2fc": -0.75}
    fields = tfbs_bed_row(row, ["c1", "c2"], "CTCF").split("\t")

    assert len(fields) == 11
    assert fields[:6] == ["chr1", "100", "110", "CTCF_.", ".", "+"]
    assert fields[6] == "1.235"        # three decimals, rounded
    assert fields[8:10] == ["1", "0"]  # bound flags
    assert fields[10] == "-0.750"


def test_observed_statistics_reduce_the_table():
    rows = [{"c1_score": 1.0, "c2_score": 0.5, "c1_c2_log2fc": 1.0},
            {"c1_score": 3.0, "c2_score": 0.5, "c1_c2_log2fc": -1.0}]
    stats = observed_statistics(rows, ["c1", "c2"], ("c1", "c2"))

    assert stats["c1_mean_score"] == pytest.approx(2.0)
    assert stats["c2_mean_score"] == pytest.approx(0.5)
    assert stats["c1_c2_mean_log2fc"] == pytest.approx(0.0)
    assert stats["n_sites"] == 2
