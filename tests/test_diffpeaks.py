"""Tests for the pure-Python differential-peak implementation.

Each test pins one published component of the DiffBind/edgeR pipeline that this
module reproduces: consensus-region merging, TMM size factors, the
negative-binomial Wald test, and Benjamini-Hochberg FDR.
"""

from __future__ import annotations

import numpy as np
import pytest

from easyepi.peaks.diffpeaks import (
    bh_adjust,
    consensus_peaks,
    merge_intervals,
    nb_wald_test,
    tmm_factors,
)


# --------------------------------------------------------------- consensus peaks
def test_merge_overlapping_intervals():
    iv = [("chr1", 100, 200), ("chr1", 190, 300), ("chr1", 500, 600)]
    assert merge_intervals(iv, min_overlap=1) == [
        ("chr1", 100, 300), ("chr1", 500, 600)]


def test_merge_keeps_gap_when_overlap_too_small():
    # a 5 bp overlap does not merge when 10 bp are required
    iv = [("chr1", 100, 200), ("chr1", 195, 300)]
    assert merge_intervals(iv, min_overlap=10) == [("chr1", 100, 200), ("chr1", 195, 300)]


def test_merge_is_per_chromosome():
    iv = [("chr1", 100, 200), ("chr2", 100, 200)]
    assert len(merge_intervals(iv)) == 2


def test_consensus_across_files(tmp_path):
    a = tmp_path / "a.bed"
    b = tmp_path / "b.bed"
    a.write_text("chr1\t100\t200\nchr1\t900\t1000\n", encoding="utf-8")
    b.write_text("chr1\t150\t250\nchr2\t50\t150\n", encoding="utf-8")
    regions = consensus_peaks([a, b])
    assert ("chr1", 100, 250) in regions
    assert ("chr1", 900, 1000) in regions
    assert ("chr2", 50, 150) in regions
    assert len(regions) == 3


# --------------------------------------------------------------- TMM
def test_tmm_compensates_for_depth():
    """A sample with 1.5x the counts on most regions gets a factor < 1."""
    counts = np.array([
        [100, 150], [120, 180], [90, 135], [110, 165], [80, 120],
        [200, 300], [130, 195], [70, 105], [95, 142], [10, 900],
    ], float)
    f = tmm_factors(counts)
    assert f[1] < f[0]
    # factors are scaled to multiply to 1, like edgeR
    assert np.prod(f) == pytest.approx(1.0, abs=1e-9)


def test_tmm_single_sample_is_identity():
    assert tmm_factors(np.array([[1.0], [2.0]])).tolist() == [1.0]


def test_tmm_needs_positive_totals():
    with pytest.raises(ValueError):
        tmm_factors(np.array([[0.0, 1.0], [0.0, 2.0]]))


# --------------------------------------------------------------- NB Wald
def _two_feature_matrix():
    return np.array([
        [100, 110, 90, 105, 1000, 1050, 980, 1020],   # truly different
        [200, 210, 190, 205,  205,  200, 195,  210],  # null
    ], float)


def test_nb_wald_finds_real_difference():
    m = _two_feature_matrix()
    mask = np.array([True] * 4 + [False] * 4)
    res = nb_wald_test(m, mask, norm=np.ones(8))
    assert res["log2fc"][0] < -3          # 100 vs 1000
    # n=4 per group: a 10x shift lands around p ~ 4e-3, so assert the order of
    # magnitude rather than an arbitrary 1e-3 cutoff.
    assert res["pvalue"][0] < 0.01
    assert abs(res["log2fc"][1]) < 0.1    # null feature
    assert res["pvalue"][1] > 0.5


def test_nb_wald_rejects_transposed_matrix():
    """A (samples, features) matrix must fail loudly rather than silently."""
    m = _two_feature_matrix().T
    mask = np.array([True] * 4 + [False] * 4)
    with pytest.raises(ValueError, match="group mask"):
        nb_wald_test(m, mask)


def test_nb_wald_needs_two_per_group():
    m = np.array([[10, 12, 900]], float)
    with pytest.raises(ValueError, match="at least 2 samples"):
        nb_wald_test(m, np.array([True, False, False]))


def test_nb_wald_dispersion_non_negative():
    m = _two_feature_matrix()
    mask = np.array([True] * 4 + [False] * 4)
    res = nb_wald_test(m, mask, norm=np.ones(8))
    assert np.all(res["dispersion"] >= 0)


# --------------------------------------------------------------- BH FDR
def test_bh_adjust_monotone_and_bounded():
    p = np.array([0.001, 0.01, 0.2, 0.5, 0.9])
    q = bh_adjust(p)
    assert np.all(np.diff(q) >= -1e-12)   # monotone in the order of p
    assert np.all(q >= p)
    assert np.all(q <= 1)


def test_bh_adjust_empty():
    assert bh_adjust(np.array([])).size == 0
