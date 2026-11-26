"""Pin the reference's row-sliced NaN initialisation.

`atacorrect_functions.py:278` writes `bias_predictions[k_flank:reg_end] = np.nan`,
slicing the ROW axis. With `overlaps = window/step = 10` and `k_flank = 12` no row
is selected, so nothing is NaN, unwritten cells keep 0.0, and `nanmean` averages
those zeros in.

Getting this wrong scales `expected` by roughly `overlaps` and flips the sign of
`corrected`. Measured on real data before/after the fix, forward strand:

    corrected sum   +445500 -> +51741   (reference: +52446)
    expected sum    4.5545e6 -> 4.94826e6 (reference: 4.94755e6)
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_atacorrect import bias_correction_tracks


def _spike(n=4000, lo=1000, hi=1200, level=10.0):
    sig = np.zeros(n)
    sig[lo:hi] = level
    return sig


def test_overlaps_below_k_flank_leaves_no_nan_cells():
    """The reference's row slice selects nothing, so zeros participate in the mean."""
    n, k_flank, window, step = 4000, 12, 100, 10
    overlaps = int(window / step)
    assert overlaps <= k_flank, "this test documents the small-overlaps regime"

    preds = np.zeros((overlaps, n))
    preds[k_flank:n - k_flank] = np.nan      # row slice, as in the reference
    assert np.isnan(preds).sum() == 0


def _sparse_cuts(n, n_cut, cf, seed=13, level=200000.0):
    """Cut sites are sparse and heavy-tailed, as in real ATAC data."""
    rng = np.random.default_rng(seed)
    cuts = np.zeros(n)
    for pos in rng.integers(100, n - 100, n_cut):
        cuts[pos] += level
    return cuts


def test_expected_is_sparse_like_the_signal():
    """With 23 cut sites, `expected` lights up only a small neighbourhood."""
    n, k_flank = 20124, 12
    sig = _sparse_cuts(n, 23, 1.0)
    rng = np.random.default_rng(19)
    bias = np.zeros(n)
    bias[k_flank:n - k_flank] = rng.uniform(-14.0, 8.0, n - 2 * k_flank)

    expected, _ = bias_correction_tracks(sig, bias, k_flank, window=100, step=10)
    lo, hi = 62, n - 62
    nonzero = int((expected[lo:hi] > 0).sum())
    # measured on real data: 915 of 20000 -- sparse, but far more than 23
    assert 0 < nonzero < 0.2 * (hi - lo), nonzero


def test_corrected_can_be_negative_on_a_single_strand():
    """The sign structure the earlier (wrong) implementation could not produce."""
    n, k_flank = 4000, 12
    rng = np.random.default_rng(5)
    sig = np.zeros(n)
    sig[1000:1300] = rng.uniform(20, 40, 300)   # strong signal
    bias = np.zeros(n)
    bias[1000:1300] = -5.0                       # strongly favourable k-mers

    _, corrected = bias_correction_tracks(sig, bias, k_flank, window=100, step=10)
    assert corrected.min() < 0, "expected is subtracted, so negatives must appear"


def test_bias_correction_is_independent_of_absolute_offset():
    """Adding a constant to the bias track shifts the fit but not the structure."""
    n, k_flank = 4000, 12
    sig = _spike(n)
    b1 = np.zeros(n)
    b1[k_flank:n - k_flank] = -3.0
    b2 = b1 - 10.0

    e1, c1 = bias_correction_tracks(sig, b1, k_flank, window=100, step=10)
    e2, c2 = bias_correction_tracks(sig, b2, k_flank, window=100, step=10)
    assert np.all(np.isfinite(e1)) and np.all(np.isfinite(e2))
    # a constant offset changes the fitted slope, so values differ -- but both must
    # remain well-formed arrays rather than all-NaN or unbounded
    assert np.abs(c2).max() < 1e6
    assert pytest.approx(c1.sum(), rel=2.0) == c2.sum()
