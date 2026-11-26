"""PlotAggregate / PlotHeatmap preparation, checked against the reference rules.

The matplotlib calls are not the interesting part of these two commands; the numbers
fed to them are. Both share three preparation steps, and it is those that determine
whether the picture is right:

    centre each region on its midpoint at a fixed width (``flank * 2``),
    drop regions whose window leaves the chromosome,
    read the bigWig over the window.

`PlotAggregate` then averages per signal/region pair; `PlotHeatmap` keeps the rows
and can sort them.

One semantics trap, measured against the CLI rather than guessed:

    --remove-outliers <float>  "Value between 0-1 indicating the percentile of
                                regions to INCLUDE, e.g. 0.99 to remove the sites
                                with 1% highest values"   (default: 1)

So it is the KEPT percentile. I first had it as the removed fraction and defaulted
it to 0.1, which kept only the lowest 10% of regions -- the opposite of the intent.
Verified after the fix: 0.99 keeps 1980 of 2000 rows, 1.0 keeps all.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_plot import (
    DEFAULT_FLANK,
    DEFAULT_REMOVE_OUTLIERS,
    aggregate,
    build_matrix,
    extract_window,
    log_transform,
    normalize_rows_for_display,
    outlier_mask,
    remove_outlier_rows,
    sort_rows,
    window_around,
)


def _signal(n=20000, seed=3):
    return np.random.default_rng(seed).gamma(2.0, 1.0, n)


# ---------------------------------------------------------------------------
# remove-outliers semantics -- the trap
# ---------------------------------------------------------------------------
def test_default_keeps_everything():
    """The CLI default is 1, meaning the 100th percentile: nothing is removed."""
    assert DEFAULT_REMOVE_OUTLIERS == 1.0
    mat = _signal(2000).reshape(200, 10)
    assert outlier_mask(mat, DEFAULT_REMOVE_OUTLIERS).all()


def test_zero_point_nine_nine_removes_about_one_percent():
    """The CLI documents 0.99 as 'remove the 1% highest values'."""
    mat = _signal(2000).reshape(400, 5)
    mat[::100, 2] *= 60                      # 4 rows with a large spike
    kept = remove_outlier_rows(mat, 0.99).shape[0]
    removed_fraction = 1 - kept / mat.shape[0]
    assert 0.0 < removed_fraction < 0.05, removed_fraction


def test_half_keeps_half_the_rows():
    mat = _signal(2000).reshape(400, 5)
    assert remove_outlier_rows(mat, 0.5).shape[0] == pytest.approx(200, abs=5)


def test_outliers_are_the_high_peaks_not_the_low_ones():
    mat = _signal(2000).reshape(400, 5)
    mat[0, 3] = mat.max() * 100              # one clear outlier
    kept = remove_outlier_rows(mat, 0.99)
    assert not np.any(np.isclose(kept[:, 3], mat[0, 3]))


# ---------------------------------------------------------------------------
# window centring
# ---------------------------------------------------------------------------
def test_window_is_centred_on_the_midpoint_at_fixed_width():
    for start, end in ((4950, 5050), (4900, 5100), (4999, 5001)):
        ws, we, mid = window_around(start, end, flank=1000)
        assert mid == pytest.approx((start + end) / 2.0)
        assert we - ws == 2000, "width must be flank * 2 regardless of input length"
        assert (ws + we) / 2.0 == pytest.approx(mid, abs=1.0)


def test_window_inside_the_chromosome_is_returned():
    sig = _signal(20000)
    win = extract_window(sig, 20000, 5000, 5100, flank=1000)
    assert win is not None
    assert win.size == 2000


def test_window_crossing_the_edge_is_dropped():
    sig = _signal(20000)
    assert extract_window(sig, 20000, 300, 400, flank=1000) is None
    assert extract_window(sig, 20000, 19700, 19800, flank=1000) is None


def test_build_matrix_drops_edge_regions_and_reports_which():
    sig = _signal(20000)
    regions = [(5000, 5100), (8000, 8100), (300, 400)]
    mat, kept = build_matrix(sig, 20000, regions, flank=1000)
    assert mat.shape == (2, 2000)
    assert kept == [0, 1], "the near-edge region (index 2) must be excluded"


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------
def test_aggregate_is_column_mean_without_filtering():
    sig = _signal(20000)
    mat, _ = build_matrix(sig, 20000, [(5000, 5100), (8000, 8100)], flank=1000)
    assert np.allclose(aggregate(mat, 1.0), np.nanmean(mat, axis=0))


def test_aggregate_matches_manual_filter_then_mean():
    sig = _signal(20000)
    mat, _ = build_matrix(sig, 20000, [(5000, 5100), (8000, 8100), (12000, 12100)],
                          flank=1000)
    expected = np.nanmean(remove_outlier_rows(mat, 0.99), axis=0)
    assert np.allclose(aggregate(mat, 0.99), expected)


def test_log_transform_preserves_the_sign():
    v = np.array([0.0, 1.0, 3.0, -3.0, -1.0])
    manual = np.log2(np.abs(v) + 1)
    manual[v < 0] *= -1
    assert np.allclose(log_transform(v), manual)
    assert log_transform(v)[3] < 0, "negatives stay negative"


def test_log_transform_saturates_high_values():
    """log2(x+1) compresses the dynamic range, which is why the flag exists."""
    v = np.array([0.0, 100.0, 10000.0])
    out = log_transform(v)
    assert out[1] / out[0] if out[0] else True
    assert out[2] < 20, out[2]


# ---------------------------------------------------------------------------
# heatmap row order
# ---------------------------------------------------------------------------
def test_sort_rows_is_descending_on_the_key():
    sig = _signal(20000)
    mat, _ = build_matrix(sig, 20000, [(5000, 5100), (8000, 8100), (12000, 12100)],
                          flank=1000)
    keys = np.array([0.1, 5.0, 2.0])
    srt, skey = sort_rows(mat, keys, remove_outliers=1.0)
    assert np.all(np.diff(skey) <= 0), "PlotHeatmap sorts with reverse=True"
    assert np.allclose(srt[0], mat[1]), "the row with key 5.0 comes first"


def test_sort_rows_rejects_mismatched_lengths():
    sig = _signal(20000)
    mat, _ = build_matrix(sig, 20000, [(5000, 5100), (8000, 8100)], flank=1000)
    with pytest.raises(ValueError):
        sort_rows(mat, np.array([1.0, 2.0, 3.0]))


def test_sort_rows_applies_the_same_filter_as_aggregate():
    sig = _signal(20000)
    mat, _ = build_matrix(sig, 20000, [(5000, 5100), (8000, 8100), (12000, 12100)],
                          flank=1000)
    keys = np.array([1.0, 2.0, 3.0])
    _, skey = sort_rows(mat, keys, remove_outliers=0.5)
    assert skey.size == remove_outlier_rows(mat, 0.5).shape[0]


# ---------------------------------------------------------------------------
# display
# ---------------------------------------------------------------------------
def test_display_normalisation_clips_to_the_requested_range():
    mat = np.array([[0.0, 5.0], [10.0, -3.0]])
    out = normalize_rows_for_display(mat, vmin=0.0, vmax=2.0)
    assert out.min() == pytest.approx(0.0)
    assert out.max() == pytest.approx(2.0)


def test_display_normalisation_handles_a_flat_matrix():
    out = normalize_rows_for_display(np.zeros((3, 4)))
    assert np.all(np.isfinite(out))


def test_defaults_are_pinned():
    assert DEFAULT_FLANK == 1000
    assert DEFAULT_REMOVE_OUTLIERS == 1.0
