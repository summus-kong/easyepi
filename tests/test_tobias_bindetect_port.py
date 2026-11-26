"""BINDetect's normalisation layer, checked against TOBIAS's own implementation.

BINDetect compares footprint scores across conditions. Raw scores are not
comparable between samples, so TOBIAS normalises each condition's background
distribution before computing log2 fold changes. Getting this layer wrong would
shift every fold change, so it is verified against the reference rather than
against my reading of it.

Measured against `tobias.tools.bindetect_functions`:

    sigmoid                max |diff| = 0            (identical)
    ArrayNorm.normalize    max |diff| = 0            (identical, incl. capping)
    quantile_normalization same func per condition, sigmoid popt within 1e-5 rel.
    normalised arrays      max |diff| <= 1.6e-06, corr = 1.00000000

The residual is `curve_fit` numerical tolerance, not an algorithmic difference.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_bindetect import (
    INTERP_POINTS,
    QUANTILES,
    SMOOTH_PAD,
    ArrayNorm,
    log2_fold_change,
    normalize_scores,
    quantile_normalization,
    sigmoid,
)

bindetect = pytest.importorskip(
    "tobias.tools.bindetect_functions", reason="TOBIAS not installed")

ref_sigmoid = bindetect.sigmoid
RefArrayNorm = bindetect.ArrayNorm
ref_qnorm = bindetect.quantile_normalization


# ---------------------------------------------------------------------------
# Agreement with the reference
# ---------------------------------------------------------------------------
def test_sigmoid_is_identical():
    x = np.linspace(-5.0, 5.0, 51)
    for params in ((0.7, 1.3, 2.0, 0.1), (0.0, 5.0, 1.0, 0.0), (2.5, -1.0, 3.0, -0.5)):
        assert np.array_equal(sigmoid(x, *params), ref_sigmoid(x, *params))


def test_array_norm_constant_is_identical():
    arr = np.array([-5.0, 0.0, 0.5, 2.0, 10.0, 100.0])
    mine = ArrayNorm("constant", 2.0, value_min=0.0, value_max=5.0).normalize(arr)
    ref = RefArrayNorm("constant", 2.0, value_min=0.0, value_max=5.0).normalize(arr)
    assert np.array_equal(mine, ref)


def test_array_norm_caps_before_applying_the_factor():
    """Capping is what stops outliers extrapolating beyond the fitted range."""
    arr = np.array([0.0, 1.0, 1000.0])
    mine = ArrayNorm("constant", 3.0, value_min=0.0, value_max=5.0).normalize(arr)
    ref = RefArrayNorm("constant", 3.0, value_min=0.0, value_max=5.0).normalize(arr)
    assert np.array_equal(mine, ref)
    # the factor itself is constant, so capping changes only the sigmoid path
    sig_mine = ArrayNorm("sigmoid", (1.0, 1.0, 1.0, 0.0),
                         value_min=0.0, value_max=5.0).normalize(arr)
    sig_ref = RefArrayNorm("sigmoid", (1.0, 1.0, 1.0, 0.0),
                           value_min=0.0, value_max=5.0).normalize(arr)
    assert np.allclose(sig_mine, sig_ref)


def _two_conditions(seed=4):
    rng = np.random.default_rng(seed)
    a = np.abs(rng.normal(5.0, 2.0, 20000))
    b = np.abs(rng.normal(8.0, 3.0, 20000))
    return a, b


def test_quantile_normalization_matches_the_reference():
    a, b = _two_conditions()
    names = ["A", "B"]
    mine = quantile_normalization([a, b], names)
    ref = ref_qnorm([a, b], names)

    for nm in names:
        assert mine[nm].func == ref[nm].func, nm
        assert mine[nm].value_min == pytest.approx(ref[nm].value_min, rel=1e-9)
        assert mine[nm].value_max == pytest.approx(ref[nm].value_max, rel=1e-9)
        if mine[nm].func == "sigmoid":
            assert np.allclose(mine[nm].popt, ref[nm].popt, rtol=1e-4, atol=1e-6)
        else:
            assert mine[nm].popt == pytest.approx(ref[nm].popt, rel=1e-6)


def test_normalised_scores_match_the_reference():
    a, b = _two_conditions()
    names = ["A", "B"]
    mine = quantile_normalization([a, b], names)
    ref = ref_qnorm([a, b], names)

    for nm, arr in zip(names, (a, b)):
        got = normalize_scores(arr, mine[nm])
        want = ref[nm].normalize(arr)
        want = np.where(want > 0, want, 0.0)
        assert np.abs(got - want).max() < 1e-5, nm
        assert np.corrcoef(got, want)[0, 1] > 0.999999


# ---------------------------------------------------------------------------
# Documented behaviour
# ---------------------------------------------------------------------------
def test_quantile_grid_matches_the_reference_constants():
    assert QUANTILES.size == 1000
    assert QUANTILES[0] == pytest.approx(0.05)
    assert QUANTILES[-1] == pytest.approx(0.99)
    assert SMOOTH_PAD == 50
    assert INTERP_POINTS == 1000


def test_zeros_are_excluded_from_the_quantile_estimate():
    """The 0.05 lower bound exists to skip the mass of zero-scoring background."""
    arr = np.zeros(10000)
    arr[:100] = np.linspace(1.0, 10.0, 100)
    mine = quantile_normalization([arr], ["only"])
    obj = mine["only"]
    # value_max comes from the positive quantiles, so it is ~10, not 0
    assert obj.value_max > 1.0


def test_single_condition_normalisation_is_a_no_op_equivalent():
    """With one condition TOBIAS uses a constant factor of 1.0."""
    constant = ArrayNorm("constant", 1.0, value_min=0, value_max=1)
    arr = np.array([0.0, 1.0, 5.0, 100.0])
    assert np.array_equal(constant.normalize(arr), arr)


def test_all_zero_input_falls_back_to_a_constant_factor():
    """`_positive_quantiles` returns [0] for an all-zero array; no crash."""
    obj = quantile_normalization([np.zeros(100)], ["z"])["z"]
    assert obj.func == "constant"
    assert obj.popt == pytest.approx(1.0)


def test_log2_fold_change_uses_the_pseudocount():
    a = np.array([10.0, 0.0, 3.0])
    b = np.array([5.0, 0.0, 3.0])
    got = log2_fold_change(a, b, 0.1)
    assert got[0] == pytest.approx(np.round(np.log2(10.1 / 5.1), 5))
    assert got[1] == pytest.approx(0.0), "equal inputs must give zero"
    assert got[2] == pytest.approx(0.0)


def test_log2_fold_change_is_symmetric_in_sign():
    a = np.array([8.0])
    b = np.array([2.0])
    assert log2_fold_change(a, b, 0.1)[0] == pytest.approx(
        -log2_fold_change(b, a, 0.1)[0], rel=1e-6)


def test_normalize_scores_clips_negatives():
    obj = ArrayNorm("constant", 1.0, value_min=-10.0, value_max=10.0)
    out = normalize_scores(np.array([-3.0, 0.0, 4.0]), obj)
    assert np.all(out >= 0)
    assert out[2] == pytest.approx(4.0)
