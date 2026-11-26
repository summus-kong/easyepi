"""BINDetect's statistics layer, verified against the reference computations.

Two things BINDetect does beyond normalising scores: decide which TFBS are bound,
and test whether binding changed between conditions. Both were transcribed from
`bindetect.py` 470-558 and `bindetect_functions.py` 500-548.

Measured against the reference formulas on shared fixtures:

    background log2fc fit   mean/std identical (max |diff| = 0)
    effect size             identical including the round to 5
    p-value                 identical (both seed np.random with n_obs)
    BH FDR                  identical to statsmodels' fdr_bh

The p-value match is not luck: the reference seeds `np.random.seed(n_obs)` before
drawing, so reproducing the seed reproduces the draws.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.stats

from easyepi.specific.atac.tobias_bindetect_stats import (
    BACKGROUND_PERCENTILE,
    GMM_COMPONENTS,
    GMM_RANDOM_STATE,
    N_SAMPLES,
    PSEUDO_FRACTION,
    background_log2fc_params,
    benjamini_hochberg,
    differential_pvalue,
    effect_size,
    estimate_background_distribution,
    estimate_threshold,
    mean_per_peak,
)


def _two_conditions(n=5000, seed=7):
    rng = np.random.default_rng(seed)
    s1 = np.concatenate([np.zeros(2000), np.abs(rng.normal(6.0, 2.0, n - 2000))])
    s2 = np.concatenate([np.zeros(2000), np.abs(rng.normal(9.0, 3.0, n - 2000))])
    return s1, s2


# ---------------------------------------------------------------------------
# Agreement with the reference computation
# ---------------------------------------------------------------------------
def test_background_log2fc_params_match_the_reference():
    s1, s2 = _two_conditions()
    pseudo = 0.31
    mine = background_log2fc_params(s1, s2, pseudo)

    # Reference: trim to the 1st-99th percentile, then norm.fit.
    included = np.logical_or(s1 > 0, s2 > 0)
    log2fcs = np.log2(np.true_divide(s1[included] + pseudo, s2[included] + pseudo))
    lower, upper = np.percentile(log2fcs, [1, 99])
    ref = scipy.stats.norm.fit(log2fcs[(log2fcs >= lower) & (log2fcs <= upper)])

    assert mine[0] == pytest.approx(ref[0], abs=1e-12)
    assert mine[1] == pytest.approx(ref[1], abs=1e-12)


def test_effect_size_matches_the_reference_formula():
    observed = np.array([0.5, -0.2, 1.1, 0.05, -0.4, 0.8])
    bg = (0.0, 1.0)
    obs_mean, obs_std = scipy.stats.norm.fit(observed)
    expected = np.round((obs_mean - bg[0]) / np.mean([obs_std, bg[1]]), 5)
    assert effect_size(observed, bg) == expected


def test_pvalue_matches_the_reference_sampling():
    """Same seed, same draws, same t-test -- so the p-value must coincide exactly."""
    observed = np.array([0.5, -0.2, 1.1, 0.05, -0.4, 0.8])
    bg = (0.0, 1.0)

    mine_p, change = differential_pvalue(observed, bg)

    np.random.seed(len(observed))
    changes = []
    for _ in range(N_SAMPLES):
        sample = scipy.stats.norm.rvs(*bg, size=len(observed))
        changes.append((np.mean(sample) - bg[0]) / np.mean([np.std(sample), bg[1]]))
    ref_p = scipy.stats.ttest_1samp(changes, change)[1]

    assert mine_p == pytest.approx(ref_p, abs=1e-15)


def test_benjamini_hochberg_matches_statsmodels():
    sm = pytest.importorskip("statsmodels.stats.multitest")
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216])
    assert np.allclose(benjamini_hochberg(p), sm.multipletests(p, method="fdr_bh")[1])


def test_mean_per_peak_matches_groupby_semantics():
    ids = ["p1", "p1", "p2", "p3", "p3", "p3"]
    vals = [1.0, 3.0, 2.0, 0.5, 1.5, 1.0]
    got_ids, got_vals = mean_per_peak(ids, vals)
    assert list(got_ids) == ["p1", "p2", "p3"]
    assert np.allclose(got_vals, [2.0, 2.0, 1.0])


# ---------------------------------------------------------------------------
# Documented behaviour
# ---------------------------------------------------------------------------
def test_reference_constants_are_pinned():
    assert BACKGROUND_PERCENTILE == 99
    assert GMM_COMPONENTS == 2
    assert GMM_RANDOM_STATE == 1
    assert PSEUDO_FRACTION == 0.5
    assert N_SAMPLES == 100


def test_zeros_are_excluded_before_fitting_the_background():
    scores = np.concatenate([np.zeros(5000), np.abs(np.random.default_rng(1)
                                                    .normal(5.0, 1.0, 500))])
    mode, log_params, pseudo = estimate_background_distribution(scores)
    assert mode > 0
    assert pseudo == pytest.approx(mode * 0.5)


def test_background_distribution_rejects_all_zero_input():
    with pytest.raises(ValueError):
        estimate_background_distribution(np.zeros(100))


def test_threshold_is_above_the_mode():
    rng = np.random.default_rng(3)
    scores = np.abs(rng.normal(5.0, 1.5, 3000))
    threshold, pseudo, mode = estimate_threshold(scores, bound_pvalue=1e-3)
    assert threshold > mode, "the bound threshold sits in the right tail"
    assert pseudo == pytest.approx(mode / 2.0)


def test_stricter_bound_pvalue_raises_the_threshold():
    rng = np.random.default_rng(5)
    scores = np.abs(rng.normal(5.0, 1.5, 3000))
    loose, _, _ = estimate_threshold(scores, bound_pvalue=1e-2)
    strict, _, _ = estimate_threshold(scores, bound_pvalue=1e-6)
    assert strict > loose


def test_identical_conditions_give_pvalue_one():
    """When both conditions have the same mean there is nothing to compare."""
    observed = np.zeros(200)
    p, change = differential_pvalue(observed, (0.0, 1.0))
    assert change == 0.0
    assert p == 1.0


def test_effect_size_sign_follows_the_shift():
    bg = (0.0, 1.0)
    up = np.concatenate([np.full(100, 2.0), np.linspace(1.8, 2.2, 100)])
    down = np.concatenate([np.full(100, -2.0), np.linspace(-2.2, -1.8, 100)])
    assert effect_size(up, bg) > 0
    assert effect_size(down, bg) < 0


def test_benjamini_hochberg_handles_empty_and_preserves_order():
    assert benjamini_hochberg(np.array([])).size == 0
    p = np.array([0.5, 0.001, 0.2])
    adj = benjamini_hochberg(p)
    assert np.all(adj >= p - 1e-12)
    assert adj[1] < adj[0], "the smallest p-value keeps the smallest adjusted value"
