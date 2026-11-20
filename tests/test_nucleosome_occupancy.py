"""Nucleosome occupancy: quantile NFR, honest repeat verdict, CLI wiring."""
import numpy as np
import pytest

from easyepi.nucleosome import (nucleosome_repeat, nucleosome_scores,
                                occupancy_profile)


def _lattice(n=20000, period_bins=18, width=1.2, noise=0.05, seed=0):
    """A phased occupancy track: one protection peak every `period_bins`."""
    x = np.arange(n)
    sig = np.zeros(n)
    for c in np.arange(50, n - 50, period_bins):
        sig += np.exp(-0.5 * ((x - c) / width) ** 2)
    sig += np.random.default_rng(seed).random(n) * noise
    return sig


# --------------------------------------------------------------------------
# Nucleosome-free fraction
# --------------------------------------------------------------------------
def test_nfr_is_not_zero_on_a_right_skewed_track():
    """The old `z < -1` rule reported 0% on every real MNase input.

    Real coverage is mostly zeros with a long tail, so z-scoring never reaches
    -1 (measured min -0.89 on GSE247785). A quantile definition cannot collapse
    to zero that way.
    """
    sig = np.zeros(10000)
    sig[:2000] = np.random.default_rng(0).random(2000)   # 80% empty
    scores = nucleosome_scores(sig)
    assert scores["nfr_score"] > 0.0
    # Ties count: the 20th percentile of a signal that is 80% zeros is 0, so
    # every zero bin is included and the fraction is 0.8, not 0.2. This is the
    # documented "depleted relative to the rest of the sample" semantics.
    assert scores["nfr_score"] == pytest.approx(0.8, abs=0.02)


def test_nfr_on_a_tie_free_signal_matches_the_requested_quantile():
    sig = np.random.default_rng(1).random(10000)   # continuous, no exact ties
    assert nucleosome_scores(sig, nfr_quantile=0.2)["nfr_score"] == pytest.approx(0.2, abs=0.01)


def test_nfr_quantile_is_honoured_and_reported():
    sig = np.random.default_rng(1).random(1000)
    a = nucleosome_scores(sig, nfr_quantile=0.1)
    b = nucleosome_scores(sig, nfr_quantile=0.3)
    assert a["nfr_score"] < b["nfr_score"]
    assert a["nfr_quantile"] == 0.1
    assert "nfr_cutoff" in a


def test_nfr_rejects_an_out_of_range_quantile():
    with pytest.raises(ValueError, match="nfr_quantile"):
        nucleosome_scores(np.ones(100), nfr_quantile=0.0)
    with pytest.raises(ValueError, match="nfr_quantile"):
        nucleosome_scores(np.ones(100), nfr_quantile=1.0)


def test_empty_signal_returns_nan_not_a_crash():
    s = nucleosome_scores(np.array([]))
    assert np.isnan(s["nfr_score"])
    assert np.isnan(s["mean_occupancy"])


def test_occupancy_profile_is_z_scored():
    sig = np.arange(100, dtype=float)
    z = occupancy_profile(sig)
    assert abs(z.mean()) < 1e-9
    assert abs(z.std() - 1.0) < 1e-6


# --------------------------------------------------------------------------
# Repeat verdict
# --------------------------------------------------------------------------
def test_repeat_is_reported_on_a_phased_lattice():
    sig = _lattice(period_bins=18)                 # 180 bp at bin_size 10
    r = nucleosome_repeat(sig, bin_size=10)
    assert r["reliable"] is True
    assert 150 <= r["dominant_period_bp"] <= 220
    assert "note" not in r


def test_repeat_is_flagged_when_no_peak_is_in_range():
    """Welch runs with nperseg=256, so a whole-track peak lands on low bins."""
    rng = np.random.default_rng(0)
    walk = np.cumsum(rng.standard_normal(50000))    # 1/f-like trend, no repeat
    r = nucleosome_repeat(walk, bin_size=10)
    assert r["reliable"] is False
    assert r["note"]


def test_repeat_reports_none_for_a_signal_too_short_to_analyse():
    r = nucleosome_repeat(np.ones(4), bin_size=10)
    assert r["dominant_period_bp"] is None
    assert r["reliable"] is False


def test_repeat_search_range_is_echoed():
    r = nucleosome_repeat(_lattice(), bin_size=10, lo_bp=140.0, hi_bp=260.0)
    assert r["search_range_bp"] == [140.0, 260.0]
