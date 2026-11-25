"""Time-course specificity: the arithmetic, and the cases it must refuse.

TAU is a descriptive statistic, so the tests pin the numbers rather than a
p-value. Three behaviours matter beyond the arithmetic:

  * a negative value makes TAU undefined -- taking absolute values would
    manufacture specificity that the data does not contain;
  * fewer than three time points cannot order anything, so the axis is refused;
  * an all-zero feature has no maximum, so TAU is reported as undefined rather
    than as zero, which would read as "measured, no specificity".
"""

from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from easyepi import timecourse as tc


def _row(name, values):
    return np.array(values, dtype=float)


def test_tau_is_one_for_a_single_spike():
    profile = _row("spike", [10.0, 0.0, 0.0, 0.0])
    assert tc.tau_specificity(profile) == pytest.approx(1.0)


def test_tau_is_zero_for_a_flat_profile():
    profile = _row("flat", [5.0, 5.0, 5.0, 5.0])
    assert tc.tau_specificity(profile) == pytest.approx(0.0)


def test_tau_matches_the_definition():
    # x/max = [1, 0.5, 0.25, 0]; sum(1 - x_hat) = 0 + 0.5 + 0.75 + 1 = 2.25
    # over n-1 = 3 -> 0.75
    profile = _row("graded", [8.0, 4.0, 2.0, 0.0])
    assert tc.tau_specificity(profile) == pytest.approx(0.75)


def test_tau_is_scale_invariant():
    """Multiplying a profile by a constant must not change specificity."""
    base = [3.0, 1.0, 0.5, 0.0]
    assert tc.tau_specificity(_row("a", base)) == pytest.approx(
        tc.tau_specificity(_row("a", [v * 17.0 for v in base])))


def test_negative_values_are_refused():
    """TAU is defined on non-negative expression only."""
    with pytest.raises(ValueError, match="negative"):
        tc.tau_specificity(_row("bad", [4.0, -1.0, 2.0, 0.0]))


def test_all_zero_is_undefined_not_zero():
    assert tc.tau_specificity(_row("z", [0.0, 0.0, 0.0, 0.0])) is None


def test_two_time_points_are_refused():
    """With two points every profile is either flat or a spike; tau says little."""
    with pytest.raises(ValueError, match="at least 3"):
        tc.tau_specificity(_row("two", [1.0, 0.0]))


def test_spearman_detects_a_monotone_decline():
    axis = np.array([0.0, 2.0, 6.0, 24.0])
    profile = np.array([9.0, 5.0, 2.0, 0.0])
    rho, p = tc._spearman(axis, profile)
    assert rho == pytest.approx(-1.0)
    assert p < 0.05


def test_spearman_needs_three_points():
    with pytest.raises(ValueError, match="at least 3"):
        tc._spearman(np.array([0.0, 1.0]), np.array([1.0, 2.0]))
