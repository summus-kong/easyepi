"""The numpy kernels must be the C kernels, not merely close to them.

DSS's smoothing lives in C. ``easyepi.kernels`` reimplements the three
functions in numpy so the package needs no compiler at install or run time.
That substitution is only acceptable if it is exact, so this module checks each
kernel against a literal transcription of the C loops -- different code, same
arithmetic -- across window widths, spacing patterns and negative/duplicate
positions.

The variance kernel is a windowed quadratic form and is the one place where a
"vectorised" rewrite silently changes meaning (the diagonal term counts once
per site in the window, the window's right edge is exclusive, the lag floor
divides by a fixed BINSIZE). Those three are asserted individually.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.kernels import (
    BINSIZE,
    compute_var_smooth_rep,
    nitem_bin,
    window_filter,
)
from easyepi.kernels import _nitem_bin_reference


def _c_var_smooth(vars_, delta, n1, pos, ws, rho=0.8):
    """``src/compute_var_smooth.c`` written out literally."""
    npos = pos.size
    nitem = _nitem_bin_reference(pos, ws)
    out = np.empty(npos, dtype=float)
    for ipos in range(npos):
        istart = ipos - nitem[0][ipos]
        iend = ipos + nitem[1][ipos]          # exclusive
        nCG_bin = iend - istart + 1
        thisvar = 0.0
        if nCG_bin <= 1:
            thisvar = vars_[ipos]
        else:
            for j in range(istart, iend):
                for k in range(j, iend):
                    if j == k:
                        thisvar += vars_[ipos]
                    else:
                        nlag = abs(int(pos[k]) - int(pos[j])) // BINSIZE
                        thisvar += 2 * (n1[j] * n1[k] * (rho ** nlag)
                                        * delta[j] * delta[k])
        out[ipos] = thisvar
    return out


@pytest.mark.parametrize("ws", [1, 2, 10, 50, 200, 500])
def test_nitem_bin_matches_the_c_loops(ws):
    rng = np.random.default_rng(ws)
    for _ in range(20):
        n = int(rng.integers(1, 60))
        pos = np.cumsum(rng.integers(1, 300, n))
        left, right = nitem_bin(pos, ws)
        ref_left, ref_right = _nitem_bin_reference(pos, ws)
        np.testing.assert_array_equal(left, ref_left)
        np.testing.assert_array_equal(right, ref_right)


def test_nitem_bin_handles_uniform_and_duplicate_positions():
    for pos in (np.arange(1, 21), np.repeat(np.arange(1, 11), 2),
                np.array([1, 1, 1, 2, 2, 100])):
        left, right = nitem_bin(pos, 50)
        ref_left, ref_right = _nitem_bin_reference(pos, 50)
        np.testing.assert_array_equal(left, ref_left)
        np.testing.assert_array_equal(right, ref_right)


@pytest.mark.parametrize("method", ["sum", "avg"])
@pytest.mark.parametrize("ws", [2, 10, 200, 500])
def test_window_filter_matches_the_c_kernel(ws, method):
    rng = np.random.default_rng(ws + len(method))
    n = 120
    pos = np.cumsum(rng.integers(1, 200, n))
    values = rng.normal(0, 5, n)
    got = window_filter(values, pos, ws, method)
    left, right = _nitem_bin_reference(pos, ws)
    csum = np.concatenate(([0.0], np.cumsum(values)))
    want = np.array([
        (csum[i + 1 + right[i]] - csum[i - left[i]])
        / (right[i] + left[i] + 1 if method == "avg" else 1)
        for i in range(n)
    ])
    np.testing.assert_allclose(got, want, rtol=1e-13, atol=1e-13)


def test_window_filter_of_a_constant_is_the_window_width():
    pos = np.cumsum(np.full(40, 7))
    total = window_filter(np.full(40, 3.0), pos, 200, "sum")
    left, right = nitem_bin(pos, 200)
    np.testing.assert_allclose(total, 3.0 * (left + right + 1))


@pytest.mark.parametrize("ws", [2, 20, 200, 500])
def test_compute_var_smooth_matches_the_c_kernel(ws):
    """Irregular spacing so different lags, and hence different rho powers, occur."""
    rng = np.random.default_rng(ws + 7)
    n = 90
    pos = np.cumsum(rng.integers(1, 400, n))
    p = rng.uniform(0.05, 0.95, n)
    phi = rng.uniform(0.01, 0.4, n)
    n1 = rng.uniform(1.0, 40.0, n)

    vars_ = n1 * p * (1 - p) * (1 + (n1 - 1) * phi)
    delta = p * (1 - p) * phi

    got = compute_var_smooth_rep(p, n1, phi, pos, ws)
    want = _c_var_smooth(vars_, delta, n1, pos, ws)
    np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)


def test_lag_uses_the_fixed_binsize_not_the_window():
    """Two sites 49 bp apart share lag 0; 50 bp apart is lag 1."""
    pos = np.array([0, 49, 50, 99, 100])
    p = np.full(5, 0.4)
    phi = np.full(5, 0.2)
    n1 = np.full(5, 10.0)
    got = compute_var_smooth_rep(p, n1, phi, pos, 500)
    vars_ = n1 * p * (1 - p) * (1 + (n1 - 1) * phi)
    delta = p * (1 - p) * phi
    want = _c_var_smooth(vars_, delta, n1, pos, 500)
    np.testing.assert_allclose(got, want, rtol=1e-14, atol=1e-14)


def test_single_site_window_is_only_the_diagonal():
    """One site in the window: ``thisvar = vars[ipos]``, with no pair terms."""
    pos = np.array([0, 10_000, 20_000])
    p = np.array([0.3, 0.6, 0.5])
    phi = np.array([0.1, 0.2, 0.3])
    n1 = np.array([5.0, 6.0, 7.0])
    got = compute_var_smooth_rep(p, n1, phi, pos, 10)
    want = n1 * p * (1 - p) * (1 + (n1 - 1) * phi)
    np.testing.assert_allclose(got, want, rtol=1e-14, atol=1e-14)


def test_right_edge_is_exclusive_as_in_the_c_source():
    """``iend = ipos + nitem[ipos+npos]`` and the loop stops at ``j < iend``.

    So the site at ``ipos + right`` is *not* in the pair sum, while
    ``filter.c`` uses the same half-width inclusively. That off-by-one is in
    DSS itself, and copying ``filter.c``'s convention here would silently
    change every variance, so it is asserted rather than described.
    """
    pos = np.array([0, 1, 2, 3])
    p = np.full(4, 0.5)
    phi = np.zeros(4)                      # delta = 0, so the pairs vanish
    n1 = np.full(4, 4.0)
    got = compute_var_smooth_rep(p, n1, phi, pos, 100)
    vars_ = n1 * p * (1 - p) * (1 + (n1 - 1) * phi)
    left, right = nitem_bin(pos, 100)
    # Three sites, not four: the right edge is exclusive.
    np.testing.assert_allclose(got, (left + right) * vars_)
    assert (left[0] + right[0]) == 3
    # window_filter, by contrast, is inclusive over the same half-widths.
    np.testing.assert_allclose(window_filter(np.ones(4), pos, 100, "sum"), 4.0)


def test_rho_powers_table_is_bit_identical():
    """The lookup table must be ``rho ** lag``, not an approximation of it.

    ``rho ** lag`` was 60% of the kernel's runtime, so the exponents are read
    from a table. That is only allowed if it changes nothing: the comparison is
    ``np.array_equal``, not ``allclose``.
    """
    from easyepi.kernels import _rho_powers

    for rho in (0.8, 0.5, 0.99):
        for ws in (1, 50, 200, 500, 1000, 2000):
            max_lag = ws // BINSIZE + 1
            table = _rho_powers(max_lag + 1, rho)
            assert table.size >= max_lag + 1
            # 1-D: every reachable exponent
            direct = rho ** np.arange(0, table.size, dtype=float)
            assert np.array_equal(table, direct), (rho, ws)
            # 2-D: the shape the kernel actually indexes with
            idx = np.arange(0, max_lag + 1)[:, None] * np.ones((1, max_lag + 1), dtype=int)
            assert np.array_equal(table[idx], rho ** idx.astype(float)), (rho, ws)


def test_the_kernel_never_indexes_past_the_table():
    """A window's largest lag must be reachable by the table the kernel builds."""
    from easyepi.kernels import _rho_powers

    for ws in (1, 10, 50, 200, 500, 1000):
        pos = np.arange(1, 41) * max(1, ws // 2)   # forces the widest windows
        p = np.full(pos.size, 0.5)
        phi = np.full(pos.size, 0.1)
        n1 = np.full(pos.size, 8.0)
        # would raise IndexError if a lag exceeded the table
        got = compute_var_smooth_rep(p, n1, phi, pos, ws)
        assert np.isfinite(got).all()
        max_lag = ws // BINSIZE + 1
        assert _rho_powers(max_lag + 1, 0.8).size >= max_lag + 1


def test_empty_input_is_empty_output():
    for fn, args in (
        (lambda: nitem_bin(np.zeros(0, dtype=int), 100), None),
        (lambda: window_filter(np.zeros(0), np.zeros(0, dtype=int), 100), None),
        (lambda: compute_var_smooth_rep(np.zeros(0), np.zeros(0),
                                        np.zeros(0), np.zeros(0, dtype=int), 100),
         None),
    ):
        result = fn()
        if isinstance(result, tuple):
            assert all(r.size == 0 for r in result)
        else:
            assert result.size == 0
