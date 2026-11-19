"""DSS's three smoothing kernels, in numpy.

DSS implements these in C (``src/nitem.c``, ``src/filter.c``,
``src/compute_var_smooth.c``). They are reproduced here in numpy instead of
being compiled and loaded through ctypes, because the package is allowed to
depend on numpy/scipy but must not assume a C toolchain on the user's machine.

``nitem_bin`` and ``window_filter`` are O(n) and transcribe the C loops exactly.
``compute_var_smooth`` is a windowed quadratic form; it is computed with a
per-window matrix product, which is the only part of this module whose cost
grows with the square of the sites per window (DSS's C pays the same
algorithmic cost). ``BINSIZE``, the lag floor and the accumulation order are
kept as in the C source.
"""
from __future__ import annotations

import functools

import numpy as np

__all__ = [
    "BINSIZE",
    "compute_var_smooth_rep",
    "kernels_available",
    "nitem_bin",
    "window_filter",
]

#: ``#define BINSIZE 50`` in DSS/src/DSS.h: lag is floor(bp distance / BINSIZE).
BINSIZE = 50


@functools.lru_cache(maxsize=32)
def _rho_powers(max_lag: int, rho: float) -> np.ndarray:
    """``rho ** k`` for every reachable integer lag, computed once.

    ``rho ** lag`` on a float array was 60% of the variance kernel's runtime
    (measured: 96.4 of 116.6 us/site at 101 sites per window). The exponents are
    integers by construction and never exceed ``span // BINSIZE + 1``, so a table
    replaces the transcendental call exactly: ``rho ** np.arange(n)`` and
    ``rho ** lag`` agree bitwise, which
    ``tests/test_dss_kernels.py::test_rho_powers_table_is_bit_identical`` pins.

    The table is not a tolerance-based approximation. DSS's own C does
    ``pow(rho, nlag)`` per term; IEEE ``pow`` on an integral ``double`` exponent
    is deterministic, and numpy uses the same libm call for both expressions.
    """
    return rho ** np.arange(0, int(max_lag))


def kernels_available() -> bool:
    """Kept for callers that used to probe the compiled library."""
    return True


def nitem_bin(pos: np.ndarray, ws: int) -> tuple:
    """``src/nitem.c``: number of sites within ``ws/2`` on each side.

    Returns ``(left, right)``; the C function writes ``[left | right]`` into one
    array of length ``2*n``. A window is the inclusive range
    ``[ipos - left[ipos], ipos + right[ipos]]``.
    """
    pos = np.ascontiguousarray(pos)
    n = pos.size
    if n == 0:
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty
    half = ws // 2
    # A Python two-pointer loop here was 20.6 ms per 20k sites per call, and the
    # kernel calls it once per replicate. Positions are sorted, so the window is
    # just two binary searches -- same answer, no interpreter loop.
    #
    #   left [i] = i - lower_bound(pos, pos[i] - half)
    #   right[i] = upper_bound(pos, pos[i] + half) - i - 1
    #
    # side="left" keeps sites exactly half away (the C test is `<= ws/2`), and
    # side="right" stops past them for the same reason.
    left = np.arange(n, dtype=np.int64) - np.searchsorted(pos, pos - half, side="left")
    right = np.searchsorted(pos, pos + half, side="right") - np.arange(n, dtype=np.int64) - 1
    return left, right


def _nitem_bin_reference(pos: np.ndarray, ws: int) -> tuple:
    """The C loops written out literally; used by the tests, not the hot path."""
    n = pos.size
    half = ws // 2
    left = np.empty(n, dtype=np.int64)
    right = np.empty(n, dtype=np.int64)
    for ipos in range(n):
        j = ipos - 1
        while j >= 0 and pos[ipos] - pos[j] <= half:
            j -= 1
        left[ipos] = ipos - j - 1
        j = ipos + 1
        while j < n and pos[j] - pos[ipos] <= half:
            j += 1
        right[ipos] = j - ipos - 1
    return left, right


def window_filter(values: np.ndarray, pos: np.ndarray, ws: int,
                  method: str = "sum") -> np.ndarray:
    """``src/filter.c``: moving sum (``"sum"``) or moving average (``"avg"``)."""
    values = np.asarray(values, dtype=float)
    pos = np.asarray(pos)
    if values.size == 0:
        return values.copy()
    left, right = nitem_bin(pos, int(ws))
    csum = np.concatenate(([0.0], np.cumsum(values)))
    idx = np.arange(values.size)
    total = csum[idx + 1 + right] - csum[idx - left]
    if method == "avg":
        return total / (left + right + 1)
    return total


def compute_var_smooth_rep(estprob: np.ndarray, n_rep: np.ndarray,
                           phi: np.ndarray, pos: np.ndarray, ws: int,
                           rho: float = 0.8) -> np.ndarray:
    """``src/compute_var_smooth.c`` for one replicate.

    The C loops are::

        istart = ipos - nitem[ipos];
        iend   = ipos + nitem[ipos + npos];        /* exclusive */
        nCG_bin = iend - istart + 1;
        if (nCG_bin <= 1) thisvar = vars[ipos];
        else for (j = istart; j < iend; j++)
                 for (k = j; k < iend; k++)
                     if (j == k) thisvar += vars[ipos];
                     else thisvar += 2 * n1[j]*n1[k]*rho^(|pos[k]-pos[j]|/BINSIZE)
                                            * delta[j]*delta[k];

    Two consequences that are easy to miss and are preserved here: the diagonal
    term contributes ``vars[ipos]`` once per ``j`` in the window (so ``m`` times,
    not once), and the site at ``iend`` is never part of the pair sum.

    The pair sum is ``u^T W u - sum(u^2)`` with ``u = n1*delta`` and
    ``W[j,k] = rho^(lag)``, which is the same number as the double loop and is
    computed here as one matrix product per site.
    """
    estprob = np.asarray(estprob, dtype=float)
    n_rep = np.asarray(n_rep, dtype=float)
    phi = np.asarray(phi, dtype=float)
    pos = np.ascontiguousarray(pos)
    n = pos.size
    if n == 0:
        return np.zeros(0, dtype=float)

    vars_ = n_rep * estprob * (1.0 - estprob) * (1.0 + (n_rep - 1.0) * phi)
    u = n_rep * estprob * (1.0 - estprob) * phi

    left, right = nitem_bin(pos, int(ws))
    # One table for the whole call: the largest lag any window can produce.
    max_lag = int(ws) // BINSIZE + 1
    table = _rho_powers(max_lag + 1, float(rho)) if rho != 0.0 else None
    out = np.empty(n, dtype=float)
    for ipos in range(n):
        istart = ipos - left[ipos]
        iend = ipos + right[ipos]
        m = iend - istart
        if m <= 1:
            out[ipos] = vars_[ipos]
            continue
        u_w = u[istart:iend]
        p_w = pos[istart:iend].astype(np.int64)
        lag = np.abs(p_w[:, None] - p_w[None, :]) // BINSIZE
        if lag.max() == 0:
            # A single lag bin: the weights all equal 1, so the quadratic form
            # collapses to (sum u)^2 - sum(u^2).
            quad = u_w.sum() ** 2 - float(np.dot(u_w, u_w))
        else:
            w = table[lag] if table is not None else rho ** lag.astype(float)
            quad = float(u_w @ w @ u_w) - float(np.dot(u_w, u_w))
        # `quad` already counts each unordered pair twice (j,k and k,j), which
        # is the factor of 2 the C loop puts in front of each off-diagonal term.
        out[ipos] = m * vars_[ipos] + quad
    return out
