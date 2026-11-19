"""Where exactly does the variance kernel spend its time, per site?

The stage timing says the kernel is 94.5% of the run. This breaks that down
further, because the three candidates have very different fixes:

  * building the lag matrix        -> O(m^2) indexing per site
  * rho ** lag                     -> transcendental, O(m^2) per site
  * u^T W u                        -> BLAS, O(m^2) per site
"""
from __future__ import annotations

import time

import numpy as np

BINSIZE = 50
RHO = 0.8


def bench(n_sites=20_000, n_span=100_000, ws=500, repeats=3):
    rng = np.random.default_rng(0)
    pos = np.sort(rng.choice(n_span, size=n_sites, replace=False)).astype(np.int64)
    u = rng.normal(0, 1, n_sites)

    # window half-widths, same rule as nitem_bin
    half = ws // 2
    left = np.searchsorted(pos, pos - half, side="left")
    right = np.searchsorted(pos, pos + half, side="right")
    m = right - left
    print(f"sites {n_sites}, window mean {m.mean():.1f}, max {m.max()}")

    def timed(fn, label):
        best = float("inf")
        for _ in range(repeats):
            t0 = time.perf_counter()
            fn()
            best = min(best, time.perf_counter() - t0)
        print(f"  {label:<38} {best * 1e3:8.1f} ms  {best / n_sites * 1e6:7.1f} us/site")
        return best

    def full():
        out = np.empty(n_sites)
        for i in range(n_sites):
            lo, hi = left[i], right[i]
            u_w = u[lo:hi]
            p_w = pos[lo:hi]
            lag = np.abs(p_w[:, None] - p_w[None, :]) // BINSIZE
            w = RHO ** lag.astype(float)
            out[i] = float(u_w @ w @ u_w) - float(np.dot(u_w, u_w))
        return out

    def lag_only():
        acc = 0.0
        for i in range(n_sites):
            lo, hi = left[i], right[i]
            p_w = pos[lo:hi]
            lag = np.abs(p_w[:, None] - p_w[None, :]) // BINSIZE
            acc += lag.sum()
        return acc

    lut = RHO ** np.arange(0, ws // BINSIZE + 3)

    def lag_plus_lut():
        acc = 0.0
        for i in range(n_sites):
            lo, hi = left[i], right[i]
            p_w = pos[lo:hi]
            lag = np.abs(p_w[:, None] - p_w[None, :]) // BINSIZE
            w = lut[lag]
            acc += w.sum()
        return acc

    def lag_plus_pow():
        acc = 0.0
        for i in range(n_sites):
            lo, hi = left[i], right[i]
            p_w = pos[lo:hi]
            lag = np.abs(p_w[:, None] - p_w[None, :]) // BINSIZE
            w = RHO ** lag.astype(float)
            acc += w.sum()
        return acc

    def matmul_only():
        acc = 0.0
        for i in range(n_sites):
            lo, hi = left[i], right[i]
            u_w = u[lo:hi]
            w = np.ones((hi - lo, hi - lo))
            acc += float(u_w @ w @ u_w)
        return acc

    print("  --- breakdown ---")
    t_full = timed(full, "whole per-site body")
    timed(lag_only, "  lag matrix only")
    timed(lag_plus_pow, "  lag + rho**lag")
    timed(lag_plus_lut, "  lag + lookup table (lut[lag])")
    timed(matmul_only, "  u^T W u with W = ones")

    # How many distinct lag values are there? A lookup table is only exact if
    # the exponents are integers, which they are by construction.
    print()
    print("  distinct lag values in a 500 bp window:",
          np.arange(0, ws // BINSIZE + 2))
    print(f"  lut covers 0..{lut.size - 1}; enough for span {ws} at BINSIZE {BINSIZE}")
    assert lut.size - 1 >= ws // BINSIZE + 1

    # Prove the table is bit-identical to the power call.
    lags = np.arange(0, ws // BINSIZE + 2, dtype=float)
    direct = RHO ** lags
    print("  lut == RHO**lag bitwise:", np.array_equal(lut[:lags.size], direct))

if __name__ == "__main__":
    bench()
