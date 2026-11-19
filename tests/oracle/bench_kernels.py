"""Measure the variance kernel at realistic scale, and find where the time goes.

The documented gap is: "on 300 sites that is milliseconds; genome-scale timing
has not been measured". A human WGBS panel has ~28M CpG sites, chr1 alone ~2.2M,
so a per-site Python loop there is the question that matters.

This times each stage separately, at increasing site counts, so a slow stage
cannot hide behind a fast one:

  smoothing (nitem_bin + window filter)   O(n) per chromosome
  variance kernel                         per-site window quadratic form
  dispersion shrinkage                    per-site bounded optimisation

The dispersion stage is expected to dominate, because each site runs
``scipy.optimize.minimize_scalar`` with tens of likelihood evaluations. If that
is what the numbers show, the honest conclusion is "site loop, not kernel".
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from easyepi.dmr_beta_binomial import (  # noqa: E402
    _chrom_groups, _variance_term, _smooth_by_chrom, beta_binomial_test,
)
from easyepi.dmr_dispersion import shrunken_dispersion  # noqa: E402
from easyepi.kernels import compute_var_smooth_rep, nitem_bin, window_filter  # noqa: E402

RNG = np.random.default_rng(0)


def fixture(n_sites: int, n_rep: int = 3, span_gap: int = 100_000):
    """One chromosome, irregular CpG spacing scaled so the site count fits."""
    pos = np.sort(RNG.choice(span_gap, size=n_sites, replace=False)).astype(np.int64)
    chrom = np.array(["chr1"] * n_sites, dtype=object)
    level = np.clip(RNG.normal(0.5, 0.12, (n_rep, n_sites)), 0.01, 0.99)
    total = RNG.integers(10, 40, (n_rep, n_sites)).astype(np.int64)
    meth = RNG.binomial(total, level).astype(np.int64)
    return meth, total, chrom, pos


def timeit(fn, *args, repeat: int = 1):
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        out = fn(*args)
        best = min(best, time.perf_counter() - t0)
    return best, out


def main() -> int:
    print(f"numpy {np.__version__}  python {sys.version.split()[0]}")
    print()
    print("=== stage timings, one chromosome, 3 replicates ===")
    header = f"{'sites':>8} {'nitem+filter':>13} {'var kernel':>12} {'dispersion':>12} {'total est':>11}"
    print(header)
    print("-" * len(header))

    rows = []
    for n in (1_000, 5_000, 20_000):
        meth, total, chrom, pos = fixture(n)
        x = meth.T.astype(float)
        nrep = total.T.astype(float)
        level = np.full(n, 0.5)
        phi = np.full(n, 0.05)
        groups = _chrom_groups(chrom, pos)

        t_smooth, _ = timeit(_smooth_by_chrom, nrep[:, 0].copy(), pos, groups, 500, "sum",
                             repeat=3)
        t_var, _ = timeit(compute_var_smooth_rep, level, nrep[:, 0].copy(), phi, pos, 500,
                          repeat=3)
        t_disp, _ = timeit(shrunken_dispersion, x.copy(), nrep.copy(), level,
                           (-3.0, 1.0), repeat=1)

        # The variance stage is called once per replicate, the others once.
        total_est = t_smooth + t_var * nrep.shape[1] + t_disp
        rows.append((n, t_smooth, t_var, t_disp, total_est))
        print(f"{n:>8} {t_smooth:>12.3f}s {t_var:>11.3f}s {t_disp:>11.2f}s {total_est:>10.2f}s")

    print()
    print("=== share of time at the largest size ===")
    n, t_smooth, t_var, t_disp, total_est = rows[-1]
    for label, t in (("dispersion (per-site optimiser)", t_disp),
                     ("variance kernel x3 reps", t_var * 3),
                     ("smoothing", t_smooth)):
        print(f"  {label:<32} {t:>7.2f}s  {100 * t / total_est:5.1f}%")

    print()
    print("=== extrapolation to chr1 (2.2M CpG sites) ===")
    per_site = total_est / n
    print(f"  measured: {per_site * 1e6:.1f} us/site at {n} sites")
    print(f"  chr1, one group: {per_site * 2.2e6 / 60:.1f} min")
    print(f"  chr1, two groups: {per_site * 2.2e6 * 2 / 60:.1f} min")
    print(f"  genome (28M), two groups: {per_site * 28e6 * 2 / 3600:.1f} h")
    print("  (linear extrapolation; the per-site stages are O(1) in n)")

    print()
    print("=== component check: is the kernel itself or the Python loop the cost? ===")
    meth, total, chrom, pos = fixture(20_000)
    x = meth.T.astype(float)
    nrep = total.T.astype(float)
    level = np.full(20_000, 0.5)
    phi = np.full(20_000, 0.05)
    left, right = nitem_bin(pos, 500)
    print(f"  sites per 500 bp window: mean {(left + right + 1).mean():.2f}, "
          f"max {(left + right + 1).max()}")
    t_nitem, _ = timeit(nitem_bin, pos, 500, repeat=3)
    t_filt, _ = timeit(window_filter, nrep[:, 0].copy(), pos, 500, "sum", repeat=3)
    t_var, _ = timeit(compute_var_smooth_rep, level, nrep[:, 0].copy(), phi, pos, 500,
                      repeat=3)
    print(f"  nitem_bin   {t_nitem * 1e3:8.1f} ms")
    print(f"  window_filter {t_filt * 1e3:6.1f} ms")
    print(f"  var kernel  {t_var * 1e3:8.1f} ms   -> "
          f"{t_var / 20_000 * 1e6:.1f} us/site")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
