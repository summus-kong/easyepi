"""Beta-binomial test with between-replicate dispersion, matching DSS 1:1.

DSS is an R package. EasyEpi does not require R at runtime: this module
reimplements ``DSS::DMLtest`` (``smoothing = TRUE``, ``equal.disp = FALSE``) on
top of numpy and scipy, which the package already depends on. The
transcription is line by line from the installed DSS 2.58.0 package; every
intermediate it produces was compared against DSS's own, not only the final
p-value (see ``tests/test_beta_binomial_matches_dss.py``).

Per group, in this order -- the order matters, each step feeds the next
--------------------------------------------------------------------------

.. code-block:: r

    estprob <- compute.mean.Smooth(x, n, allchr, allpos, span)
    prior   <- est.prior.BSseq.logN(x, n)
    phi     <- dispersion.shrinkage.BSseq(x, n, prior, estprob, ncores)
    wt      <- 1/(1 + (n-1)*phi);  wt <- wt/mean(wt)
    estprob <- compute.mean.Smooth(x*wt, n*wt, allchr, allpos, span)
    var     <- compute.var.Smooth(estprob[,1], n*wt, phi, span, allchr, allpos)

then ``diff = estprob1 - estprob2``, ``vv = var1 + var2`` with DSS's two-sided
guard, ``se = sqrt(vv)``, ``pval = 2*pnorm(-|diff/se|)``.

Details that a plausible-looking reimplementation gets wrong
------------------------------------------------------------

1. ``compute.mean.Smooth`` is not a plain smoothed ratio. It adds a
   pseudocount *before* smoothing, and the pseudocount is
   ``const = mean(X/N)`` over every finite entry, not 0.5:
   ``N <- rowSums(N) + 1;  X <- rowSums(X) + const`` (``N + 0.4`` and
   ``X + 0.4*const`` when a group has a single replicate).
2. ``compute.var.Smooth`` receives ``n*wt``, the *weighted* depth, and adds
   ``const = 1/nreps`` to it. Feeding it the raw depth changes every variance.
3. Windows are genomic: ``nitem_bin`` compares positions against ``span/2``, so
   two CpG sites 5 bp apart and two 400 bp apart see different windows. An
   index-based window is wrong unless the sites are uniformly spaced.
4. ``est.prior.BSseq.logN`` uses ``IQR/1.39`` and has no minimum floor.
5. In ``compute.var.Smooth`` the diagonal contributes ``vars[ipos]`` once per
   site in the window, and the window is ``[ipos-left, ipos+right)`` -- the
   right edge is exclusive here but inclusive in ``filter.c``.
6. ``DMLtest`` drops sites with too little coverage (``hasCoverage``) *before*
   anything is computed, and that drop changes every later window. The filter
   is applied inside :func:`beta_binomial_test` for that reason, and the
   dropped sites are reported as ``NaN`` with ``tested == False``.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from scipy.stats import norm  # type: ignore

from .dmr_dispersion import (
    estimate_prior_log_dispersion,
    shrunken_dispersion,
)
from .kernels import compute_var_smooth_rep, nitem_bin, window_filter

__all__ = [
    "beta_binomial_test",
    "compute_mean_no_smooth",
    "compute_mean_smooth",
    "dss_has_coverage",
    "estimate_prior_log_dispersion",
    "shrunken_dispersion",
]

#: DSS's rho, hard-coded inside compute.var.Smooth.
_RHO = 0.8
#: DSS's window width for the coverage filter (``hasCoverage``: ``ws <- 200``).
_COVERAGE_SPAN = 200
#: DSS truncates the summed variance outside [1e-5, 1-1e-5] to [1e-4, 1-1e-4].
_VAR_LOW_TRIGGER = 1e-5
_VAR_LOW_SET = 1e-4
_VAR_HIGH_TRIGGER = 1 - 1e-5
_VAR_HIGH_SET = 1 - 1e-4
#: compute.waldStat.noSmooth floors to 1e-05 (not 1e-04) and has no upper cap.
_NO_SMOOTH_VAR_FLOOR = 1e-5


# --------------------------------------------------------------------------
# chromosome cohorts
# --------------------------------------------------------------------------

def _as_chrom(chrom, n_sites: int) -> np.ndarray:
    if chrom is None:
        return np.zeros(n_sites, dtype=np.int64)
    arr = np.asarray(chrom)
    if arr.shape != (n_sites,):
        raise ValueError(f"chrom must have shape ({n_sites},)")
    return arr


def _chrom_groups(chrom: np.ndarray, pos: np.ndarray) -> list:
    """Index arrays, one per chromosome, each in ascending position order.

    ``smooth.chr`` and the variance kernel are both applied per chromosome in
    DSS (``idx <- split(1:length(allchr), allchr)``), so a site list spanning
    several chromosomes must be split before either is called.
    """
    if pos.size == 0:
        return []
    labels = np.asarray([str(c) for c in chrom])
    order = np.lexsort((pos, labels))
    sorted_labels = labels[order]
    breaks = np.flatnonzero(sorted_labels[1:] != sorted_labels[:-1]) + 1
    return [g for g in np.split(order, breaks) if g.size]


def _smooth_by_chrom(values: np.ndarray, pos: np.ndarray, groups: list,
                     span: int, method: str) -> np.ndarray:
    """``DSS:::smooth.chr``: the window filter, per chromosome."""
    out = np.zeros(values.size, dtype=float)
    for idx in groups:
        out[idx] = window_filter(values[idx], pos[idx], int(span), method)
    return out


# --------------------------------------------------------------------------
# DSS's own DMLtest entry steps
# --------------------------------------------------------------------------

def dss_has_coverage(n_rep: np.ndarray, pos: np.ndarray,
                     chrom=None, thresh: int = 2) -> np.ndarray:
    """``DSS:::hasCoverage``: smoothed total depth above ``thresh * n_samples``.

    ``n_rep`` is ``(n_samples, n_sites)``. ``DMLtest`` requires a site to pass
    this in **both** groups; sites that fail are dropped before any statistic
    is computed.
    """
    n_rep = np.atleast_2d(np.asarray(n_rep, dtype=float))
    n_sites = n_rep.shape[1]
    chrom = _as_chrom(chrom, n_sites)
    pos = np.asarray(pos)
    groups = _chrom_groups(chrom, pos)
    smoothed = _smooth_by_chrom(n_rep.sum(axis=0), pos, groups,
                                _COVERAGE_SPAN, "sum")
    return smoothed > thresh * n_rep.shape[0]


def compute_mean_smooth(x: np.ndarray, n: np.ndarray, span: int,
                        chrom=None, pos=None) -> np.ndarray:
    """``DSS:::compute.mean.Smooth``: one smoothed level per sample.

    ``x`` and ``n`` are ``(n_sites, n_reps)``, the layout DSS passes. The
    pseudocount goes on the *summed* counts before smoothing, and the result is
    repeated down the columns (``matrix(rep(p, nreps), ncol = nreps)``),
    because the smoothed level does not distinguish replicates.
    """
    x = np.asarray(x, dtype=float)
    n = np.asarray(n, dtype=float)
    if x.ndim != 2:
        raise ValueError("x must be (n_sites, n_reps)")
    n_sites, nreps = x.shape
    pos = np.arange(n_sites) if pos is None else np.asarray(pos)
    chrom = _as_chrom(chrom, n_sites)
    groups = _chrom_groups(chrom, pos)

    with np.errstate(invalid="ignore", divide="ignore"):
        p1 = x / n
    finite = np.isfinite(p1)
    const = float(p1[finite].mean()) if finite.any() else 0.0

    if nreps > 1:
        n_sum = n.sum(axis=1) + 1.0
        x_sum = x.sum(axis=1) + const
    else:
        n_sum = n[:, 0] + 0.4
        x_sum = x[:, 0] + 0.4 * const

    x_sm = _smooth_by_chrom(np.asarray(x_sum, dtype=float), pos, groups, span, "sum")
    n_sm = _smooth_by_chrom(np.asarray(n_sum, dtype=float), pos, groups, span, "sum")
    with np.errstate(invalid="ignore", divide="ignore"):
        p = x_sm / n_sm
    # R would carry NA here; an empty window cannot occur after the coverage
    # filter, and 0 keeps the arithmetic finite if a caller skips it.
    p = np.where(np.isfinite(p), p, 0.0)
    return np.repeat(p[:, None], nreps, axis=1)


def compute_mean_no_smooth(x: np.ndarray, n: np.ndarray) -> np.ndarray:
    """``DSS:::compute.mean.noSmooth``: the per-site level, no window.

    .. code-block:: r

        p <- X/N
        const <- mean(p, na.rm = TRUE)
        p <- (rowSums(X) + const)/(rowSums(N) + 1)
        res <- matrix(rep(p, nreps), ncol = nreps)

    The pseudocount is the same ``mean(X/N)`` as the smoothed path, but there is
    no window smoothing: the level is the pooled ratio of the site itself.
    """
    x = np.asarray(x, dtype=float)
    n = np.asarray(n, dtype=float)
    if x.ndim != 2:
        raise ValueError("x must be (n_sites, n_reps)")
    with np.errstate(invalid="ignore", divide="ignore"):
        p = x / n
    finite = np.isfinite(p)
    const = float(p[finite].mean()) if finite.any() else 0.0
    with np.errstate(invalid="ignore", divide="ignore"):
        level = (x.sum(axis=1) + const) / (n.sum(axis=1) + 1.0)
    level = np.where(np.isfinite(level), level, 0.0)
    return np.repeat(level[:, None], x.shape[1], axis=1)


def _variance_term_no_smooth(estprob: np.ndarray, n_rep: np.ndarray,
                             phi: np.ndarray) -> np.ndarray:
    """``DSS:::compute.waldStat.noSmooth``'s variance, one group.

    .. code-block:: r

        n1m  <- rowSums(n1)
        var1 <- rowSums(n1 * estprob1 * (1-estprob1) * (1 + (n1-1)*phi1)) / n1m^2

    ``n1`` and ``phi1`` recycle down the columns, so each replicate keeps its own
    depth and the same per-site dispersion.
    """
    n_rep = np.atleast_2d(np.asarray(n_rep, dtype=float))
    n_sum = n_rep.sum(axis=0)
    per_rep = n_rep * estprob[None, :] * (1.0 - estprob[None, :]) \
        * (1.0 + (n_rep - 1.0) * phi[None, :])
    with np.errstate(invalid="ignore", divide="ignore"):
        out = per_rep.sum(axis=0) / (n_sum ** 2)
    return np.where(np.isfinite(out), out, 0.0)


def _variance_term(estprob: np.ndarray, n_rep: np.ndarray, phi: np.ndarray,
                   span: int, groups: list, pos: np.ndarray) -> np.ndarray:
    """``DSS:::compute.var.Smooth`` -- one group's variance of the level.

    ``n_rep`` is the **weighted** depth ``n*wt``, shape ``(n_reps, n_sites)``.
    DSS adds ``const = 1/nreps`` to it, sums each replicate's variance-
    covariance window sum, then divides by ``rowSums(n1.sm)^2`` where ``n1.sm``
    is the per-replicate smoothing of the same adjusted depth.
    """
    n_rep = np.atleast_2d(np.asarray(n_rep, dtype=float))
    nreps = n_rep.shape[0]
    n1 = n_rep + 1.0 / nreps

    acc = np.zeros(estprob.size, dtype=float)
    for irep in range(nreps):
        acc = acc + _var_one_rep(estprob, n1[irep], phi, span, groups, pos)

    n1_sm = np.zeros(estprob.size, dtype=float)
    for irep in range(nreps):
        n1_sm = n1_sm + _smooth_by_chrom(n1[irep], pos, groups, span, "sum")
    denom = n1_sm ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        out = acc / denom
    return np.where(np.isfinite(out), out, 0.0)


def _var_one_rep(estprob: np.ndarray, n_rep: np.ndarray, phi: np.ndarray,
                 span: int, groups: list, pos: np.ndarray) -> np.ndarray:
    """One replicate of ``compute.var.Smooth``, chromosome by chromosome."""
    out = np.zeros(estprob.size, dtype=float)
    for idx in groups:
        out[idx] = compute_var_smooth_rep(estprob[idx], n_rep[idx], phi[idx],
                                          pos[idx], int(span), _RHO)
    return out


def _guard_no_smooth(var: np.ndarray) -> np.ndarray:
    """``compute.waldStat.noSmooth``'s guard: one-sided, and to 1e-05.

    .. code-block:: r

        vv[vv < 1e-05] <- 1e-05

    That is not a typo on DSS's side, and it is not the rule the smoothed path
    uses -- that one floors to 1e-04 *and* caps the top. Reusing the smoothed
    guard here would change every no-smoothing p-value at low coverage.
    """
    var = np.array(var, dtype=float, copy=True)
    var[var < _NO_SMOOTH_VAR_FLOOR] = _NO_SMOOTH_VAR_FLOOR
    return var


def _guard(var: np.ndarray) -> np.ndarray:
    """``compute.waldStat.Smooth``'s two-sided guard on the *summed* variance."""
    var = np.array(var, dtype=float, copy=True)
    var[var < _VAR_LOW_TRIGGER] = _VAR_LOW_SET
    var[var > _VAR_HIGH_TRIGGER] = _VAR_HIGH_SET
    return var


# --------------------------------------------------------------------------
# the test
# --------------------------------------------------------------------------

def beta_binomial_test(
    meth: np.ndarray,
    total: np.ndarray,
    design: np.ndarray,
    smoothing_span: int = 500,
    chrom: Optional[Sequence[str]] = None,
    pos: Optional[Sequence[int]] = None,
    equal_disp: bool = False,
    smoothing: bool = True,
) -> dict:
    """``DSS::DMLtest`` on CpG sites, two groups, with replicates.

    ``meth``/``total`` are ``(n_samples, n_sites)``; ``design`` is 0 for group
    1 and 1 for group 2. ``chrom``/``pos`` are the real genomic coordinates --
    they decide the smoothing windows, so omitting them puts every site on a
    one-base grid and changes the result.

    ``smoothing=False`` reproduces ``DMLtest(smoothing = FALSE)``, i.e.
    ``compute.mean.noSmooth`` plus ``compute.waldStat.noSmooth``. That path uses
    a *different* variance floor (``1e-05``, one-sided) and no window at all, so
    it is a separate code path, not a span of zero. With ``smoothing=False``
    DSS also shares one dispersion whenever either group has a single
    replicate; this function then reports ``equal_disp`` as True.

    Returns arrays in the caller's site order: ``diff``, ``se``, ``stat``,
    ``pvalue``, ``phi1``/``phi2``, ``var1``/``var2``, the first-pass levels
    ``level1``/``level2``, the weighted levels ``level1_weighted``/
    ``level2_weighted``, both priors, and ``tested``, a boolean mask of the
    sites DSS keeps. Sites outside that mask are ``NaN``.
    """
    meth = np.asarray(meth, dtype=float)
    total = np.asarray(total, dtype=float)
    design = np.asarray(design).astype(int)
    if meth.ndim != 2 or total.shape != meth.shape:
        raise ValueError("meth and total must both be (n_samples, n_sites)")
    n_samples, n_sites = meth.shape
    if design.shape != (n_samples,):
        raise ValueError(f"design must have shape ({n_samples},)")

    g1 = design == 0
    g2 = design == 1
    if g1.sum() < 1 or g2.sum() < 1:
        raise ValueError("both groups need at least one sample")

    chrom_arr = _as_chrom(chrom, n_sites)
    pos_arr = np.arange(n_sites) if pos is None else np.asarray(pos)

    # DSS works in (chr, pos) order; the caller gets its own order back.
    order = np.lexsort((pos_arr, np.asarray([str(c) for c in chrom_arr])))

    # DMLtest drops sites that fail hasCoverage in either group *before*
    # computing anything, and that drop moves every later window, so the filter
    # has to happen here and not only in the caller.
    keep_sorted = (dss_has_coverage(total[np.ix_(g1, order)], pos_arr[order],
                                    chrom_arr[order])
                   & dss_has_coverage(total[np.ix_(g2, order)], pos_arr[order],
                                      chrom_arr[order]))
    kept = order[keep_sorted]
    kept_mask = np.zeros(n_sites, dtype=bool)
    kept_mask[kept] = True

    x1 = meth[np.ix_(g1, kept)].T          # (sites, reps)
    n1 = total[np.ix_(g1, kept)].T
    x2 = meth[np.ix_(g2, kept)].T
    n2 = total[np.ix_(g2, kept)].T
    chr_sorted = chrom_arr[kept]
    pos_sorted = pos_arr[kept]
    groups = _chrom_groups(chr_sorted, pos_sorted)
    nreps1, nreps2 = x1.shape[1], x2.shape[1]

    # DSS switches the whole first/recompute pair on `smoothing`:
    #   Smooth   -> compute.mean.Smooth   (window over allchr/allpos, pseudocount)
    #   noSmooth -> compute.mean.noSmooth (per-site pooled ratio, same pseudocount)
    mean_fn = compute_mean_smooth if smoothing else compute_mean_no_smooth
    mean_args = (smoothing_span, chr_sorted, pos_sorted) if smoothing else ()

    level1_raw = mean_fn(x1, n1, *mean_args)
    level2_raw = mean_fn(x2, n2, *mean_args)

    prior1 = estimate_prior_log_dispersion(x1, n1)
    prior2 = estimate_prior_log_dispersion(x2, n2)

    # DMLtest.noSmooth: `if (equal.disp | nreps1 == 1 | nreps2 == 1)` -- the
    # no-smoothing path forces a shared dispersion when a group has no replicate,
    # because a per-group estimate would not exist. DMLtest.Smooth does not.
    if not smoothing and (nreps1 == 1 or nreps2 == 1):
        equal_disp = True

    if equal_disp:
        # DSS: est.dispersion.BSseq(cbind(x1, x2), cbind(n1, n2),
        #                           cbind(estprob1, estprob2))
        xs = np.hstack([x1, x2])
        ns = np.hstack([n1, n2])
        lv = np.hstack([level1_raw, level2_raw])
        prior_eq = estimate_prior_log_dispersion(xs, ns)
        phi1 = phi2 = shrunken_dispersion(xs.T, ns.T, lv, prior_eq)
        prior1 = prior2 = prior_eq
    else:
        phi1 = shrunken_dispersion(x1.T, n1.T, level1_raw, prior1)
        phi2 = shrunken_dispersion(x2.T, n2.T, level2_raw, prior2)

    # wt <- 1/(1 + (n-1)*phi); wt <- wt/mean(wt)
    # R recycles phi down the columns of the (sites, reps) matrix, so each site
    # gets its own phi, and mean(wt) averages over every entry.
    wt1 = 1.0 / (1.0 + (n1 - 1.0) * phi1[:, None])
    wt1 = wt1 / np.mean(wt1)
    wt2 = 1.0 / (1.0 + (n2 - 1.0) * phi2[:, None])
    wt2 = wt2 / np.mean(wt2)

    n1_wt = n1 * wt1
    n2_wt = n2 * wt2
    level1 = mean_fn(x1 * wt1, n1_wt, *mean_args)
    level2 = mean_fn(x2 * wt2, n2_wt, *mean_args)

    e1 = level1[:, 0]
    e2 = level2[:, 0]
    if smoothing:
        var1 = _variance_term(e1, n1_wt.T, phi1, smoothing_span, groups, pos_sorted)
        var2 = _variance_term(e2, n2_wt.T, phi2, smoothing_span, groups, pos_sorted)
        vv = _guard(var1 + var2)
    else:
        # compute.waldStat.noSmooth: per-replicate sums, no window, and a
        # different floor (1e-05, one-sided) -- see _guard_no_smooth.
        var1 = _variance_term_no_smooth(e1, n1_wt.T, phi1)
        var2 = _variance_term_no_smooth(e2, n2_wt.T, phi2)
        vv = _guard_no_smooth(var1 + var2)
    se = np.sqrt(vv)
    stat = (e1 - e2) / se

    def back(a):
        out = np.full(n_sites, np.nan, dtype=float)
        out[kept] = np.asarray(a, dtype=float)
        return out

    return {
        "diff": back(e1 - e2),
        "se": back(se),
        "stat": back(stat),
        "pvalue": back(2.0 * norm.sf(np.abs(stat))),
        "phi1": back(phi1),
        "phi2": back(phi2),
        "var1": back(var1),
        "var2": back(var2),
        "level1": back(level1_raw[:, 0]),
        "level2": back(level2_raw[:, 0]),
        "level1_weighted": back(e1),
        "level2_weighted": back(e2),
        "tested": kept_mask,
        "prior1": prior1,
        "prior2": prior2,
        "smoothing_span": int(smoothing_span),
        "equal_disp": bool(equal_disp),
        "smoothing": bool(smoothing),
        "dispersion_method": "dss_logN_shrinkage",
        "reference": (
            "DSS 2.58.0 DMLtest(smoothing=TRUE)" if smoothing
            else "DSS 2.58.0 DMLtest(smoothing=FALSE)"
        ),
    }
