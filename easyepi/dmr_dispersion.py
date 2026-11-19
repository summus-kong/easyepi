"""DSS's shrunken dispersion estimator, in numpy/scipy.

Transcribed from DSS 2.58.0 so the pure-Python beta-binomial test reports the
same dispersion DSS would, instead of the method-of-moments approximation that
made ``diff`` systematically wrong (measured: 0.0727 vs 0.0431 on the 300-site
fixture).

Both functions below are taken line by line from the installed package, which
is the only authority here; the R source is reproduced in each docstring.

``est.prior.BSseq.logN(X, N)``

.. code-block:: r

    if (ncol(X) == 1) return(c(-3, 1))
    ix = rowMeans(N > 10) == 1 & rowSums(N == 0) == 0
    if (sum(ix) < 50) return(c(-3, 1))
    X = X[ix, ]; N = N[ix, ]
    p  = X/N
    mm = rowMeans(p); mm[mm == 0] = 1e-05; mm[mm == 1] = 1 - 1e-05
    vv = rowVars(p)                    # divides by n-1
    phi = vv/mm/(1 - mm)
    phi = phi[vv > 0]
    lphi = log(phi[phi > 0])
    prior.mean = median(lphi, na.rm = TRUE)
    prior.sd   = IQR(lphi, na.rm = TRUE)/1.39
    c(prior.mean, prior.sd)

Two details that were wrong in the first attempt and are easy to get wrong
again: the prior sd is ``IQR/1.39``, **not** ``1.4826 * IQR`` and there is no
``sqrt(2/n)`` floor; and ``phi`` is filtered on ``vv > 0`` *before* the log, so
the prior is fitted to the log of the positive dispersions only.

``dispersion.shrinkage.BSseq(X, N, prior, estprob, ncores)``

.. code-block:: r

    plik.logN <- function(size, X, mu, m0, tau, phi)
        -(sum(dbb(size, X, mu, exp(phi))) + dnorm(phi, mean = m0, sd = tau, log = TRUE))
    shrk.phi = exp(rep(prior[1], nrow(N)))          # prior[1] is the mean
    ix <- rowSums(N > 0) > 0
    for (i in 1:nrow(X2))
        shrk.phi2[i] = exp(optimize(plik.logN, size = N2[i, ], X = X2[i, ],
                                    mu = estprob2[i, ], m0 = prior[1], tau = prior[2],
                                    interval = c(-5, log(0.99)), tol = 0.001)$minimum)

``estprob2[i, ]`` is the whole row of per-sample smoothed levels. ``dbb`` takes
the level as a vector, so a site whose per-sample levels differ would see a
vector ``mu``; ``compute.mean.Smooth`` returns every column equal, which is why
a scalar reproduces DSS here. That is a property of the caller, not of this
function, so the row is passed through as a vector.

``ncores`` only selects the loop (serial vs ``mclapply``); it never enters the
arithmetic. ``const <- 1/ncores`` lives in ``compute.var.Smooth`` and is a
different quantity -- see :mod:`easyepi.dmr_beta_binomial`.
"""
from __future__ import annotations

import numpy as np
from scipy import optimize, special, stats

#: Prior used when the coverage is too thin to fit one (R: ``c(-3, 1)``).
FALLBACK_PRIOR = (-3.0, 1.0)

#: Optimisation interval on ``log(phi)``, from the R source.
LOG_PHI_BOUNDS = (-5.0, float(np.log(0.99)))

#: R's ``optimize(..., tol = 0.001)``.
OPT_TOL = 1e-3


def dbb_logpmf(size: np.ndarray, x: np.ndarray, mu, phi: float) -> np.ndarray:
    """``DSS:::dbb(size, x, mu, phi, log = TRUE)``, vectorised over samples.

    ``mu`` may be a scalar or an array broadcastable against ``size``.
    """
    tmp = 1.0 / phi - 1.0
    mu = np.asarray(mu, dtype=float)
    alpha = mu * tmp
    beta = tmp - alpha
    return (special.gammaln(size + 1.0) - special.gammaln(x + 1.0)
            - special.gammaln(size - x + 1.0)
            - special.betaln(beta, alpha)
            + special.betaln(size - x + beta, x + alpha))


def estimate_prior_log_dispersion(x: np.ndarray, n: np.ndarray) -> tuple:
    """``DSS:::est.prior.BSseq.logN``: ``(prior.mean, prior.sd)`` on the log scale.

    ``x`` and ``n`` are ``(n_sites, n_samples)``, which is the layout
    ``est.dispersion.BSseq`` receives from ``DMLtest.Smooth``.
    """
    x = np.asarray(x, dtype=float)
    n = np.asarray(n, dtype=float)
    if x.ndim != 2:
        raise ValueError("x must be (n_sites, n_samples)")
    if x.shape[1] == 1:
        return FALLBACK_PRIOR

    # rowMeans(N > 10) == 1 & rowSums(N == 0) == 0
    ix = (n > 10).all(axis=1) & ((n == 0).sum(axis=1) == 0)
    if int(ix.sum()) < 50:
        return FALLBACK_PRIOR

    xs = x[ix]
    ns = n[ix]
    with np.errstate(invalid="ignore", divide="ignore"):
        p = xs / ns
    mm = p.mean(axis=1)
    mm = np.where(mm == 0, 1e-5, mm)
    mm = np.where(mm == 1, 1 - 1e-5, mm)
    vv = p.var(axis=1, ddof=1)          # R rowVars divides by n-1
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = vv / mm / (1.0 - mm)

    # phi = phi[vv > 0]; lphi = log(phi[phi > 0])
    lphi = np.log(phi[(vv > 0) & (phi > 0)])
    if lphi.size == 0:
        return FALLBACK_PRIOR

    prior_mean = float(np.median(lphi))
    q25, q75 = np.percentile(lphi, [25.0, 75.0])   # R IQR uses the same rule
    prior_sd = float((q75 - q25) / 1.39)
    return prior_mean, prior_sd


def _shrink_one(x_row: np.ndarray, n_row: np.ndarray, mu_row,
                prior_mean: float, prior_sd: float) -> float:
    """``optimize(plik.logN, ..., interval = c(-5, log(0.99)))`` for one site.

    ``mu_row`` is recycled against ``x_row`` the way R recycles ``mu`` inside
    ``dbb``: under ``equal.disp`` DSS passes a two-column ``estprob`` against
    six count columns, and R silently repeats it.
    """
    x_row = np.asarray(x_row, dtype=float)
    n_row = np.asarray(n_row, dtype=float)
    mu = np.resize(np.asarray(mu_row, dtype=float), x_row.shape)

    def neg_log_posterior(log_phi: float) -> float:
        phi = float(np.exp(log_phi))
        tmp = 1.0 / phi - 1.0
        alpha = mu * tmp
        beta = tmp - alpha
        if (np.any(alpha <= 0) or np.any(beta <= 0)
                or not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(beta))):
            return np.inf
        ll = dbb_logpmf(n_row, x_row, mu, phi).sum()
        ll += stats.norm.logpdf(log_phi, loc=prior_mean, scale=prior_sd)
        return -float(ll) if np.isfinite(ll) else np.inf

    res = optimize.minimize_scalar(neg_log_posterior, bounds=LOG_PHI_BOUNDS,
                                   method="bounded", options={"xatol": OPT_TOL})
    if not res.success or not np.isfinite(res.fun):
        return float(np.exp(prior_mean))
    return float(np.exp(float(res.x)))


def shrunken_dispersion(x: np.ndarray, n: np.ndarray, level: np.ndarray,
                        prior: tuple) -> np.ndarray:
    """``DSS:::dispersion.shrinkage.BSseq``: one dispersion per site.

    ``x``, ``n`` are ``(n_reps, n_sites)`` and ``level`` is the per-site
    smoothed level used as ``mu``. ``level`` may be ``(n_sites,)`` or
    ``(n_sites, n_reps)``; DSS passes ``estprob2[i, ]``, the whole row, so the
    two-dimensional form is the faithful one. It only differs from the scalar
    form where the columns disagree, which cannot happen within one group
    (``compute.mean.Smooth`` repeats the same column) but can under
    ``equal.disp``, which ``cbind``s two groups.

    Sites with no coverage in any sample keep ``exp(prior.mean)``, exactly as
    ``shrk.phi`` does in R.
    """
    x = np.asarray(x, dtype=float)
    n = np.asarray(n, dtype=float)
    prior_mean, prior_sd = float(prior[0]), float(prior[1])
    if not np.isfinite(prior_sd) or prior_sd <= 0:
        prior_sd = 1.0

    x_sites = x.T
    n_sites_arr = n.T
    level = np.asarray(level, dtype=float)
    if level.ndim == 1:
        level = level[:, None]

    phi = np.full(x_sites.shape[0], float(np.exp(prior_mean)), dtype=float)
    covered = (n_sites_arr > 0).any(axis=1)
    for i in np.flatnonzero(covered):
        phi[i] = _shrink_one(x_sites[i], n_sites_arr[i], level[i],
                             prior_mean, prior_sd)
    return phi
