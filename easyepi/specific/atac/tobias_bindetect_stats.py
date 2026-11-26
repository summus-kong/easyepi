"""1:1 reimplementation of TOBIAS BINDetect's statistics layer.

Two pieces beyond the normalisation already ported:

**Bound/unbound threshold** (`bindetect.py` 470-558)
    bg = background scores, drop zeros, drop values above the 99th percentile
    fit a 2-component GaussianMixture on log(bg), keep the component with the
    larger mean
    lognorm.fit(bg, f0=std, fscale=exp(mean)) -> mode via fmin of -pdf
    pseudo = mode / 2
    mirror the left side of that lognorm around the mode, fit a normal to the
    mirrored curve, then threshold = norm.ppf(1 - bound_pvalue, mode, std)

**Differential binding** (`bindetect_functions.py` 500-548)
    per comparison: background log2fc = log2((s1+p)/(s2+p)) on sites where either
    score > 0, trimmed to the 1st-99th percentile, then a normal fit gives
    (bg_mean, bg_std)
    per TF: observed log2fc averaged per peak (`groupby("peak_id")`)
    effect size = (obs_mean - bg_mean) / mean([obs_std, bg_std]), rounded to 5
    p-value = ttest_1samp over 100 draws from the background normal, seeded with
    `np.random.seed(n_obs)` where `n_obs = len(observed)`

Both depend on scipy/sklearn in the same way the reference does, so the numeric
behaviour should match rather than merely resemble.
"""
from __future__ import annotations

import numpy as np
import scipy.optimize
import scipy.stats
from sklearn.mixture import GaussianMixture

#: Reference caps the background filter at this percentile.
BACKGROUND_PERCENTILE = 99

#: Reference uses a 2-component mixture on log(scores): zeros vs signal.
GMM_COMPONENTS = 2
GMM_RANDOM_STATE = 1

#: `pseudo = mode / 2`.
PSEUDO_FRACTION = 0.5

#: Number of background samples drawn for the t-test.
N_SAMPLES = 100


def estimate_background_distribution(scores):
    """Fit the lognormal mode of the non-zero background scores.

    Returns ``(mode, lognorm_params, pseudo)``.
    """
    bg = np.asarray(scores, dtype=float)
    bg = bg[~np.isclose(bg, 0.0)]
    if bg.size == 0:
        raise ValueError("no non-zero background scores")

    x_max = np.percentile(bg, [BACKGROUND_PERCENTILE])
    bg = bg[bg < x_max]
    if bg.size == 0:
        raise ValueError("no background scores below the 99th percentile")

    log_vals = np.log(bg).reshape(-1, 1)
    gmm = GaussianMixture(n_components=GMM_COMPONENTS, random_state=GMM_RANDOM_STATE)
    gmm.fit(log_vals)

    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_).flatten()

    chosen = int(np.argmax(means))
    log_params = scipy.stats.lognorm.fit(bg, f0=stds[chosen],
                                         fscale=np.exp(means[chosen]))
    mode = scipy.optimize.fmin(
        lambda x: -scipy.stats.lognorm.pdf(x, *log_params), 0, disp=False)[0]
    return float(mode), log_params, float(mode * PSEUDO_FRACTION)


def estimate_threshold(scores, bound_pvalue=1e-3):
    """Bound/unbound threshold, from the mirrored left side of the lognormal."""
    mode, log_params, pseudo = estimate_background_distribution(scores)

    left_x = np.linspace(scipy.stats.lognorm(*log_params).ppf([0.01]), mode, 100)
    left_pdf = scipy.stats.lognorm.pdf(left_x, *log_params)

    scaled = left_x - np.min(left_x)
    mirrored_x = np.concatenate([left_x, np.max(left_x) + scaled]).flatten()
    mirrored_pdf = np.concatenate([left_pdf, left_pdf[::-1]]).flatten()

    popt, _ = scipy.optimize.curve_fit(
        lambda x, std, sc: sc * scipy.stats.norm.pdf(x, mode, std),
        mirrored_x, mirrored_pdf)
    norm_params = (mode, popt[0])

    threshold = round(float(scipy.stats.norm.ppf(1 - bound_pvalue, *norm_params)), 5)
    return threshold, pseudo, mode


def background_log2fc_params(scores_1, scores_2, pseudo):
    """Fit the background log2 fold-change distribution for one comparison."""
    s1 = np.asarray(scores_1, dtype=float)
    s2 = np.asarray(scores_2, dtype=float)
    included = np.logical_or(s1 > 0, s2 > 0)
    s1, s2 = s1[included], s2[included]
    if s1.size == 0:
        raise ValueError("no sites with signal in either condition")

    log2fcs = np.log2(np.true_divide(s1 + pseudo, s2 + pseudo))
    lower, upper = np.percentile(log2fcs, [1, 99])
    trimmed = log2fcs[(log2fcs >= lower) & (log2fcs <= upper)]
    return scipy.stats.norm.fit(trimmed)


def mean_per_peak(peak_ids, log2fcs):
    """``groupby("peak_id").mean()`` -- multiple TFBS in one peak are averaged."""
    peak_ids = np.asarray(peak_ids)
    vals = np.asarray(log2fcs, dtype=float)
    order = np.argsort(peak_ids, kind="stable")
    ids_sorted = peak_ids[order]
    vals_sorted = vals[order]

    out_ids, out_vals = [], []
    start = 0
    for i in range(1, len(ids_sorted) + 1):
        if i == len(ids_sorted) or ids_sorted[i] != ids_sorted[start]:
            out_ids.append(ids_sorted[start])
            out_vals.append(float(np.mean(vals_sorted[start:i])))
            start = i
    return np.asarray(out_ids), np.asarray(out_vals)


def effect_size(observed, bg_params):
    """``(obs_mean - bg_mean) / mean([obs_std, bg_std])``, rounded to 5."""
    obs_mean, obs_std = scipy.stats.norm.fit(np.asarray(observed, dtype=float))
    bg_mean, bg_std = bg_params
    if obs_mean == bg_mean:
        return 0.0
    return float(np.round((obs_mean - bg_mean) / np.mean([obs_std, bg_std]), 5))


def differential_pvalue(observed, bg_params, n_samples=N_SAMPLES):
    """Sampling t-test, seeded exactly as the reference seeds it."""
    observed = np.asarray(observed, dtype=float)
    n_obs = len(observed)
    bg_mean, bg_std = bg_params

    change = effect_size(observed, bg_params)
    if change == 0.0:
        return 1.0, change

    np.random.seed(n_obs)
    sample_changes = []
    for _ in range(n_samples):
        sample = scipy.stats.norm.rvs(*bg_params, size=n_obs)
        s_mean, s_std = np.mean(sample), np.std(sample)
        sample_changes.append((s_mean - bg_mean) / np.mean([s_std, bg_std]))

    pvalue = scipy.stats.ttest_1samp(sample_changes, change)[1]
    return float(pvalue), change


def benjamini_hochberg(pvalues):
    """BH FDR, returned in the input order."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return out
