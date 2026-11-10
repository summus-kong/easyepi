"""Core statistics primitives: descriptive stats, significance tests, correction.

Implemented directly with numpy/scipy so the whole statistical layer is
dependency-free beyond the hard requirement set. These primitives are shared by
the qc / normalize / differential modules.
"""

from __future__ import annotations

from typing import Iterable, List, Optional

import numpy as np
from scipy import stats as _st


# ---------------------------------------------------------------------------
# Descriptive statistics
# ---------------------------------------------------------------------------
def describe(values: Iterable[float]) -> dict:
    """Return descriptive statistics for a 1-D array."""
    a = np.asarray(values, dtype=np.float64)
    a = a[~np.isnan(a)]
    if a.size == 0:
        return {"n": 0, "mean": np.nan, "median": np.nan, "std": np.nan,
                "min": np.nan, "max": np.nan, "q25": np.nan, "q75": np.nan}
    q25, q75 = np.percentile(a, [25, 75])
    return {
        "n": int(a.size),
        "mean": float(np.mean(a)),
        "median": float(np.median(a)),
        "std": float(np.std(a, ddof=1)) if a.size > 1 else 0.0,
        "min": float(np.min(a)),
        "max": float(np.max(a)),
        "q25": float(q25),
        "q75": float(q75),
        "iqr": float(q75 - q25),
    }


def correlation(x: Iterable[float], y: Iterable[float], method: str = "pearson") -> float:
    """Pearson and Spearman correlation with an explicit method choice."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if method == "spearman":
        return float(_st.spearmanr(x, y)[0])
    return float(_st.pearsonr(x, y)[0])


# ---------------------------------------------------------------------------
# Significance tests
# ---------------------------------------------------------------------------
def t_test(a: Iterable[float], b: Iterable[float], equal_var: bool = True) -> dict:
    """Welch or Student t-test between two groups; returns statistic and p-value."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return {"statistic": np.nan, "pvalue": np.nan, "test": "t-test",
                "mean_a": float(np.mean(a)) if a.size else np.nan,
                "mean_b": float(np.mean(b)) if b.size else np.nan,
                "n_a": int(a.size), "n_b": int(b.size)}
    stat, p = _st.ttest_ind(a, b, equal_var=equal_var, nan_policy="omit")
    return {"statistic": float(stat), "pvalue": float(p), "test": "t-test",
            "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "n_a": int(a.size), "n_b": int(b.size),
            "equal_var": bool(equal_var)}


def wilcoxon(a: Iterable[float], b: Iterable[float]) -> dict:
    """Paired Wilcoxon signed-rank test; falls back to nan when undefined."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size != b.size:
        raise ValueError("paired Wilcoxon requires equal-length vectors")
    finite = np.isfinite(a) & np.isfinite(b)
    a = a[finite]
    b = b[finite]
    if a.size < 2:
        return {"statistic": np.nan, "pvalue": np.nan, "test": "wilcoxon",
                "mean_a": float(np.mean(a)) if a.size else np.nan,
                "mean_b": float(np.mean(b)) if b.size else np.nan,
                "n_pairs": int(a.size)}
    try:
        stat, p = _st.wilcoxon(a, b, zero_method="wilcox", correction=False)
    except ValueError:
        return {"statistic": np.nan, "pvalue": np.nan, "test": "wilcoxon",
                "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                "n_pairs": int(a.size)}
    return {"statistic": float(stat), "pvalue": float(p), "test": "wilcoxon",
            "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "n_pairs": int(a.size)}


def mann_whitney_u(a: Iterable[float], b: Iterable[float]) -> dict:
    """Two-sided Mann-Whitney U test for unpaired groups."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return {"statistic": np.nan, "pvalue": np.nan, "test": "mann-whitney-u",
                "mean_a": float(np.mean(a)) if a.size else np.nan,
                "mean_b": float(np.mean(b)) if b.size else np.nan,
                "n_a": int(a.size), "n_b": int(b.size)}
    stat, p = _st.mannwhitneyu(a, b, alternative="two-sided")
    return {"statistic": float(stat), "pvalue": float(p), "test": "mann-whitney-u",
            "mean_a": float(a.mean()), "mean_b": float(b.mean()),
            "n_a": int(a.size), "n_b": int(b.size)}


def fisher_exact(table: Iterable[int]) -> dict:
    """2x2 Fisher exact test: [[a,b],[c,d]]."""
    t = np.asarray(table, dtype=np.int64)
    if t.shape != (2, 2):
        raise ValueError("fisher_exact expects a 2x2 table")
    odds, p = _st.fisher_exact(t)
    return {"odds_ratio": float(odds), "pvalue": float(p), "test": "fisher-exact"}


def hypergeom_survival(k: int, N: int, n: int, K: int) -> float:
    """P-value of observing >=k overlaps under the hypergeometric null.

    Args:
        k: number of overlapping items.
        N: population size.
        n: number of draws (e.g. query peak count).
        K: number of success states in population (e.g. annotated gene count).
    """
    dist = _st.hypergeom(M=int(N), n=int(n), N=int(K))
    return float(dist.sf(k - 1))


def cohens_d(a: Iterable[float], b: Iterable[float]) -> float:
    """Standardised mean difference (pooled SD) as an effect size."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return np.nan
    sp = np.sqrt(((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1))
                 / (a.size + b.size - 2))
    if sp == 0:
        return 0.0
    return float((a.mean() - b.mean()) / sp)


# ---------------------------------------------------------------------------
# Multiple testing correction
# ---------------------------------------------------------------------------
def benjamini_hochberg(pvalues: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR. Returns adjusted p-values (same length/order)."""
    p = np.asarray(pvalues, dtype=np.float64)
    n = p.size
    if n == 0:
        return p
    if np.any((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)):
        raise ValueError("p-values must lie in [0, 1] or be NaN")
    finite = np.isfinite(p)
    adj = np.full(n, np.nan, dtype=np.float64)
    if not finite.any():
        return adj
    # Missing/non-finite p-values are excluded from the multiplicity burden;
    # their positions remain NaN so a failed test cannot look significant.
    idx = np.flatnonzero(finite)
    order_local = np.argsort(p[idx])
    order = idx[order_local]
    m = order.size
    ranks = np.arange(1, m + 1)
    # BH: q_i = p_(i) * n / rank_i, in sorted order
    q_sorted = p[order] * m / ranks
    # enforce monotone non-decreasing from the largest p downwards
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    # place back into original order
    adj[order] = q_sorted
    return adj


def log2fc(a: Iterable[float], b: Iterable[float]) -> float:
    """log2 fold change between two group means, with a pseudocount floor."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    m_a, m_b = a.mean(), b.mean()
    if m_a == 0 and m_b == 0:
        return 0.0
    return float(np.log2((m_a + 1.0) / (m_b + 1.0)))
