"""Replicate-aware differential methylation via DSS's beta-binomial GLM.

Why this module exists
----------------------
`dmr.group_counts` sums the per-site counts of a group and `dmr.differential_sites`
then treats that sum as ONE binomial observation. Biological replicates are not
technical replicates, so between-replicate variance is discarded and the
effective sample size is overstated.

Measured on simulated data (2000 sites, 30x, 3 replicates/group, BH<0.05):

    regime                              pooled calls
    A  replicates agree, real effect            19
    B  replicates scatter, real effect         516
    C  replicates scatter, NO effect           330

Regime C is decisive: with no true difference anywhere, pooling still calls
16.5% of sites. That is a false-positive rate, not a power difference, and no
tuning of the pooled Wald test removes it, because the missing quantity --
between-replicate variance -- is not in that model at all.

This module fits the model that does contain it: DSS's beta-binomial GLM with
dispersion estimated from the replicates. It is used only when both groups have
two or more replicates; with one replicate per group the dispersion is not
identifiable and the caller must say so instead of guessing.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

#: Below this many replicates per group a dispersion cannot be estimated.
MIN_REPLICATES = 2


def can_estimate_dispersion(n_a: int, n_b: int) -> bool:
    return n_a >= MIN_REPLICATES and n_b >= MIN_REPLICATES


def _stack(frame: pd.DataFrame, samples: list[int], field: str) -> np.ndarray:
    """(n_samples, n_sites) int matrix for one count field."""
    return np.vstack([frame[(s, field)].to_numpy(np.int64) for s in samples])


def differential_sites_replicate_aware(
    per_sample_a: pd.DataFrame,
    per_sample_b: pd.DataFrame,
    *,
    min_coverage: int = 5,
    min_diff: float = 0.0,
    alpha: float = 0.05,
    smoothing: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """DSS beta-binomial test on per-sample counts.

    ``per_sample_a`` / ``per_sample_b`` come from
    ``group_counts(..., per_sample=True)``: one column per sample under a
    ``n_meth``/``n_total`` second level, indexed by (chrom, pos, strand).
    """
    from easyepi.core.stats import benjamini_hochberg
    from easyepi.dmr_beta_binomial import beta_binomial_test

    n_a = int(per_sample_a.shape[1] // 2)
    n_b = int(per_sample_b.shape[1] // 2)
    if not can_estimate_dispersion(n_a, n_b):
        raise ValueError(
            "the replicate-aware test needs >=2 samples per group to estimate "
            f"dispersion; got {n_a} and {n_b}. Use the pooled test and state "
            "that between-replicate variance is not modelled.")

    joined = pd.concat({"a": per_sample_a, "b": per_sample_b}, axis=1,
                       join="outer").fillna(0.0)
    n_sites = len(joined)
    if n_sites == 0:
        return pd.DataFrame(), {"method": "beta_binomial", "sites_tested": 0, "significant": 0}

    a_samples = list(range(n_a))
    b_samples = list(range(n_b))
    meth = np.vstack([_stack(joined["a"], a_samples, "n_meth"),
                      _stack(joined["b"], b_samples, "n_meth")])
    total = np.vstack([_stack(joined["a"], a_samples, "n_total"),
                       _stack(joined["b"], b_samples, "n_total")])

    # Per-group coverage filter: a group whose mean depth is below min_coverage
    # cannot support a level estimate, and letting it through would let a single
    # read drive a call.
    cov_a = total[:n_a].sum(axis=0) / n_a
    cov_b = total[n_a:].sum(axis=0) / n_b
    keep = (cov_a >= min_coverage) & (cov_b >= min_coverage)
    kept = int(keep.sum())
    if kept == 0:
        return pd.DataFrame(), {"method": "beta_binomial", "sites_tested": 0, "significant": 0,
                                "dropped_low_coverage": n_sites}

    design = np.concatenate([np.zeros(n_a, dtype=int), np.ones(n_b, dtype=int)])
    labels = list(joined.index[keep])
    # Real coordinates: the R bridge used to fabricate chr1 with positions
    # 1..n_sites, which made every smoothing window fictional (F04). Windows
    # are genomic (span/2 in bp), so they must be the real ones.
    site_chrom = np.array([lbl[0] for lbl in labels], dtype=object)
    site_pos = np.array([int(lbl[1]) for lbl in labels], dtype=np.int64)
    result = beta_binomial_test(meth[:, keep], total[:, keep], design,
                                chrom=site_chrom, pos=site_pos,
                                smoothing=bool(smoothing))

    # beta_binomial_test applies DSS's own hasCoverage filter, because the
    # sites it drops change every later smoothing window. Those sites come back
    # as NaN and must not enter the table: a NaN p-value would poison the BH
    # correction and the region merge.
    tested = np.asarray(result["tested"], dtype=bool)
    dropped_by_dss = int((~tested).sum())
    pvalue = np.asarray(result["pvalue"], dtype=float)[tested]
    diff = np.asarray(result["diff"], dtype=float)[tested]
    labels = [labels[i] for i in np.flatnonzero(tested)]
    keep_idx = np.flatnonzero(keep)[tested]

    summed_a = total[:n_a][:, keep_idx].sum(axis=0)
    summed_b = total[n_a:][:, keep_idx].sum(axis=0)
    # Integer methylated counts, taken straight from the matrix the test was
    # run on. merge_dmr_regions sums these per region; the ratios below are
    # derived from them for display, never the other way round.
    meth_a = meth[:n_a][:, keep_idx].sum(axis=0)
    meth_b = meth[n_a:][:, keep_idx].sum(axis=0)

    frame = pd.DataFrame({
        "chrom": [lbl[0] for lbl in labels],
        "pos": [int(lbl[1]) for lbl in labels],
        "strand": [lbl[2] for lbl in labels],
        "diff": diff,
        "pvalue": pvalue,
        "padj": benjamini_hochberg(pvalue),
        "mean_a": meth_a / np.maximum(1, summed_a),
        "mean_b": meth_b / np.maximum(1, summed_b),
        "m_a": meth_a,
        "n_a": summed_a,
        "m_b": meth_b,
        "n_b": summed_b,
    })
    frame = frame[np.abs(frame["diff"]) >= min_diff].copy()
    frame["significant"] = frame["padj"] < alpha

    stats = {
        "method": "beta_binomial",
        # DSS reports phi1/phi2 per site; the mean over tested sites is the
        # single number worth putting in a run summary.
        "dispersion": (
            None if not tested.any()
            else float(np.mean(np.concatenate([
                np.asarray(result["phi1"], dtype=float)[tested],
                np.asarray(result["phi2"], dtype=float)[tested]])))
        ),
        "sites_tested": int(len(frame)),
        "significant": int(frame["significant"].sum()),
        "dropped_low_coverage": int(n_sites - kept),
        "dropped_by_dss_coverage_filter": dropped_by_dss,
        "dropped_below_min_diff": int(kept - dropped_by_dss - len(frame)),
        "replicates_a": n_a,
        "replicates_b": n_b,
        "between_replicate_variance_modelled": True,
        "smoothing": ("dss_compute_mean_smooth" if smoothing
                      else "dss_compute_mean_noSmooth"),
        "dispersion_method": "dss_logN_shrinkage",
    }
    return frame.sort_values(["chrom", "pos", "strand"]).reset_index(drop=True), stats
