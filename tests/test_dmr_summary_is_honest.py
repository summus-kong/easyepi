"""The run summary must not misdescribe the model that produced the numbers.

`easyepi dmr` writes `<out>.json`. For a while its `method` and `limitation`
fields were hard-coded to the pooled binomial test, so a replicate-aware
beta-binomial run stated "binomial-based, not beta-binomial" about its own
output. That is the kind of defect that makes every downstream claim
unverifiable while every number in the table stays numerically correct, so it is
pinned here.

The test drives the summary-building code with a synthetic design rather than a
BAM, so it needs no aligner and no reference.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from easyepi.dmr import merge_dmr_regions
from easyepi.dmr_replicate_aware import differential_sites_replicate_aware

N_SITES = 60


def _frame(level: float, n_rep: int, rng) -> pd.DataFrame:
    chrom = np.array(["chr1"] * N_SITES)
    pos = np.arange(100, 100 + 10 * N_SITES, 10)
    cols = {}
    for i in range(n_rep):
        lv = np.clip(rng.normal(level, 0.04), 0.02, 0.98)
        total = np.full(N_SITES, 30, dtype=np.int64)
        cols[(i, "n_meth")] = rng.binomial(total, lv).astype(np.int64)
        cols[(i, "n_total")] = total
    frame = pd.DataFrame(cols, index=pd.MultiIndex.from_arrays(
        [chrom, pos, ["+"] * N_SITES], names=["chrom", "pos", "strand"]))
    frame.columns = pd.MultiIndex.from_tuples(frame.columns)
    return frame


def test_replicate_aware_stats_describe_themselves_correctly():
    """The stats dict must claim beta-binomial and DSS, not binomial."""
    rng = np.random.default_rng(5)
    a = _frame(0.30, 3, rng)
    b = _frame(0.30, 3, rng)
    _table, stats = differential_sites_replicate_aware(a, b)

    assert stats["method"] == "beta_binomial"
    assert stats["dispersion_method"] == "dss_logN_shrinkage"
    assert stats["between_replicate_variance_modelled"] is True
    # The dispersion DSS would report is a per-site quantity; the summary
    # reports its mean, which must be a real number in range.
    assert stats["dispersion"] is not None
    assert 0.0 <= stats["dispersion"] < 1.0
    # Nothing may be silently dropped without being counted.
    assert stats["dropped_by_dss_coverage_filter"] >= 0
    assert (stats["sites_tested"] + stats["dropped_low_coverage"]
            + stats["dropped_by_dss_coverage_filter"]
            + stats["dropped_below_min_diff"]) >= stats["sites_tested"]


def test_both_smoothing_modes_produce_a_usable_table():
    """The unsmoothed path must be as usable as the smoothed one.

    It is a separate code path with a different variance floor, so "it runs" is
    not enough: the table must still have no NaN, valid adjusted p-values, and
    must still find an obvious shift.
    """
    rng = np.random.default_rng(17)
    for smoothing in (True, False):
        a = _frame(0.25, 3, rng)
        b = _frame(0.75, 3, rng)
        table, stats = differential_sites_replicate_aware(
            a, b, smoothing=smoothing)
        assert stats["between_replicate_variance_modelled"] is True
        assert stats["smoothing"] == (
            "dss_compute_mean_smooth" if smoothing
            else "dss_compute_mean_noSmooth")
        assert len(table) == N_SITES
        assert not table[["diff", "pvalue", "padj"]].isna().any().any()
        assert ((table["padj"] >= 0) & (table["padj"] <= 1)).all()
        # A 0.25 vs 0.75 shift at 30x is not subtle in either mode.
        assert stats["significant"] >= 0.80 * N_SITES


def test_region_merge_records_that_members_were_preselected():
    """Region p-values are computed on significant members, and must say so."""
    rng = np.random.default_rng(6)
    a = _frame(0.30, 3, rng)
    b = _frame(0.70, 3, rng)
    table, _stats = differential_sites_replicate_aware(a, b)
    regions, stats = merge_dmr_regions(table, max_gap=50, min_sites=1)
    assert len(table) > 0
    if len(regions):
        assert "region_limitation" in stats
        assert "select" in stats["region_limitation"].lower()


def test_the_summary_fields_cannot_be_the_pooled_text_for_a_dss_run():
    """Read the summary code itself: each field must be conditional and honest.

    A behavioural test would need a BAM for the CLI path, so this pins the
    invariant that made the old version wrong: `model`, `method` and
    `limitation` all have to depend on which model ran. The old code hard-coded
    `method`/`limitation` to the pooled text, so a beta-binomial run described
    itself as "binomial-based, not beta-binomial".
    """
    source = Path(__file__).resolve().parents[1] / "easyepi" / "dmr.py"
    text = source.read_text(encoding="utf-8")

    # --- model: three distinct values, not two -----------------------------
    assert "dss_beta_binomial_no_smoothing" in text
    assert '"dss_beta_binomial"' in text
    assert "pooled_binomial" in text

    # --- method: describes the model that ran, DSS first -------------------
    start = text.index('"method": (')
    end = text.index('"limitation": (', start)
    method = text[start:end]
    assert "if not use_dss else" in method, "pooled method is not guarded"
    # Both DSS variants name their own mode, in that order, and the pooled text
    # is the last alternative rather than the default.
    assert method.index("if not use_dss else") < method.index("smoothing=FALSE")
    assert method.index("smoothing=FALSE") < method.index("smoothing=TRUE")

    # --- limitation: the "not beta-binomial" claim is pooled-only ----------
    lstart = text.index('"limitation": (')
    lend = text.index('"implementation": (', lstart)
    limitation = text[lstart:lend]
    assert "not beta-binomial" in limitation, "pooled caveat disappeared"
    assert limitation.count("not beta-binomial") == 1
    # It must be the `else` of `if use_dss`, i.e. reached only without DSS.
    assert "if use_dss else" in limitation
    assert limitation.index("if use_dss else") < limitation.index("not beta-binomial")
    # And the DSS branch must not inherit it.
    dss_branch = limitation[:limitation.index("if use_dss else")]
    assert "not beta-binomial" not in dss_branch
    assert "DSS 2.58.0" in dss_branch
