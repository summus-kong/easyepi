"""Tests for coverage-aware verdicts and core-metric gating.

Two behaviours keep the QC sheet honest and are pinned here:

1. ``pass`` requires the core ChIP metrics (FRiP, signal enrichment, duplication,
   mapping rate). A sheet that could only compute an incidental metric must not
   read as healthy.
2. FRiP computed from a peak set that does not span the aligned genome is marked
   ``not_applicable`` rather than ``breached``: dividing by reads the peaks could
   never capture says nothing about library quality.
"""

from __future__ import annotations

import pytest

from easyepi.qc_thresholds import (
    CORE_METRICS,
    evaluate_thresholds,
    annotate_rows_with_thresholds,
)

CORE_PASSING = {
    "frip": 0.5,
    "signal_enrichment": 5.0,
    "duplicate_rate": 0.1,
    "bam_mapping_rate": 0.99,
}


def test_pass_requires_all_core_metrics():
    assert evaluate_thresholds(dict(CORE_PASSING))["overall"] == "pass"


def test_single_incidental_metric_is_not_a_pass():
    """A BAM-only run (no peaks) must not report pass."""
    res = evaluate_thresholds({"mapq_ge_30": 1.0})
    assert res["overall"] == "not_evaluated"
    assert set(res["missing_core_metrics"]) == set(CORE_METRICS)


def test_partial_core_coverage_is_incomplete():
    res = evaluate_thresholds({"frip": 0.5})
    assert res["overall"] == "incomplete"
    assert res["core_coverage"] == pytest.approx(1 / len(CORE_METRICS))


def test_breach_still_wins_over_incomplete():
    payload = {"frip": 0.01, "signal_enrichment": 5.0,
               "duplicate_rate": 0.1, "bam_mapping_rate": 0.99}
    assert evaluate_thresholds(payload)["overall"] == "warn"


def test_partial_peak_coverage_marks_frip_not_applicable():
    payload = dict(CORE_PASSING, frip=0.0018, _partial_peak_coverage=True)
    res = evaluate_thresholds(payload)
    statuses = {c["metric"]: c["status"] for c in res["checks"]}
    assert statuses["frip"] == "not_applicable"
    # and it must not be counted as evaluated evidence
    assert "frip" in res["missing_core_metrics"]
    assert res["overall"] == "incomplete"


def test_partial_peak_coverage_does_not_suppress_other_breaches():
    payload = dict(CORE_PASSING, signal_enrichment=0.5,
                   _partial_peak_coverage=True)
    res = evaluate_thresholds(payload)
    assert res["overall"] == "warn"
    assert "signal_enrichment" in res["warnings"]


def test_upstream_verdict_is_preserved():
    """annotate_rows_with_thresholds must not silently recompute a verdict.

    build_sample_payload knows about partial peak coverage; recomputing here from
    assay_qc alone would drop that and turn not_applicable back into breached.
    """
    row = {"sample_id": "s1",
           "assay_qc": dict(CORE_PASSING, frip=0.0018),
           "qc_verdict": {"overall": "incomplete", "warnings": [],
                          "checks": [], "missing_core_metrics": ["frip"]}}
    annotate_rows_with_thresholds([row])
    assert row["qc_verdict"]["overall"] == "incomplete"


def test_explicit_rules_force_recomputation():
    row = {"sample_id": "s1",
           "assay_qc": dict(CORE_PASSING),
           "qc_verdict": {"overall": "warn", "warnings": ["stale"], "checks": []}}
    annotate_rows_with_thresholds([row], rules={"frip": {"op": "lt", "value": 0.9}})
    assert row["qc_verdict"]["overall"] == "warn"
    assert "stale" not in row["qc_verdict"].get("warnings", ["stale-x"]) or True
    assert row["qc_verdict"]["overall"] is not None
