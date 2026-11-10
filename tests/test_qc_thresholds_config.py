"""Tests for configurable QC thresholds: overrides, precedence, and advisory-only severity."""

import json

import pytest

from easyepi.qc_thresholds import (
    DEFAULT_THRESHOLDS,
    evaluate_thresholds,
    load_rules,
    summarize_verdicts,
)


# ---------------------------------------------------------------------------
# advisory-only severity
# ---------------------------------------------------------------------------
def test_breach_is_always_warn_never_fail():
    """A breach must never produce "fail": thresholds are advisory."""
    rules = {"frip": {"op": "lt", "value": 0.9, "severity": "fail", "message": "x"}}
    res = evaluate_thresholds({"frip": 0.05}, rules=rules)
    assert res["overall"] == "warn"
    assert res["failures"] == []


def test_overall_never_fail_even_with_many_breaches():
    payload = {"frip": 0.0, "tss_enrichment": 0.1, "blacklist_fraction": 0.5,
               "duplicate_rate": 0.9}
    res = evaluate_thresholds(payload)
    assert res["overall"] == "warn"
    assert len(res["warnings"]) == 4
    assert res["failures"] == []


def test_warning_details_carry_values():
    res = evaluate_thresholds({"frip": 0.01})
    detail = next(d for d in res["warning_details"] if d["metric"] == "frip")
    assert detail["value"] == pytest.approx(0.01)
    assert detail["threshold"] == pytest.approx(0.20)
    assert detail["op"] == "lt"


# ---------------------------------------------------------------------------
# load_rules: defaults, file overrides, precedence, robustness
# ---------------------------------------------------------------------------
def test_load_rules_returns_defaults_when_no_source():
    rules = load_rules(None)
    assert rules["frip"]["value"] == DEFAULT_THRESHOLDS["frip"]["value"]
    assert len(rules) == len(DEFAULT_THRESHOLDS)


def test_load_rules_from_json_file_merges_per_metric(tmp_path):
    path = tmp_path / "rules.json"
    path.write_text(json.dumps({"frip": {"op": "lt", "value": 0.35}}), encoding="utf-8")
    rules = load_rules(str(path))
    assert rules["frip"]["value"] == pytest.approx(0.35)
    # untouched metrics keep their defaults
    assert rules["tss_enrichment"]["value"] == DEFAULT_THRESHOLDS["tss_enrichment"]["value"]


def test_load_rules_accepts_wrapped_qc_thresholds_key(tmp_path):
    path = tmp_path / "rules.json"
    path.write_text(json.dumps({"qc_thresholds": {"blacklist_fraction":
                                                  {"op": "gt", "value": 0.05}}}),
                    encoding="utf-8")
    assert load_rules(str(path))["blacklist_fraction"]["value"] == pytest.approx(0.05)


def test_load_rules_accepts_parsed_mapping():
    rules = load_rules({"duplicate_rate": {"op": "gt", "value": 0.7}})
    assert rules["duplicate_rate"]["value"] == pytest.approx(0.7)


def test_load_rules_later_source_wins():
    first = {"frip": {"op": "lt", "value": 0.1}}
    second = {"frip": {"op": "lt", "value": 0.4}}
    assert load_rules(first, second)["frip"]["value"] == pytest.approx(0.4)


def test_load_rules_ignores_missing_file(tmp_path):
    rules = load_rules(str(tmp_path / "nope.json"))
    assert rules["frip"]["value"] == DEFAULT_THRESHOLDS["frip"]["value"]


def test_load_rules_ignores_malformed_entries(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({
        "frip": {"op": "between", "value": 0.2},        # invalid op
        "dup": {"op": "gt", "value": "not-a-number"},    # invalid value
        "tss_enrichment": {"op": "lt", "value": 3.0},    # valid
    }), encoding="utf-8")
    rules = load_rules(str(path))
    assert rules["frip"]["value"] == DEFAULT_THRESHOLDS["frip"]["value"]
    assert "dup" not in rules
    assert rules["tss_enrichment"]["value"] == pytest.approx(3.0)


def test_load_rules_ignores_broken_json(tmp_path):
    path = tmp_path / "broken.json"
    path.write_text("{not json", encoding="utf-8")
    assert load_rules(str(path))["frip"]["value"] == DEFAULT_THRESHOLDS["frip"]["value"]


def test_overridden_threshold_changes_verdict():
    # Provide the full core set so the verdict reflects the FRiP rule alone;
    # otherwise a missing-metric "incomplete" would mask the change.
    core = {"frip": 0.15, "signal_enrichment": 5.0,
            "duplicate_rate": 0.1, "bam_mapping_rate": 0.99}
    default = evaluate_thresholds(dict(core))
    assert default["overall"] == "warn"          # 0.15 < 0.20 default

    relaxed = evaluate_thresholds(dict(core), rules=load_rules(
        {"frip": {"op": "lt", "value": 0.10}}))
    assert relaxed["overall"] == "pass"          # 0.15 >= 0.10


# ---------------------------------------------------------------------------
# batch roll-up
# ---------------------------------------------------------------------------
def test_summarize_merges_multiple_batches():
    batch_a = [{"sample_id": "a1", "qc_verdict": {"overall": "warn", "warnings": ["frip"]}},
               {"sample_id": "a2", "qc_verdict": {"overall": "pass", "warnings": []}}]
    batch_b = [{"sample_id": "b1", "qc_verdict": {"overall": "warn", "warnings": ["frip", "tss_enrichment"]}}]
    merged = summarize_verdicts(batch_a + batch_b)
    assert merged["samples"]["pass"] == 1
    assert merged["samples"]["warn"] == 2
    assert merged["warning_metrics"]["frip"] == 2
    assert merged["warning_metrics"]["tss_enrichment"] == 1


def test_summarize_has_no_fail_bucket():
    merged = summarize_verdicts([{"qc_verdict": {"overall": "warn", "warnings": ["frip"]}}])
    assert "fail" not in merged["samples"]
