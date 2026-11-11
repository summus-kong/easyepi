"""Tests for merging multiple qc_report.json files into a batch-level summary."""

import json

import pytest

from easyepi.qc_report_merge import collect_reports, merge_reports


def _write_report(path, samples, module="qc-report", project_id=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"module": module, "samples": samples}
    if project_id:
        payload["project_id"] = project_id
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _sample(overall, warnings=()):
    return {"assay": "ChIP-seq",
            "qc_verdict": {"overall": overall, "warnings": list(warnings),
                           "checks": []}}


# ---------------------------------------------------------------------------
# collect_reports: files, directories, globs
# ---------------------------------------------------------------------------
def test_collect_direct_file(tmp_path):
    report = _write_report(tmp_path / "qc_report.json", {"s1": _sample("pass")})
    found, unresolved = collect_reports([str(report)])
    assert found == [report]
    assert unresolved == []


def test_collect_from_directory(tmp_path):
    _write_report(tmp_path / "batch_a" / "qc_report.json", {"s1": _sample("pass")})
    _write_report(tmp_path / "batch_b" / "qc_report.json", {"s2": _sample("warn", ["frip"])})
    found, unresolved = collect_reports([str(tmp_path)])
    assert len(found) == 2
    assert unresolved == []


def test_collect_from_glob(tmp_path):
    _write_report(tmp_path / "batch_a" / "qc_report.json", {"s1": _sample("pass")})
    _write_report(tmp_path / "batch_b" / "qc_report.json", {"s2": _sample("pass")})
    found, unresolved = collect_reports([str(tmp_path / "batch_*" / "qc_report.json")])
    assert len(found) == 2
    assert unresolved == []


def test_collect_deduplicates_same_file(tmp_path):
    report = _write_report(tmp_path / "qc_report.json", {"s1": _sample("pass")})
    found, _ = collect_reports([str(report), str(report), str(tmp_path)])
    assert len(found) == 1


def test_collect_reports_missing_pattern():
    found, unresolved = collect_reports(["/nonexistent/path/qc_report.json"])
    assert found == []
    assert len(unresolved) == 1


# ---------------------------------------------------------------------------
# merge_reports
# ---------------------------------------------------------------------------
def test_merge_single_report(tmp_path):
    report = _write_report(tmp_path / "qc_report.json", {
        "s1": _sample("pass"),
        "s2": _sample("warn", ["frip"]),
    })
    merged = merge_reports([report])
    assert merged["n_samples"] == 2
    assert merged["n_reports"] == 1
    assert merged["verdict_summary"]["samples"] == {"pass": 1, "warn": 1, "not_evaluated": 0}
    assert merged["verdict_summary"]["warning_metrics"] == {"frip": 1}
    assert merged["problems"] == []


def test_merge_multiple_batches_aggregates(tmp_path):
    a = _write_report(tmp_path / "a" / "qc_report.json",
                      {"s1": _sample("pass"), "s2": _sample("warn", ["tss_enrichment"])},
                      project_id="batch_a")
    b = _write_report(tmp_path / "b" / "qc_report.json",
                      {"s1": _sample("warn", ["frip"]), "s3": _sample("pass")},
                      project_id="batch_b")
    merged = merge_reports([a, b])
    assert merged["n_reports"] == 2
    assert merged["n_samples"] == 4
    assert merged["batches"] == {"batch_a": 2, "batch_b": 2}
    summary = merged["verdict_summary"]
    assert summary["samples"]["pass"] == 2
    assert summary["samples"]["warn"] == 2
    assert summary["warning_metrics"] == {"tss_enrichment": 1, "frip": 1}


def test_duplicate_sample_ids_are_namespaced(tmp_path):
    """Two batches both containing 's1' must not collapse into one sample."""
    a = _write_report(tmp_path / "a" / "qc_report.json", {"s1": _sample("pass")},
                      project_id="batch_a")
    b = _write_report(tmp_path / "b" / "qc_report.json", {"s1": _sample("warn", ["frip"])},
                      project_id="batch_b")
    merged = merge_reports([a, b])
    assert merged["n_samples"] == 2
    assert set(merged["samples"]) == {"s1", "batch_b::s1"}
    assert merged["samples"]["s1"]["batch"] == "batch_a"
    assert merged["samples"]["batch_b::s1"]["batch"] == "batch_b"


def test_malformed_report_is_reported_not_fatal(tmp_path):
    good = _write_report(tmp_path / "good" / "qc_report.json", {"s1": _sample("pass")})
    bad = tmp_path / "bad" / "qc_report.json"
    bad.parent.mkdir(parents=True)
    bad.write_text("{not json", encoding="utf-8")
    merged = merge_reports([good, bad])
    assert merged["n_samples"] == 1
    assert len(merged["problems"]) == 1
    assert "bad" in merged["problems"][0]


def test_report_without_samples_key_is_reported(tmp_path):
    path = tmp_path / "odd.json"
    path.write_text(json.dumps({"module": "x"}), encoding="utf-8")
    merged = merge_reports([path])
    assert merged["n_samples"] == 0
    assert any("samples" in p for p in merged["problems"])


def test_sample_without_verdict_counts_as_not_evaluated(tmp_path):
    report = _write_report(tmp_path / "qc_report.json", {"s1": {"assay": "ChIP-seq"}})
    merged = merge_reports([report])
    assert merged["verdict_summary"]["samples"]["not_evaluated"] == 1


def test_merged_output_has_no_fail_bucket(tmp_path):
    report = _write_report(tmp_path / "qc_report.json", {"s1": _sample("warn", ["frip"])})
    merged = merge_reports([report])
    assert "fail" not in merged["verdict_summary"]["samples"]


def test_merge_preserves_source_report_paths(tmp_path):
    report = _write_report(tmp_path / "qc_report.json", {"s1": _sample("pass")})
    merged = merge_reports([report])
    assert merged["reports"] == [str(report)]


def test_report_with_list_samples_shape_is_accepted(tmp_path):
    """Assay commands write {"samples": [ {...}, ... ]}; that must merge too."""
    path = tmp_path / "chip_qc.json"
    path.write_text(json.dumps({
        "assay": "ChIP-seq",
        "samples": [
            {"sample_id": "s1", "qc_verdict": {"overall": "warn", "warnings": ["frip"]}},
            {"sample_id": "s2", "qc_verdict": {"overall": "pass", "warnings": []}},
        ],
    }), encoding="utf-8")
    merged = merge_reports([path])
    assert merged["n_samples"] == 2
    assert merged["problems"] == []
    assert merged["verdict_summary"]["samples"] == {"pass": 1, "warn": 1, "not_evaluated": 0}


def test_mixed_shapes_merge_together(tmp_path):
    """A dict-shaped report and a list-shaped report must combine cleanly."""
    a = tmp_path / "a" / "qc_report.json"
    a.parent.mkdir(parents=True)
    a.write_text(json.dumps({"project_id": "proj_a",
                             "samples": {"s1": _sample("pass")}}), encoding="utf-8")
    b = tmp_path / "b" / "qc_report.json"
    b.parent.mkdir(parents=True)
    b.write_text(json.dumps({"project_id": "proj_b",
                             "samples": [{"sample_id": "s1",
                                          "qc_verdict": {"overall": "warn",
                                                         "warnings": ["tss_enrichment"]}}]}),
                 encoding="utf-8")
    merged = merge_reports([a, b])
    assert merged["n_reports"] == 2
    assert merged["n_samples"] == 2
    assert merged["problems"] == []
    assert merged["batches"] == {"proj_a": 1, "proj_b": 1}
