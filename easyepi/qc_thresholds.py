"""Pass / warn thresholds for QC metrics.

Turns the raw metric dictionary into an auditable verdict per metric so a JSON
report can be machine-read (CI gates, batch triage) instead of relying on a human
eyeballing the TSV.

IMPORTANT: every rule is **advisory**. A breached rule always yields ``warn``;
this module never returns ``fail``. QC thresholds encode conventions, not
biology, so a single metric must not be able to reject a library on its own — the
report surfaces the concern, and a human decides.

Rules are overridable at two levels (highest precedence last):
  1. ``DEFAULT_THRESHOLDS`` (this module),
  2. name/value overrides from a JSON file (``--qc-thresholds rules.json``),
  3. ``[qc_thresholds]`` tables in the project TOML.

Overrides are *merged per metric*: specifying ``frip`` only changes ``frip`` and
leaves every other default rule intact.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

# direction: "lt" means a LOW value is bad, "gt" means a HIGH value is bad.
# Metrics that must be evaluated before a sample may be called "pass". These
# are the read- and peak-level signals a reviewer expects on a ChIP/CUT&Tag QC
# sheet; if none of them could be computed the sample is unassessed, not healthy.
# Ratios whose value depends on the peak set spanning the aligned genome.
# Reporting them as breached for a partial peak set would be a false alarm.
_COVERAGE_SENSITIVE: frozenset[str] = frozenset({"frip"})

# Repressive and enhancer marks are genuinely diffuse: their peak-vs-flank
# contrast and FRiP sit far below what a promoter mark reaches. Judging them with
# one global cutoff produces a permanent warning that trains the reader to ignore
# the column. These per-marker overrides adjust only the contrast-sensitive rules;
# every other rule still applies.
MARKER_THRESHOLDS: dict[str, dict[str, dict[str, Any]]] = {
    "h3k27me3": {
        "signal_enrichment": {"op": "lt", "value": 1.2,
                              "message": "weak peak-vs-flank enrichment (<1.2x, repressive mark)"},
        "frip": {"op": "lt", "value": 0.01,
                 "message": "low FRiP (<1%, broad repressive domain)"},
    },
    "h3k9me3": {
        "signal_enrichment": {"op": "lt", "value": 1.2,
                              "message": "weak peak-vs-flank enrichment (<1.2x, heterochromatin)"},
        "frip": {"op": "lt", "value": 0.01,
                 "message": "low FRiP (<1%, heterochromatin domain)"},
    },
    "h3k36me3": {
        "signal_enrichment": {"op": "lt", "value": 1.3,
                              "message": "weak peak-vs-flank enrichment (<1.3x, gene-body mark)"},
    },
    "h3k4me1": {
        "signal_enrichment": {"op": "lt", "value": 1.5,
                              "message": "weak peak-vs-flank enrichment (<1.5x, enhancer mark)"},
    },
    "h3k27ac": {
        "signal_enrichment": {"op": "lt", "value": 1.8,
                              "message": "weak peak-vs-flank enrichment (<1.8x, enhancer mark)"},
    },
}


def _marker_key(marker: str | None) -> str | None:
    """Normalise a mark label such as 'H3K27me3' / 'h3k27me3' to a lookup key."""
    if not marker:
        return None
    return str(marker).strip().lower().replace("-", "").replace("_", "").replace(" ", "")


def resolve_rules(marker: str | None,
                  *sources: Any) -> dict[str, dict[str, Any]]:
    """Merge defaults, per-marker overrides, then explicit user sources.

    Precedence (later wins, per metric):
        DEFAULT_THRESHOLDS < MARKER_THRESHOLDS[marker] < sources in order
    """
    merged = {k: dict(v) for k, v in DEFAULT_THRESHOLDS.items()}
    key = _marker_key(marker)
    if key and key in MARKER_THRESHOLDS:
        for metric, rule in MARKER_THRESHOLDS[key].items():
            merged[metric] = dict(rule)
    # load_rules starts from the defaults, so only the metrics the user actually
    # provided should override the marker-specific values.
    overridden = _collect_overrides(*sources)
    for metric, rule in overridden.items():
        merged[metric] = rule
    return merged


def _collect_overrides(*sources: Any) -> dict[str, dict[str, Any]]:
    """Return only the metrics explicitly provided by the user sources."""
    found: dict[str, dict[str, Any]] = {}
    for source in sources:
        if source is None:
            continue
        payload: Any = source
        if isinstance(source, (str, Path)):
            path = Path(source)
            if not path.exists():
                continue
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except Exception:
                continue
        if not isinstance(payload, Mapping):
            continue
        candidate = payload.get("qc_thresholds", payload)
        if not isinstance(candidate, Mapping):
            continue
        for metric, rule in candidate.items():
            normalised = _normalise_rule(str(metric), rule)
            if normalised is not None:
                found[str(metric)] = normalised
    return found


CORE_METRICS: tuple[str, ...] = (
    "frip",
    "signal_enrichment",
    "duplicate_rate",
    "bam_mapping_rate",
)

DEFAULT_THRESHOLDS: dict[str, dict[str, Any]] = {
    # --- alignment / library ---
    "bam_mapping_rate": {"op": "lt", "value": 0.80,
                         "message": "low mapping rate (<80%)"},
    "duplicate_rate": {"op": "gt", "value": 0.50,
                       "message": "high duplication (>50%)"},
    "mitochondrial_fraction": {"op": "gt", "value": 0.20,
                               "message": "high mitochondrial fraction (>20%), possible lysis"},
    "mapq_lt_10": {"op": "gt", "value": 0.20,
                   "message": "many low-MAPQ reads (>20% below MAPQ 10)"},
    # --- signal quality ---
    "frip": {"op": "lt", "value": 0.20,
             "message": "low FRiP (<0.20)"},
    "signal_enrichment": {"op": "lt", "value": 2.0,
                          "message": "weak peak-vs-flank enrichment (<2x)"},
    "tss_enrichment": {"op": "lt", "value": 2.0,
                       "message": "low TSS enrichment (<2)"},
    # --- peak support distribution ---
    "fraction_peaks_zero": {"op": "gt", "value": 0.50,
                            "message": "over half of peaks have zero reads"},
    "top_decile_read_share": {"op": "gt", "value": 0.90,
                              "message": "reads concentrated in top 10% of peaks (>90%)"},
    # --- contamination ---
    "blacklist_fraction": {"op": "gt", "value": 0.01,
                           "message": "blacklist fraction above 1%"},
    # --- ATAC fragment structure ---
    "fraction_nfr": {"op": "lt", "value": 0.05,
                     "message": "few nucleosome-free fragments (<5%), poor tagmentation"},
    "mono_to_nfr_ratio": {"op": "gt", "value": 20.0,
                          "message": "mono/NFR ratio very high, under-digested"},
}

# Backwards-compatible alias used by earlier call sites.
THRESHOLDS = DEFAULT_THRESHOLDS

_VALID_OPS = {"lt", "gt"}


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and not (
        isinstance(value, float) and math.isnan(value))


def _normalise_rule(metric: str, rule: Any) -> dict[str, Any] | None:
    """Validate one rule mapping; return None when malformed."""
    if not isinstance(rule, Mapping):
        return None
    op = str(rule.get("op", "")).lower()
    if op not in _VALID_OPS:
        return None
    try:
        value = float(rule["value"])
    except (KeyError, TypeError, ValueError):
        return None
    return {"op": op, "value": value,
            "message": str(rule.get("message") or f"{metric} threshold")}


def load_rules(*sources: Any) -> dict[str, dict[str, Any]]:
    """Merge threshold overrides over :data:`DEFAULT_THRESHOLDS`.

    Each source may be:
      * ``None`` (ignored),
      * a path to a JSON object ``{"frip": {"op": "lt", "value": 0.3}}``,
      * an already-parsed mapping.

    Unknown or malformed entries are skipped rather than raising, so one typo in
    a config file cannot break a whole batch. The returned table always contains
    every default metric.
    """
    merged: dict[str, dict[str, Any]] = {k: dict(v) for k, v in DEFAULT_THRESHOLDS.items()}
    for source in sources:
        if source is None:
            continue
        payload: Any = source
        if isinstance(source, (str, Path)):
            path = Path(source)
            if not path.exists():
                continue
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
                # bad region or file cannot abort a batch QC run. The caller reports
                # the metric as not_evaluated rather than failing.
                continue
        if not isinstance(payload, Mapping):
            continue
        # Accept both a flat rule map and {"qc_thresholds": {...}}.
        candidate = payload.get("qc_thresholds", payload)
        if not isinstance(candidate, Mapping):
            continue
        for metric, rule in candidate.items():
            normalised = _normalise_rule(str(metric), rule)
            if normalised is not None:
                merged[str(metric)] = normalised
    return merged


# Metrics whose stored name differs from the rule name. ``bam_qc`` emits
# ``mapping_rate``; without this alias the ``bam_mapping_rate`` rule never fires
# and every sample reports it as not_evaluated.
_METRIC_ALIASES: dict[str, str] = {
    "mapping_rate": "bam_mapping_rate",
}


def _collect(metrics: dict[str, Any]) -> dict[str, float]:
    """Flatten nested QC payloads into a single scalar metric map."""
    flat: dict[str, float] = {}
    for key, value in (metrics or {}).items():
        if _is_number(value):
            flat.setdefault(key, float(value))
            alias = _METRIC_ALIASES.get(key)
            if alias:
                flat.setdefault(alias, float(value))
        elif isinstance(value, dict):
            for sub_key, sub_value in value.items():
                if not _is_number(sub_value):
                    continue
                flat.setdefault(sub_key, float(sub_value))
                alias = _METRIC_ALIASES.get(sub_key)
                if alias:
                    flat.setdefault(alias, float(sub_value))
    return flat


def evaluate_thresholds(metrics: dict[str, Any],
                        rules: dict[str, dict[str, Any]] | None = None) -> dict[str, Any]:
    """Evaluate QC metrics against the threshold table.

    Returns ``{"overall": "pass"|"warn"|"not_evaluated", "checks": [...],
    "warnings": [...]}``. Every breach is advisory, so ``overall`` is never
    ``"fail"``; ``failures`` is kept as an always-empty list for callers that
    still read it.
    """
    table = rules if rules is not None else DEFAULT_THRESHOLDS
    flat = _collect(metrics)
    checks: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []

    # A peak set that only covers part of the genome makes FRiP meaningless: the
    # denominator includes reads from regions the peaks could never overlap. The
    # caller signals this by passing ``partial_peak_coverage=True``.
    partial_peak_coverage = bool(metrics.get("_partial_peak_coverage"))

    for metric, rule in table.items():
        if metric in _COVERAGE_SENSITIVE and partial_peak_coverage:
            checks.append({"metric": metric, "status": "not_applicable",
                           "value": flat.get(metric), "threshold": rule["value"],
                           "op": rule["op"],
                           "message": "peak set does not span the aligned genome; "
                                      "this ratio is not comparable to genome-wide "
                                      "conventions"})
            continue
        if metric not in flat:
            checks.append({"metric": metric, "status": "not_evaluated",
                           "value": None, "threshold": rule["value"],
                           "op": rule["op"], "message": rule.get("message", "")})
            continue
        value = flat[metric]
        op = rule["op"]
        threshold = float(rule["value"])
        breached = value < threshold if op == "lt" else value > threshold
        entry = {"metric": metric, "status": "breached" if breached else "pass",
                 "value": value, "threshold": threshold, "op": op,
                 "severity": "warn",  # advisory only; there is no "fail" severity
                 "message": rule.get("message", "")}
        checks.append(entry)
        if breached:
            warnings.append(entry)

    # A verdict of "pass" must mean the assay was actually checked, not merely that
    # one incidental metric happened to be in range. Peak/read-level signals are the
    # core of ChIP QC: if none of them could be computed (e.g. a BAM with no peak
    # file), reporting "pass" misleads the reader into thinking the library passed.
    evaluated = {c["metric"] for c in checks
                 if c["status"] in {"pass", "breached"}}
    missing_core = [m for m in CORE_METRICS if m not in evaluated]
    coverage = (len(CORE_METRICS) - len(missing_core)) / len(CORE_METRICS)

    if warnings:
        overall = "warn"
    elif not missing_core:
        overall = "pass"
    elif coverage > 0:
        # Some core evidence exists but not all: do not claim a clean pass.
        overall = "incomplete"
    else:
        overall = "not_evaluated"
    return {"overall": overall, "checks": checks,
            "core_coverage": round(coverage, 3),
            "missing_core_metrics": missing_core,
            "warnings": [w["metric"] for w in warnings],
            "failures": [],  # deprecated: breaches are never fatal
            "warning_details": warnings}


def annotate_rows_with_thresholds(rows: Iterable[dict[str, Any]],
                                  rules: dict[str, dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    """Attach a ``qc_verdict`` block to each sample row (in place) and return them.

    A verdict already computed upstream is preserved: ``build_sample_payload``
    knows things this function cannot see (it marked FRiP not_applicable for a
    partial peak set), and recomputing here from ``assay_qc`` alone silently
    dropped that context. Callers that genuinely want a fresh verdict pass
    ``rules`` explicitly, which forces recomputation.
    """
    out: list[dict[str, Any]] = []
    for row in rows:
        if row.get("qc_verdict") is not None and rules is None:
            out.append(row)
            continue
        payload = dict(row.get("assay_qc") or {})
        if row.get("bam_qc"):
            payload.setdefault("bam_qc", row["bam_qc"])
        row["qc_verdict"] = evaluate_thresholds(payload, rules)
        out.append(row)
    return out


def summarize_verdicts(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate per-sample verdicts into a run-level summary."""
    counts = {"pass": 0, "warn": 0, "not_evaluated": 0}
    warned: dict[str, int] = {}
    for row in rows:
        verdict = row.get("qc_verdict") or {}
        overall = verdict.get("overall", "not_evaluated")
        counts[overall] = counts.get(overall, 0) + 1
        for metric in verdict.get("warnings", []):
            warned[metric] = warned.get(metric, 0) + 1
    return {"samples": counts, "warning_metrics": warned}
