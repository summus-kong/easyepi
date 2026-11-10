"""Merge multiple ``qc_report.json`` files into one batch-level verdict summary.

``make-qc-report`` produces one report per project. Real studies routinely span
several projects / sequencing batches, and the question "is the whole study
healthy?" then needs an answer that no single report can give. This module
collects any number of reports (files, directories or shell globs) and rolls their
per-sample verdicts up with :func:`easyepi.qc_thresholds.summarize_verdicts`.

Design notes
------------
* A sample id colliding across batches is namespaced as ``batch::sample`` so two
  projects that both contain ``sample_1`` stay distinguishable.
* Unreadable or malformed reports are collected in ``problems`` instead of
  aborting the merge: a partial answer is more useful than none.
* Only verdicts are aggregated; raw metrics stay in the source reports.
"""

from __future__ import annotations

import glob as _glob
import json
from pathlib import Path
from typing import Any, Iterable

from .qc_thresholds import summarize_verdicts

DEFAULT_REPORT_NAME = "qc_report.json"


def _expand_one(pattern: str) -> list[Path]:
    """Expand a single file / directory / glob argument into report paths."""
    text = str(pattern)
    if any(ch in text for ch in "*?["):
        return sorted(Path(x) for x in _glob.glob(text, recursive=True) if Path(x).is_file())
    path = Path(text)
    if path.is_dir():
        # Prefer the canonical name at any depth; fall back to any *qc_report.json.
        found = sorted(path.rglob(DEFAULT_REPORT_NAME))
        if found:
            return found
        return sorted(p for p in path.rglob("*qc_report.json") if p.is_file())
    return [path] if path.is_file() else []


def collect_reports(patterns: Iterable[str]) -> tuple[list[Path], list[str]]:
    """Resolve patterns to unique report paths; return ``(found, unresolved)``."""
    found: list[Path] = []
    unresolved: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        paths = _expand_one(pattern)
        if not paths:
            unresolved.append(pattern)
            continue
        for path in paths:
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            found.append(path)
    return found, unresolved


def merge_reports(paths: Iterable[Path]) -> dict[str, Any]:
    """Merge report files into one batch-level structure.

    Returns a dict with ``reports``, ``batches`` (batch -> sample count),
    ``samples`` (namespaced sample -> verdict payload), ``verdict_summary`` and
    ``problems`` (unreadable / malformed inputs).
    """
    rows: list[dict[str, Any]] = []
    batches: dict[str, int] = {}
    problems: list[str] = []
    used_ids: set[str] = set()

    for path in paths:
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except Exception as exc:
            problems.append(f"{path}: {type(exc).__name__}: {exc}")
            continue
        if not isinstance(data, dict):
            problems.append(f"{path}: top level is not an object")
            continue
        samples = data.get("samples")
        # Two shapes reach us: make-qc-report writes {"samples": {id: {...}}} while
        # the assay commands write {"samples": [{sample_id: ..., ...}, ...]}. Accept
        # both so any report in the tree can be merged.
        if isinstance(samples, list):
            as_dict: dict[str, Any] = {}
            for entry in samples:
                if not isinstance(entry, dict):
                    continue
                sid = entry.get("sample_id")
                if sid is None:
                    continue
                as_dict[str(sid)] = entry
            samples = as_dict
        if not isinstance(samples, dict):
            problems.append(f"{path}: missing or unusable 'samples'")
            continue
        batch = str(data.get("project_id")
                    or data.get("module")
                    or Path(path).parent.name)
        batches[batch] = batches.get(batch, 0) + len(samples)
        for sample_id, payload in samples.items():
            payload = payload if isinstance(payload, dict) else {}
            verdict = payload.get("qc_verdict")
            if not isinstance(verdict, dict):
                verdict = {}
            sample_key = str(sample_id)
            if sample_key in used_ids:
                sample_key = f"{batch}::{sample_key}"
            used_ids.add(sample_key)
            rows.append({"batch": batch, "sample_id": sample_key,
                         "assay": payload.get("assay"), "qc_verdict": verdict})

    return {
        "module": "merged-qc-report",
        "reports": [str(p) for p in paths],
        "batches": batches,
        "n_reports": len(batches),
        "n_samples": len(rows),
        "samples": {r["sample_id"]: {"batch": r["batch"], "assay": r["assay"],
                                     "qc_verdict": r["qc_verdict"]} for r in rows},
        "verdict_summary": summarize_verdicts(rows),
        "problems": problems,
    }
