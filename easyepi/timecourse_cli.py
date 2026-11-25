"""CLI entry point for `easyepi timecourse`."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from easyepi import timecourse as tc
from easyepi.differential import read_count_matrix


def run_timecourse(args) -> int:
    """Time-course specificity (TAU) and monotonic trend per feature."""
    counts_path = Path(args.counts)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    region_ids, mat, sample_ids = read_count_matrix(counts_path)
    if mat.size == 0:
        out.write_text("# easyepi timecourse - no data\n", encoding="utf-8")
        print("[easyepi] timecourse: empty input, nothing to test")
        return 0

    design_path = getattr(args, "design", None)
    if not design_path:
        raise ValueError("--design is required; time points are never inferred")
    raw_time, levels, time_values = tc.read_time_design(design_path, sample_ids)

    rows = tc.timecourse_summary(
        mat, region_ids, sample_ids, raw_time, levels, time_values)
    ordered, padj = tc.write_timecourse(rows, out, levels, time_values)

    measured = [r for r in ordered if r["tau"] is not None]
    undefined = len(ordered) - len(measured)
    top = measured[0] if measured else None

    summary = {
        "features": len(ordered),
        "time_points": len(levels),
        "replicates_per_time": _replicate_counts(sample_ids, raw_time, levels),
        "tau_defined": len(measured),
        "tau_undefined_all_zero": undefined,
        "tau_median": (float(np.median([r["tau"] for r in measured]))
                       if measured else None),
        "trend_test": ("spearman" if time_values is not None else None),
        "trend_note": (
            None if time_values is not None else
            "time labels are not numeric, so no monotonic trend was tested"),
        "most_specific": None if top is None else {
            "feature": top["feature"],
            "tau": top["tau"],
            "peak_time": top["peak_time"],
        },
        "tau_is_descriptive": (
            "tau measures how concentrated a profile is on one time point; it is "
            "not a test statistic and carries no p-value"),
    }
    (out.parent / (out.name + ".json")).write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print("[easyepi] timecourse: %d features over %d time points -> %s"
          % (len(ordered), len(levels), out))
    if top is not None:
        print("           most specific: %s tau=%.3f at %s"
              % (top["feature"], top["tau"], top["peak_time"]))
    if undefined:
        print("           %d feature(s) all-zero at every time point; tau undefined"
              % undefined)
    return 0


def _replicate_counts(sample_ids, raw_time, levels):
    counts = {lvl: 0 for lvl in levels}
    for sid in sample_ids:
        counts[raw_time[sample_ids.index(sid)]] += 1
    return counts
