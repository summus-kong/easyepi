"""`easyepi concordance`: do the samples that should agree actually agree?

Why this is a command and not a line in the QC report
-----------------------------------------------------
Per-sample QC answers "is this library good on its own". This answers a question
no single sample can: whether the samples *assigned to the same group* resemble
each other more than they resemble the other group. A cohort where every library
looks fine but the groups do not separate is measuring batch, not biology, and
nothing in a per-sample verdict can see that. It also needs a design file, which
a per-sample metric does not, so it stays a separate step that says what it used.

Two inputs, because they answer the question at two resolutions
--------------------------------------------------------------
* ``--matrix`` — a features x samples table (the same layout ``diff`` and
  ``normalize`` take). Gives correlation only.
* ``--peaks`` — one peak file per sample. Gives correlation *and* region-set
  Jaccard, because only peak files carry genomic coordinates.

What is deliberately not produced
---------------------------------
No ``qc_verdict`` block and no pass/fail threshold. ``concordance`` reports the
numbers and a within-group versus between-group summary; what counts as
concordant depends on the assay, the depth and the question, so it is a project
decision rather than a package default. The module docstring says the same thing
and this command keeps that promise instead of quietly inventing a cutoff.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from easyepi.concordance import (json_safe, pairwise_concordance,
                                 read_region_sets, replicate_concordance)


def _correlation_matrix(sample_ids, matrix, *, method: str, region_sets=None):
    """Square correlation matrix in the caller's sample order.

    ``pairwise_concordance`` is the single source of the numbers, so the matrix
    cannot drift from the pair table it is drawn beside. The diagonal is 1.0 by
    definition rather than by self-correlation, which would be NaN for any
    sample with too few finite features.
    """
    n = len(sample_ids)
    grid = np.full((n, n), np.nan)
    np.fill_diagonal(grid, 1.0)
    for pair in pairwise_concordance(sample_ids, matrix, method=method,
                                     region_sets=region_sets):
        i = sample_ids.index(pair["sample_a"])
        j = sample_ids.index(pair["sample_b"])
        grid[i, j] = grid[j, i] = pair["correlation"]
    return grid


def _read_design_groups(path: str) -> dict[str, str]:
    """Read ``sample_id`` + ``group``/``condition`` as a plain mapping.

    Deliberately looser than ``differential.read_design``: that one requires
    exactly two levels because a two-group test needs it, while concordance is
    meaningful with any number of groups. Sample-to-matrix matching is still
    exact and still refuses a mismatch by name.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"design file {path} has no rows")
    if "sample_id" not in rows[0]:
        raise ValueError(
            f"design file {path} must have a sample_id column; found "
            f"{sorted(rows[0])}")
    group_key = "group" if "group" in rows[0] else (
        "condition" if "condition" in rows[0] else None)
    if group_key is None:
        raise ValueError(
            f"design file {path} must have a group or condition column; found "
            f"{sorted(rows[0])}")
    groups: dict[str, str] = {}
    for row in rows:
        sample = (row.get("sample_id") or "").strip()
        group = (row.get(group_key) or "").strip()
        if not sample or not group:
            raise ValueError(
                f"design file {path}: every row needs a nonempty sample_id and "
                f"{group_key}")
        if sample in groups:
            raise ValueError(f"design file {path}: duplicate sample_id {sample!r}")
        groups[sample] = group
    return groups


def _write_pairs_tsv(path: Path, pairs) -> None:
    """Long-form pair table, for pandas/R/awk.

    The same records as the JSON, in the layout a spreadsheet or a one-liner can
    read. A value that could not be computed is written as an empty field, not as
    0 and not as "nan": an empty field reads as missing in every table tool,
    whereas "nan" is a string that silently compares as zero in some of them.
    """

    def fmt(value):
        if value is None:
            return ""
        if isinstance(value, float) and not np.isfinite(value):
            return ""
        return value

    columns = ["sample_a", "sample_b", "comparison", "group_a", "group_b",
               "correlation", "jaccard", "method"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(columns)
        for pair in pairs:
            writer.writerow([fmt(pair.get(c)) for c in columns])


def _console_report(result: dict[str, Any], out: Path) -> None:
    """Print the summary, and the two things a reader would otherwise miss."""
    summary = result["summary"]
    print(f"[easyepi] concordance: {summary['pairs']} pair(s) from "
          f"{result['n_samples']} sample(s), level={result['level']}")

    def fmt(value):
        return "n/a" if value is None else f"{value:.4f}"

    # Each mean is formatted defensively: with one sample per group the
    # between-group mean is legitimately absent while the within-group one
    # exists, and a bare ``{value:.4f}`` on that None raises instead of
    # reporting the absence.
    within = summary.get("within_group_mean_correlation")
    between = summary.get("between_group_mean_correlation")
    if within is not None or between is not None:
        print(f"[easyepi] mean correlation: within-group {fmt(within)}, "
              f"between-group {fmt(between)}")
    wj = summary.get("within_group_mean_jaccard")
    bj = summary.get("between_group_mean_jaccard")
    if wj is not None or bj is not None:
        print(f"[easyepi] mean Jaccard: within-group {fmt(wj)}, "
              f"between-group {fmt(bj)}")

    supported = summary.get("grouping_supported")
    if supported is not None:
        print(f"[easyepi] grouping_supported: {supported}")
        # This caveat is printed on every run, not only odd ones. The value is a
        # comparison of means, so one mislabelled sample among six still reads
        # True -- measured on synthetic designs with real group structure. A
        # bare "True" would therefore be read as "no sample is swapped", which
        # the number does not support.
        print("[easyepi]   note: a group-level comparison of means; it does not "
              "identify a mislabelled sample")

    for note in result.get("notes") or ():
        print(f"[easyepi] note: {note}")
    print(f"[easyepi] written: {out / 'concordance.json'}, "
          f"{out / 'concordance_pairs.tsv'}")


def run_concordance(args) -> int:
    """Entry point for ``easyepi concordance``."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    design_path = getattr(args, "design", None)
    if not design_path:
        raise ValueError(
            "--design is required: sample groups are never inferred from column "
            "order or from the filename")
    groups = _read_design_groups(design_path)

    peaks = getattr(args, "peaks", None)
    matrix_path = getattr(args, "matrix", None)
    if peaks and matrix_path:
        raise ValueError("--peaks and --matrix are mutually exclusive; pass one")
    if not peaks and not matrix_path:
        raise ValueError(
            "no input: pass --peaks (peak files, gives correlation and Jaccard) "
            "or --matrix (a features x samples table, gives correlation only)")

    method = getattr(args, "method", "spearman")
    region_sets = None

    if peaks:
        region_sets = read_region_sets(peaks)
        sample_ids = list(region_sets)
        level = "peaks"
        # Correlation across regions: presence in a sample's peak set becomes a
        # 0/1 signal, which is exactly what "do these peak sets agree" means when
        # there is no count matrix to correlate.
        universe = set().union(*region_sets.values()) if region_sets else set()
        ordered = sorted(universe)
        if not ordered:
            raise ValueError(
                "every --peaks file is empty; there is nothing to compare")
        matrix = np.array([[1.0 if region in region_sets[s] else 0.0
                            for s in sample_ids] for region in ordered],
                          dtype=float)
    else:
        from easyepi.differential import read_count_matrix

        _regions, matrix, sample_ids = read_count_matrix(matrix_path)
        if matrix.size == 0:
            raise ValueError(f"{matrix_path} has no data rows")
        level = "matrix"

    # Refuse before doing arithmetic rather than after producing a table of
    # NaN: a "concordance" of one sample is not a weak result, it is not a
    # result. Two is the floor because a pair needs two samples; --min-samples
    # can raise it, never lower it.
    floor = max(2, int(getattr(args, "min_samples", 2) or 2))
    if len(sample_ids) < floor:
        raise ValueError(
            f"concordance needs at least {floor} samples, "
            f"got {len(sample_ids)}: {sample_ids}")

    missing = [s for s in sample_ids if s not in groups]
    extra = [s for s in groups if s not in sample_ids]
    if missing or extra:
        # Same wording as differential.read_design, so a user who has hit one
        # message recognises the other.
        raise ValueError(
            f"design/input sample mismatch; missing={missing}, extra={extra}")

    design_order = list(groups)
    ordered_ids = [s for s in design_order if s in sample_ids]
    if ordered_ids != sample_ids:
        # Reorder to the design so the heatmap reads group by group, and remap
        # the matrix columns to match.
        index = {s: i for i, s in enumerate(sample_ids)}
        matrix = matrix[:, [index[s] for s in ordered_ids]]
        sample_ids = ordered_ids

    result = replicate_concordance(sample_ids, matrix, groups=groups,
                                   region_sets=region_sets, method=method)

    corr = _correlation_matrix(sample_ids, matrix, method=method,
                               region_sets=region_sets)

    payload: dict[str, Any] = {
        "level": level,
        "method": method,
        "n_samples": len(sample_ids),
        "samples": sample_ids,
        "groups": {s: groups[s] for s in sample_ids},
        "design_file": str(design_path),
        "summary": result["summary"],
        "notes": [],
        "correlation_matrix": corr.tolist(),
        "pairs": result["pairs"],
    }

    note = result["summary"].get("note")
    if note:
        payload["notes"].append(note)
    if level == "matrix":
        payload["notes"].append(
            "Jaccard was not computed: a count matrix has no genomic "
            "coordinates. Pass --peaks for region-set overlap.")
    if payload["summary"].get("grouping_supported") is not None:
        payload["notes"].append(
            "grouping_supported compares mean within-group against mean "
            "between-group agreement; it is a group-level statement and does "
            "not identify which sample is mislabelled.")

    # NaN -> null before serialising: json.dumps would otherwise emit a bare
    # `NaN`, which is not JSON and makes the whole document unreadable to a
    # strict parser.
    payload = json_safe(payload)
    (out / "concordance.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_pairs_tsv(out / "concordance_pairs.tsv", result["pairs"])

    _console_report(payload, out)
    return 0
