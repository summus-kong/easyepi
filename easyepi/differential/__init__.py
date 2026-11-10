"""Differential analysis and statistical modeling (native).

Built directly on ``easyepi.core.stats``, so the whole statistical layer is
dependency-free beyond numpy/scipy. Supports the descriptive summaries and
frequentist tests requested (t-test / Wilcoxon / Mann-Whitney), BH-FDR
correction, effect sizes, and a peak-overlap enrichment test.

Two-layer design:
  * Semantic functions (``differential_features`` / ``differential_test``)
    take plain arrays / matrix and return structured dicts / results. These are
    the functions a Python user calls via ``import easyepi``.
  * ``run_diff`` parses CLI args and wraps the semantic functions to write a
    TSV report.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from easyepi.core import stats as S


def differential_test(a: np.ndarray, b: np.ndarray, method: str = "mannwhitney",
                      equal_var: bool = False) -> dict:
    """Run the requested two-group test on one feature."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("test inputs must be one-dimensional")
    if method == "wilcoxon":
        if a.size != b.size:
            raise ValueError("paired tests require equal sample counts")
        keep = np.isfinite(a) & np.isfinite(b)
        a, b = a[keep], b[keep]
    else:
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        raise ValueError("at least two finite biological replicates per group are required")
    if (a < 0).any() or (b < 0).any():
        raise ValueError("log2 fold change requires nonnegative count/signal values")
    if method == "ttest":
        res = S.t_test(a, b, equal_var=equal_var)
    elif method == "wilcoxon":
        res = S.wilcoxon(a, b)
    elif method == "mannwhitney":
        res = S.mann_whitney_u(a, b)
    else:
        raise ValueError(f"unknown test: {method}")
    res["log2fc"] = S.log2fc(a, b)
    res["cohens_d"] = S.cohens_d(a, b)
    res["mean_a"] = float(a.mean()) if a.size else np.nan
    res["mean_b"] = float(b.mean()) if b.size else np.nan
    return res


def differential_features(matrix: np.ndarray, group_a_mask: np.ndarray,
                          method: str = "mannwhitney", equal_var: bool = False) -> dict:
    """Run a per-row differential test between two sample groups.

    Args:
        matrix: 2-D count/signal matrix (rows = regions, cols = samples).
        group_a_mask: boolean mask (length cols) selecting group A samples.
        method: one of 'mannwhitney' (default, independent groups),
            'wilcoxon' (paired observations), or 'ttest'.
        equal_var: for 'ttest' only, whether to assume equal variance.

    Returns:
        dict with per-region arrays: mean_a, mean_b, log2fc, cohens_d,
        pvalue, padj (BH-FDR corrected) plus the test name.
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError("matrix must be 2-D (regions x samples)")
    mask_a = np.asarray(group_a_mask, dtype=bool)
    if mask_a.ndim != 1 or mask_a.shape[0] != matrix.shape[1]:
        raise ValueError("group_a_mask must match number of samples")
    mask_b = ~mask_a
    if mask_a.sum() == 0 or mask_b.sum() == 0:
        raise ValueError("both groups must contain at least one sample")

    if method not in {"ttest", "wilcoxon", "mannwhitney"}:
        raise ValueError(f"unknown test: {method}")
    if min(mask_a.sum(), mask_b.sum()) < 2:
        raise ValueError("at least two biological replicates per group are required")
    if method == "wilcoxon" and mask_a.sum() != mask_b.sum():
        raise ValueError("paired Wilcoxon requires equal group sizes")
    n_regions = matrix.shape[0]
    mean_a = np.full(n_regions, np.nan)
    mean_b = np.full(n_regions, np.nan)
    log2fc = np.full(n_regions, np.nan)
    cohens_d = np.full(n_regions, np.nan)
    pvalue = np.full(n_regions, np.nan)

    for i in range(n_regions):
        av = matrix[i, mask_a]
        bv = matrix[i, mask_b]
        if method == "wilcoxon":
            keep = np.isfinite(av) & np.isfinite(bv)
            av, bv = av[keep], bv[keep]
        else:
            av, bv = av[np.isfinite(av)], bv[np.isfinite(bv)]
        if av.size < 2 or bv.size < 2:
            continue
        res = differential_test(av, bv, method=method, equal_var=equal_var)
        mean_a[i] = res["mean_a"]
        mean_b[i] = res["mean_b"]
        log2fc[i] = res["log2fc"]
        cohens_d[i] = res["cohens_d"]
        pvalue[i] = res["pvalue"]

    padj = S.benjamini_hochberg(pvalue)
    return {
        "mean_a": mean_a, "mean_b": mean_b, "log2fc": log2fc,
        "cohens_d": cohens_d, "pvalue": pvalue, "padj": padj,
        "test": method, "n_regions": int(n_regions),
    }


def read_count_matrix(path) -> tuple:
    """Read a TSV count matrix into (region_ids, matrix, sample_ids)."""
    import csv
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if row:
                rows.append(row)
    if not rows:
        return [], np.empty((0, 0)), []
    if header is None or len(header) < 2 or len(set(header[1:])) != len(header[1:]) or any(not v for v in header):
        raise ValueError("Matrix header requires unique, nonempty sample IDs")
    if any(len(row) != len(header) for row in rows):
        raise ValueError("Matrix rows must match the header width")
    region_ids = [r[0] for r in rows]
    if len(region_ids) != len(set(region_ids)) or any(not v for v in region_ids):
        raise ValueError("Region IDs must be nonempty and unique")
    mat = np.array([[float(x) for x in r[1:]] for r in rows], dtype=np.float64)
    sample_ids = list(header[1:]) if header else [f"s{i}" for i in range(mat.shape[1])]
    return region_ids, mat, sample_ids


def read_design(path, sample_ids):
    """Read a design TSV and return group mask plus validated metadata.

    Required columns are ``sample_id`` and either ``group`` or ``condition``.
    The order of the design rows is deliberately irrelevant; sample IDs are
    matched to the count-matrix header to prevent accidental column splitting.
    """
    import csv
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or "sample_id" not in rows[0]:
        raise ValueError("design TSV must contain a sample_id column")
    group_key = "group" if "group" in rows[0] else ("condition" if "condition" in rows[0] else None)
    if group_key is None:
        raise ValueError("design TSV must contain group or condition")
    if len(sample_ids) != len(set(sample_ids)) or len(rows) != len({r["sample_id"] for r in rows}):
        raise ValueError("duplicate sample IDs in matrix or design")
    if any(not r["sample_id"] or not r.get(group_key) for r in rows):
        raise ValueError("design sample IDs and groups must be non-empty")
    by_id = {r["sample_id"]: r for r in rows}
    missing = [s for s in sample_ids if s not in by_id]
    extra = [s for s in by_id if s not in sample_ids]
    if missing or extra:
        raise ValueError(f"design/count sample mismatch; missing={missing}, extra={extra}")
    groups = [by_id[s][group_key] for s in sample_ids]
    unique = list(dict.fromkeys(groups))
    if len(unique) != 2:
        raise ValueError(f"two-group differential requires exactly 2 levels, got {unique}")
    return np.array([g == unique[0] for g in groups], dtype=bool), by_id, unique


def run_diff(args) -> int:
    """Design-aware differential analysis across all features."""
    import csv

    counts_path = Path(args.counts)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    region_ids, mat, sample_ids = read_count_matrix(counts_path)
    if mat.size == 0:
        out.write_text("# easyepi differential - no data\n", encoding="utf-8")
        print("[easyepi] differential: empty input")
        return 0
    if mat.shape[1] < 2:
        raise ValueError("need at least two sample columns for differential")

    design_path = getattr(args, "design", None)
    if not design_path:
        raise ValueError("--design is required; groups are never inferred from column order")
    mask_a, design, levels = read_design(design_path, sample_ids)
    method = getattr(args, "method", "mannwhitney")
    if method == "wilcoxon":
        pair_key = getattr(args, "pair", None)
        if not pair_key:
            raise ValueError("paired Wilcoxon requires --pair column in the design")
        a_ids = [sid for sid, selected in zip(sample_ids, mask_a) if selected]
        b_ids = [sid for sid, selected in zip(sample_ids, mask_a) if not selected]
        a_pairs = [design[sid].get(pair_key, "") for sid in a_ids]
        b_pairs = [design[sid].get(pair_key, "") for sid in b_ids]
        if (not all(a_pairs + b_pairs) or len(set(a_pairs)) != len(a_pairs)
                or len(set(b_pairs)) != len(b_pairs) or set(a_pairs) != set(b_pairs)):
            raise ValueError("paired design requires one sample per group for every pair")
        b_by_pair = dict(zip(b_pairs, b_ids))
        ordered = a_ids + [b_by_pair[pair] for pair in a_pairs]
        mat = mat[:, [sample_ids.index(sid) for sid in ordered]]
        mask_a = np.array([True] * len(a_ids) + [False] * len(a_ids))
    res = differential_features(mat, mask_a, method=method,
                                equal_var=bool(getattr(args, "equal_var", False)))
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["region", "group_a", "group_b", "mean_a", "mean_b", "log2FC", "cohens_d", "pvalue", "padj", "test"])
        for i, rid in enumerate(region_ids):
            w.writerow([rid, levels[0], levels[1], res["mean_a"][i], res["mean_b"][i],
                        res["log2fc"][i], res["cohens_d"][i], res["pvalue"][i],
                        res["padj"][i], method])
    print(f"[easyepi] differential results written to {out}")
    return 0
