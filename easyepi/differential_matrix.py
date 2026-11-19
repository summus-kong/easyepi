"""Descriptive differential analysis of integrated BigWig signal.

Builds a shared-peak matrix from BigWig signal across samples, then runs
a per-region independent-group differential test (Mann-Whitney by default)
and returns log2FC, p-values, BH-adjusted q-values, and a volcano figure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

# Last non-strict region read error, for diagnostics.
_LAST_REGION_ERROR: list[str | None] = [None]

from easyepi.core.stats import benjamini_hochberg


def read_narrowpeak_union(peak_files: Sequence[Optional[str]],
                          bin_size: int = 200) -> List[tuple]:
    """Merge peaks from several narrowPeak files into a non-overlapping union of
    regions (chrom, start, end). Regions are sorted and merged when overlapping."""
    intervals: List[tuple] = []
    for pf in peak_files:
        if not pf:
            continue
        if not Path(pf).is_file():
            raise FileNotFoundError(pf)
        with open(pf) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                p = line.split("\t")
                if len(p) >= 3:
                    try:
                        start, end = int(p[1]), int(p[2])
                        if start < 0 or end <= start:
                            raise ValueError("invalid interval")
                        intervals.append((p[0], start, end))
                    except ValueError as exc:
                        raise ValueError(f"Invalid BED row in {pf}: {line.strip()}") from exc
                else:
                    raise ValueError(f"Invalid BED row in {pf}: {line.strip()}")
    if not intervals:
        return []
    # sort and merge
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    merged: List[tuple] = []
    for chrom, s, e in intervals:
        if not merged or chrom != merged[-1][0]:
            merged.append((chrom, s, e))
            continue
        last_c, last_s, last_e = merged[-1]
        if chrom == last_c and s <= last_e:
            merged[-1] = (chrom, last_s, max(last_e, e))
        else:
            merged.append((chrom, s, e))
    return merged


def bigwig_region_sum(bw_handle, chrom: str, start: int, end: int, strict: bool = True) -> float:
    """Sum the BigWig signal over a region using an already-open handle."""
    try:
        vals = bw_handle.values(chrom, int(start), int(end))
        if vals is None:
            return 0.0
        arr = np.array([0.0 if v is None or np.isnan(v) else float(v) for v in vals])
        return float(arr.sum())
    except Exception as exc:
        if strict:
            raise
        # Non-strict mode: an unreadable region contributes zero signal, but the
        # reason is recorded so the caller can surface it instead of guessing.
        _LAST_REGION_ERROR[0] = f"{type(exc).__name__}: {exc}"
        return 0.0


def build_count_matrix(bigwigs: Sequence[str], peak_files: Sequence[Optional[str]],
                       sample_names: Optional[Sequence[str]] = None,
                       max_regions: int | None = None, seed: int = 0,
                       strict: bool = True) -> dict:
    """Build a (n_regions, n_samples) integrated-signal matrix from union peak regions.

    Reuses a single pyBigWig handle per sample (the naive per-region open was a
    performance bottleneck with 100k+ regions). If the union has more than
    an explicitly supplied ``max_regions`` limit, a seeded subset is used.
    All regions are retained by default. These values are not raw read counts.
    """
    try:
        import pyBigWig
    except ImportError:
        raise ImportError("pyBigWig is required to build a count matrix (pip install pyBigWig)")
    regions = read_narrowpeak_union(peak_files)
    if not regions:
        raise ValueError("No valid peak regions were supplied")
    if max_regions is not None and max_regions < 1:
        raise ValueError("max_regions must be positive or None")
    if max_regions is not None and len(regions) > max_regions:
        rng = np.random.RandomState(seed)
        idx = np.sort(rng.choice(len(regions), size=max_regions, replace=False))
        regions = [regions[i] for i in idx]
    names = list(sample_names) if sample_names else [Path(b).stem for b in bigwigs]
    if len(names) != len(bigwigs) or len(set(names)) != len(names):
        raise ValueError("sample names must be unique and match BigWig files")
    n_regions = len(regions)
    matrix = np.zeros((n_regions, len(bigwigs)), dtype=np.float64)
    for j, bw in enumerate(bigwigs):
        try:
            handle = pyBigWig.open(bw)
        except Exception as exc:
            if strict:
                raise RuntimeError(f"cannot open BigWig: {bw}") from exc
            continue
        try:
            for i, (chrom, s, e) in enumerate(regions):
                matrix[i, j] = bigwig_region_sum(handle, chrom, s, e, strict=strict)
        finally:
            handle.close()
    return {"regions": regions, "matrix": matrix, "sample_names": names}


def differential_matrix(matrix: np.ndarray, group_a: np.ndarray,
                         method: str = "mannwhitney") -> dict:
    """Per-region differential test between two groups.

    Returns log2FC, pvalues, padj (BH), and a significance flag per region.
    """
    from easyepi.differential import differential_features
    result = differential_features(matrix, group_a, method=method)
    result["significant"] = result["padj"] < 0.05
    return result


def volcano_figure(log2fc: np.ndarray, pvalue: np.ndarray,
                   title: str = "Differential volcano") -> object:
    """Render a real volcano figure (matplotlib Figure) from log2FC & p-values."""
    import pandas as pd
    from easyepi.plot import plot_volcano
    return plot_volcano(pd.DataFrame({"log2fc": log2fc, "pvalue": pvalue,
                                      "padj": benjamini_hochberg(pvalue)}), title=title)
