"""Strict zero-based methylation tables, coverage summaries and exploratory two-group site comparisons. This module does not call DMRs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .core.stats import benjamini_hochberg


class MethylationError(ValueError):
    """Invalid methylation table input."""


def read_methylation_table(path: str | Path) -> pd.DataFrame:
    """Read a position/beta table; raises when required columns are absent."""
    frame = pd.read_csv(path, sep="\t", comment="#")
    lower = {str(c).lower(): c for c in frame.columns}
    def pick(*names: str):
        """Choose the first present column from a list of candidate names."""
        for name in names:
            if name in lower:
                return lower[name]
        return None
    chrom, start, end = pick("chrom", "chr", "chromosome"), pick("start", "pos", "position"), pick("end")
    if not chrom or not start:
        raise MethylationError("methylation table requires chrom/chr and start/pos/position columns")
    if frame.empty or frame[chrom].isna().any() or frame[chrom].astype(str).str.strip().eq("").any():
        raise MethylationError("methylation rows require a non-empty chromosome")
    out = pd.DataFrame({"chrom": frame[chrom].astype(str), "start": pd.to_numeric(frame[start], errors="raise")})
    out["end"] = pd.to_numeric(frame[end], errors="coerce") if end else out["start"] + 1
    beta = pick("beta", "methylation", "methylation_fraction", "fraction")
    meth = pick("methylated", "methylated_count", "m")
    unmeth = pick("unmethylated", "unmethylated_count", "u")
    cov = pick("coverage", "depth", "total_count", "n")
    if meth and unmeth:
        m = pd.to_numeric(frame[meth], errors="raise")
        u = pd.to_numeric(frame[unmeth], errors="raise")
        if not np.isfinite(m).all() or not np.isfinite(u).all() or (m < 0).any() or (u < 0).any():
            raise MethylationError("methylated and unmethylated counts must be finite and nonnegative")
        if ((m % 1 != 0) | (u % 1 != 0)).any() or ((m + u) <= 0).any():
            raise MethylationError("counts must be integers with positive total coverage")
        out["methylated"], out["unmethylated"] = m, u
    if beta:
        out["beta"] = pd.to_numeric(frame[beta], errors="coerce")
    elif meth and unmeth:
        m = pd.to_numeric(frame[meth], errors="coerce")
        u = pd.to_numeric(frame[unmeth], errors="coerce")
        out["beta"] = m / (m + u).replace(0, np.nan)
        out["methylated"] = m
        out["unmethylated"] = u
    else:
        raise MethylationError("methylation table requires beta or methylated+unmethylated columns")
    if cov:
        out["coverage"] = pd.to_numeric(frame[cov], errors="coerce")
    else:
        out["coverage"] = m + u if meth and unmeth else np.nan
    numeric = out[["start", "end", "beta"]].to_numpy(float)
    if not np.isfinite(numeric).all():
        raise MethylationError("coordinates and beta must be finite numbers")
    if (out["start"] < 0).any() or (out["end"] <= out["start"]).any() or ((out["start"] % 1 != 0) | (out["end"] % 1 != 0)).any():
        raise MethylationError("coordinates must be zero-based half-open integer intervals")
    if ((out["beta"] < 0) | (out["beta"] > 1)).any():
        raise MethylationError("beta must lie in [0, 1]; convert percentages explicitly")
    if cov and (not np.isfinite(out["coverage"]).all() or (out["coverage"] < 0).any()):
        raise MethylationError("coverage must be finite and nonnegative")
    if cov and meth and unmeth and not np.allclose(out["coverage"], m + u):
        raise MethylationError("coverage must equal methylated + unmethylated counts")
    if beta and meth and unmeth and not np.allclose(out["beta"], m / (m + u), atol=1e-6):
        raise MethylationError("beta is inconsistent with methylated/unmethylated counts")
    if out.duplicated(["chrom", "start", "end"]).any():
        raise MethylationError("duplicate methylation intervals are not supported")
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    if out.empty:
        raise MethylationError(f"no valid methylation rows: {path}")
    return out


def summarize_methylation(frame: pd.DataFrame, min_coverage: float = 0) -> dict[str, float | int]:
    """Coverage and beta-value summary, optionally filtered by minimum coverage."""
    keep = _filter_coverage(frame, min_coverage)
    beta = keep["beta"].to_numpy(float)
    return {"n_rows": int(len(frame)), "n_retained": int(len(keep)), "n_chromosomes": int(keep["chrom"].nunique()),
            "mean_beta": float(np.mean(beta)) if len(beta) else float("nan"),
            "median_beta": float(np.median(beta)) if len(beta) else float("nan"),
            "mean_coverage": float(keep["coverage"].mean()) if len(keep) else float("nan"),
            "min_coverage": float(min_coverage)}


def profile_methylation(frame: pd.DataFrame, out: str | Path, min_coverage: float = 0) -> dict[str, str | int]:
    """Methylation profile across the supplied regions."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    keep = _filter_coverage(frame, min_coverage).copy()
    keep.to_csv(out, sep="\t", index=False)
    return {"output": str(out), "n_rows": int(len(keep))}


def differential_methylation(frames: Sequence[pd.DataFrame], sample_ids: Sequence[str], groups: Sequence[str], min_coverage: float = 5) -> pd.DataFrame:
    """Exploratory two-group methylation site comparison; not DMR calling."""
    if len(frames) != len(groups) or len(frames) != len(sample_ids):
        raise MethylationError("frames, sample_ids and groups must have equal length")
    levels = list(dict.fromkeys(groups))
    if len(levels) != 2:
        raise MethylationError("differential methylation requires exactly two groups")
    if len(set(sample_ids)) != len(sample_ids) or any(not s for s in sample_ids):
        raise MethylationError("sample IDs must be non-empty and unique")
    if any(not g for g in groups) or any(groups.count(g) < 2 for g in levels):
        raise MethylationError("at least two biological replicates per group are required")
    keyed = []
    for frame, sid in zip(frames, sample_ids):
        x = _filter_coverage(frame, min_coverage)[["chrom", "start", "end", "beta"]].copy()
        x = x.rename(columns={"beta": sid}).set_index(["chrom", "start", "end"])
        keyed.append(x)
    matrix = pd.concat(keyed, axis=1, join="outer")
    a = [sid for sid, g in zip(sample_ids, groups) if g == levels[0]]
    b = [sid for sid, g in zip(sample_ids, groups) if g == levels[1]]
    from scipy.stats import mannwhitneyu
    rows = []
    for idx, row in matrix.iterrows():
        va, vb = row[a].dropna().to_numpy(float), row[b].dropna().to_numpy(float)
        if len(va) < 2 or len(vb) < 2:
            p = np.nan
        else:
            p = float(mannwhitneyu(va, vb, alternative="two-sided").pvalue)
        rows.append({"chrom": idx[0], "start": idx[1], "end": idx[2], "mean_beta_group1": float(np.mean(va)) if len(va) else np.nan,
                     "mean_beta_group2": float(np.mean(vb)) if len(vb) else np.nan,
                     "delta_beta": float(np.mean(vb) - np.mean(va)) if len(va) and len(vb) else np.nan, "pvalue": p})
    result = pd.DataFrame(rows, columns=["chrom", "start", "end", "mean_beta_group1", "mean_beta_group2", "delta_beta", "pvalue"])
    result["padj"] = benjamini_hochberg(result["pvalue"].to_numpy(float))
    result["group1"] = levels[0]
    result["group2"] = levels[1]
    return result


def _filter_coverage(frame: pd.DataFrame, minimum: float) -> pd.DataFrame:
    """Do not fabricate coverage for beta-only input tables."""
    if minimum < 0:
        raise MethylationError("minimum coverage must be nonnegative")
    if minimum == 0:
        return frame
    if frame["coverage"].isna().any():
        raise MethylationError("coverage filtering requires measured coverage; use min_coverage=0 for beta-only tables")
    return frame[frame["coverage"] >= minimum]

