"""Additional analysis figures: MA plot, peak-set overlap, TSS heatmap, locus track.

Each function returns a Matplotlib Figure so it composes with the existing
``easyepi.plot`` / ``easyepi.pubfig`` helpers, and writes nothing unless ``save``
is given. That keeps the package's "no implicit figure export" contract.

Why these four
--------------
* **MA plot** — the standard sanity view for a differential table. A volcano shows
  significance vs effect; MA shows effect vs abundance, which is where library
  composition and low-count bias become visible.
* **Peak overlap** — replicate concordance and cross-mark sharing are argued from
  an overlap matrix, not from a single union count.
* **TSS heatmap** — the per-region view of a metagene profile. A mean profile can
  look correct while half the loci are flat; the heatmap exposes that.
* **Locus track** — the only way to show that a called position, peak or state
  actually lines up with the signal it was called from.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _finish(fig, save, formats, dpi=300):
    if not save:
        return fig
    from easyepi.pubfig import finalize_figure
    return finalize_figure(fig, save,
                           formats=["png"] if formats is None else list(formats),
                           dpi=dpi)


def _clean(ax, spines=("top", "right")):
    for spine in spines:
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# MA plot
# ---------------------------------------------------------------------------
def plot_ma(frame, log2fc_col: str = "log2fc", mean_col: Optional[str] = None,
            padj_col: str = "padj", alpha: float = 0.05,
            title: str = "MA plot", pseudo: float = 0.5,
            save: Optional[str] = None, formats: Optional[Sequence[str]] = None):
    """Mean-abundance vs fold-change, with the significance cut-off shaded.

    ``mean_col`` is used when the table carries a mean/logCPM column. Otherwise the
    abundance proxy is recovered as ``log2`` of the mean count when a ``baseMean``
    column exists; if neither is present the function raises instead of inventing
    an x axis, because a MA plot without abundance is just a strip plot.
    """
    plt = _mpl()
    import pandas as pd

    if log2fc_col not in frame.columns:
        raise ValueError(f"table has no '{log2fc_col}' column")
    lfc = pd.to_numeric(frame[log2fc_col], errors="coerce")
    if mean_col is not None:
        if mean_col not in frame.columns:
            raise ValueError(f"table has no '{mean_col}' column")
        abundance = pd.to_numeric(frame[mean_col], errors="coerce")
        xlabel = mean_col
    elif "baseMean" in frame.columns:
        abundance = np.log2(pd.to_numeric(frame["baseMean"], errors="coerce") + pseudo)
        xlabel = f"log2(baseMean + {pseudo})"
    else:
        raise ValueError(
            "need an abundance column to draw an MA plot: pass mean_col= or "
            "provide a 'baseMean' column")
    padj = (pd.to_numeric(frame[padj_col], errors="coerce")
            if padj_col in frame.columns else pd.Series(np.nan, index=frame.index))

    ok = lfc.notna() & abundance.notna()
    if not ok.any():
        raise ValueError("no rows with both fold change and abundance")

    sig = ok & padj.notna() & (padj < alpha)
    ns = ok & ~sig
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.scatter(abundance[ns], lfc[ns], s=5, color="#9a9a9a", lw=0, alpha=0.6,
               label="not significant")
    if sig.any():
        ax.scatter(abundance[sig], lfc[sig], s=6, color="#B64342", lw=0,
                   label=f"padj < {alpha}")
    ax.axhline(0, color="#444444", lw=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("log2 fold change")
    ax.set_title(title, fontsize=8)
    ax.legend(frameon=False, fontsize=6, markerscale=2)
    _clean(ax)
    fig.tight_layout()
    return _finish(fig, save, formats)


# ---------------------------------------------------------------------------
# Peak-set overlap
# ---------------------------------------------------------------------------
def read_peak_sets(paths: Sequence[str | Path]) -> dict[str, list[tuple[str, int, int]]]:
    """Read each BED/narrowPeak into ``{label: [(chrom, start, end), ...]}``."""
    sets: dict[str, list[tuple[str, int, int]]] = {}
    for path in paths:
        p = Path(path)
        opener = __import__("gzip").open if p.suffix == ".gz" else open
        rows: list[tuple[str, int, int]] = []
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                try:
                    s, e = int(f[1]), int(f[2])
                except ValueError:
                    continue
                if e > s:
                    rows.append((f[0], s, e))
        if not rows:
            raise ValueError(f"no intervals read from {path}")
        sets[p.stem] = rows
    if not sets:
        raise ValueError("no peak sets were read")
    return sets


def _sorted_by_chrom(rows):
    by: dict[str, np.ndarray] = {}
    for chrom, s, e in rows:
        by.setdefault(chrom, []).append((s, e))
    return {c: np.asarray(sorted(v), dtype=np.int64) for c, v in by.items()}


def overlap_fraction(a: list[tuple[str, int, int]],
                     b: list[tuple[str, int, int]],
                     min_overlap: int = 1) -> float:
    """Fraction of ``a`` intervals that overlap any interval of ``b``.

    Uses a binary search per query interval, so it stays O(n log n) on real peak
    sets. ``min_overlap`` is the number of shared bp required (1 = "any overlap"),
    matching how union/consensus sets are usually built.
    """
    idx = _sorted_by_chrom(b)
    total = len(a)
    if total == 0:
        return 0.0
    hits = 0
    for chrom, s, e in a:
        arr = idx.get(chrom)
        if arr is None:
            continue
        starts = arr[:, 0]
        # candidates whose start is before this interval ends
        hi = int(np.searchsorted(starts, e, side="left"))
        if hi == 0:
            continue
        cand_starts = arr[:hi, 0]
        cand_ends = arr[:hi, 1]
        # shared bp is bounded on BOTH sides: clipping only one end reports the
        # whole length of the candidate for a partial overlap.
        shared = np.minimum(cand_ends, e) - np.maximum(cand_starts, s)
        if np.any(shared >= min_overlap):
            hits += 1
    return hits / total


def overlap_count(a: list[tuple[str, int, int]], b: list[tuple[str, int, int]],
                  min_overlap: int = 1) -> int:
    """Number of intervals in ``a`` that overlap any interval of ``b``."""
    idx = _sorted_by_chrom(b)
    hits = 0
    for chrom, s, e in a:
        arr = idx.get(chrom)
        if arr is None:
            continue
        hi = int(np.searchsorted(arr[:, 0], e, side="left"))
        if hi == 0:
            continue
        shared = np.minimum(arr[:hi, 1], e) - np.maximum(arr[:hi, 0], s)
        if np.any(shared >= min_overlap):
            hits += 1
    return hits


def _merge(rows: list[tuple[str, int, int]]) -> dict[str, np.ndarray]:
    """Merge overlapping/adjacent intervals per chromosome."""
    by: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e in rows:
        if e > s:
            by.setdefault(chrom, []).append((s, e))
    merged: dict[str, np.ndarray] = {}
    for chrom, spans in by.items():
        spans.sort()
        out = []
        cs, ce = spans[0]
        for s, e in spans[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                out.append((cs, ce))
                cs, ce = s, e
        out.append((cs, ce))
        merged[chrom] = np.asarray(out, dtype=np.int64)
    return merged


def basepair_overlap(a: list[tuple[str, int, int]],
                     b: list[tuple[str, int, int]]) -> dict[str, int]:
    """Total bp shared, bp in a, bp in b — computed on merged intervals.

    Interval *counts* mislead when one set is fragmented or broad: a broadPeak
    domain covering many narrowPeaks can touch dozens of intervals and inflate a
    count-based overlap while sharing little sequence. Merging first and summing
    bp removes both the fragmentation and the width bias.
    """
    ma, mb = _merge(a), _merge(b)
    shared = 0
    for chrom, spans_a in ma.items():
        spans_b = mb.get(chrom)
        if spans_b is None:
            continue
        i = j = 0
        while i < len(spans_a) and j < len(spans_b):
            lo = max(spans_a[i, 0], spans_b[j, 0])
            hi = min(spans_a[i, 1], spans_b[j, 1])
            if hi > lo:
                shared += int(hi - lo)
            if spans_a[i, 1] <= spans_b[j, 1]:
                i += 1
            else:
                j += 1
    bp_a = int(sum(int(e - s) for arr in ma.values() for s, e in arr))
    bp_b = int(sum(int(e - s) for arr in mb.values() for s, e in arr))
    return {"shared_bp": shared, "bp_a": bp_a, "bp_b": bp_b}


def jaccard_bp(a: list[tuple[str, int, int]],
               b: list[tuple[str, int, int]]) -> float:
    """Base-pair Jaccard: shared bp / (bp_a + bp_b - shared)."""
    o = basepair_overlap(a, b)
    denom = o["bp_a"] + o["bp_b"] - o["shared_bp"]
    return float(o["shared_bp"] / denom) if denom > 0 else 0.0


def jaccard_index(a: list[tuple[str, int, int]], b: list[tuple[str, int, int]],
                  min_overlap: int = 1) -> float:
    """Symmetric set similarity: shared / (|a| + |b| - shared).

    Preferred over a directional "fraction of a covered by b" for a matrix,
    because the directional version is intentionally asymmetric and reads like a
    bug: a subset set is fully covered by a superset, not the other way round.
    """
    shared = overlap_count(a, b, min_overlap)
    denom = len(a) + len(b) - shared
    return float(shared / denom) if denom > 0 else 0.0


def plot_peak_overlap(paths: Sequence[str | Path], min_overlap: int = 1,
                      title: str = "Peak-set overlap",
                      save: Optional[str] = None,
                      formats: Optional[Sequence[str]] = None,
                      metric: str = "bp_jaccard"):
    """Symmetric overlap matrix between peak sets.

    ``metric="bp_jaccard"`` (default) is symmetric, base-pair weighted, and is the
    right choice when sets differ in width (broad vs narrow peaks) or
    fragmentation; ``metric="jaccard"`` counts intervals;
    ``metric="covered"`` is the directional fraction of the row set that the
    column set covers, which is occasionally what a replicate-concordance argument
    needs but is deliberately asymmetric.
    """
    plt = _mpl()
    sets = read_peak_sets(paths)
    if len(sets) < 2:
        raise ValueError("need at least two peak sets to compare")
    labels = list(sets)
    n = len(labels)
    if metric not in {"bp_jaccard", "jaccard", "covered"}:
        raise ValueError("metric must be 'bp_jaccard', 'jaccard' or 'covered'")
    mat = np.zeros((n, n))
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            if i == j:
                mat[i, j] = 1.0
            elif metric == "bp_jaccard":
                mat[i, j] = jaccard_bp(sets[a], sets[b])
            elif metric == "jaccard":
                mat[i, j] = jaccard_index(sets[a], sets[b], min_overlap)
            else:
                mat[i, j] = overlap_fraction(sets[a], sets[b], min_overlap)
    fig, ax = plt.subplots(figsize=(0.62 * n + 2.0, 0.55 * n + 1.8))
    im = ax.imshow(mat, vmin=0, vmax=1, cmap="Blues")
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=6)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if mat[i, j] > 0.6 else "#222222")
    label = {"bp_jaccard": "base-pair Jaccard",
             "jaccard": "interval Jaccard",
             "covered": "row covered by column"}[metric]
    ax.set_title(f"{title} ({label})", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=6)
    fig.tight_layout()
    return _finish(fig, save, formats)


# ---------------------------------------------------------------------------
# TSS heatmap
# ---------------------------------------------------------------------------
def loci_matrix(signal: np.ndarray, sites: Sequence[int], flank: int,
                offset: int = 0) -> np.ndarray:
    """Rows = loci, columns = ``-flank..+flank`` bins around each site."""
    sig = np.asarray(signal, dtype=np.float64)
    width = 2 * flank + 1
    rows = []
    for pos in sites:
        p = int(pos) - int(offset)
        lo, hi = p - flank, p + flank + 1
        if lo < 0 or hi > sig.size:
            continue
        rows.append(sig[lo:hi])
    if not rows:
        return np.empty((0, width))
    return np.vstack(rows)


def plot_tss_heatmap(matrix, flank: int, title: str = "Signal around TSS",
                     vmax_percentile: float = 98.0, cluster: bool = True,
                     max_rows: int = 4000, save: Optional[str] = None,
                     formats: Optional[Sequence[str]] = None):
    """Heatmap of per-locus profiles, sorted by signal, with the mean on top.

    Rows are normalised per locus only for ordering, so a few very strong loci do
    not flatten everything else. ``vmax_percentile`` clips the colour scale at a
    robust quantile instead of the maximum for the same reason.
    """
    plt = _mpl()
    from easyepi.pubfig import make_heatmap

    mat = np.asarray(matrix, dtype=np.float64)
    if mat.ndim != 2 or mat.shape[0] == 0:
        raise ValueError("matrix must be a non-empty (loci, bins) array")
    if mat.shape[1] != 2 * flank + 1:
        raise ValueError(f"matrix has {mat.shape[1]} columns; flank={flank} needs "
                         f"{2 * flank + 1}")
    if mat.shape[0] > max_rows:
        step = int(np.ceil(mat.shape[0] / max_rows))
        mat = mat[::step]

    centre = flank
    strength = mat[:, centre]
    order = np.argsort(-strength) if cluster else np.arange(mat.shape[0])
    ordered = mat[order]
    vmax = float(np.percentile(ordered, vmax_percentile)) or 1.0

    fig, (ax_top, ax_heat) = plt.subplots(
        2, 1, figsize=(3.6, 4.2), gridspec_kw={"height_ratios": [1, 3]},
        sharex=True)
    mean = ordered.mean(axis=0)
    x = np.arange(-flank, flank + 1)
    ax_top.plot(x, mean, color="#0F4D92", lw=1.2)
    ax_top.set_ylabel("mean", fontsize=6)
    ax_top.set_title(title, fontsize=8)
    ax_top.tick_params(labelsize=6)
    _clean(ax_top)

    im = ax_heat.imshow(ordered, aspect="auto", cmap="viridis", vmin=0, vmax=vmax,
                        extent=(-flank, flank, ordered.shape[0], 0),
                        interpolation="nearest")
    ax_heat.axvline(0, color="white", lw=0.6, ls="--")
    ax_heat.set_xlabel("distance from site (bp)")
    ax_heat.set_ylabel(f"{ordered.shape[0]} loci", fontsize=6)
    ax_heat.tick_params(labelsize=6)
    fig.colorbar(im, ax=ax_heat, fraction=0.046, pad=0.02).ax.tick_params(labelsize=6)
    fig.tight_layout()
    return _finish(fig, save, formats)


# ---------------------------------------------------------------------------
# Locus track
# ---------------------------------------------------------------------------
def plot_locus_track(tracks: Sequence[dict], chrom: str, start: int, end: int,
                     annotations: Optional[Sequence[dict]] = None,
                     title: Optional[str] = None, save: Optional[str] = None,
                     formats: Optional[Sequence[str]] = None):
    """Stacked genome-browser view of one region.

    ``tracks`` is a sequence of ``{"name": str, "signal": array, "bin_size": int}``
    (or ``{"name":..., "intervals": [(s, e), ...]}`` for feature tracks).
    ``annotations`` adds labelled spans drawn under the tracks. Coordinates are
    interpreted as genomic bp; arrays are placed by ``bin_size`` from ``start``.
    """
    plt = _mpl()
    tracks = list(tracks or [])
    if not tracks:
        raise ValueError("at least one track is required")
    if end <= start:
        raise ValueError("end must exceed start")

    n = len(tracks) + (1 if annotations else 0)
    fig, axes = plt.subplots(n, 1, figsize=(6.4, 1.0 * n + 0.9), sharex=True,
                             squeeze=False)
    axes = axes[:, 0]
    for ax, tr in zip(axes, tracks):
        name = tr.get("name", "track")
        if "signal" in tr and tr["signal"] is not None:
            sig = np.asarray(tr["signal"], dtype=np.float64)
            bin_size = int(tr.get("bin_size", 1))
            x = start + np.arange(sig.size) * bin_size
            ax.fill_between(x, sig, color="#0F4D92", alpha=0.3, lw=0)
            ax.plot(x, sig, color="#0F4D92", lw=0.7)
            ax.set_ylabel("signal", fontsize=6)
        elif "intervals" in tr:
            for (s, e) in tr["intervals"]:
                if e <= start or s >= end:
                    continue
                ax.add_patch(plt.Rectangle((max(s, start), 0.25),
                                           min(e, end) - max(s, start), 0.5,
                                           color="#B64342", lw=0))
            ax.set_ylim(0, 1)
            ax.set_yticks([])
        ax.set_title(name, fontsize=7, loc="left")
        ax.tick_params(labelsize=6)
        _clean(ax)

    if annotations:
        ax = axes[-1]
        for ann in annotations:
            s, e = int(ann.get("start", start)), int(ann.get("end", end))
            if e <= start or s >= end:
                continue
            ax.add_patch(plt.Rectangle((max(s, start), 0.3), min(e, end) - max(s, start),
                                       0.4, color=ann.get("color", "#4d4d4d"), lw=0))
            ax.text(max(s, start), 0.78, ann.get("label", ""), fontsize=5,
                    ha="left", va="bottom")
        ax.set_ylim(0, 1); ax.set_yticks([])
        ax.set_title("annotations", fontsize=7, loc="left")
        _clean(ax)

    axes[-1].set_xlim(start, end)
    axes[-1].set_xlabel(f"{chrom} position (bp)")
    fig.suptitle(title or f"{chrom}:{start:,}-{end:,}", fontsize=8)
    fig.tight_layout()
    return _finish(fig, save, formats)
