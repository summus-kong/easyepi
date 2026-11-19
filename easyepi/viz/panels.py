"""Composite multi-panel figures: gene track panels, UpSet, sample correlation.

`plot_gene_panel` is the "several genes x several marks" figure that opens most
epigenomics papers: one column per locus, one row per mark, with the gene model
underneath. `plot_upset` replaces a pair-wise overlap matrix when more than two
peak sets are compared. `plot_sample_correlation` gives the replicate/batch view.

All functions return a Figure and write nothing unless ``save`` is given, matching
the package-wide "no implicit figure export" contract.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence

import numpy as np


def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _finish(fig, save, formats, dpi=300):
    if not save:
        return fig
    from easyepi._plots import DEFAULT_FIGURE_FORMAT
    from easyepi.pubfig import finalize_figure
    return finalize_figure(fig, save,
                           formats=[DEFAULT_FIGURE_FORMAT] if formats is None
                           else list(formats),
                           dpi=dpi)


def _clean(ax, spines=("top", "right")):
    for spine in spines:
        ax.spines[spine].set_visible(False)


# ---------------------------------------------------------------------------
# bigWig region reader
# ---------------------------------------------------------------------------
def bigwig_bins(path, chrom: str, start: int, end: int, n_bins: int) -> np.ndarray:
    """Mean signal per bin over ``[start, end)`` from a bigWig.

    Uses ``stats(..., nBins=)`` rather than fetching per-base values: a 20 kb window
    at 25 bp bins would otherwise pull 20,000 values per locus per track just to
    average them in Python.
    """
    import pyBigWig
    if end <= start:
        raise ValueError("end must exceed start")
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")
    bw = pyBigWig.open(str(path))
    try:
        sizes = bw.chroms()
        if chrom not in sizes:
            raise ValueError(
                f"chromosome '{chrom}' absent from {path} "
                f"(has e.g. {list(sizes)[:5]})")
        vals = bw.stats(chrom, int(start), int(end), type="mean", nBins=int(n_bins))
    finally:
        bw.close()
    out = np.array([0.0 if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)
                    for v in (vals or [])], dtype=np.float64)
    if out.size != n_bins:
        # a region beyond the contig end can return fewer bins
        out = np.pad(out, (0, max(0, n_bins - out.size)))[:n_bins]
    return np.clip(out, 0, None)


def gene_window(gene: Mapping, flank: int = 2000) -> tuple[int, int]:
    """Locus window for a gene, padded by ``flank`` and clamped at 0.

    The span check is on the *gene* coordinates: padding a zero-length gene with a
    positive flank would still produce a positive-width window and hide the fact
    that the input has no gene body.
    """
    gs, ge = int(gene["start"]), int(gene["end"])
    if ge <= gs:
        raise ValueError(f"gene {gene.get('name', '?')} has an empty span")
    return max(0, gs - int(flank)), ge + int(flank)


# ---------------------------------------------------------------------------
# Multi-gene track panel
# ---------------------------------------------------------------------------
def plot_gene_panel(track_specs: Sequence[Mapping],
                    genes: Sequence[Mapping],
                    flank: int = 2000,
                    bin_size: int = 25,
                    exons: Optional[Mapping[str, Sequence[tuple[int, int]]]] = None,
                    shared_ylim: bool = True,
                    width_per_gene: float = 2.4,
                    height_per_track: float = 0.9,
                    title: Optional[str] = None,
                    save: Optional[str] = None,
                    formats: Optional[Sequence[str]] = None):
    """Grid of gene tracks: one column per locus, one row per mark.

    ``track_specs`` entries are ``{"name": str, "bigwig": path}``.
    ``genes`` entries are ``{"chrom", "start", "end", "name", "strand"}``; a
    ``"label"`` key overrides the displayed name. ``exons`` maps a gene name to
    ``[(start, end), ...]`` for a real gene model; without it the gene body is drawn
    as a single box, which is honest but coarse.

    y-limits are shared per track row by default so marks are comparable across
    loci. The shared limit uses a high quantile rather than the maximum, because a
    single saturated locus would otherwise flatten every other panel.
    """
    plt = _mpl()
    track_specs = list(track_specs or [])
    genes = list(genes or [])
    if not track_specs:
        raise ValueError("at least one track spec is required")
    if not genes:
        raise ValueError("at least one gene is required")
    for spec in track_specs:
        if "bigwig" not in spec:
            raise ValueError(f"track {spec.get('name', '?')} needs a 'bigwig' path")
    for gene in genes:
        for key in ("chrom", "start", "end"):
            if key not in gene:
                raise ValueError(f"gene entry needs '{key}': {gene}")
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")

    n_track, n_gene = len(track_specs), len(genes)
    windows = [gene_window(g, flank) for g in genes]

    # Pre-read every panel so a per-row y-limit can be chosen before drawing.
    data: dict[tuple[int, int], np.ndarray] = {}
    for ti, spec in enumerate(track_specs):
        for gi, (gene, (ws, we)) in enumerate(zip(genes, windows)):
            n_bins = max(2, int(np.ceil((we - ws) / bin_size)))
            data[(ti, gi)] = bigwig_bins(spec["bigwig"], gene["chrom"], ws, we, n_bins)

    # A single shared limit per track row can be dominated by one saturated locus
    # and flatten the others. The default is therefore a robust per-panel limit,
    # and "shared" means one limit per row computed from per-locus quantiles so
    # rows stay comparable without hiding low-signal loci.
    ylims: dict[int, tuple[float, float]] = {}
    for ti in range(n_track):
        per_locus = [float(np.percentile(data[(ti, gi)], 99))
                     for gi in range(n_gene) if data[(ti, gi)].size]
        per_locus = [v for v in per_locus if np.isfinite(v)]
        if shared_ylim and per_locus:
            hi = max(per_locus)                       # one limit, comparable rows
        else:
            hi = float(np.median(per_locus)) if per_locus else 1.0
        ylims[ti] = (0.0, hi if hi > 0 else 1.0)

    n_rows = n_track + 1  # +1 for the gene model
    fig, axes = plt.subplots(
        n_rows, n_gene,
        figsize=(width_per_gene * n_gene + 1.4, height_per_track * n_rows + 1.2),
        squeeze=False, sharex="col", sharey="row")

    for ti, spec in enumerate(track_specs):
        for gi, (gene, (ws, we)) in enumerate(zip(genes, windows)):
            ax = axes[ti, gi]
            sig = data[(ti, gi)]
            x = np.linspace(ws, we, sig.size)
            ax.fill_between(x, sig, color="#0F4D92", alpha=0.32, lw=0)
            ax.plot(x, sig, color="#0F4D92", lw=0.6)
            ax.set_xlim(ws, we)
            if shared_ylim:
                ax.set_ylim(*ylims[ti])
            else:
                local = float(np.percentile(sig, 99)) if sig.size else 1.0
                ax.set_ylim(0.0, local if local > 0 else 1.0)
            ax.tick_params(labelsize=5)
            _clean(ax)
            if gi == 0:
                ax.set_ylabel(str(spec.get("name", "")), fontsize=6)
            if ti == 0:
                ax.set_title(str(gene.get("label") or gene.get("name") or f"gene {gi+1}"),
                             fontsize=7)

    for gi, (gene, (ws, we)) in enumerate(zip(genes, windows)):
        ax = axes[-1, gi]
        name = gene.get("name")
        spans = (exons or {}).get(name) or [(int(gene["start"]), int(gene["end"]))]
        ax.hlines(0.5, ws, we, color="#4d4d4d", lw=0.6)
        # Draw the tall exon boxes first, then thin connector lines for the parts
        # of the gene that fall between exons. Without the connector a gene with
        # many short exons reads as a row of unrelated bars.
        span_w = max(1, we - ws)
        min_visible = max(1, int(0.0015 * span_w))
        for (s, e) in spans:
            left = max(s, ws)
            width = max(min_visible, min(e, we) - left)
            ax.add_patch(plt.Rectangle((left, 0.28), width, 0.44,
                                       color="#B64342", lw=0))
        # strand arrow gives direction without an extra legend entry
        if str(gene.get("strand", "+")) == "-":
            ax.annotate("", xy=(we - (we - ws) * 0.02, 0.5), xytext=(ws + (we - ws) * 0.12, 0.5),
                        arrowprops=dict(arrowstyle="->", color="#4d4d4d", lw=0.6))
        else:
            ax.annotate("", xy=(we - (we - ws) * 0.12, 0.5), xytext=(ws + (we - ws) * 0.02, 0.5),
                        arrowprops=dict(arrowstyle="->", color="#4d4d4d", lw=0.6))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        # Genomic coordinates carry the axis: without ticks a reader cannot tell
        # where a peak sits relative to the gene, which is the point of the figure.
        span = we - ws
        step = max(1_000, int(round(span / 4 / 1_000)) * 1_000)
        ticks = np.arange(ws, we + 1, step)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{int(t) / 1000:.0f}" for t in ticks], fontsize=5)
        ax.set_xlabel(f"{gene.get('chrom', '')} (kb)", fontsize=5)
        ax.tick_params(labelsize=5, labelbottom=True)
        _clean(ax)

    fig.suptitle(title or f"{n_gene} loci x {n_track} tracks", fontsize=8)
    fig.tight_layout()
    return _finish(fig, save, formats)



# ---------------------------------------------------------------------------
# Gene models from GTF
# ---------------------------------------------------------------------------
def read_gene_models(path, genes: Sequence[str] | None = None,
                     max_genes: int | None = None) -> tuple[list[dict], dict[str, list[tuple[int, int]]]]:
    """Read gene bodies and exons from a GTF/GFF.

    Returns ``(gene_records, exons_by_name)``. Records carry
    ``chrom/start/end/name/strand`` in BED-like 0-based half-open coordinates.
    ``genes`` restricts the output to those symbols; a requested symbol that is
    absent raises, because silently drawing fewer panels than asked would hide a
    name mismatch between the GTF and the user's list.
    """
    wanted = {g.strip() for g in genes} if genes else None
    bodies: dict[str, dict] = {}
    exons: dict[str, list[tuple[int, int]]] = {}
    opener = __import__("gzip").open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _src, feature = f[0], f[1], f[2]
            try:
                start, end = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            attrs = f[8]
            name = None
            for key in ("gene_name", "gene", "gene_id", "transcript_id"):
                marker = key + ' "'
                idx = attrs.find(marker)
                if idx >= 0:
                    rest = attrs[idx + len(marker):]
                    name = rest.split('"')[0]
                    break
            if not name:
                continue
            if wanted is not None and name not in wanted:
                continue
            if feature == "gene" or (feature in {"transcript", "mRNA"} and name not in bodies):
                bodies.setdefault(name, {"chrom": chrom, "start": start, "end": end,
                                         "name": name, "strand": f[6]})
            elif feature == "exon":
                exons.setdefault(name, []).append((start, end))
                bodies.setdefault(name, {"chrom": chrom, "start": start, "end": end,
                                         "name": name, "strand": f[6]})
    if wanted is not None:
        missing = sorted(wanted - set(bodies))
        if missing:
            raise ValueError(
                f"gene symbol(s) not found in {path}: {', '.join(missing[:10])}")
    records = list(bodies.values())
    records.sort(key=lambda g: (g["chrom"], g["start"]))
    if max_genes is not None:
        records = records[: int(max_genes)]
    for rec in records:
        spans = exons.get(rec["name"])
        if spans:
            merged = sorted(spans)
            out = [merged[0]]
            for s0, e0 in merged[1:]:
                if s0 <= out[-1][1]:
                    out[-1] = (out[-1][0], max(out[-1][1], e0))
                else:
                    out.append((s0, e0))
            exons[rec["name"]] = out
    return records, {k: v for k, v in exons.items() if k in {r["name"] for r in records}}


def parse_track_specs(items):
    """Resolve ``--bigwig`` entries into ``{"name", "bigwig"}`` records.

    ``--bigwig`` takes a path or ``NAME=PATH``; a directory flag was removed because
    it made the common case (a path you already have) into two arguments, and it
    could not express a path outside that one folder.

    * ``tracks/H3K4me3.bw`` -> name ``H3K4me3`` (stem), absolute or relative
    * ``H3K4me3=tracks/H3K4me3.bw`` -> explicit label, useful when the filename is
      not what the figure should say

    Relative paths are resolved against the current working directory, so a run from
    a project root behaves as expected. Every path is checked here: a missing track
    must fail before any figure work starts, not part-way through.
    """
    out: list[dict] = []
    for raw in items or []:
        for piece in str(raw).split(","):
            entry = piece.strip()
            if not entry:
                continue
            if "=" in entry:
                name, path = entry.split("=", 1)
                name, path = name.strip(), path.strip()
            else:
                path = entry
                name = ""
            if not path:
                raise ValueError(f"malformed --bigwig entry {entry!r}: no path")
            from easyepi.errors import require_file

            resolved = require_file(path, purpose="BigWig track",
                                    suffixes=(".bw", ".bigWig", ".bigwig"))
            if not name:
                name = resolved.stem
            out.append({"name": name, "bigwig": str(resolved)})
    if not out:
        raise ValueError("at least one --bigwig is required")
    names = [t["name"] for t in out]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        raise ValueError(
            f"duplicate track name(s): {', '.join(dupes)}; "
            "label one explicitly with NAME=PATH")
    return out


def expand_csv_values(items):
    """Split repeated and comma-separated values into one flat list.

    ``--signal H3K4me3=a.bg,H3K27ac=b.bg`` behaves like two flags; this is the one
    place that rule is expressed so every command agrees on it.
    """
    out = []
    for raw in items or []:
        for piece in str(raw).split(","):
            value = piece.strip()
            if value:
                out.append(value)
    return out

def parse_gene_names(items):
    """Flatten repeatable or comma-separated gene symbols into one ordered list.

    Both ``--gene A --gene B`` and ``--gene A,B`` are accepted, and the order the
    user wrote is preserved so the panels appear in that order.
    """
    names = []
    for raw in items or []:
        for piece in str(raw).split(","):
            name = piece.strip()
            if name and name not in names:
                names.append(name)
    return names


def run_gene_panel(args) -> int:
    """CLI driver: build the multi-locus x multi-mark panel figure."""
    import json
    from pathlib import Path

    from easyepi._plots import selected_formats

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    specs = parse_track_specs(args.bigwig)
    gene_names = parse_gene_names(getattr(args, "gene", None))
    if args.gtf:
        genes, exons = read_gene_models(args.gtf, genes=gene_names or None,
                                        max_genes=args.max_genes)
        if not genes:
            raise ValueError(f"no genes selected from {args.gtf}")
    else:
        if not args.bed:
            raise ValueError("provide --gtf or --bed with chrom/start/end/name columns")
        genes, exons = [], {}
        opener = __import__("gzip").open if str(args.bed).endswith(".gz") else open
        with opener(args.bed, "rt", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                genes.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                              "name": f[3], "strand": f[5] if len(f) > 5 else "+"})
        if gene_names:
            keep = set(gene_names)
            genes = [g for g in genes if g["name"] in keep]
        if not genes:
            raise ValueError("no genes to draw")

    formats = selected_formats(args)
    fig = plot_gene_panel(specs, genes, flank=int(args.flank),
                          bin_size=int(args.bin_size), exons=exons or None,
                          shared_ylim=not args.independent_ylim,
                          title=args.title,
                          save=str(out / "gene_panel") if formats else None,
                          formats=list(formats) if formats else None)

    manifest = {
        "tracks": specs,
        "genes": [{"name": g["name"], "chrom": g["chrom"],
                   "start": int(g["start"]), "end": int(g["end"]),
                   "strand": g.get("strand", "+")} for g in genes],
        "flank": int(args.flank),
        "bin_size": int(args.bin_size),
        "shared_ylim": not args.independent_ylim,
        "figures": sorted(p.name for p in out.glob("gene_panel.*")),
    }
    (out / "gene_panel.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[easyepi] gene-panel: {len(specs)} track(s) x {len(genes)} gene(s) "
          f"-> {out}")
    for g in genes:
        print(f"    {g['name']:12s} {g['chrom']}:{int(g['start'])}-{int(g['end'])} "
              f"({g.get('strand', '+')})")
    if formats:
        print(f"[easyepi] wrote gene_panel.{formats[0]}")
    return 0

# ---------------------------------------------------------------------------
# UpSet
# ---------------------------------------------------------------------------
def _segments(per_chrom: dict[str, list[tuple[int, int]]]) -> list[tuple[str, int, int]]:
    """Split intervals into elementary segments where set membership is constant.

    Overlaps of two sets can partially intersect, so testing whole merged
    intervals against the other set mislabels the remainder: A=[100,200],
    B=[150,250] contains an A-only part, a shared part and a B-only part. Emitting
    the three elementary segments and testing each one is the only way an
    intersection figure can be correct.
    """
    out: list[tuple[str, int, int]] = []
    for chrom, spans in per_chrom.items():
        if not spans:
            continue
        edges = sorted({p for span in spans for p in span})
        for left, right in zip(edges, edges[1:]):
            if right > left:
                out.append((chrom, left, right))
    return out


def _membership(sets: Mapping[str, np.ndarray]
                ) -> tuple[list[str], list[tuple[int, ...]], np.ndarray]:
    """Membership vector and size of every elementary union segment.

    Segments (not whole intervals) are the unit, so a region shared by A and B and
    a region unique to either set are all counted exactly once.
    """
    labels = list(sets)
    all_spans: dict[str, list[tuple[int, int]]] = {}
    for arr in sets.values():
        for chrom_i, arr_c in arr.items():
            for s, e in arr_c:
                all_spans.setdefault(chrom_i, []).append((int(s), int(e)))
    if not all_spans:
        raise ValueError("peak sets contain no intervals")

    segments = _segments(all_spans)
    if not segments:
        raise ValueError("peak sets contain no intervals")

    grouped: dict[str, np.ndarray] = {}
    for name, arr in sets.items():
        for chrom, sub in arr.items():
            grouped.setdefault(name, {})[chrom] = sub  # type: ignore[index]

    patterns: dict[tuple[int, ...], int] = {}
    for chrom, s, e in segments:
        key: list[int] = []
        for name in labels:
            sub = grouped[name].get(chrom)  # type: ignore[union-attr]
            hit = 0
            if sub is not None and len(sub):
                hi = int(np.searchsorted(sub[:, 0], e, side="left"))
                if hi:
                    shared = np.minimum(sub[:hi, 1], e) - np.maximum(sub[:hi, 0], s)
                    hit = int(np.any(shared >= 1))
            key.append(hit)
        if not any(key):
            # a gap between intervals belongs to no set and must not be counted
            continue
        tup = tuple(key)
        patterns[tup] = patterns.get(tup, 0) + 1
    if not patterns:
        raise ValueError("no region belongs to any peak set")
    return labels, list(patterns), np.array(list(patterns.values()), dtype=np.int64)


def plot_upset(peak_sets: Mapping[str, Sequence[tuple[str, int, int]]],
               max_intersections: int = 15, min_size: int = 1,
               title: str = "Intersection of peak sets",
               save: Optional[str] = None,
               formats: Optional[Sequence[str]] = None):
    """UpSet plot: bar chart of intersection sizes with a membership dot matrix.

    Better than a pairwise matrix once there are three or more sets: it shows the
    exact combination each region belongs to, which a heatmap of pair-wise
    overlaps cannot express.
    """
    plt = _mpl()
    from easyepi.viz.genomeviews import _merge

    if len(peak_sets) < 2:
        raise ValueError("UpSet needs at least two sets")
    indexed = {name: _merge(rows) for name, rows in peak_sets.items()}
    labels, patterns, counts = _membership(indexed)

    order = np.argsort(-counts)
    kept = [(patterns[i], int(counts[i])) for i in order
            if counts[i] >= min_size and sum(patterns[i]) > 0][:max_intersections]
    if not kept:
        raise ValueError(f"no intersection has at least {min_size} region(s)")

    set_sizes = []
    for name in labels:
        arr = indexed[name]
        set_sizes.append(int(sum(int(e - s) for a in arr.values() for s, e in a)))

    n_col = len(kept)
    fig = plt.figure(figsize=(max(3.2, 0.42 * n_col + 2.2), 3.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, max(2.0, 0.42 * n_col)],
                          height_ratios=[2.2, 1.3], wspace=0.08, hspace=0.06)
    ax_sets = fig.add_subplot(gs[:, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1])

    y = np.arange(len(labels))
    ax_sets.barh(y, set_sizes, color="#4d4d4d", height=0.6)
    ax_sets.set_yticks(y)
    ax_sets.set_yticklabels(labels, fontsize=6)
    ax_sets.invert_yaxis()
    ax_sets.invert_xaxis()
    ax_sets.set_xlabel("regions (bp total)", fontsize=6)
    ax_sets.tick_params(labelsize=5)
    _clean(ax_sets)

    x = np.arange(n_col)
    ax_top.bar(x, [c for _, c in kept], color="#B64342", width=0.62)
    ax_top.set_ylabel("intersections", fontsize=6)
    ax_top.tick_params(labelsize=5, labelbottom=False)
    for xi, (pattern, c) in enumerate(kept):
        if c > 0:
            ax_top.text(xi, c, str(c), ha="center", va="bottom", fontsize=5)
    _clean(ax_top)

    for xi, (pattern, _c) in enumerate(kept):
        members = [i for i, v in enumerate(pattern) if v]
        if len(members) > 1:
            ax_mat.plot([xi, xi], [min(members), max(members)], color="#4d4d4d", lw=0.8, zorder=1)
        for i in range(len(labels)):
            on = bool(pattern[i])
            ax_mat.scatter([xi], [i], s=16, zorder=2,
                           color="#B64342" if on else "#d9d9d9")
    ax_mat.set_xlim(-0.6, n_col - 0.4)
    ax_mat.set_ylim(len(labels) - 0.4, -0.6)
    ax_mat.set_yticks(y)
    ax_mat.set_yticklabels(labels, fontsize=6)
    ax_mat.set_xticks(x)
    ax_mat.set_xticklabels([])
    ax_mat.tick_params(labelsize=5)
    ax_mat.set_ylabel("", fontsize=6)
    _clean(ax_mat, spines=("top", "right", "bottom", "left"))
    ax_mat.grid(axis="x", color="#e6e6e6", lw=0.5, zorder=0)

    fig.suptitle(title, fontsize=8)
    return _finish(fig, save, formats)


# ---------------------------------------------------------------------------
# Sample correlation
# ---------------------------------------------------------------------------
def correlation_matrix(frame, method: str = "pearson",
                       min_features: int = 10) -> tuple[np.ndarray, list[str]]:
    """Pair-wise sample correlation, dropping features with no variance.

    Constant rows carry no information and make Pearson undefined, so they are
    removed first instead of being allowed to emit NaN for a whole row.
    """
    import pandas as pd

    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame (features x samples)")
    if frame.shape[1] < 2:
        raise ValueError("need at least two samples")
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.loc[:, numeric.notna().sum(axis=0) > 0]
    numeric = numeric.fillna(0.0)
    # Drop features with no variance: they cannot correlate and would emit NaN for
    # a whole row AND column.
    varying = numeric.loc[numeric.std(axis=1) > 0]
    # Drop samples with no variance for the same reason: a perfectly flat sample
    # leaves a NaN row/column and makes the whole matrix non-finite.
    varying = varying.loc[:, varying.std(axis=0) > 0]
    if varying.shape[0] < min_features:
        raise ValueError(
            f"only {varying.shape[0]} features vary across samples; "
            f"at least {min_features} are needed for a meaningful correlation")
    if varying.shape[1] < 2:
        raise ValueError(
            "fewer than two samples retain any variance; every pair would be "
            "undefined")
    corr = varying.corr(method=method).to_numpy(dtype=np.float64)
    if not np.isfinite(corr).all():
        raise ValueError(
            "correlation is undefined for this matrix; check for samples with "
            "constant values across all features")
    return corr, [str(c) for c in varying.columns]


def plot_sample_correlation(frame, method: str = "pearson", cluster: bool = True,
                            title: str = "Sample correlation",
                            cmap: str = "RdBu_r",
                            save: Optional[str] = None,
                            formats: Optional[Sequence[str]] = None):
    """Correlation heatmap with an optional hierarchical-clustering dendrogram.

    A heatmap alone shows which samples are alike but not the grouping; the
    dendrogram makes the replicate structure explicit, which is the question the
    figure is usually included to answer.
    """
    plt = _mpl()
    corr, labels = correlation_matrix(frame, method=method)
    n = len(labels)

    if cluster and n >= 3:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform
        # 1 - r is a proper distance for correlation
        dist = 1.0 - corr
        np.fill_diagonal(dist, 0.0)
        dist = (dist + dist.T) / 2.0
        condensed = squareform(dist, checks=False)
        link = linkage(condensed, method="average")
        fig = plt.figure(figsize=(0.42 * n + 2.4, 0.42 * n + 2.0))
        gs = fig.add_gridspec(2, 2, width_ratios=[0.28, 1], height_ratios=[0.28, 1],
                              wspace=0.04, hspace=0.04)
        ax_d = fig.add_subplot(gs[0, 1])
        ax_h = fig.add_subplot(gs[1, 1])
        ax_l = fig.add_subplot(gs[1, 0])

        d = dendrogram(link, ax=ax_d, no_labels=True, color_threshold=0,
                       above_threshold_color="#4d4d4d")
        ax_d.set_axis_off()
        order = d["leaves"]

        _draw_corr(plt, ax_h, corr, labels, order, cmap)
        ax_l.axis("off")
        for i, idx in enumerate(order):
            ax_l.text(0.98, 1 - (i + 0.5) / n, labels[idx], ha="right", va="center",
                      fontsize=6, transform=ax_l.transAxes)
    else:
        fig, ax_h = plt.subplots(figsize=(0.42 * n + 2.0, 0.42 * n + 1.4))
        _draw_corr(plt, ax_h, corr, labels, list(range(n)), cmap)

    fig.suptitle(f"{title} ({method})", fontsize=8)
    return _finish(fig, save, formats)


def _draw_corr(plt, ax, corr, labels, order, cmap):
    sub = corr[np.ix_(order, order)]
    im = ax.imshow(sub, vmin=-1, vmax=1, cmap=cmap)
    ticks = np.arange(len(order))
    ax.set_xticks(ticks)
    ax.set_xticklabels([labels[i] for i in order], rotation=45, ha="right", fontsize=6)
    ax.set_yticks(ticks)
    ax.set_yticklabels([labels[i] for i in order], fontsize=6)
    for i in range(len(order)):
        for j in range(len(order)):
            ax.text(j, i, f"{sub[i, j]:.2f}", ha="center", va="center", fontsize=5,
                    color="white" if abs(sub[i, j]) > 0.6 else "#222222")
    fig = ax.figure
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(labelsize=5)
