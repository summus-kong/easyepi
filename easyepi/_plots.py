"""Shared plotting functions and optional module-output rendering. Plotting reads completed analysis tables without altering statistical results."""
from __future__ import annotations

import logging

import json
from pathlib import Path
from typing import Any


#: The format every figure is written in unless the caller asks for another one.
#: PDF is a vector format, so a figure stays sharp at any zoom and can be placed
#: in a manuscript without re-exporting; PNG was the old default and is raster, so
#: text and thin lines degraded as soon as a reader zoomed in. This is the single
#: place the default lives: every plotting entry point reads it from here, so
#: changing it changes all of them at once.
DEFAULT_FIGURE_FORMAT = "pdf"

#: Accepted --figure-format values. Raster options stay available because a
#: raster copy is sometimes wanted for a slide or a web page.
FIGURE_FORMATS = ("pdf", "svg", "png", "tiff", "jpg")


def _mpl():
    import matplotlib.pyplot as plt
    return plt


def save_figure(fig: Any, fig_dir: str | Path, stem: str,
                formats: tuple[str, ...] | list[str] = (DEFAULT_FIGURE_FORMAT,),
                dpi: int = 300) -> list[str]:
    """Export a Figure in the selected formats and return written paths.

    ``formats`` defaults to :data:`DEFAULT_FIGURE_FORMAT` (PDF). A vector target
    carries no pixel grid, so ``dpi`` only affects raster outputs in the list;
    passing it is harmless and keeps one call signature for every format.
    """
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    allowed = {"png", "svg", "pdf", "tiff", "jpg", "jpeg"}
    selected = tuple(str(x).lower().lstrip(".") for x in formats)
    invalid = set(selected) - allowed
    if invalid:
        raise ValueError(f"unsupported figure format(s): {', '.join(sorted(invalid))}")
    paths = [fig_dir / f"{stem}.{fmt}" for fmt in selected]
    for path in paths:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    _mpl().close(fig)
    return [str(p) for p in paths]


def _save(fig: Any, fig_dir: Path, stem: str,
          formats: tuple[str, ...] = (DEFAULT_FIGURE_FORMAT,)) -> list[str]:
    """Internal alias for save_figure."""
    return save_figure(fig, fig_dir, stem, formats=formats)


def figure_dir(out: str | Path) -> Path:
    """Resolve the shared figures directory for file or directory outputs."""
    path = Path(out)
    if path.suffix.lower() in {".tsv", ".csv", ".json", ".txt"}:
        result = path.parent / "figures"
    else:
        result = path / "figures"
    result.mkdir(parents=True, exist_ok=True)
    return result


def add_figure_args(parser) -> None:
    """Add optional figure export formats; no files are exported by default."""
    parser.add_argument("--figure-format", action="append", dest="figure_formats",
                        choices=FIGURE_FORMATS,
                        default=None,
                        help="Optional export format; repeat to export multiple "
                             f"formats. Defaults to {DEFAULT_FIGURE_FORMAT} (vector) "
                             "when a figure is produced")


def selected_formats(args) -> tuple[str, ...]:
    """Figure formats requested on the command line (empty means no figures)."""
    return tuple(getattr(args, "figure_formats", None) or ())


def _numeric(frame):
    import pandas as pd
    cols = []
    for col in frame.columns:
        values = pd.to_numeric(frame[col], errors="coerce")
        if values.notna().sum() > 0:
            cols.append((col, values))
    return cols


def _read_tables(out_dir: Path):
    import pandas as pd
    tables = []
    for path in sorted(out_dir.rglob("*.tsv")) + sorted(out_dir.rglob("*.csv")):
        try:
            frame = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",")
            if not frame.empty:
                tables.append((path, frame))
        except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
            # bad region or file cannot abort a batch QC run. The caller reports
            # the metric as not_evaluated rather than failing.
            continue
    return tables


def plot_volcano(frame, title: str = "Differential volcano"):
    """Return a volcano Figure while preserving missing statistical values."""
    import numpy as np
    plt = _mpl()
    lfc_col = next((c for c in ("log2fc", "log2_fold_change", "delta_beta", "delta") if c in frame), None)
    p_col = next((c for c in ("padj", "qvalue", "pvalue", "p_value") if c in frame), None)
    if not lfc_col or not p_col:
        return None
    x = pd_numeric(frame[lfc_col])
    y = -np.log10(np.clip(pd_numeric(frame[p_col]), 1e-300, 1))
    sig = pd_numeric(frame[p_col]) < (0.05 if p_col != "pvalue" else 0.05)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x[~sig], y[~sig], s=8, alpha=.45, color="#8c8c8c", label="not significant")
    ax.scatter(x[sig], y[sig], s=10, alpha=.8, color="#d95f02", label=f"{p_col} < 0.05")
    ax.axvline(0, color="black", lw=.6)
    ax.set_xlabel(lfc_col); ax.set_ylabel("-log10(p/q)"); ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return fig


def _volcano(frame, fig_dir: Path, stem: str = "volcano",
             formats: tuple[str, ...] = (DEFAULT_FIGURE_FORMAT,)) -> list[str]:
    fig = plot_volcano(frame)
    return [] if fig is None else _save(fig, fig_dir, stem, formats=formats)


def pd_numeric(values):
    """Coerce a column to numeric, dropping unparseable values."""
    import pandas as pd
    return pd.to_numeric(values, errors="coerce").to_numpy(float)


def plot_heatmap(frame, title: str = ""):
    """Return a heatmap Figure for a numeric matrix."""
    import numpy as np
    plt = _mpl()
    numeric = frame.select_dtypes(include="number")
    if numeric.shape[0] < 2 or numeric.shape[1] < 2:
        return None
    # Large matrices are sampled deterministically for a readable overview.
    data = numeric.iloc[: min(len(numeric), 200), : min(numeric.shape[1], 80)]
    fig, ax = plt.subplots(figsize=(max(5, data.shape[1] * .18), max(4, data.shape[0] * .08)))
    im = ax.imshow(data.to_numpy(float), aspect="auto", interpolation="nearest", cmap="RdBu_r")
    ax.set_xlabel("numeric columns"); ax.set_ylabel("rows")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="value")
    return fig


def _heatmap(frame, fig_dir: Path, stem: str = "sample_heatmap",
             formats: tuple[str, ...] = (DEFAULT_FIGURE_FORMAT,)) -> list[str]:
    fig = plot_heatmap(frame)
    return [] if fig is None else _save(fig, fig_dir, stem, formats=formats)


def plot_bar(frame, title: str = ""):
    """Return a bar chart Figure for a sample metric."""
    plt = _mpl()
    nums = _numeric(frame)
    if not nums:
        return None
    label = next((c for c in ("sample_id", "sample", "motif", "tf", "state", "group") if c in frame), None)
    value_name, values = nums[0]
    labels = frame[label].astype(str).tolist() if label else [str(i + 1) for i in range(len(values))]
    fig, ax = plt.subplots(figsize=(max(5, min(14, len(labels) * .35)), 4))
    ax.bar(range(len(values)), values, color="#2c7fb8")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel(value_name); ax.set_title(title or value_name)
    return fig


def _bar(frame, fig_dir: Path, stem: str = "metrics",
         formats: tuple[str, ...] = (DEFAULT_FIGURE_FORMAT,)) -> list[str]:
    fig = plot_bar(frame)
    return [] if fig is None else _save(fig, fig_dir, stem, formats=formats)


def plot_distribution(frame, title: str = ""):
    """Return a distribution Figure."""
    plt = _mpl()
    nums = _numeric(frame)
    if not nums:
        return None
    value_name, values = nums[0]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist(values.dropna(), bins=40, color="#41ab5d", alpha=.85)
    ax.set_xlabel(value_name); ax.set_ylabel("count"); ax.set_title(title or value_name)
    return fig


def _distribution(frame, fig_dir: Path, stem: str = "distribution",
                  formats: tuple[str, ...] = (DEFAULT_FIGURE_FORMAT,)) -> list[str]:
    fig = plot_distribution(frame)
    return [] if fig is None else _save(fig, fig_dir, stem, formats=formats)


def _motif_json(path: Path, fig_dir: Path,
                formats: tuple[str, ...] = (DEFAULT_FIGURE_FORMAT,)) -> list[str]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
        # bad region or file cannot abort a batch QC run. The caller reports
        # the metric as not_evaluated rather than failing.
        return []
    rows = payload.get("motifs") if isinstance(payload, dict) else None
    if not isinstance(rows, list) or not rows:
        return []
    import pandas as pd
    frame = pd.DataFrame(rows)
    value = "padj" if "padj" in frame else ("odds_ratio" if "odds_ratio" in frame else None)
    if not value:
        return []
    plt = _mpl()
    vals = pd.to_numeric(frame[value], errors="coerce").fillna(0)
    if value == "padj":
        vals = np_log10(vals)
        ylabel = "-log10 adjusted p-value"
    else:
        ylabel = value
    names = frame["motif"].astype(str) if "motif" in frame else frame.index.to_series().astype(str)
    order = vals.argsort().to_numpy()[::-1][:30]
    fig, ax = plt.subplots(figsize=(7, max(3.5, len(order) * .22)))
    ax.barh(range(len(order)), vals.iloc[order].to_numpy()[::-1], color="#756bb1")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(names.iloc[order].to_numpy()[::-1], fontsize=8)
    ax.set_xlabel(ylabel); ax.set_title("Motif enrichment")
    return _save(fig, fig_dir, "motif_enrichment", formats=formats)


def _json_sample_metrics(path: Path, fig_dir: Path,
                         formats: tuple[str, ...]) -> list[str]:
    """ json sample metrics helper."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
        # bad region or file cannot abort a batch QC run. The caller reports
        # the metric as not_evaluated rather than failing.
        return []
    records = payload.get("samples") if isinstance(payload, dict) else None
    if isinstance(records, dict):
        records = [{"sample_id": k, **(v if isinstance(v, dict) else {"value": v})}
                   for k, v in records.items()]
    if not isinstance(records, list) or not records:
        # normalization factors are sample -> factor dictionaries
        if isinstance(payload, dict) and payload and all(isinstance(v, (int, float)) for v in payload.values()):
            import pandas as pd
            return _bar(pd.DataFrame({"sample_id": list(payload), "factor": list(payload.values())}),
                        fig_dir, "normalization_factors", formats=formats)
        return []
    import pandas as pd
    rows = []
    for record in records:
        if not isinstance(record, dict):
            continue
        row = {"sample_id": record.get("sample_id", record.get("name", str(len(rows) + 1)))}
        for key in ("mapping_rate", "total_reads", "n_intervals", "mean_beta", "size_bytes"):
            if key in record and isinstance(record[key], (int, float)):
                row[key] = record[key]
            elif isinstance(record.get("bam_qc"), dict) and key in record["bam_qc"]:
                row[key] = record["bam_qc"][key]
        if len(row) > 1:
            rows.append(row)
    if not rows:
        return []
    frame = pd.DataFrame(rows)
    return _bar(frame, fig_dir, "qc_metrics", formats=formats)


def np_log10(values):
    """log10 transform that keeps zeros finite via a small offset."""
    import numpy as np
    return -np.log10(np.clip(values.to_numpy(float), 1e-300, 1))


def plot_fragment_length(frame, title: str = "Fragment length distribution",
                         raw_fragments=None):
    """Fragment-length histogram with NFR / mononucleosome / dinucleosome bands.

    Needs the per-fragment lengths (not a per-sample median) to draw a real
    distribution; ``raw_fragments`` carries those values when available.
    """
    plt = _mpl()
    import numpy as np
    values = np.asarray([], dtype=float)
    if raw_fragments is not None and len(raw_fragments):
        values = np.asarray(raw_fragments, dtype=float)
    else:
        import pandas as pd
        col = next((c for c in ("fragment_length", "median_fragment_length", "width", "length")
                    if c in frame.columns), None)
        if col is not None and not frame.empty:
            values = pd.to_numeric(frame[col], errors="coerce").dropna().to_numpy(dtype=float)
    if values.size == 0:
        return None
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values, bins=60, color="#2c7fb8", alpha=.85)
    for lo, hi, colour, label in ((0, 100, "#41ab5d", "NFR (<100)"),
                                  (100, 250, "#fec44f", "mono (100-250)"),
                                  (250, 500, "#d95f02", "di (250-500)")):
        ax.axvspan(lo, hi, color=colour, alpha=.12, lw=0, label=label)
    ax.set_xlabel("fragment length (bp)")
    ax.set_ylabel("fragments")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return fig


def plot_frip(frame, title: str = "FRiP per sample"):
    """Bar chart of Fraction of Reads in Peaks."""
    plt = _mpl()
    if "frip" not in frame.columns or frame.empty:
        return None
    import pandas as pd
    values = pd.to_numeric(frame["frip"], errors="coerce")
    if values.notna().sum() == 0:
        return None
    labels = (frame["sample_id"].astype(str) if "sample_id" in frame.columns
              else pd.Series([f"s{i+1}" for i in range(len(frame))]))
    fig, ax = plt.subplots(figsize=(max(4, len(labels) * .6), 3.6))
    ax.bar(range(len(values)), values.fillna(0), color="#d95f02")
    ax.axhline(0.2, color="#666666", ls="--", lw=.8, label="typical good FRiP (0.2)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("FRiP")
    ax.set_ylim(0, max(0.35, float(values.max() or 0) * 1.15))
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return fig


def plot_beta_distribution(beta_values, title: str = "Methylation beta distribution"):
    """Bimodal beta-value histogram with the hypo/hyper regions shaded."""
    plt = _mpl()
    import numpy as np
    values = np.asarray([v for v in beta_values if v is not None], dtype=float)
    values = values[(values >= 0) & (values <= 1)]
    if values.size == 0:
        return None
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values, bins=40, range=(0, 1), color="#9ecae1", edgecolor="white")
    ax.axvspan(0, .2, color="#41ab5d", alpha=.12, lw=0, label="hypo (<0.2)")
    ax.axvspan(.8, 1, color="#d95f02", alpha=.12, lw=0, label="hyper (>0.8)")
    ax.set_xlabel("beta value")
    ax.set_ylabel("CpG sites")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return fig


def plot_tss_enrichment(frame, title: str = "TSS enrichment"):
    """Bar chart of the ENCODE TSS enrichment score with the pass threshold."""
    plt = _mpl()
    if "tss_enrichment" not in frame.columns or frame.empty:
        return None
    import pandas as pd
    values = pd.to_numeric(frame["tss_enrichment"], errors="coerce")
    if values.notna().sum() == 0:
        return None
    labels = (frame["sample_id"].astype(str) if "sample_id" in frame.columns
              else pd.Series([f"s{i+1}" for i in range(len(frame))]))
    fig, ax = plt.subplots(figsize=(max(4, len(labels) * .6), 3.6))
    ax.bar(range(len(values)), values.fillna(0), color="#41ab5d")
    ax.axhline(2.0, color="#666666", ls="--", lw=.8, label="pass threshold (2.0)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("TSS enrichment")
    ax.set_ylim(0, max(3.0, float(values.max() or 0) * 1.2))
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return fig


def plot_read_quality(frame, title: str = "Duplication and mitochondrial fraction"):
    """Grouped bars for duplicate rate and mitochondrial fraction.

    Both are common failure signals (over-amplified library / lysed nuclei) and
    belong side by side on the QC sheet.
    """
    plt = _mpl()
    cols = [c for c in ("duplicate_rate", "mitochondrial_fraction") if c in frame.columns]
    if not cols or frame.empty:
        return None
    import numpy as np
    import pandas as pd
    labels = (frame["sample_id"].astype(str) if "sample_id" in frame.columns
              else pd.Series([f"s{i+1}" for i in range(len(frame))]))
    fig, ax = plt.subplots(figsize=(max(4.5, len(labels) * .8), 3.6))
    width = 0.8 / max(1, len(cols))
    x = np.arange(len(labels))
    for k, col in enumerate(cols):
        vals = pd.to_numeric(frame[col], errors="coerce").fillna(0).to_numpy()
        ax.bar(x + k * width, vals, width=width, label=col.replace("_", " "))
    ax.set_xticks(x + width * (len(cols) - 1) / 2)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("fraction of reads")
    ax.set_ylim(0, 1.0)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return fig


def plot_peak_read_distribution(out_dir, title: str = "Reads per peak"):
    """Histogram of log10(reads+1) per peak, showing how evenly peaks are supported.

    A single FRiP number hides whether signal is spread across peaks or carried by
    a handful. The x axis is log-scale because per-peak counts are heavy-tailed.
    """
    plt = _mpl()
    import numpy as np
    payload = None
    for path in sorted(Path(out_dir).rglob("*_qc.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
            # bad region or file cannot abort a batch QC run. The caller reports
            # the metric as not_evaluated rather than failing.
            continue
        for record in (data.get("samples") or []) if isinstance(data, dict) else []:
            specific = (record or {}).get("assay_qc") or {}
            prd = specific.get("peak_reads") or {}
            if prd.get("hist_counts") and prd.get("hist_edges"):
                payload = prd
                break
        if payload:
            break
    if not payload:
        return None
    counts = np.asarray(payload["hist_counts"], dtype=float)
    edges = np.asarray(payload["hist_edges"], dtype=float)
    if counts.size == 0 or edges.size != counts.size + 1:
        return None
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(edges[:-1], counts, width=np.diff(edges), align="edge",
           color="#2c7fb8", edgecolor="white")
    zero_frac = payload.get("fraction_peaks_zero")
    subtitle = f"peaks with 0 reads: {zero_frac:.1%}" if isinstance(zero_frac, float) else ""
    ax.set_xlabel("log10(1 + reads per peak)")
    ax.set_ylabel("peaks")
    ax.set_title(f"{title}\n{subtitle}".strip())
    return fig


def plot_blacklist(frame, title: str = "Blacklist fraction"):
    """Bar chart of the fraction of reads in blacklist regions."""
    plt = _mpl()
    if "blacklist_fraction" not in frame.columns or frame.empty:
        return None
    import pandas as pd
    values = pd.to_numeric(frame["blacklist_fraction"], errors="coerce")
    if values.notna().sum() == 0:
        return None
    labels = (frame["sample_id"].astype(str) if "sample_id" in frame.columns
              else pd.Series([f"s{i+1}" for i in range(len(frame))]))
    fig, ax = plt.subplots(figsize=(max(4, len(labels) * .6), 3.6))
    ax.bar(range(len(values)), values.fillna(0), color="#756bb1")
    ax.axhline(0.01, color="#666666", ls="--", lw=.8, label="flag threshold (1%)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("fraction of reads")
    ax.set_ylim(0, max(0.02, float(values.max() or 0) * 1.2))
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    return fig


def _assay_metric_frame(out_dir: Path):
    """Flatten the assay_qc JSON payloads into one tidy frame.

    Returns ``(frame, beta_lists)`` where ``beta_lists`` carries any stored beta
    vectors, which the JSON keeps as nested ``bin_fractions`` rather than raw
    values.
    """
    import pandas as pd
    rows: list[dict] = []
    beta_lists: list[list[float]] = []
    fragment_samples: list[list[float]] = []
    for path in sorted(out_dir.rglob("*_qc.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
            # bad region or file cannot abort a batch QC run. The caller reports
            # the metric as not_evaluated rather than failing.
            continue
        for record in payload.get("samples", []) if isinstance(payload, dict) else []:
            if not isinstance(record, dict):
                continue
            specific = record.get("assay_qc") or {}
            row = {"sample_id": record.get("sample_id"), "assay": record.get("assay")}
            for key, value in specific.items():
                if isinstance(value, (int, float, str, bool)):
                    row[key] = value
                elif isinstance(value, list) and key == "fragment_length_sample":
                    fragment_samples.append([float(v) for v in value if isinstance(v, (int, float))])
                elif isinstance(value, dict) and key == "beta_detail":
                    bins = value.get("bin_fractions")
                    if isinstance(bins, list) and bins:
                        # Reconstruct a coarse beta sample from the histogram so the
                        # diagnostic figure can be drawn without the raw table.
                        expanded: list[float] = []
                        for i, frac in enumerate(bins):
                            centre = (i + 0.5) / len(bins)
                            expanded.extend([centre] * max(0, int(round(frac * 1000))))
                        if expanded:
                            beta_lists.append(expanded)
            if len(row) > 2:
                rows.append(row)
    return (pd.DataFrame(rows) if rows else None), beta_lists, fragment_samples


_LOGGER = logging.getLogger(__name__)


def render_module_figures(module: str, out: str | Path,
                          formats: tuple[str, ...] | list[str] = (DEFAULT_FIGURE_FORMAT,)) -> list[str]:
    """Render requested formats from completed module output tables."""
    if not formats:
        return []
    try:
        out_path = Path(out)
        out_dir = out_path.parent if out_path.suffix.lower() in {".tsv", ".csv", ".json", ".txt"} else out_path
        fig_dir = figure_dir(out)
        module_key = module.lower().replace("_", "-")
        tables = _read_tables(out_dir)
        made: list[str] = []
        # motif enrichment/scan may intentionally write JSON only.
        for path in sorted(out_dir.rglob("*.json")):
            if "motif" in path.name.lower() or "scan" in path.name.lower():
                made.extend(_motif_json(path, fig_dir, formats=tuple(formats)))
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                    scores = payload.get("scores") if isinstance(payload, dict) else None
                    if isinstance(scores, list) and scores:
                        import pandas as pd
                        made.extend(_distribution(pd.DataFrame({"motif_score": scores}), fig_dir,
                                                  "motif_score_distribution", formats=tuple(formats)))
                except Exception as exc:  # noqa: BLE001
                    # One unreadable/unusable table must not stop the other
                    # figures; note the reason instead of discarding it.
                    _LOGGER.debug("figure input skipped: %s: %s",
                                       type(exc).__name__, exc)
                    pass
            elif "qc" in module_key or "report" in module_key or "normal" in module_key:
                made.extend(_json_sample_metrics(path, fig_dir, formats=tuple(formats)))
        # Peak callers/consensus modules commonly write BED and a JSON manifest.
        bed_rows = []
        for bed in sorted(out_dir.rglob("*.bed")) + sorted(out_dir.rglob("*.narrowPeak")):
            try:
                with bed.open(encoding="utf-8", errors="replace") as handle:
                    for line in handle:
                        if line.strip() and not line.startswith(("#", "track", "browser")):
                            fields = line.rstrip("\n").split("\t")
                            if len(fields) >= 3:
                                start, end = int(fields[1]), int(fields[2])
                                if end > start:
                                    bed_rows.append({"width": end - start, "chrom": fields[0]})
            except Exception as exc:  # noqa: BLE001
                # One unreadable/unusable table must not stop the other
                # figures; note the reason instead of discarding it.
                _LOGGER.debug("figure input skipped: %s: %s",
                                   type(exc).__name__, exc)
                continue
        if bed_rows and any(token in module_key for token in ("peak", "consensus", "annotation", "qc")):
            import pandas as pd
            made.extend(_distribution(pd.DataFrame(bed_rows)[["width"]], fig_dir,
                                      "peak_width_distribution", formats=tuple(formats)))
        if bed_rows and ("tn5" in module_key or "shift" in module_key):
            import pandas as pd
            chrom_counts = pd.DataFrame(bed_rows)["chrom"].value_counts().rename_axis("chrom").reset_index(name="insertions")
            made.extend(_bar(chrom_counts, fig_dir, "insertion_counts", formats=tuple(formats)))
        for path in sorted(out_dir.glob("*.json")):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                records = payload.get("samples") if isinstance(payload, dict) else None
                if isinstance(records, list) and records and any("n_peaks" in x for x in records if isinstance(x, dict)):
                    import pandas as pd
                    made.extend(_bar(pd.DataFrame(records), fig_dir, "peak_counts", formats=tuple(formats)))
            except Exception as exc:  # noqa: BLE001
                # One unreadable/unusable table must not stop the other
                # figures; note the reason instead of discarding it.
                _LOGGER.debug("figure input skipped: %s: %s",
                                   type(exc).__name__, exc)
                continue
        diff = next(((p, f) for p, f in tables if {"pvalue", "padj"}.intersection(f.columns)
                     and ({"log2fc", "delta_beta", "delta"}.intersection(f.columns))), None)
        if diff:
            made.extend(_volcano(diff[1], fig_dir, "volcano", formats=tuple(formats)))
        # Count matrices and TF/sample activity use a heatmap as the primary view.
        for path, frame in tables:
            lower = path.name.lower()
            if any(token in module_key for token in ("count", "activity", "integrat", "normalize")) or "matrix" in lower:
                made.extend(_heatmap(frame, fig_dir, "sample_heatmap", formats=tuple(formats)))
                break
        # Assay-specific QC figures: the generic bar chart cannot show fragment
        # periodicity, FRiP quality or methylation bimodality.
        if "qc" in module_key:
            assay_frame, beta_lists, fragment_samples = _assay_metric_frame(out_dir)
            if assay_frame is not None:
                if "median_fragment_length" in assay_frame.columns:
                    fig = plot_fragment_length(assay_frame,
                                               raw_fragments=(fragment_samples[0] if fragment_samples else None))
                    if fig:
                        made.extend(_save(fig, fig_dir, "fragment_length_distribution",
                                          formats=tuple(formats)))
                if "frip" in assay_frame.columns:
                    fig = plot_frip(assay_frame)
                    if fig:
                        made.extend(_save(fig, fig_dir, "frip", formats=tuple(formats)))
                if "signal_enrichment" in assay_frame.columns:
                    renamed = assay_frame.rename(columns={"signal_enrichment": "peak_vs_flank"})
                    fig = plot_bar(renamed, title="Peak signal enrichment")
                    if fig:
                        made.extend(_save(fig, fig_dir, "signal_enrichment",
                                          formats=tuple(formats)))
                if "tss_enrichment" in assay_frame.columns:
                    fig = plot_tss_enrichment(assay_frame)
                    if fig:
                        made.extend(_save(fig, fig_dir, "tss_enrichment",
                                          formats=tuple(formats)))
                if {"duplicate_rate", "mitochondrial_fraction"} & set(assay_frame.columns):
                    fig = plot_read_quality(assay_frame)
                    if fig:
                        made.extend(_save(fig, fig_dir, "read_quality",
                                          formats=tuple(formats)))
            if "blacklist_fraction" in getattr(assay_frame, "columns", []):
                fig = plot_blacklist(assay_frame)
                if fig:
                    made.extend(_save(fig, fig_dir, "blacklist_fraction",
                                      formats=tuple(formats)))
            fig = plot_peak_read_distribution(out_dir)
            if fig:
                made.extend(_save(fig, fig_dir, "peak_read_distribution",
                                  formats=tuple(formats)))
            for values in beta_lists:
                fig = plot_beta_distribution(values)
                if fig:
                    made.extend(_save(fig, fig_dir, "beta_distribution",
                                      formats=tuple(formats)))
        if "qc" in module_key or "report" in module_key:
            qc = next(((p, f) for p, f in tables if any(c in f.columns for c in ("sample_id", "mapping_rate", "n_intervals", "mean_beta"))), None)
            if qc:
                made.extend(_bar(qc[1], fig_dir, "qc_metrics", formats=tuple(formats)))
        if "methylation" in module_key:
            profile = next(((p, f) for p, f in tables if "beta" in f.columns), None)
            if profile:
                made.extend(_distribution(profile[1], fig_dir, "beta_distribution", formats=tuple(formats)))
            annotated = next(((p, f) for p, f in tables if "signed_distance_to_tss" in f.columns), None)
            if annotated:
                import pandas as pd
                made.extend(_distribution(annotated[1][["signed_distance_to_tss"]], fig_dir,
                                          "tss_distance_distribution", formats=tuple(formats)))
        if "chromhmm" in module_key or "histone" in module_key:
            state = next(((p, f) for p, f in tables if any(str(c).lower().endswith("state") or str(c).lower() == "state" for c in f.columns)), None)
            if state:
                import pandas as pd
                col = next(c for c in state[1].columns if str(c).lower() == "state" or str(c).lower().endswith("state"))
                counts = state[1][col].astype(str).value_counts().rename_axis("state").reset_index(name="n_regions")
                made.extend(_bar(counts, fig_dir, "state_coverage", formats=tuple(formats)))
        if "mnase" in module_key or "nucleosome" in module_key:
            for path in sorted(out_dir.rglob("*.json")):
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                    values = payload.get("occupancy") or payload.get("profile") or payload.get("signal")
                    if isinstance(values, list) and values:
                        plt = _mpl(); fig, ax = plt.subplots(figsize=(6, 3.5)); ax.plot(values, color="#1f78b4")
                        ax.set_xlabel("position/bin"); ax.set_ylabel("signal"); made.extend(_save(fig, fig_dir, "nucleosome_profile", formats=tuple(formats))); break
                except Exception as exc:  # noqa: BLE001
                    # One unreadable/unusable table must not stop the other
                    # figures; note the reason instead of discarding it.
                    _LOGGER.debug("figure input skipped: %s: %s",
                                       type(exc).__name__, exc)
                    continue
        if "peak" in module_key or "consensus" in module_key or "annotation" in module_key:
            peak = next(((p, f) for p, f in tables if any(c in f.columns for c in ("start", "end", "n_intervals", "width"))), None)
            if peak:
                frame = peak[1].copy()
                if "start" in frame and "end" in frame:
                    import pandas as pd
                    frame["width"] = pd.to_numeric(frame["end"], errors="coerce") - pd.to_numeric(frame["start"], errors="coerce")
                made.extend(_distribution(frame[["width"]] if "width" in frame else frame, fig_dir, "peak_width_distribution", formats=tuple(formats)))
        # Every successful tabular module gets at least one useful overview figure.
        if not made and bed_rows:
            import pandas as pd
            made.extend(_distribution(pd.DataFrame(bed_rows)[["width"]], fig_dir,
                                      "interval_width_distribution", formats=tuple(formats)))
        if not made and tables:
            made.extend(_bar(tables[0][1], fig_dir, "module_overview", formats=tuple(formats)))
        return sorted(set(made))
    except Exception as exc:  # noqa: BLE001
        # Figure rendering is best-effort; the analysis itself
        # already succeeded, so record why and return what we have.
        _LOGGER.warning("figure rendering skipped: %s: %s",
                            type(exc).__name__, exc)
        return []
