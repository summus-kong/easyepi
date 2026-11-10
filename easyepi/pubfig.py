"""Easyepi publication figure helpers (Nature-style, Python/matplotlib).

Wraps the publication-figure rules (editable SVG text, restrained palette,
export to SVG/PDF/TIFF) into one reusable set of helpers so every Easyepi plot
uses the same journal-quality conventions. Based on the configured figure style.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np

# Nature-style palette (semantic: neutral family, signal family, accent family).
PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_3": "#8BCF8B",
    "green_2": "#AADCA9",
    "red_strong": "#B64342",
    "red_soft": "#E9A6A1",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
}
# Default ordered set for series/panels.
DEFAULT_COLORS = [PALETTE["blue_main"], PALETTE["green_3"], PALETTE["red_strong"],
                  PALETTE["teal"], PALETTE["violet"], PALETTE["neutral_mid"]]


def apply_publication_style(font_size: float = 7, axes_linewidth: float = 0.8,
                            use_tex: bool = False) -> None:
    """Apply Nature-style rcParams. Call once before creating figures."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
    plt.rcParams['svg.fonttype'] = 'none'   # editable text in SVG
    plt.rcParams['pdf.fonttype'] = 42       # editable TrueType text in PDF
    plt.rcParams['font.size'] = font_size
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.linewidth'] = axes_linewidth
    plt.rcParams['legend.frameon'] = False
    if use_tex:
        plt.rcParams['text.usetex'] = True


def finalize_figure(fig, out_path: str, formats: Optional[List[str]] = None,
                    dpi: int = 300, pad: float = 1.5,
                    bbox_inches: str = "tight", close: bool = True) -> List[str]:
    """Tight-layout and export a figure to SVG/PDF/TIFF/PNG.

    ``formats`` defaults to the extension of ``out_path``, else ['pdf','svg','png'].
    Returns the list of saved paths.
    """
    import os
    from pathlib import Path
    fig.tight_layout(pad=pad)
    base = Path(out_path)
    if formats is None:
        formats = [base.suffix.lstrip('.') or "png"]
        base = base.with_suffix("")
    os.makedirs(base.parent, exist_ok=True)
    saved = []
    for fmt in formats:
        p = str(base) + f".{fmt}"
        kw = {}
        if bbox_inches is not None:
            kw["bbox_inches"] = bbox_inches
        if fmt == "svg":
            kw["edgecolor"] = "none"  # avoid svg artifacts
        fig.savefig(p, dpi=dpi, **kw)
        saved.append(p)
    if close:
        import matplotlib.pyplot as plt
        plt.close(fig)
    return saved


def add_panel_label(ax, label: str, x: float = -0.06, y: float = 1.02,
                    fontsize: float = 8, fontweight: str = "bold") -> None:
    """Nature-style bold panel label near the top-left edge."""
    ax.text(x, y, label, transform=ax.transAxes, fontsize=fontsize,
            fontweight=fontweight, ha="left", va="bottom")


def make_heatmap(ax, matrix, x_labels=None, y_labels=None, cmap: str = "RdBu_r",
                 cbar_label: Optional[str] = None, annotate: bool = False,
                 fmt: str = "{:.2f}", fontsize: float = 6, vmin=None, vmax=None):
    """Journal-style heatmap with optional colorbar + annotated cells."""
    import matplotlib as mpl
    im = ax.imshow(np.asarray(matrix), cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    if cbar_label:
        cbar = ax.figure.colorbar(im, ax=ax)
        cbar.set_label(cbar_label)
    if x_labels:
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=6)
    if y_labels:
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=6)
    if annotate:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cm_obj = mpl.colormaps.get_cmap(cmap)
        arr = np.asarray(matrix)
        for (i, j), val in np.ndenumerate(arr):
            col = cm_obj(norm(val)) if vmin is not None else None
            lum = 0.299*col[0]+0.587*col[1]+0.114*col[2] if col else 0.5
            ax.text(j, i, fmt.format(val), ha="center", va="center",
                    fontsize=fontsize, color="white" if lum < 0.5 else "black")
    ax.set_frame_on(False)


def make_trend(ax, x, y_series, labels, colors=None, ylabel: Optional[str] = None,
               xlabel: Optional[str] = None, lw: float = 1.2, marker: str = "o",
               markersize: float = 3.5, show_shadow: bool = False):
    """Multi-line trend plot with optional ± std band."""
    if colors is None:
        colors = DEFAULT_COLORS
    for y, label, color in zip(y_series, labels, colors):
        y = np.asarray(y)
        if y.ndim == 2:
            mean, std = y.mean(0), y.std(0)
        else:
            mean, std = y, None
        ax.plot(x, mean, color=color, lw=lw, marker=marker, markersize=markersize, label=label)
        if show_shadow and std is not None:
            ax.fill_between(x, mean - std, mean + std, color=color, alpha=0.15)
    if ylabel:
        ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.legend(frameon=False, fontsize=6)


def make_grouped_bar(ax, categories, series, labels, colors=None,
                     ylabel: Optional[str] = None, bar_width: float = 0.8,
                     annotate: bool = False):
    """Grouped bar chart (journal style)."""
    if colors is None:
        colors = DEFAULT_COLORS
    n_groups = len(series)
    n_cats = len(categories)
    w = bar_width / n_groups
    x = np.arange(n_cats)
    for i, (vals, label, color) in enumerate(zip(series, labels, colors)):
        offset = (i - (n_groups - 1) / 2) * w
        bars = ax.bar(x + offset, np.asarray(vals), width=w, label=label,
                      color=color, edgecolor="black", linewidth=0.6)
        if annotate:
            for bar, val in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()+0.01,
                        f"{val:.2f}", ha="center", va="bottom", fontsize=5)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=6)
    ax.set_ylabel(ylabel or "")
    ax.legend(frameon=False, fontsize=6)
