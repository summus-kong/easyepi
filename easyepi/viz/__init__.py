"""Visualization (native matplotlib/seaborn).

Shared signal profiles and heatmaps. The plotting functions are
dependency-free beyond the matplotlib hard requirement.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _plt():
    import matplotlib.pyplot as plt
    return plt


def _save_path(fig, save) -> None:
    """Write ``fig`` to ``save``.

    A path with a suffix is honoured exactly (``plot.png`` writes a PNG); a path
    without one gets :data:`easyepi._plots.DEFAULT_FIGURE_FORMAT`, so the caller
    does not have to know the default to get vector output.
    """
    from easyepi._plots import DEFAULT_FIGURE_FORMAT

    path = Path(save)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.suffix:
        path = path.with_suffix("." + DEFAULT_FIGURE_FORMAT)
    fig.savefig(path, dpi=150)


def plot_profile(signal: np.ndarray, title: str = "", save: Path = None,
                 xlabel: str = "position (bp)", ylabel: str = "signal"):
    """Plot a 1-D signal profile and return the matplotlib Figure."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(signal, lw=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    if save:
        _save_path(fig, save)
    return fig


def plot_matrix_heatmap(matrix: np.ndarray, title: str = "", save: Path = None,
                        xlabel: str = "title", cmap: str = "viridis"):
    """Plot a 2-D array as a heatmap and return the matplotlib Figure.

    Named ``plot_matrix_heatmap`` rather than ``plot_heatmap`` so it does not
    collide with the table-based :func:`easyepi.plot.plot_heatmap`.
    """
    import pandas as pd
    from easyepi._plots import plot_heatmap as shared_heatmap
    fig = shared_heatmap(pd.DataFrame(np.asarray(matrix)),title=title)
    if fig is None:
        raise ValueError("Heatmap requires at least two rows and two columns")
    ax = fig.axes[0]
    ax.set_xlabel(xlabel)
    ax.images[0].set_cmap(cmap)
    if save:
        _save_path(fig, save)
    return fig


def make_profile(signal, **kwargs):
    """Alias for :func:`plot_profile` returning a Figure (no save)."""
    return plot_profile(np.asarray(signal), **{k: v for k, v in kwargs.items() if k != "save"})


def make_heatmap(matrix, **kwargs):
    """Alias for :func:`plot_matrix_heatmap` returning a Figure (no save)."""
    return plot_matrix_heatmap(np.asarray(matrix), **{k: v for k, v in kwargs.items() if k != "save"})


def plot_heatmap(matrix: np.ndarray, *args, **kwargs):
    """Deprecated: use :func:`plot_matrix_heatmap`.

    Retained because this name previously resolved through ``easyepi`` to the
    array-based function while ``easyepi.plot.plot_heatmap`` was table-based,
    which made the same name mean two different things.
    """
    import warnings
    warnings.warn(
        "easyepi.viz.plot_heatmap is deprecated; use plot_matrix_heatmap",
        DeprecationWarning, stacklevel=2)
    return plot_matrix_heatmap(matrix, *args, **kwargs)


def run_viz(args) -> int:
    """Render the requested visualisation from prepared tables."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    regions = Path(args.regions)
    if not regions.exists():
        raise FileNotFoundError(f"regions input not found: {regions}")
    rows = []
    with regions.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                try:
                    rows.append((parts[0], int(parts[1]), int(parts[2]), float(parts[3])))
                except ValueError:
                    continue
    if not rows:
        raise ValueError("regions input must contain chrom, start, end and numeric signal columns")
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    sig = np.asarray([r[3] for r in rows], dtype=float)
    from easyepi._plots import selected_formats, save_figure
    formats = selected_formats(args)
    if not formats:
        # Only reachable with --no-figures: `viz` exists to draw, so declining to
        # draw is a contradiction the user should hear about, not a silent no-op.
        raise ValueError("viz writes figures only; drop --no-figures "
                         "(or call easyepi.plot_profile to get a Figure back)")
    fig = plot_profile(sig, title="Easyepi input signal profile", xlabel="region order")
    save_figure(fig, out, "profile", formats=formats)
    print(f"[easyepi] visualization written to {out}")
    return 0
