"""Shared public plotting API. Functions return Matplotlib Figures; callers choose whether and how to export them."""
from __future__ import annotations

from easyepi._plots import (figure_dir, plot_bar, plot_distribution,
                                   plot_heatmap, plot_volcano, render_module_figures,
                                   save_figure)

__all__ = [
    "figure_dir", "plot_bar", "plot_distribution", "plot_heatmap",
    "plot_volcano", "render_module_figures", "save_figure",
]

from easyepi.viz import plot_matrix_heatmap, plot_profile  # noqa: F401
__all__ += ["plot_profile", "plot_matrix_heatmap"]
