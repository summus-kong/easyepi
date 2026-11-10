"""Motif figure helpers: sequence logo + enrichment bar chart.

Renders a Position-Weight-Matrix as a sequence logo and a known-motif
enrichment bar chart, ready to embed in the Easyepi report.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def sequence_logo(pwm: np.ndarray, title: str = "", save: Optional[str] = None):
    """Render a sequence logo from a (4, L) PWM (rows ACGT) as a matplotlib Figure.

    Letter heights are scaled by base probability (information-content-weighted
    heights), the standard sequence-logo convention.
    """
    pwm = np.asarray(pwm, dtype=np.float64)
    if pwm.shape[0] != 4:
        raise ValueError("pwm must be (4, L), rows ACGT")
    freq = pwm / pwm.sum(axis=0, keepdims=True) if pwm.sum() else pwm
    L = pwm.shape[1]
    # information content
    freq_safe = np.clip(freq, 1e-9, 1.0)
    ic = np.sum(freq_safe * np.log2(freq_safe / 0.25), axis=0)
    ic = np.clip(ic, 0, 2.0)
    plt = _plt()
    fig, ax = plt.subplots(figsize=(0.45 * L + 2, 3.2))
    bases = "ACGT"
    colors = {"A": "#2E86AB", "C": "#A23B72", "G": "#F18F01", "T": "#C73E1D"}
    for j in range(L):
        # order letters by frequency (descending)
        order = np.argsort(freq[:, j])[::-1]
        y = 0.0
        for bi in order:
            h = freq[bi, j] * ic[j]
            if h > 0:
                ax.text(j, y + h / 2, bases[bi], ha="center", va="center",
                        fontsize=14 if h > 0.4 else 9,
                        color=colors[bases[bi]],
                        fontweight="bold")
                y += h
    ax.set_xlim(-0.5, L - 0.5)
    ax.set_ylim(0, 2.05)
    ax.set_xticks(range(L))
    ax.set_title(title)
    ax.set_ylabel("bits")
    fig.tight_layout()
    if save:
        Path(save).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save, dpi=300)
    return fig


def enrichment_bar(rows: List[Dict], top: int = 10, save: Optional[str] = None):
    """Render a known-motif enrichment bar chart.

    Filters out sequence-bias / low-information entries (names starting with
    'SeqBias', or generic low-complexity tags) and sorts by enrichment
    significance (``log_p``, i.e. -log10 p) so genuine TF motifs rise to the
    top. Falls back to ``log_odds`` when ``log_p`` is unavailable.
    """
    plt = _plt()
    # Drop sequence-bias / non-TF entries.
    cleaned = [r for r in rows
               if not str(r.get("name", "")).startswith("SeqBias")
               and not str(r.get("name", "")).startswith("AT")]  # orphan/plant IDs are noise here
    # Sort by enrichment significance: highest -log_p first.
    def _sig(r):
        lp = r.get("log_p")
        lo = r.get("log_odds")
        if lp is not None:
            return float(lp) if float(lp) >= 0 else float(-lp)  # -log10 p is positive
        return float(lo) if lo is not None else 0.0
    cleaned.sort(key=lambda r: -_sig(r))
    rows_sorted = cleaned[:top]
    if not rows_sorted:
        rows_sorted = sorted(rows, key=lambda r: r.get("log_odds", 0), reverse=True)[:top]
        vals = [r.get("log_odds", 0) for r in rows_sorted][::-1]
        xlabel = "log2 enrichment (target vs background)"
    else:
        vals = [abs(_sig(r)) for r in rows_sorted][::-1]
        xlabel = "-log10 p (enrichment significance)"
    names = [r.get("name", "")[:26] for r in rows_sorted][::-1]
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(names) + 2))
    ax.barh(range(len(names)), vals, color="#1f78b4")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel(xlabel)
    ax.set_title("Known motif enrichment")
    fig.tight_layout()
    if save:
        Path(save).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save, dpi=300)
    return fig
