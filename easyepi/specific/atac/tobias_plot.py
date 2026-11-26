"""1:1 reimplementation of TOBIAS ``PlotAggregate`` and ``PlotHeatmap``.

Both commands do the same three things to a set of TFBS/BED regions before any
matplotlib call, and it is that shared preparation which determines the numbers:

    1. centre each region on its midpoint and give it a fixed width
       (``set_width(flank * 2)``, so the output spans ``flank`` either side)
    2. drop any region whose extended window leaves the chromosome
    3. read the bigWig signal over the window

``PlotAggregate`` then averages the resulting matrix per signal/region pair, while
``PlotHeatmap`` keeps the rows and optionally sorts them. Reference:
``tobias/tools/plot_aggregate.py`` and ``tobias/tools/plot_heatmap.py``.

Outlier removal is the one step that changes the aggregate value, and the
reference applies it per row on the row maximum: a region is dropped when its
maximum exceeds the ``100 * remove_outliers`` percentile of all row maxima
(``--remove-outliers`` is a fraction, default 0.1, i.e. the top 10% by peak).
"""
from __future__ import annotations

import numpy as np

#: Percentile of regions to KEEP, per the CLI's own description
#: ("Value between 0-1 indicating the percentile of regions to include, e.g. 0.99
#: to remove the sites with 1% highest values", default 1). 1.0 keeps everything.
DEFAULT_REMOVE_OUTLIERS = 1.0

#: Default flank in bp; the window is ``flank * 2`` wide, centred on the midpoint.
DEFAULT_FLANK = 1000


def window_around(start: float, end: float, flank: int = DEFAULT_FLANK):
    """Return ``(win_start, win_end, midpoint)`` for one region.

    The reference uses ``OneRegion.set_width(width)``, which recentres on the
    region midpoint and forces the span to ``flank * 2``.
    """
    midpoint = (start + end) / 2.0
    width = flank * 2
    win_start = int(midpoint - width / 2.0)
    win_end = win_start + width
    return win_start, win_end, midpoint


def extract_window(signal, chrom_len: int, start: float, end: float,
                   flank: int = DEFAULT_FLANK):
    """Signal over the centred window, or ``None`` if it leaves the chromosome.

    ``signal`` is a 1-D array already covering ``[0, chrom_len)``, as a bigWig
    reader would return for the whole chromosome. Returning ``None`` mirrors
    ``check_boundary(action="remove")``.
    """
    win_start, win_end, _ = window_around(start, end, flank)
    if win_start < 0 or win_end > chrom_len:
        return None
    return np.asarray(signal[win_start:win_end], dtype=float)


def build_matrix(signal, chrom_len: int, regions, flank: int = DEFAULT_FLANK):
    """Matrix of shape ``(n_regions, flank*2)`` plus the indices kept.

    Regions whose window leaves the chromosome are dropped, matching the
    reference's boundary removal before reading signal.
    """
    rows, kept = [], []
    for i, (start, end) in enumerate(regions):
        win = extract_window(signal, chrom_len, start, end, flank)
        if win is None:
            continue
        rows.append(win)
        kept.append(i)
    if not rows:
        return np.zeros((0, flank * 2)), []
    return np.vstack(rows), kept


def remove_outlier_rows(signalmat, remove_outliers: float = DEFAULT_REMOVE_OUTLIERS):
    """Drop regions whose peak exceeds the chosen row-maximum percentile.

    ``remove_outliers`` is the percentile to INCLUDE, not the fraction to remove.
    """
    mat = np.asarray(signalmat, dtype=float)
    if mat.shape[0] == 0:
        return mat
    return mat[outlier_mask(mat, remove_outliers)]


def log_transform(signalmat):
    """``log2(|x| + 1)`` with the sign of the original values restored."""
    signalmat = np.asarray(signalmat, dtype=float)
    log = np.log2(np.abs(signalmat) + 1)
    log[signalmat < 0] *= -1
    return log


def aggregate(signalmat, remove_outliers: float = DEFAULT_REMOVE_OUTLIERS,
              do_log: bool = False):
    """Column-wise mean of the filtered matrix, as ``PlotAggregate`` computes it."""
    mat = remove_outlier_rows(signalmat, remove_outliers)
    if mat.shape[0] == 0:
        return np.zeros(signalmat.shape[1] if signalmat.ndim == 2 else 0)
    if do_log:
        mat = log_transform(mat)
    with np.errstate(invalid="ignore"):
        return np.nanmean(mat, axis=0)


def outlier_mask(signalmat, remove_outliers: float = DEFAULT_REMOVE_OUTLIERS):
    """Boolean mask of rows to KEEP, by the row-maximum percentile rule.

    ``remove_outliers`` is the percentile to include (1.0 keeps all, 0.99 removes the
    highest 1%), matching the CLI help and the reference's
    ``np.percentile(row_maxima, 100 * remove_outliers)`` followed by ``<=``.
    """
    mat = np.asarray(signalmat, dtype=float)
    if mat.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    if remove_outliers is None or remove_outliers >= 1:
        return np.ones(mat.shape[0], dtype=bool)
    max_values = np.max(mat, axis=1)
    upper_limit = np.percentile(max_values, [100 * remove_outliers])[0]
    return max_values <= upper_limit


def sort_rows(signalmat, key_values, remove_outliers: float = DEFAULT_REMOVE_OUTLIERS):
    """Filter as ``aggregate`` does, then sort rows by ``key_values`` descending.

    ``PlotHeatmap`` sorts with ``reverse=True`` on the requested column, so the
    largest key comes first.
    """
    mat = np.asarray(signalmat, dtype=float)
    keys = np.asarray(key_values, dtype=float)
    if mat.shape[0] != keys.shape[0]:
        raise ValueError("signalmat rows (%d) != key_values (%d)"
                         % (mat.shape[0], keys.shape[0]))
    mask = outlier_mask(mat, remove_outliers)
    mat, keys = mat[mask], keys[mask]
    if mat.shape[0] == 0:
        return mat, keys
    order = np.argsort(keys)[::-1]
    return mat[order], keys[order]


def normalize_rows_for_display(signalmat, vmin: float = None, vmax: float = None):
    """Clip to ``[vmin, vmax]``; defaults come from the data, as the heatmap does."""
    mat = np.asarray(signalmat, dtype=float)
    lo = np.nanmin(mat) if vmin is None else vmin
    hi = np.nanmax(mat) if vmax is None else vmax
    if hi <= lo:
        hi = lo + 1e-9
    return np.clip(mat, lo, hi)
