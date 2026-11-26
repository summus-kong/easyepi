"""1:1 reimplementation of TOBIAS BINDetect's normalisation layer.

BINDetect compares footprint scores between conditions. Raw scores are not
comparable across samples, so TOBIAS normalises each condition's background score
distribution onto a common one before computing log2 fold changes.

Reference: `tobias/tools/bindetect_functions.py`

    sigmoid(x, a, b, L, shift) = L / (1 + exp(-b*(x-a))) + shift
    ArrayNorm.normalize(arr):
        cap arr to [value_min, value_max]
        factor = sigmoid(arr_capped, *popt)   or  constant popt
        return arr * factor

`quantile_normalization` builds each condition's factor curve from the quantiles of
its positive background scores, averaged across conditions:

    quantiles = linspace(0.05, 0.99, 1000)
    xdata_i   = quantile(arr_i[arr_i > 0], quantiles)
    ydata     = mean_i(xdata_i) / xdata_i          # the correction needed
    ydata     = centred rolling mean(ydata, 50)    # smoothed, edge-padded
    y_inter   = interp1d(xdata, ydata_smooth, linear, extrapolate)(linspace(min, max, 1000))
    fit sigmoid with bounds ((min(x), -inf, 0, 0), (max(x), inf, L_range, inf))
        where L_range = diff(percentile(y_inter, [5, 95]))
    on failure: constant factor = mean(y_inter)
    ArrayNorm(func, popt, value_min=min(xdata), value_max=max(xdata))

The smoothing step reuses `fast_rolling_math(..., "mean")`, which this package
already ports for ATACorrect.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

from ..atac.tobias_atacorrect import fast_rolling_math

#: TOBIAS quantile grid: 0.05..0.99 inclusive, 1000 points. The lower bound
#: deliberately excludes most zeros.
QUANTILES = np.linspace(0.05, 0.99, 1000, endpoint=True)

#: Smoothing window applied to the correction curve.
SMOOTH_PAD = 50

#: Number of points used for the interpolation grid before the sigmoid fit.
INTERP_POINTS = 1000


def sigmoid(x, a, b, L, shift):
    """Reference sigmoid: `a` midpoint, `b` slope, `L` y-range, `shift` offset."""
    return L / (1 + np.exp(-b * (x - a))) + shift


class ArrayNorm:
    """Save a normalisation function and apply it to new arrays."""

    def __init__(self, function: str, popt, value_min: float, value_max: float):
        self.func = function
        self.popt = popt
        self.value_min = value_min
        self.value_max = value_max

    def get_norm_factor(self, arr):
        arr = np.asarray(arr, dtype=float)
        if self.func == "sigmoid":
            return sigmoid(arr, *self.popt)
        if self.func == "constant":
            return arr * 0 + self.popt
        raise ValueError("unknown norm function: %s" % self.func)

    def normalize(self, arr):
        arr = np.asarray(arr, dtype=float)
        # Cap before evaluating the factor, as the reference does, so outliers do
        # not extrapolate beyond the range seen during estimation.
        capped = arr * (arr <= self.value_max) + self.value_max * (arr > self.value_max)
        capped = (capped * (capped >= self.value_min)
                  + self.value_min * (capped < self.value_min))
        return arr * self.get_norm_factor(capped)


def _positive_quantiles(arr):
    """Quantiles over positive values only; all-zero input yields [0]."""
    arr = np.asarray(arr, dtype=float)
    if np.sum(arr > 0) == 0:
        return np.quantile([0], QUANTILES)
    return np.quantile(arr[arr > 0], QUANTILES)


def quantile_normalization(list_of_arrays, names, logger=None):
    """Build one `ArrayNorm` per condition, matching the reference procedure."""
    arrays = [np.asarray(a, dtype=float) for a in list_of_arrays]
    n = len(arrays)
    if n != len(names):
        raise ValueError("list_of_arrays and names must have the same length")

    array_quantiles = [_positive_quantiles(a) for a in arrays]
    mean_quantiles = np.mean(np.vstack(array_quantiles), axis=0)

    norm_objects = {}
    for i, name in enumerate(names):
        xdata = array_quantiles[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            ydata = mean_quantiles / xdata      # nan where xdata == 0

        # Smooth with a centred rolling mean over an edge-padded copy.
        padded = np.pad(ydata, SMOOTH_PAD, "edge")
        smoothed = fast_rolling_math(padded, SMOOTH_PAD, "mean")
        smoothed = smoothed[SMOOTH_PAD:-SMOOTH_PAD]

        finite = ~np.isnan(smoothed)
        if finite.sum() < 2:
            # Degenerate input: no usable correction curve.
            factor = 1.0
            norm_objects[name] = ArrayNorm("constant", factor,
                                           value_min=float(np.min(xdata)),
                                           value_max=float(np.max(xdata)))
            continue

        interp = interp1d(xdata[finite], smoothed[finite], kind="linear",
                          fill_value="extrapolate")
        xvals = np.linspace(np.min(xdata), np.max(xdata), INTERP_POINTS)
        y_inter = interp(xvals)

        keep = ~np.isnan(y_inter)
        xvals, y_inter = xvals[keep], y_inter[keep]

        try:
            a_range = np.min(xdata), np.max(xdata)
            l_range = np.diff(np.percentile(y_inter, [5, 95]))[0]
            bounds = ((a_range[0], -np.inf, 0, 0), (a_range[1], np.inf, l_range, np.inf))
            popt, _ = curve_fit(sigmoid, xvals, y_inter, bounds=bounds)
            func = "sigmoid"
        except Exception:                       # noqa: BLE001 - reference behaviour
            popt = float(np.mean(y_inter))
            func = "constant"

        norm_objects[name] = ArrayNorm(func, popt,
                                       value_min=float(np.min(xdata)),
                                       value_max=float(np.max(xdata)))
    return norm_objects


def normalize_scores(values, norm_object, clip_negative: bool = True):
    """Apply one `ArrayNorm` to a score array the way `process_tfbs` does."""
    values = np.asarray(values, dtype=float)
    out = norm_object.normalize(values)
    if clip_negative:
        out = np.where(out > 0, out, 0.0)
    return out


def log2_fold_change(score_a, score_b, pseudo):
    """`log2((a + pseudo) / (b + pseudo))`, as used per TFBS."""
    a = np.asarray(score_a, dtype=float) + pseudo
    b = np.asarray(score_b, dtype=float) + pseudo
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.round(np.log2(a / b), 5)
