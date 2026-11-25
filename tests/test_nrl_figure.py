"""The NRL autocorrelation figure is produced, and draws the lag it reports.

The computation is covered by tests/test_nucleosome_repeat_length.py, but nothing
pinned the *figure*: whether `mnase-qc` renders it, and whether the drawn curve
peaks where the reported NRL says. Measured on a synthetic ladder (fragment
centres 80/185/370/555 bp) `mnase-qc` wrote four figures including
`nrl_autocorrelation.png` and reported nrl_bp = 184.7 +/- 12.1.

The axes always carry two Line2D objects: the zero reference (`axhline`) and the
curve (`ax.plot`). The NRL marker is a third, added only when nrl_bp is given.
"""

from __future__ import annotations

import pytest

pytest.importorskip("matplotlib")

from easyepi import _plots


def _ladder_lags():
    """Lags on a 1 bp grid so the peak can be located exactly."""
    lags = list(range(0, 400))
    values = []
    for lag in lags:
        v = 1.0 * pow(2.718, -((lag - 185) ** 2) / (2 * 15.0 ** 2))
        v += 0.45 * pow(2.718, -((lag - 370) ** 2) / (2 * 15.0 ** 2))
        values.append(v)
    return lags, values


def test_the_curve_and_the_marker_agree():
    lags, values = _ladder_lags()
    fig = _plots.plot_nrl_autocorrelation(lags, values, nrl_bp=185.0,
                                          uncertainty=12.0)
    assert fig is not None
    ax = fig.axes[0]

    # lines[0] is the zero reference, lines[1] the curve, lines[2] the marker.
    curve = ax.lines[1]
    xs = [float(v) for v in curve.get_xdata()]
    ys = [float(v) for v in curve.get_ydata()]
    peak_lag = xs[ys.index(max(ys))]
    assert peak_lag == pytest.approx(185.0, abs=1.0)

    marker = float(ax.lines[2].get_xdata()[0])
    assert marker == pytest.approx(185.0)
    # The marker is drawn at the peak rather than at some fixed offset from it.
    assert abs(marker - peak_lag) <= 1.0

    assert "NRL" in ax.get_legend().get_texts()[0].get_text()
    assert "185" in ax.get_legend().get_texts()[0].get_text()


def test_without_an_nrl_only_the_zero_line_and_curve_are_drawn():
    """No ladder means no reported repeat: no marker, no legend, no invented value."""
    lags = list(range(0, 400, 10))
    fig = _plots.plot_nrl_autocorrelation(lags, [0.0] * len(lags), nrl_bp=None)
    assert fig is not None
    ax = fig.axes[0]
    assert len(ax.lines) == 2
    assert ax.get_legend() is None


def test_mismatched_lengths_are_refused():
    assert _plots.plot_nrl_autocorrelation([1, 2, 3], [1.0], nrl_bp=185.0) is None
    assert _plots.plot_nrl_autocorrelation([1.0], [1.0], nrl_bp=185.0) is None
