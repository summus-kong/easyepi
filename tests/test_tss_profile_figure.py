"""The TSS metagene must reach the figure path, and describe itself honestly.

BACKLOG item 11 asked for "ATAC: fragment-size and TSS-enrichment profiles". The
fragment-size half already worked; the TSS half was missing because
``_assay_metric_frame`` keeps only scalars and two special-cased lists, so the
nested ``tss_enrichment_detail`` -- and the pileup inside it -- never reached a
plotter. The reader had been computing ``profile`` all along.

These tests pin the three things that were actually wrong or could regress:

  * the detail dict, and the profile in it, survive the flattening;
  * ``profile`` keeps its documented length (2*flank+1) rather than being
    silently resampled for the sake of a smaller JSON;
  * the plotter builds its x axis from ``profile_center_index``, so a peak at
    the TSS lands at 0 instead of wherever the midpoint happens to be.
"""

from __future__ import annotations

import json

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("matplotlib")

from easyepi import _plots


def _payload(tmp_path, profile, centre):
    record = {
        "sample_id": "s1",
        "assay": "atac",
        "assay_qc": {
            "tss_enrichment": 3.4,
            "tss_enrichment_detail": {
                "tss_enrichment_status": "measured",
                "profile": profile,
                "profile_center_index": centre,
            },
        },
    }
    path = tmp_path / "qc.json"
    path.write_text(json.dumps({"samples": {"s1": record}}), encoding="utf-8")
    return tmp_path


def test_the_profile_survives_the_flattening(tmp_path):
    _payload(tmp_path, [0, 0, 9, 0, 0], 2)
    _frame, _beta, _frag, profiles = _plots._assay_metric_frame(tmp_path)
    assert len(profiles) == 1
    label, detail = profiles[0]
    assert label == "s1"
    assert detail["profile"] == [0, 0, 9, 0, 0]


def test_a_record_without_a_profile_is_skipped(tmp_path):
    record = {"sample_id": "s1", "assay": "atac",
              "assay_qc": {"tss_enrichment": 1.0,
                           "tss_enrichment_detail": {"tss_enrichment_status": "no_tss"}}}
    (tmp_path / "qc.json").write_text(json.dumps({"samples": {"s1": record}}),
                                      encoding="utf-8")
    _frame, _beta, _frag, profiles = _plots._assay_metric_frame(tmp_path)
    assert profiles == []


def test_the_plot_puts_the_tss_at_zero(tmp_path):
    """A peak at profile_center_index must be drawn at x == 0."""
    profile = [0] * 2001
    profile[1400] = 50          # deliberately off-centre
    fig = _plots.plot_tss_profile([("s1", {"profile": profile,
                                           "profile_center_index": 1400})])
    assert fig is not None
    ax = fig.axes[0]
    xdata, ydata = ax.lines[0].get_xdata(), ax.lines[0].get_ydata()
    peak_at = float(np.asarray(xdata)[int(np.argmax(ydata))])
    assert peak_at == pytest.approx(0.0, abs=1e-9)
    # The axis spans the whole window: +/-1000 for a 2001-value profile.
    assert float(np.min(xdata)) == pytest.approx(-1400.0)
    assert float(np.max(xdata)) == pytest.approx(600.0)


def test_no_profiles_means_no_figure():
    assert _plots.plot_tss_profile([]) is None
    assert _plots.plot_tss_profile([("s1", {"profile": [1]})]) is None
