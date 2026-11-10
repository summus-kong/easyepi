"""Unit tests for the Easyepi importable API.

Run in the ``easyepi`` conda environment:

    python -m pytest tests/ -q
"""

import numpy as np
import pytest

import easyepi as epi


# --- core stats ------------------------------------------------------------
def test_describe():
    d = epi.describe([1, 2, 3, 4, 5])
    assert d["n"] == 5
    assert d["mean"] == 3.0
    assert d["median"] == 3.0
    assert d["min"] == 1.0 and d["max"] == 5.0


def test_describe_empty():
    d = epi.describe([])
    assert d["n"] == 0
    assert np.isnan(d["mean"])


def test_significance_tests():
    a = np.array([1.0, 2, 3, 4, 5])
    b = np.array([6.0, 7, 8, 9, 10])
    assert epi.t_test(a, b)["pvalue"] < 0.01
    assert epi.mann_whitney_u(a, b)["pvalue"] < 0.05
    assert epi.wilcoxon(a, b)["test"] == "wilcoxon"


def test_effect_sizes():
    a = np.array([1.0, 2, 3])
    b = np.array([4.0, 5, 6])
    assert abs(epi.cohens_d(a, b)) > 0
    assert abs(epi.log2fc(a, b)) > 0


def test_bh_correction_monotone():
    p = np.array([0.01, 0.04, 0.031, 0.6, 0.2])
    adj = epi.benjamini_hochberg(p)
    # BH never reduces a p-value (adjusted >= original) and is bounded by 1.
    assert np.all(adj >= p - 1e-9)
    assert np.all(adj <= 1.0 + 1e-9)
    # Adjusted q-values, when mapped back to the sorted p-order, are monotone
    # non-decreasing (the key defining property of BH).
    order = np.argsort(p)
    sorted_adj = adj[order]
    assert np.all(np.diff(sorted_adj) >= -1e-9)


def test_bh_preserves_missing_values_without_global_contamination():
    adj = epi.benjamini_hochberg([0.01, np.nan, 0.2])
    assert np.isfinite(adj[[0, 2]]).all()
    assert np.isnan(adj[1])


def test_independent_default_and_paired_validation():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([4.0, 5.0, 6.0])
    assert epi.differential_test(a, b)["test"] == "mann-whitney-u"
    with pytest.raises(ValueError):
        epi.differential_test(a, b[:2], method="wilcoxon")


# --- data model ------------------------------------------------------------
def test_peakset():
    ps = epi.PeakSet(
        chrom=np.array(["chr1", "chr1"]),
        start=np.array([100, 900]),
        end=np.array([300, 1100]),
    )
    assert len(ps) == 2
    assert ps.width.tolist() == [200, 200]


def test_epigenomedata_groups():
    ed = epi.EpigenomeData(name="x")
    ed.add_sample("S1", {"group": "ctl"})
    ed.add_sample("S2", {"group": "ctl"})
    ed.add_sample("S3", {"group": "trt"})
    assert ed.sample_groups() == {"ctl": ["S1", "S2"], "trt": ["S3"]}


# --- differential ----------------------------------------------------------
def test_differential_features_significance():
    matrix = np.array([
        [12, 14, 13, 55, 60, 58],
        [5, 4, 6, 3, 2, 4],
        [20, 21, 18, 2, 1, 3],
    ], dtype=float)
    mask_a = np.array([True, True, True, False, False, False])
    res = epi.differential_features(matrix, mask_a, method="wilcoxon")
    assert res["n_regions"] == 3
    assert res["log2fc"].shape == (3,)
    assert res["padj"].shape == (3,)
    # g3 strongly up in A
    assert res["log2fc"][2] > 0


def test_differential_test_methods():
    a = np.array([1, 2, 3], dtype=float)
    b = np.array([5, 6, 7], dtype=float)
    for method in ("ttest", "wilcoxon", "mannwhitney"):
        r = epi.differential_test(a, b, method=method)
        assert "pvalue" in r and "log2fc" in r


# --- normalization ---------------------------------------------------------
def test_cpm():
    m = np.array([[1000, 2000], [500, 900]], dtype=float)
    out = epi.cpm_normalize(m, np.array([1e6, 2e6]))
    assert out.shape == (2, 2)
    # cpm of first sample = counts / lib * 1e6
    assert abs(out[0, 0] - 1000.0) < 1e-6


def test_normalization_factors():
    f = epi.compute_normalization_factors(
        ["A", "B"], lib_sizes={"A": 1e6, "B": 2e6}, spikeins={"B": 0.5}
    )
    assert f.lib_sizes["A"] == 1e6
    assert f.spikein["B"] == 0.5
    assert "scale_factors" in f.as_dict()


# --- nucleosome ------------------------------------------------------------
def test_occupancy_and_periodicity():
    sig = np.random.RandomState(0).randn(100).cumsum()
    occ = epi.occupancy_profile(sig)
    assert occ.shape == sig.shape
    assert abs(occ.mean()) < 1e-6  # ~z-scored
    per = epi.periodicity(sig)
    assert "dominant_period_bins" in per


# --- mnase profile CLI -----------------------------------------------------
def test_mnase_profile_cli_serialises_occupancy_json(tmp_path):
    """Regression: nucleosome_scores returns an ndarray occupancy, which is not
    JSON serialisable. The CLI must convert it instead of crashing."""
    import json
    from easyepi.specific.mnase.profile import main as mnase_main

    signal = tmp_path / "signal.tsv"
    signal.write_text("".join(f"{1.0 + (0.5 if i % 10 < 5 else -0.2):.3f}\n" for i in range(200)),
                      encoding="utf-8")
    out = tmp_path / "profile.json"
    assert mnase_main(["--signal", str(signal), "--out", str(out)]) == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert isinstance(payload["occupancy"], list) and len(payload["occupancy"]) == 200
    assert "nfr_score" in payload and "periodicity" in payload


def test_mnase_profile_cli_rejects_nonfinite_signal(tmp_path):
    from easyepi.specific.mnase.profile import main as mnase_main

    signal = tmp_path / "bad.tsv"
    signal.write_text("1.0\nnan\n2.0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="finite"):
        mnase_main(["--signal", str(signal), "--out", str(tmp_path / "out.json")])


# --- visualization ---------------------------------------------------------
def test_make_profile_returns_figure():
    fig = epi.make_profile(np.random.RandomState(0).randn(50).cumsum())
    assert fig is not None


def test_utr_is_strand_and_end_aware():
    from easyepi.annotate import _extract_utr
    exons = [(100, 200), (300, 400)]
    cds = [(150, 350)]
    assert _extract_utr(exons, cds, "+", "5") == [(100, 150)]
    assert _extract_utr(exons, cds, "+", "3") == [(350, 400)]
    assert _extract_utr(exons, cds, "-", "5") == [(350, 400)]
    assert _extract_utr(exons, cds, "-", "3") == [(100, 150)]
