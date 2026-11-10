"""Composite panels: gene tracks, UpSet, sample correlation."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from easyepi.viz.panels import (correlation_matrix, expand_csv_values,
                                gene_window,
                                parse_gene_names, parse_track_specs,
                                plot_gene_panel, plot_sample_correlation,
                                plot_upset)


GENES = [
    {"chrom": "chr1", "start": 1000, "end": 2000, "name": "GeneA", "strand": "+"},
    {"chrom": "chr1", "start": 5000, "end": 6200, "name": "GeneB", "strand": "-"},
]


# --------------------------------------------------------------------------
# gene_window
# --------------------------------------------------------------------------
def test_gene_window_pads_and_clamps():
    assert gene_window({"start": 1000, "end": 2000}, flank=500) == (500, 2500)
    # must not go negative
    assert gene_window({"start": 100, "end": 200}, flank=500)[0] == 0


def test_gene_window_rejects_empty_span():
    with pytest.raises(ValueError, match="empty span"):
        gene_window({"start": 500, "end": 500, "name": "X"})


# --------------------------------------------------------------------------
# plot_gene_panel validation (no pyBigWig needed for these)
# --------------------------------------------------------------------------
def test_gene_panel_requires_tracks():
    with pytest.raises(ValueError, match="at least one track"):
        plot_gene_panel([], GENES)


def test_gene_panel_requires_genes():
    with pytest.raises(ValueError, match="at least one gene"):
        plot_gene_panel([{"name": "t", "bigwig": "x.bw"}], [])


def test_gene_panel_requires_bigwig_key():
    with pytest.raises(ValueError, match="needs a 'bigwig' path"):
        plot_gene_panel([{"name": "t"}], GENES)


def test_gene_panel_requires_gene_keys():
    with pytest.raises(ValueError, match="needs 'start'"):
        plot_gene_panel([{"name": "t", "bigwig": "x.bw"}],
                        [{"chrom": "chr1", "end": 100}])


def test_gene_panel_rejects_bad_bin_size():
    with pytest.raises(ValueError, match="bin_size"):
        plot_gene_panel([{"name": "t", "bigwig": "x.bw"}], GENES, bin_size=0)


def test_gene_panel_uses_stubbed_reader(monkeypatch):
    """Draw a 2 loci x 2 marks grid with a synthetic signal provider."""
    from easyepi.viz import panels

    def fake_bins(path, chrom, start, end, n_bins):
        x = np.linspace(0, 1, n_bins)
        return np.exp(-0.5 * ((x - 0.5) / 0.12) ** 2) * (5.0 if "H3K27ac" in str(path) else 2.0)

    monkeypatch.setattr(panels, "bigwig_bins", fake_bins)
    fig = plot_gene_panel(
        [{"name": "H3K4me3", "bigwig": "a.bw"}, {"name": "H3K27ac", "bigwig": "b.bw"}],
        GENES, flank=500, bin_size=25)
    # 2 marks + 1 gene-model row, 2 gene columns
    assert len(fig.axes) == 6
    titles = [ax.get_title() for ax in fig.axes[:2]]
    assert titles == ["GeneA", "GeneB"]


def test_gene_panel_shared_ylim_makes_rows_comparable(monkeypatch):
    from easyepi.viz import panels

    def fake_bins(path, chrom, start, end, n_bins):
        # locus 1 is 10x stronger than locus 2
        amp = 10.0 if start < 3000 else 1.0
        return np.full(n_bins, amp)

    monkeypatch.setattr(panels, "bigwig_bins", fake_bins)
    fig = plot_gene_panel([{"name": "mark", "bigwig": "a.bw"}], GENES,
                          flank=100, bin_size=25, shared_ylim=True)
    top, right = fig.axes[0], fig.axes[1]
    assert top.get_ylim() == right.get_ylim()


def test_gene_panel_draws_exons_when_given(monkeypatch):
    from easyepi.viz import panels
    monkeypatch.setattr(panels, "bigwig_bins",
                        lambda *a, **k: np.ones(10))
    fig = plot_gene_panel([{"name": "m", "bigwig": "a.bw"}], GENES, flank=100,
                          exons={"GeneA": [(1000, 1200), (1500, 2000)]})
    # rows: 1 track + 1 gene model; axes are row-major, so the model row for GeneA
    # is the first of the last two axes (2 gene columns).
    model_ax = fig.axes[-2]
    assert len(model_ax.patches) == 2  # two exons


def test_gene_panel_reports_missing_chromosome(monkeypatch):
    from easyepi.viz import panels

    def boom(path, chrom, start, end, n_bins):
        raise ValueError(f"chromosome '{chrom}' absent from {path}")

    monkeypatch.setattr(panels, "bigwig_bins", boom)
    with pytest.raises(ValueError, match="absent from"):
        plot_gene_panel([{"name": "m", "bigwig": "a.bw"}], GENES)


# --------------------------------------------------------------------------
# UpSet
# --------------------------------------------------------------------------
def test_upset_needs_two_sets():
    with pytest.raises(ValueError, match="at least two"):
        plot_upset({"A": [("chr1", 1, 10)]})


def test_upset_counts_exclusive_and_shared_patterns():
    sets = {
        "A": [("chr1", 100, 200), ("chr1", 1000, 1100)],
        "B": [("chr1", 150, 250)],
    }
    fig = plot_upset(sets)
    assert fig is not None
    ax_top = fig.axes[1]
    heights = sorted(int(round(p.get_height())) for p in ax_top.patches)
    # Three distinct membership patterns: {A only}, {A,B}, {B only}. A contributes
    # two elementary segments to the "A only" pattern, so that bar has height 2 and
    # the total region count is 4.
    assert len(heights) == 3
    assert sum(heights) == 4
    assert max(heights) == 2


def test_upset_merges_fragments_before_membership():
    """A fragmented set must not double count a single region."""
    sets = {
        "A": [("chr1", 100, 150), ("chr1", 150, 200)],
        "B": [("chr1", 100, 200)],
    }
    fig = plot_upset(sets)
    ax_top = fig.axes[1]
    total = sum(int(round(p.get_height())) for p in ax_top.patches)
    assert total == 1  # one merged region, shared by both sets


def test_upset_honours_min_size():
    sets = {
        "A": [("chr1", 100, 200), ("chr1", 5000, 5100)],
        "B": [("chr1", 150, 250)],
    }
    with pytest.raises(ValueError, match="at least"):
        plot_upset(sets, min_size=99)


def test_upset_rejects_empty_sets():
    with pytest.raises(ValueError, match="no intervals"):
        plot_upset({"A": [], "B": []})


# --------------------------------------------------------------------------
# Sample correlation
# --------------------------------------------------------------------------
def _frame():
    rng = np.random.default_rng(0)
    base = rng.random(200) * 100
    return pd.DataFrame({
        "s1": base + rng.random(200) * 3,
        "s2": base + rng.random(200) * 3,        # a duplicate of s1
        "s3": rng.random(200) * 100,             # unrelated
    })


def test_correlation_matrix_is_symmetric_with_unit_diagonal():
    corr, labels = correlation_matrix(_frame())
    assert labels == ["s1", "s2", "s3"]
    assert np.allclose(corr, corr.T)
    assert np.allclose(np.diag(corr), 1.0)


def test_correlation_recovers_the_planted_duplicate():
    corr, labels = correlation_matrix(_frame())
    i = labels.index("s1")
    j = labels.index("s2")
    k = labels.index("s3")
    assert corr[i, j] > 0.99
    assert corr[i, k] < 0.5


def test_correlation_drops_constant_rows():
    frame = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0],
                          "b": [2.0, 4.0, 6.0, 8.0],
                          "c": [5.0, 5.0, 5.0, 5.0]})
    # 'c' is constant, so it cannot correlate; it must not produce NaN
    corr, labels = correlation_matrix(frame, min_features=1)
    assert np.isfinite(corr).all()
    assert "c" not in labels


def test_correlation_rejects_too_few_varying_features():
    frame = pd.DataFrame({"a": [1.0, 1.0, 1.0], "b": [2.0, 2.0, 2.0]})
    with pytest.raises(ValueError, match="features vary"):
        correlation_matrix(frame, min_features=10)


def test_correlation_needs_two_samples():
    with pytest.raises(ValueError, match="at least two samples"):
        correlation_matrix(pd.DataFrame({"only": [1.0, 2.0]}))


def test_correlation_requires_a_dataframe():
    with pytest.raises(TypeError, match="pandas DataFrame"):
        correlation_matrix(np.zeros((5, 3)))


def test_spearman_is_available():
    corr, _ = correlation_matrix(_frame(), method="spearman")
    assert np.isfinite(corr).all()


def test_plot_sample_correlation_clusters_without_error():
    fig = plot_sample_correlation(_frame(), cluster=True)
    assert fig is not None
    assert len(fig.axes) >= 3


def test_plot_sample_correlation_without_clustering():
    fig = plot_sample_correlation(_frame(), cluster=False)
    assert fig is not None


def test_save_writes_png(tmp_path):
    plot_sample_correlation(_frame(), cluster=False,
                            save=str(tmp_path / "corr"), formats=["png"])
    assert (tmp_path / "corr.png").exists()
    assert (tmp_path / "corr.png").stat().st_size > 1000

# --------------------------------------------------------------------------
# Compact CLI syntax: --bigwig NAME[,NAME] and --gene A,B
# --------------------------------------------------------------------------
def test_parse_gene_names_accepts_repeats_and_commas():
    assert parse_gene_names(["AGRN,TP73", "SKI"]) == ["AGRN", "TP73", "SKI"]
    assert parse_gene_names(["AGRN", "AGRN"]) == ["AGRN"]      # de-duplicated
    assert parse_gene_names([" AGRN , PLCH2 "]) == ["AGRN", "PLCH2"]
    assert parse_gene_names(None) == []


def test_parse_gene_names_preserves_user_order():
    assert parse_gene_names(["TP73,AGRN,SKI"]) == ["TP73", "AGRN", "SKI"]


def test_parse_track_specs_takes_a_plain_path(tmp_path):
    """A single path is the common case and needs no label."""
    bw = tmp_path / "H3K4me3.bw"
    bw.write_bytes(b"x")
    specs = parse_track_specs([str(bw)])
    assert specs == [{"name": "H3K4me3", "bigwig": str(bw)}]


def test_parse_track_specs_accepts_a_relative_path(tmp_path, monkeypatch):
    """Relative paths resolve against the working directory."""
    (tmp_path / "tracks").mkdir()
    (tmp_path / "tracks" / "H3K27ac.bw").write_bytes(b"x")
    monkeypatch.chdir(tmp_path)
    specs = parse_track_specs(["tracks/H3K27ac.bw"])
    assert specs[0]["name"] == "H3K27ac"
    assert specs[0]["bigwig"] == "tracks/H3K27ac.bw"


def test_parse_track_specs_explicit_label_overrides_the_stem(tmp_path):
    bw = tmp_path / "ENCFF123ABC.bw"
    bw.write_bytes(b"x")
    specs = parse_track_specs([f"H3K4me3={bw}"])
    assert specs[0]["name"] == "H3K4me3"


def test_parse_track_specs_accepts_comma_lists(tmp_path):
    a = tmp_path / "a.bw"; a.write_bytes(b"x")
    b = tmp_path / "b.bw"; b.write_bytes(b"x")
    specs = parse_track_specs([f"{a},{b}"])
    assert [s["name"] for s in specs] == ["a", "b"]


def test_parse_track_specs_reports_a_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="input not found"):
        parse_track_specs([str(tmp_path / "absent.bw")])


def test_parse_track_specs_rejects_a_duplicate_label(tmp_path):
    a = tmp_path / "sub1" / "x.bw"; a.parent.mkdir(); a.write_bytes(b"x")
    b = tmp_path / "sub2" / "x.bw"; b.parent.mkdir(); b.write_bytes(b"x")
    with pytest.raises(ValueError, match="duplicate track name"):
        parse_track_specs([str(a), str(b)])


def test_parse_track_specs_requires_a_path():
    with pytest.raises(ValueError, match="no path"):
        parse_track_specs(["NAME="])


def test_parse_track_specs_requires_at_least_one():
    with pytest.raises(ValueError, match="at least one"):
        parse_track_specs([])


# --------------------------------------------------------------------------
# Shared comma-list expansion
# --------------------------------------------------------------------------
def test_expand_csv_values_splits_and_strips():
    assert expand_csv_values(["a,b", " c "]) == ["a", "b", "c"]


def test_expand_csv_values_keeps_colons_and_equals():
    """Values like NAME=PATH and sample=BAM:PEAK must survive intact."""
    assert expand_csv_values(["S=a.bam:p.bed"]) == ["S=a.bam:p.bed"]
    assert expand_csv_values(["M=x.bg"]) == ["M=x.bg"]


def test_expand_csv_values_empty():
    assert expand_csv_values([]) == []
    assert expand_csv_values(None) == []


def test_histone_states_accepts_comma_separated_signals():
    from easyepi import histone_states_cli
    assert expand_csv_values(["H3K4me3=a.bg,H3K27ac=b.bg"]) == ["H3K4me3=a.bg", "H3K27ac=b.bg"]
