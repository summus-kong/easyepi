"""Additional figures: MA, peak overlap, TSS heatmap, locus track."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from easyepi.viz.genomeviews import (basepair_overlap, jaccard_bp, jaccard_index,
                                     loci_matrix, overlap_count, overlap_fraction,
                                     plot_locus_track, plot_ma, plot_peak_overlap,
                                     plot_tss_heatmap, read_peak_sets)


# --------------------------------------------------------------------------
# MA plot
# --------------------------------------------------------------------------
def test_ma_requires_an_abundance_column():
    frame = pd.DataFrame({"log2fc": [1.0, -1.0], "padj": [0.01, 0.5]})
    with pytest.raises(ValueError, match="abundance column"):
        plot_ma(frame)


def test_ma_requires_fold_change_column():
    with pytest.raises(ValueError, match="log2fc"):
        plot_ma(pd.DataFrame({"baseMean": [10.0], "padj": [0.1]}))


def test_ma_uses_base_mean_when_present():
    frame = pd.DataFrame({"log2fc": [2.0, -2.0, 0.1],
                          "baseMean": [100.0, 50.0, 5.0],
                          "padj": [0.001, 0.02, 0.9]})
    fig = plot_ma(frame)
    assert fig is not None
    ax = fig.axes[0]
    assert "baseMean" in ax.get_xlabel()
    assert ax.get_ylabel() == "log2 fold change"


def test_ma_accepts_an_explicit_mean_column():
    frame = pd.DataFrame({"log2FC": [1.0, -1.0], "mean_log2CPM": [3.0, 4.0],
                          "padj": [0.01, 0.5]})
    fig = plot_ma(frame, log2fc_col="log2FC", mean_col="mean_log2CPM")
    assert fig.axes[0].get_xlabel() == "mean_log2CPM"


def test_ma_rejects_a_table_with_no_complete_rows():
    frame = pd.DataFrame({"log2fc": [np.nan, np.nan], "baseMean": [np.nan, 1.0]})
    with pytest.raises(ValueError, match="no rows"):
        plot_ma(frame)


# --------------------------------------------------------------------------
# Peak overlap
# --------------------------------------------------------------------------
def _write(tmp_path, name, rows):
    p = tmp_path / name
    p.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n",
                 encoding="utf-8")
    return p


def test_overlap_fraction_identical_sets_is_one():
    rows = [("chr1", 100, 200), ("chr1", 500, 600)]
    assert overlap_fraction(rows, rows) == 1.0


def test_overlap_fraction_disjoint_sets_is_zero():
    a = [("chr1", 100, 200)]
    b = [("chr1", 5000, 6000)]
    assert overlap_fraction(a, b) == 0.0


def test_overlap_fraction_counts_partial_overlap_by_min_overlap():
    a = [("chr1", 100, 200)]
    b = [("chr1", 195, 400)]  # 5 bp shared
    assert overlap_fraction(a, b, min_overlap=1) == 1.0
    assert overlap_fraction(a, b, min_overlap=10) == 0.0


def test_overlap_fraction_ignores_other_chromosomes():
    a = [("chr1", 100, 200)]
    b = [("chr2", 100, 200)]
    assert overlap_fraction(a, b) == 0.0


def test_overlap_fraction_empty_query():
    assert overlap_fraction([], [("chr1", 1, 2)]) == 0.0


def test_read_peak_sets_skips_headers_and_rejects_single_set(tmp_path):
    a = _write(tmp_path, "A_peaks.bed",
               [("track", "name=x"), ("chr1", 1, 100), ("chr1", 200, 300)])
    sets = read_peak_sets([a])
    assert "A_peaks" in sets and len(sets["A_peaks"]) == 2

    b = _write(tmp_path, "B_peaks.bed", [("chr1", 50, 150)])
    two = read_peak_sets([a, b])
    assert set(two) == {"A_peaks", "B_peaks"}
    assert overlap_fraction(two["A_peaks"], two["B_peaks"]) == 0.5


def test_jaccard_is_symmetric_and_between_zero_and_one():
    a = [("chr1", 100, 200), ("chr1", 300, 400)]
    b = [("chr1", 150, 250)]
    assert jaccard_index(a, b) == jaccard_index(b, a)
    assert 0.0 < jaccard_index(a, b) < 1.0
    assert jaccard_index(a, a) == 1.0
    assert jaccard_index(a, [("chr2", 1, 2)]) == 0.0


def test_jaccard_of_a_subset_is_not_one():
    """A subset must not score 1.0 against its superset: that is the whole point."""
    small = [("chr1", 100, 200)]
    big = [("chr1", 100, 200), ("chr1", 300, 400)]
    assert overlap_fraction(small, big) == 1.0
    assert jaccard_index(small, big) == pytest.approx(1 / 2)


def test_overlap_count_matches_fraction_times_size():
    a = [("chr1", 100, 200), ("chr1", 300, 400), ("chr1", 900, 950)]
    b = [("chr1", 150, 250)]
    assert overlap_count(a, b) == 1
    assert overlap_fraction(a, b) == pytest.approx(1 / 3)


def test_basepair_overlap_counts_shared_bp_not_intervals():
    a = [("chr1", 100, 200)]          # 100 bp
    b = [("chr1", 150, 400)]          # 250 bp, shares 50
    o = basepair_overlap(a, b)
    assert o == {"shared_bp": 50, "bp_a": 100, "bp_b": 250}
    assert jaccard_bp(a, b) == pytest.approx(50 / (100 + 250 - 50))


def test_basepair_overlap_merges_adjacent_fragments():
    """Two touching intervals in one set must not be double counted."""
    fragmented = [("chr1", 100, 200), ("chr1", 200, 300)]
    contiguous = [("chr1", 100, 300)]
    o = basepair_overlap(fragmented, contiguous)
    assert o["bp_a"] == 200 and o["shared_bp"] == 200
    assert jaccard_bp(fragmented, contiguous) == 1.0


def test_bp_jaccard_does_not_inflate_a_broad_domain():
    """One wide domain must not look highly similar to many narrow peaks."""
    narrow = [("chr1", 1000 * i, 1000 * i + 200) for i in range(1, 6)]
    broad = [("chr1", 1000, 6000)]   # covers all five by interval count
    assert overlap_fraction(narrow, broad) == 1.0     # interval metric says 1.0
    assert jaccard_bp(narrow, broad) < 0.3            # bp metric is honest


def test_plot_peak_overlap_rejects_unknown_metric(tmp_path):
    a = _write(tmp_path, "A.bed", [("chr1", 1, 100)])
    b = _write(tmp_path, "B.bed", [("chr1", 50, 150)])
    with pytest.raises(ValueError, match="metric"):
        plot_peak_overlap([a, b], metric="cosine")


def test_plot_peak_overlap_is_symmetric_on_the_diagonal(tmp_path):
    a = _write(tmp_path, "A.bed", [("chr1", 100, 200), ("chr1", 300, 400)])
    b = _write(tmp_path, "B.bed", [("chr1", 150, 250)])
    fig = plot_peak_overlap([a, b])
    assert fig is not None
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert labels == ["A", "B"]


def test_plot_peak_overlap_needs_two_sets(tmp_path):
    a = _write(tmp_path, "only.bed", [("chr1", 1, 10)])
    with pytest.raises(ValueError, match="at least two"):
        plot_peak_overlap([a])


# --------------------------------------------------------------------------
# TSS heatmap
# --------------------------------------------------------------------------
def test_loci_matrix_extracts_symmetric_windows():
    signal = np.arange(1000, dtype=float)
    mat = loci_matrix(signal, [500, 700], flank=10)
    assert mat.shape == (2, 21)
    assert mat[0, 10] == signal[500]


def test_loci_matrix_drops_loci_that_would_be_truncated():
    signal = np.zeros(100)
    mat = loci_matrix(signal, [2, 50, 98], flank=10)
    assert mat.shape == (1, 21)  # only the centred locus survives


def test_loci_matrix_empty_when_no_locus_fits():
    mat = loci_matrix(np.zeros(10), [1], flank=50)
    assert mat.shape[0] == 0


def test_tss_heatmap_validates_shape():
    with pytest.raises(ValueError, match="non-empty"):
        plot_tss_heatmap(np.empty((0, 21)), flank=10)
    with pytest.raises(ValueError, match="columns"):
        plot_tss_heatmap(np.zeros((5, 10)), flank=10)


def test_tss_heatmap_sorts_rows_by_centre_signal():
    mat = np.zeros((4, 21))
    mat[:, 10] = [1.0, 4.0, 2.0, 3.0]
    fig = plot_tss_heatmap(mat, flank=10, max_rows=10)
    assert fig is not None
    assert len(fig.axes) >= 2


def test_tss_heatmap_downsampled_to_max_rows():
    mat = np.ones((100, 21))
    fig = plot_tss_heatmap(mat, flank=10, max_rows=10)
    assert fig is not None


# --------------------------------------------------------------------------
# Locus track
# --------------------------------------------------------------------------
def test_locus_track_requires_a_track():
    with pytest.raises(ValueError, match="at least one track"):
        plot_locus_track([], "chr1", 0, 100)


def test_locus_track_rejects_inverted_region():
    with pytest.raises(ValueError, match="end must exceed"):
        plot_locus_track([{"name": "t", "signal": np.ones(5)}], "chr1", 100, 100)


def test_locus_track_draws_signal_and_interval_tracks():
    tracks = [
        {"name": "signal", "signal": np.linspace(0, 1, 50), "bin_size": 20},
        {"name": "peaks", "intervals": [(100, 300), (700, 900)]},
    ]
    fig = plot_locus_track(tracks, "chr1", 0, 1000,
                           annotations=[{"start": 400, "end": 600, "label": "gene A"}])
    assert fig is not None
    titles = [ax.get_title(loc="left") for ax in fig.axes]
    assert "signal" in titles and "peaks" in titles and "annotations" in titles


def test_locus_track_keeps_out_of_view_intervals_harmless():
    fig = plot_locus_track([{"name": "p", "intervals": [(-500, -100), (5000, 6000)]}],
                           "chr1", 0, 1000)
    assert fig is not None


# --------------------------------------------------------------------------
# Files are actually written when asked
# --------------------------------------------------------------------------
def test_save_writes_a_png(tmp_path):
    frame = pd.DataFrame({"log2fc": [1.0, -1.0], "baseMean": [10.0, 20.0],
                          "padj": [0.01, 0.3]})
    out = tmp_path / "ma"
    plot_ma(frame, save=str(out), formats=["png"])
    assert (tmp_path / "ma.png").exists()
    assert (tmp_path / "ma.png").stat().st_size > 1000
