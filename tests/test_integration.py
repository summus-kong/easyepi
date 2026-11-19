from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import easyepi as epi
from easyepi.integration import IntegrationError, GeneRegion, cluster_tracks, link_peaks_to_genes, read_gene_bed, read_rna_results


def test_cluster_tracks_is_deterministic_and_keeps_track_contract():
    index = ["g1", "g2", "g3", "g4"]
    track_a = pd.DataFrame([[1, 2], [1, 3], [9, 8], [10, 9]], index=index)
    track_b = pd.DataFrame([[2, 1], [3, 1], [8, 9], [9, 10]], index=index)
    labels_a, features_a, inertia_a = cluster_tracks({"ATAC": track_a, "H3K27ac": track_b}, 2, seed=7)
    labels_b, features_b, inertia_b = cluster_tracks({"ATAC": track_a, "H3K27ac": track_b}, 2, seed=7)
    assert labels_a.equals(labels_b)
    assert features_a.equals(features_b)
    assert inertia_a == inertia_b
    assert all(column.startswith(("ATAC|", "H3K27ac|")) for column in features_a.columns)


def test_peak_gene_links_are_explicit_nearest_tss_candidates(tmp_path: Path):
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("chr1\t90\t110\nchr1\t850\t870\n", encoding="utf-8")
    genes = [GeneRegion("plus", "chr1", 100, 300, "+"), GeneRegion("minus", "chr1", 500, 860, "-")]
    links = link_peaks_to_genes(peaks, genes, max_distance=50)
    assert links["gene_id"].tolist() == ["plus", "minus"]
    assert links["link_method"].unique().tolist() == ["nearest_tss"]
    assert links["signed_distance_to_tss"].tolist() == [0, 0]


def test_gene_bed_and_rna_id_validation(tmp_path: Path):
    genes = tmp_path / "genes.bed"
    genes.write_text("chr1\t1\t100\tG1\t0\t+\n", encoding="utf-8")
    assert read_gene_bed(genes)[0].gene_id == "G1"
    rna = tmp_path / "rna.tsv"
    rna.write_text("gene_id\tlog2FoldChange\nENSG1.4\t1.5\n", encoding="utf-8")
    assert read_rna_results(rna)["gene_id"].tolist() == ["ENSG1"]
    with pytest.raises(IntegrationError, match="gene column"):
        read_rna_results(rna, "missing")


def test_integration_is_exported_from_package():
    assert epi.run_integration is not None
