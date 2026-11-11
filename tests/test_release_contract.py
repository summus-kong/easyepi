"""Regression tests for migration, method boundaries and data integrity."""
import argparse
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from easyepi._peak_calling import macs3_command
from easyepi._plots import add_figure_args, selected_formats, plot_volcano
from easyepi.core.stats import benjamini_hochberg
from easyepi.differential import differential_features, run_diff, read_design
from easyepi.differential_matrix import read_narrowpeak_union
from easyepi.methylation import read_methylation_table, MethylationError
from easyepi.motifs import read_pwm, scan_sequences
from easyepi.motifs import external_tools
from easyepi.specific.histone.chromhmm import main as chromhmm_main


def test_peak_commands_preserve_assay_and_fragment_semantics():
    sample = SimpleNamespace(bam="sample.bam", sample_id="S1")
    atac = macs3_command("macs3", sample, "out", genome_size="hs", assay="ATAC")
    chip = macs3_command("macs3", sample, "out", genome_size="hs", assay="CHIP",
                         broad=True, paired_end=False, control="input.bam")
    assert "BAMPE" in atac and "--shift" not in atac
    assert "--broad" in chip and "input.bam" in chip and "--shift" not in chip
    single = macs3_command("macs3", sample, "out", genome_size="hs", assay="ATAC", paired_end=False)
    assert single[single.index("--shift")+1] == "-100"


def test_chromhmm_headerless_table_and_official_commands(tmp_path):
    project = tmp_path / "project.toml"
    project.write_text("[[samples]]\nsample_id='S1'\nassay='histone'\npeak='p.bed'\ncell='C1'\nmark='H3K27ac'\n")
    output = tmp_path / "chromhmm"
    assert chromhmm_main(["--project", str(project), "--out", str(output), "--no-check-files",
                         "--dry-run", "--run-binarize", "--run-learn", "--chromhmm-jar", "ChromHMM.jar",
                         "--genome", "chrom.sizes", "--genome-build", "hg38", "--threads", "4"]) == 0
    assert (output / "cellmarkfiletable.tsv").read_text().splitlines() == ["C1\tH3K27ac\ttrack_1.bed"]
    commands = json.loads((output / "chromhmm_manifest.json").read_text())["commands"]
    assert "-peaks" in commands[0] and commands[0][commands[0].index("-b")+1] == "200"
    assert str(output / "input_beds") in commands[0]
    assert commands[1][commands[1].index("-p")+1] == "4" and "-s" in commands[1]
    assert not (output / "input_beds" / "track_1.bed").exists()


@pytest.mark.parametrize("row", ["chr1\t0\t1.2\t10", "chr1\t-1\t0.2\t10", "chr1\t0.5\t0.2\t10"])
def test_invalid_methylation_is_rejected(tmp_path, row):
    path = tmp_path / "meth.tsv"
    path.write_text("chrom\tstart\tbeta\tcoverage\n" + row + "\n")
    with pytest.raises(MethylationError):
        read_methylation_table(path)


def test_methylation_counts_supply_coverage(tmp_path):
    path = tmp_path / "meth.tsv"
    path.write_text("chrom\tstart\tmethylated\tunmethylated\nchr1\t10\t3\t7\n")
    frame = read_methylation_table(path)
    assert frame.coverage.tolist() == [10] and frame.beta.tolist() == [0.3]


def test_peak_union_does_not_drop_substring_chromosome_names(tmp_path):
    path = tmp_path / "p.bed"
    path.write_text("chr10\t0\t10\nchr1\t0\t10\nchr1\t5\t15\n")
    assert read_narrowpeak_union([str(path)]) == [("chr1",0,15),("chr10",0,10)]


def test_bh_uses_only_tested_features():
    actual = benjamini_hochberg([0.01, np.nan, 0.2])
    assert np.allclose(actual[[0,2]], [0.02,0.2]) and np.isnan(actual[1])


def test_paired_differential_reorders_by_subject(tmp_path, monkeypatch):
    matrix = tmp_path / "matrix.tsv"
    matrix.write_text("region\ta1\ta2\tb2\tb1\nr1\t1\t2\t20\t10\n")
    design = tmp_path / "design.tsv"
    design.write_text("sample_id\tgroup\tsubject\na1\tA\t1\na2\tA\t2\nb2\tB\t2\nb1\tB\t1\n")
    captured = {}
    def fake_test(values, mask, **kwargs):
        captured["matrix"] = values
        return {key: np.array([1.]) for key in ("mean_a","mean_b","log2fc","cohens_d","pvalue","padj")}
    monkeypatch.setattr("easyepi.differential.differential_features",fake_test)
    args = SimpleNamespace(counts=str(matrix),out=str(tmp_path/"out.tsv"),design=str(design),
                           method="wilcoxon",pair="subject",equal_var=False)
    assert run_diff(args) == 0
    assert captured["matrix"].tolist() == [[1.,2.,10.,20.]]


def test_missing_values_keep_paired_observations_aligned():
    from scipy.stats import wilcoxon
    matrix = np.array([[1.,2.,3.,5.,np.nan,7.]])
    actual = differential_features(matrix,[True,True,True,False,False,False],method="wilcoxon")
    assert actual["pvalue"][0] == wilcoxon([1.,3.],[5.,7.]).pvalue
    with pytest.raises(ValueError,match="two biological"):
        differential_features([[1.,2.]],[True,False])


def test_duplicate_design_is_rejected(tmp_path):
    path = tmp_path / "design.tsv"
    path.write_text("sample_id\tgroup\na\tA\na\tB\n")
    with pytest.raises(ValueError,match="duplicate"):
        read_design(path,["a"])


def test_figures_export_as_pdf_by_default_and_nan_is_not_significant():
    """A run writes figures without being asked; --no-figures turns that off.

    The previous contract was "no flag, no figure", which made a run that should
    have produced a plot indistinguishable from one that did. The default is now
    the vector format, with an explicit opt-out for table-only runs.
    """
    parser = argparse.ArgumentParser()
    add_figure_args(parser)
    assert selected_formats(parser.parse_args([])) == ("pdf",)
    assert selected_formats(parser.parse_args(["--figure-format","pdf"])) == ("pdf",)
    assert selected_formats(parser.parse_args(["--no-figures"])) == ()
    assert selected_formats(
        parser.parse_args(["--no-figures","--figure-format","png"])) == ("png",)
    fig = plot_volcano(pd.DataFrame({"log2fc":[1.,2.],"padj":[np.nan,.01]}))
    assert len(fig.axes[0].collections[1].get_offsets()) == 1


def test_reverse_complement_flag_and_forward_coordinates():
    pwm = np.array([[3,-5],[-5,3],[-5,-5],[-5,-5]],dtype=float)
    assert scan_sequences(pwm,["TTGT"],reverse_complement=False)["scores"][0] < 6
    reverse = scan_sequences(pwm,["TTGT"])
    assert reverse["strands"] == ["-"] and reverse["positions"] == [2]


def test_labelled_jaspar_parser(tmp_path):
    path=tmp_path/"motif.pfm"
    path.write_text(">M1 demo\nA [ 10 0 ]\nC [ 0 10 ]\nG [ 0 0 ]\nT [ 0 0 ]\n")
    assert read_pwm(path).shape == (4,2)


def test_external_failure_propagates(monkeypatch):
    def fail(*args,**kwargs):
        assert kwargs["check"] is True
        raise subprocess.CalledProcessError(2,args[0],stderr="invalid input")
    monkeypatch.setattr(external_tools.subprocess,"run",fail)
    with pytest.raises(subprocess.CalledProcessError):
        external_tools._run(["fake-tool"])


def test_fimo_saves_hits_and_excludes_header(tmp_path,monkeypatch):
    monkeypatch.setattr(external_tools,"_require",lambda *a,**k:"fimo")
    monkeypatch.setattr(external_tools,"_run",lambda *a,**k:SimpleNamespace(
        returncode=0,stderr="",stdout="motif_id\tsequence_name\nM1\tS1\n"))
    result=external_tools.run_fimo_scan("motifs.meme","target.fa",str(tmp_path))
    assert result["n_hits"] == 1
    assert Path(result["hits_file"]).is_file()


def test_known_motif_enrichment_uses_score_threshold(tmp_path):
    pytest.importorskip("Bio")
    from easyepi.universal.motif_enrichment import main
    target=tmp_path/"target.fa";background=tmp_path/"background.fa";pwm=tmp_path/"m.pfm"
    target.write_text(">t1\nAAAAAA\n>t2\nAAAAAA\n")
    background.write_text(">b1\nCCCCCC\n>b2\nCCCCCC\n")
    pwm.write_text("100 100 100 100 100 100\n0 0 0 0 0 0\n0 0 0 0 0 0\n0 0 0 0 0 0\n")
    output=tmp_path/"enrich"
    assert main(["--target",str(target),"--background",str(background),"--pwm",str(pwm),"--out",str(output),"--fpr","0.01"]) == 0
    row=json.loads((output/"motif_enrichment.json").read_text())["motifs"][0]
    assert row["target_hits"]==2 and row["background_hits"]==0
    # A run writes its figure without being asked; the vector format is the default.
    written = sorted(p.name for p in (output / "figures").glob("*"))
    assert any(n.endswith(".pdf") for n in written), written
    assert not any(n.endswith(".png") for n in written), written


def test_consensus_counts_overlapping_samples(tmp_path):
    from easyepi.specific.atac.consensus import main
    (tmp_path/"a.bed").write_text("chr1\t0\t20\nchr1\t5\t15\n")
    (tmp_path/"b.bed").write_text("chr1\t10\t30\n")
    project=tmp_path/"project.toml"
    project.write_text("[[samples]]\nsample_id='A'\nassay='ATAC-seq'\npeak='a.bed'\n[[samples]]\nsample_id='B'\nassay='ATAC-seq'\npeak='b.bed'\n")
    out=tmp_path/"consensus"
    assert main(["--project",str(project),"--out",str(out),"--min-samples","2"]) == 0
    assert (out/"atac_consensus.bed").read_text().split("\t")[:3] == ["chr1","10","20"]


def test_normalization_rejects_invented_quantities():
    from easyepi.normalize import cpm_normalize, rpgc_normalize, compute_normalization_factors
    with pytest.raises(ValueError):
        cpm_normalize(np.ones((2,2)),[0,10])
    with pytest.raises(ValueError,match="Missing measured"):
        compute_normalization_factors(["S1"])
    with pytest.raises(ValueError,match="fragment length"):
        rpgc_normalize(np.ones((1,2)),1000,[10,20])
    assert np.allclose(rpgc_normalize(np.ones((1,2)),1000,[10,20],100),[[1.,.5]])
    factors=compute_normalization_factors(["S1"],lib_sizes={"S1":2e6},spikeins={"S1":.5})
    assert factors.factor_for("S1") == factors.scale_factors["S1"] == .25


def test_activity_rejects_sequence_only_and_constant_networks():
    from easyepi.activity import estimate_tf_activity
    matrix=pd.DataFrame({"S1":[1,2,3],"S2":[2,3,4]},index=["p1","p2","p3"])
    network=pd.DataFrame({"source":["TF"]*3,"target":matrix.index,"weight":[1.]*3})
    with pytest.raises(ValueError,match="constant network weights"):
        estimate_tf_activity(matrix,network,min_targets=2)


def test_nearest_tss_selects_correct_gene_across_chromosomes():
    """Regression: the per-chromosome TSS list was indexed with a global index,
    so on multi-chromosome input the reported gene did not match the distance."""
    from easyepi.annotate import _nearest_tss_distance

    tss_chrom = np.array(['chr2', 'chr3', 'chr1', 'chr1', 'chr1'], dtype=object)
    tss_pos = np.array([100, 200, 300, 500, 700])
    tss_id = np.array(['g_chr2', 'g_chr3', 'g_chr1_a', 'g_chr1_b', 'g_chr1_c'], dtype=object)
    tss_strand = np.array(['+', '+', '+', '-', '+'], dtype=object)
    ids, dists = _nearest_tss_distance(
        np.array([490]), np.array([510]), np.array(['chr1'], dtype=object),
        tss_pos, tss_chrom, tss_id, tss_strand)
    assert ids == ['g_chr1_b']
    assert dists[0] == 0


def test_heatmap_names_are_unambiguous():
    """Regression: ``plot_heatmap`` resolved to an array-based function via
    ``easyepi`` but to a table-based one via ``easyepi.plot``."""
    import easyepi
    from easyepi import plot as plot_module
    from easyepi.viz import plot_matrix_heatmap

    assert easyepi.plot_heatmap is plot_module.plot_heatmap
    assert easyepi.plot_heatmap is not plot_matrix_heatmap
    assert easyepi.__all__.count('plot_heatmap') == 1
    assert easyepi.__all__.count('plot_matrix_heatmap') == 1
    assert callable(plot_matrix_heatmap)


def test_deprecated_viz_heatmap_alias_still_works():
    import warnings
    from easyepi import viz

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fig = viz.plot_heatmap(np.random.RandomState(0).rand(4, 4))
    assert fig is not None
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)


def test_no_placeholder_positioning_results(tmp_path):
    from easyepi.specific.mnase.positioning import main
    with pytest.raises(SystemExit) as exc:
        main([])
    assert exc.value.code == 2
