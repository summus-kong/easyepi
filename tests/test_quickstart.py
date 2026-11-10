"""Guided entry points: `easyepi init` and `easyepi quick`."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from easyepi.quickflow import make_plan
from easyepi.quickstart import (build_project, guess_assay, guess_genome,
                                project_to_toml, scan_samples)


def _touch(path: Path, text: str = "x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# Assay and genome inference
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name,expected", [
    ("ATAC_T0", "ATAC-seq"),
    ("atac_rep1", "ATAC-seq"),
    ("CUTandTAG_H3K27ac", "CUT&Tag"),
    ("CUTTag_K4", "CUT&Tag"),
    ("CUTandRUN_K27", "CUT&RUN"),
    ("MNase_30min", "MNase-seq"),
    ("H3K4me3_rep1", "histone"),
    ("H3K27me3_rep2", "histone"),
    ("H3K27ac", "histone"),
    ("WGBS_liver", "methylation"),
    ("ChIP_ER", "ChIP-seq"),
    ("random_sample_1", "unknown"),
])
def test_guess_assay(name, expected):
    assert guess_assay(name) == expected


@pytest.mark.parametrize("text,expected", [
    ("/ref/GRCh38/genome.fa", "hg38"),
    ("hg19", "hg19"),
    ("GRCh37.primary.fa", "hg19"),
    ("mm10_ref", "mm10"),
    ("/data/GRCm39/", "mm39"),
    ("nothing-here", None),
])
def test_guess_genome(text, expected):
    assert guess_genome(text) == expected


# --------------------------------------------------------------------------
# Folder scanning
# --------------------------------------------------------------------------
def test_scan_groups_tracks_and_ignores_indexes(tmp_path):
    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s1.bam.bai")
    _touch(tmp_path / "s1.bw")
    _touch(tmp_path / "s1.narrowPeak")
    got = scan_samples(tmp_path)
    assert list(got) == ["s1"]
    assert set(got["s1"]) == {"bam", "bigwig", "peak"}


def test_scan_rejects_two_files_of_the_same_kind(tmp_path):
    _touch(tmp_path / "s1.bw")
    _touch(tmp_path / "s1.bigWig")
    with pytest.raises(ValueError, match="two bigwig files"):
        scan_samples(tmp_path)


def test_scan_rejects_empty_folder(tmp_path):
    with pytest.raises(ValueError, match="no BAM/BigWig/peak"):
        scan_samples(tmp_path)


def test_scan_rejects_missing_folder(tmp_path):
    with pytest.raises(NotADirectoryError, match="not a directory"):
        scan_samples(tmp_path / "nope")


# --------------------------------------------------------------------------
# Project generation and TOML output
# --------------------------------------------------------------------------
def test_build_project_records_guesses_in_notes(tmp_path):
    _touch(tmp_path / "ATAC_rep1.bam")
    _touch(tmp_path / "ATAC_rep1.narrowPeak")
    project, notes = build_project(tmp_path)
    assert project["samples"][0]["assay"] == "ATAC-seq"
    # every inference must be surfaced, never silent
    assert any("assay" in n or "genome" in n for n in notes)


def test_build_project_flags_unknown_assay(tmp_path):
    _touch(tmp_path / "sampleA.bam")
    _project, notes = build_project(tmp_path)
    assert any("sampleA" in n and "not detected" in n for n in notes)


def test_build_project_records_references(tmp_path):
    _touch(tmp_path / "s1.bam")
    project, _notes = build_project(tmp_path, genome_build="hg38",
                                    annotation="/ref/g.gtf",
                                    motif_database="/ref/j.pfm")
    assert project["project"]["genome_build"] == "hg38"
    assert project["references"]["annotation"] == "/ref/g.gtf"
    assert project["references"]["motif_database"] == "/ref/j.pfm"


def test_project_to_toml_roundtrips_through_tomllib(tmp_path):
    import tomllib

    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s1.bw")
    project, _ = build_project(tmp_path, genome_build="hg38")
    text = project_to_toml(project)
    parsed = tomllib.loads(text)
    assert parsed["project"]["genome_build"] == "hg38"
    assert parsed["samples"][0]["sample_id"] == "s1"
    assert parsed["samples"][0]["bam"].endswith("s1.bam")


def test_toml_escapes_windows_style_paths(tmp_path):
    import tomllib

    _touch(tmp_path / "s1.bam")
    project, _ = build_project(tmp_path, annotation=r"C:\ref\genes.gtf")
    parsed = tomllib.loads(project_to_toml(project))
    assert parsed["references"]["annotation"] == r"C:\ref\genes.gtf"


# --------------------------------------------------------------------------
# Plan: what can run, given the inputs
# --------------------------------------------------------------------------
def test_plan_with_bam_and_annotation_runs_qc_peaks_annotate(tmp_path):
    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s1.narrowPeak")
    project, _ = build_project(tmp_path, genome_build="hg38",
                               annotation="/ref/g.gtf")
    pfile = tmp_path / "p.toml"
    pfile.write_text(project_to_toml(project), encoding="utf-8")
    plan = make_plan(pfile)
    assert "qc" in plan.steps
    assert "annotate" in plan.steps
    assert "motifs" not in plan.steps        # no motif database
    assert plan.steps[-1] == "report"


def test_plan_from_folder_source(tmp_path):
    _touch(tmp_path / "s1.bam")
    plan = make_plan(tmp_path)
    assert "qc" in plan.steps
    assert any("annotation" in n for n in plan.notes)


def test_plan_rejects_unavailable_explicit_step(tmp_path):
    _touch(tmp_path / "s1.bam")
    with pytest.raises(ValueError, match="not available"):
        make_plan(tmp_path, steps=["motifs"])


def test_plan_errors_on_missing_source(tmp_path):
    with pytest.raises(FileNotFoundError, match="no such project"):
        make_plan(tmp_path / "absent")


def test_plan_errors_on_project_without_samples(tmp_path):
    pfile = tmp_path / "p.toml"
    pfile.write_text("[project]\nproject_id = \"x\"\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no samples"):
        make_plan(pfile)


# --------------------------------------------------------------------------
# CLI wiring
# --------------------------------------------------------------------------
def test_cli_exposes_init():
    """`init` remains: writing a project file is a genuine convenience."""
    from easyepi import cli
    parser = cli.build_parser()
    a = parser.parse_args(["init", "/data/tracks", "--out", "p.toml"])
    assert a.module == "init" and a.source == "/data/tracks"


def test_quick_command_is_gone():
    """The one-shot pipeline was removed in favour of step-by-step tools."""
    from easyepi import cli
    parser = cli.build_parser()
    assert "quick" not in set(parser._subparsers._group_actions[0].choices)

def test_init_then_quick_plan_are_consistent(tmp_path):
    import tomllib

    _touch(tmp_path / "H3K4me3_rep1.bam")
    _touch(tmp_path / "H3K4me3_rep1.narrowPeak")
    from easyepi.quickflow import run_init

    class Args:
        source = str(tmp_path)
        out = str(tmp_path / "project.toml")
        project_id = "demo"
        genome_build = "hg38"
        annotation = None
        motif_database = None

    assert run_init(Args) == 0
    parsed = tomllib.loads((tmp_path / "project.toml").read_text(encoding="utf-8"))
    assert parsed["project"]["project_id"] == "demo"
    assert parsed["samples"][0]["assay"] == "histone"

    plan = make_plan(tmp_path / "project.toml")
    assert "qc" in plan.steps
    assert "report" in plan.steps


def test_quick_records_a_failed_step_without_aborting(tmp_path, monkeypatch):
    """One broken step must not stop the chain; the summary reports it."""
    from easyepi import quickflow

    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s1.narrowPeak")
    project, _ = build_project(tmp_path, genome_build="hg38")
    pfile = tmp_path / "p.toml"
    pfile.write_text(project_to_toml(project), encoding="utf-8")

    def boom(*_a, **_k):
        raise RuntimeError("simulated step failure")

    monkeypatch.setattr(quickflow, "_step_qc", boom)

    class Args:
        source = str(pfile)
        out = str(tmp_path / "run")
        steps = ["qc", "report"]
        verbose = False

    assert quickflow.run_quick(Args) == 1        # failure is signalled
    summary = json.loads((tmp_path / "run" / "easyepi_report.json").read_text())
    statuses = {s["step"]: s["status"] for s in summary["steps"]}
    assert statuses["qc"] == "failed"
    # the report step writes the file that describes itself and the other steps
    assert (tmp_path / "run" / "SUMMARY.md").exists()


# --------------------------------------------------------------------------
# Small studies: no differential analysis, few samples
# --------------------------------------------------------------------------
def test_plan_omits_differential_without_groups(tmp_path):
    """A study with no group column must not be given a differential step."""
    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s1.narrowPeak")
    project, _ = build_project(tmp_path, genome_build="hg38")
    pfile = tmp_path / "p.toml"
    pfile.write_text(project_to_toml(project), encoding="utf-8")
    plan = make_plan(pfile)
    assert "differential" not in plan.steps
    assert any("differential" in n for n in plan.notes)


def test_plan_offers_differential_with_two_groups(tmp_path):
    _touch(tmp_path / "ctrl1.bam")
    _touch(tmp_path / "ctrl1.narrowPeak")
    _touch(tmp_path / "trt1.bam")
    _touch(tmp_path / "trt1.narrowPeak")
    pfile = tmp_path / "p.toml"
    pfile.write_text(_project_with_groups(tmp_path, ["ctrl", "treat"]), encoding="utf-8")
    plan = make_plan(pfile)
    assert "differential" in plan.steps


def test_single_replicate_per_group_still_plans_differential(tmp_path):
    """n=1 per group is legitimate: the implementation reports fold change only."""
    _touch(tmp_path / "a1.bam")
    _touch(tmp_path / "a1.narrowPeak")
    _touch(tmp_path / "b1.bam")
    _touch(tmp_path / "b1.narrowPeak")
    pfile = tmp_path / "p.toml"
    pfile.write_text(_project_with_groups(tmp_path, ["A", "B"]), encoding="utf-8")
    assert "differential" in make_plan(pfile).steps


def test_differential_step_skips_when_a_peak_file_is_missing(tmp_path):
    from easyepi import quickflow

    _touch(tmp_path / "a1.bam")
    _touch(tmp_path / "a1.narrowPeak")
    _touch(tmp_path / "b1.bam")            # no peak file for b1
    pfile = tmp_path / "p.toml"
    pfile.write_text(_project_with_groups(tmp_path, ["A", "B"]), encoding="utf-8")
    plan = quickflow.make_plan(pfile)
    result = quickflow._step_differential(plan.project, pfile, tmp_path / "out")
    assert result.status == "skipped"
    assert "b1" in result.detail


def test_two_or_three_samples_are_not_rejected(tmp_path):
    """The chain must run for a two-sample study, not demand a cohort."""
    _touch(tmp_path / "s1.bam")
    _touch(tmp_path / "s2.bam")
    plan = make_plan(tmp_path)
    assert "qc" in plan.steps
    assert "report" in plan.steps


def _project_with_groups(folder: Path, labels, genome_build: str | None = "hg38") -> str:
    """Build a project TOML from scanned tracks, then add group labels.

    ``build_project`` never invents groups, so a test that needs a two-group
    design has to write them explicitly, exactly as a user would edit the file.
    """
    project, _ = build_project(folder, genome_build=genome_build)
    for row, label in zip(project["samples"], labels):
        row["group"] = label
    return project_to_toml(project)
