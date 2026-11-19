"""Guided entry points: `easyepi init` and `easyepi quick`."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

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

