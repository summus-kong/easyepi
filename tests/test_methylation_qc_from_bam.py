"""Methylation QC must read beta values from the BAM, as the WGBS page promises.

The documented command is `easyepi qc --assay methylation --bam ctrl1.bam`: the
bisulfite BAM carries MM/ML or Bismark XM tags, so no table is needed. Before this
test existed, that command reported generic read flags and an `incomplete` verdict
with no beta distribution at all, because the payload only looked for a table path.

These tests build a small BAM with XM tags and assert the beta metrics appear, and
that a BAM without either tag is reported as an error rather than as 0% methylated.
"""
from __future__ import annotations

import json

import pytest

pysam = pytest.importorskip("pysam")

from easyepi.qc_payload import build_sample_payload  # noqa: E402

SITES = 6
DEPTH = 10


def _write_bisulfite_bam(path, xm: bool = True, sites: int = SITES, depth: int = DEPTH):
    """Site i is methylated in i of its `depth` reads, giving betas 0 .. (depth-1)/depth."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 100000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for site in range(sites):
            for read in range(depth):
                a = pysam.AlignedSegment()
                a.query_name = f"s{site}_r{read}"
                a.query_sequence = "C"
                a.flag = 0
                a.reference_id = 0
                a.reference_start = 1000 + site * 50
                a.mapping_quality = 60
                a.cigarstring = "1M"
                a.query_qualities = pysam.qualitystring_to_array("I")
                if xm:
                    a.set_tag("XM", "Z" if read < site else "z")
                out.write(a)
    return path


def test_beta_metrics_come_from_the_bam(tmp_path):
    bam = _write_bisulfite_bam(tmp_path / "wgbs.bam")
    built = build_sample_payload("methylation", bam=str(bam))

    assay_qc = built["assay_qc"]
    assert "methylation_error" not in assay_qc, assay_qc.get("methylation_error")
    assert assay_qc["mean_beta"] == pytest.approx(0.25, abs=0.05)
    assert assay_qc["beta_detail"]["n_sites"] == SITES
    # the tag reader's own accounting travels with the result
    assert assay_qc["methylation_bam"]["encoding"] == "XM"
    assert assay_qc["methylation_bam"]["sites"] == SITES


def test_a_bam_without_tags_is_an_explained_error(tmp_path):
    bam = _write_bisulfite_bam(tmp_path / "plain.bam", xm=False)
    built = build_sample_payload("methylation", bam=str(bam))

    assay_qc = built["assay_qc"]
    # Not "0% methylated": the absence of the tag is reported as a problem.
    assert "mean_beta" not in assay_qc
    assert "MM/ML" in assay_qc["methylation_error"]
    assert "XM" in assay_qc["methylation_error"]


def test_a_table_still_works_and_wins_when_both_are_given(tmp_path):
    # The table reader accepts `beta`, or `methylated` + `unmethylated` counts.
    table = tmp_path / "calls.tsv"
    table.write_text("chrom\tpos\tstrand\tcontext\tmethylated\tunmethylated\n"
                     "chr1\t10\t+\tCpG\t8\t2\nchr1\t20\t+\tCpG\t2\t8\n",
                     encoding="utf-8")
    bam = _write_bisulfite_bam(tmp_path / "wgbs.bam")

    from_table = build_sample_payload("methylation", methylation=str(table))["assay_qc"]
    both = build_sample_payload("methylation", bam=str(bam),
                                methylation=str(table))["assay_qc"]

    # the explicit table is what the caller asked for, so it is what is used
    assert from_table["beta_detail"] == both["beta_detail"]
    assert "methylation_bam" not in both


def test_low_coverage_sites_are_dropped_not_averaged(tmp_path):
    bam = _write_bisulfite_bam(tmp_path / "thin.bam", sites=3, depth=1)
    assay_qc = build_sample_payload("methylation", bam=str(bam))["assay_qc"]

    # One read per site cannot support a beta, and the count dropped is reported.
    assert assay_qc["methylation_bam"]["dropped_below_coverage"] == 3
    assert "mean_beta" not in assay_qc
