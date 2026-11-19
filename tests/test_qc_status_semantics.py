"""QC must say *why* a metric is absent, not just that it is absent.

Before this, an unreadable annotation and an annotation that genuinely has no
TSS both produced an empty list, so the report could not distinguish "we failed
to read your file" from "your file has nothing to read". These tests pin the
four outcomes apart.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from easyepi.qc_extended import read_tss, tss_enrichment

GTF_GENE = 'chr1\thavana\tgene\t100\t200\t.\t+\t.\tgene_id "g1";\n'
GTF_EXON_ONLY = (
    'chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
)


def test_missing_annotation_is_reported_as_input_missing(tmp_path: Path):
    status: dict[str, str] = {}
    assert read_tss(str(tmp_path / "nope.gtf"), status=status) == []
    assert status["state"] == "input_missing"
    assert "does not exist" in status["detail"]


def test_annotation_without_gene_lines_is_not_a_failure(tmp_path: Path):
    path = tmp_path / "exon_only.gtf"
    path.write_text(GTF_EXON_ONLY, encoding="utf-8")
    status: dict[str, str] = {}
    assert read_tss(str(path), status=status) == []
    # A readable file that simply has no TSS is a different finding from a
    # crash, and must not be reported as one.
    assert status["state"] == "no_tss"


def test_unreadable_annotation_is_reported_as_calculation_failed(tmp_path: Path):
    status: dict[str, str] = {}
    # A directory exists, so the existence check passes, but opening it raises.
    assert read_tss(str(tmp_path), status=status) == []
    assert status["state"] == "calculation_failed"
    assert "Error" in status["detail"]


def test_readable_annotation_reports_ok(tmp_path: Path):
    path = tmp_path / "ok.gtf"
    path.write_text(GTF_GENE, encoding="utf-8")
    status: dict[str, str] = {}
    assert read_tss(str(path), status=status) == [("chr1", 99)]
    assert status["state"] == "ok"


def test_status_is_optional_so_existing_callers_are_unaffected(tmp_path: Path):
    path = tmp_path / "ok.gtf"
    path.write_text(GTF_GENE, encoding="utf-8")
    assert read_tss(str(path)) == [("chr1", 99)]


def test_tss_enrichment_carries_the_reason_into_its_result(tmp_path: Path):
    missing = tss_enrichment(str(tmp_path / "x.bw"), str(tmp_path / "x.gtf"))
    assert missing["status"] == "input_missing"
    assert "does not exist" in missing["error"]
