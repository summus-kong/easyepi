"""GFF3 children must attach to every parent they name.

A shared exon in GFF3 is written once with several parents:

    chr1  x  exon  66346  66499  .  -  .  Parent=TxID:1,TxID:2,TxID:3

Keeping only the first id silently dropped the exon from transcripts 2 and 3. On
a real TxDb export that cost 36% of transcripts their middle exons, which then
surfaced as missing introns and empty UTR bins -- with no error anywhere.

These tests pin the one-to-many relation, and cross-check that a GTF written from
the same model yields the same per-transcript exon sets.

Note on the reader: this file targets ``tx_features``, the module
``annotate/__init__.py`` actually imports. Its ``read_transcripts`` returns
``TranscriptRecord`` objects that expose attributes, not dicts, and its
intervals are 0-based half-open, so GFF3 ``1000-1500`` becomes ``(999, 1500)``.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from easyepi.annotate.tx_features import detect_format, read_transcripts

GFF3_SHARED_EXON = """\
##gff-version 3
chr1\tx\tgene\t1000\t9000\t.\t+\t.\tID=GeneID:g1
chr1\tx\tmRNA\t1000\t3000\t.\t+\t.\tParent=GeneID:g1;ID=TxID:t1
chr1\tx\tmRNA\t1000\t9000\t.\t+\t.\tParent=GeneID:g1;ID=TxID:t2
chr1\tx\texon\t1000\t1500\t.\t+\t.\tParent=TxID:t1,TxID:t2
chr1\tx\texon\t2000\t2500\t.\t+\t.\tParent=TxID:t1,TxID:t2
chr1\tx\texon\t3000\t3200\t.\t+\t.\tParent=TxID:t1
chr1\tx\texon\t8000\t9000\t.\t+\t.\tParent=TxID:t2
chr1\tx\tCDS\t1000\t1500\t.\t+\t.\tParent=TxID:t1,TxID:t2
chr1\tx\tCDS\t2000\t2500\t.\t+\t.\tParent=TxID:t1,TxID:t2
"""


def _write(tmp_path: Path, name: str, text: str) -> str:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def _by_id(path: str) -> dict:
    return {t.transcript_id: t for t in read_transcripts(path)}


def test_a_shared_exon_reaches_every_parent(tmp_path):
    txs = _by_id(_write(tmp_path, "shared.gff3", GFF3_SHARED_EXON))

    assert set(txs) == {"TxID:t1", "TxID:t2"}
    # The two shared exons must be present in BOTH transcripts.
    assert len(txs["TxID:t1"].exons) == 3   # 1000-1500, 2000-2500, 3000-3200
    assert len(txs["TxID:t2"].exons) == 3   # 1000-1500, 2000-2500, 8000-9000
    assert len(txs["TxID:t1"].cds) == 2
    assert len(txs["TxID:t2"].cds) == 2


def test_shared_exons_are_the_same_intervals_in_each_parent(tmp_path):
    txs = _by_id(_write(tmp_path, "shared2.gff3", GFF3_SHARED_EXON))
    # 0-based half-open, so GFF3 1000-1500 (1-based, inclusive) is (999, 1500).
    shared = {(999, 1500), (1999, 2500)}
    assert shared <= set(txs["TxID:t1"].exons)
    assert shared <= set(txs["TxID:t2"].exons)


def test_transcript_start_and_end_span_all_its_exons(tmp_path):
    txs = _by_id(_write(tmp_path, "span.gff3", GFF3_SHARED_EXON))
    # The parser stores 0-based half-open intervals, so a GFF3 column-1 coordinate
    # of 1000 becomes 999 and the end stays exclusive.
    assert txs["TxID:t1"].start == 999 and txs["TxID:t1"].end == 3200
    assert txs["TxID:t2"].start == 999 and txs["TxID:t2"].end == 9000


def test_gtf_and_gff3_of_the_same_model_agree(tmp_path):
    """The two input paths must not disagree once the gene model is identical."""
    gtf = """\
chr1\tx\ttranscript\t1000\t3000\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t1";
chr1\tx\ttranscript\t1000\t9000\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t2";
chr1\tx\texon\t1000\t1500\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t1";
chr1\tx\texon\t2000\t2500\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t1";
chr1\tx\texon\t3000\t3200\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t1";
chr1\tx\texon\t1000\t1500\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t2";
chr1\tx\texon\t2000\t2500\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t2";
chr1\tx\texon\t8000\t9000\t.\t+\t.\tgene_id "g1"; transcript_id "TxID:t2";
"""
    a = _by_id(_write(tmp_path, "same.gff3", GFF3_SHARED_EXON))
    b = _by_id(_write(tmp_path, "same.gtf", gtf))

    assert set(a) == set(b)
    for tx in a:
        assert sorted(a[tx].exons) == sorted(b[tx].exons), tx
        assert a[tx].tss == b[tx].tss, tx


def test_format_detection_on_both(tmp_path):
    assert detect_format(_write(tmp_path, "a.gff3", GFF3_SHARED_EXON)) == "gff3"
    gtf = 'chr1\tx\ttranscript\t1\t9\t.\t+\t.\tgene_id "g";\n'
    assert detect_format(_write(tmp_path, "a.gtf", gtf)) == "gtf"
