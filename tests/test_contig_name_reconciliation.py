"""Align the mitochondrial test with what the code deliberately does.

`_SKIP_PREFIXES` contains `MT`: mitochondrial contigs are excluded from TSS
enrichment on purpose, because a TSS profile dominated by mitochondrial reads
measures contamination rather than regulation. The reconciliation added for
`chr1` vs `1` must therefore *not* resurrect them.

The first version of this test asserted the opposite and failed, which is the
right outcome: the skip is intentional and the test was wrong. It is replaced by
one that pins the intention, plus a check that a non-mitochondrial rename still
works so the skip cannot hide a broken reconciler.
"""

from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

from easyepi.tss_enrichment import tss_enrichment  # noqa: E402

TSS = (20000, 60000, 100000, 140000)


def _build(path, contig="chr1"):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": 200000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for t in TSS:
            for i in range(60):
                a = pysam.AlignedSegment()
                a.query_name = "t%d_%d" % (t, i)
                a.query_sequence = "A" * 50
                a.flag = 0
                a.reference_id = 0
                a.reference_start = t - 25 + (i % 50)
                a.mapping_quality = 60
                a.cigarstring = "50M"
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(a)
            for off in (20, 45, 70, 95, 3905, 3930, 3955, 3980):
                a = pysam.AlignedSegment()
                a.query_name = "b%d_%d" % (t, off)
                a.query_sequence = "A" * 50
                a.flag = 0
                a.reference_id = 0
                a.reference_start = t - 2000 + off
                a.mapping_quality = 60
                a.cigarstring = "50M"
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(a)
    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


@pytest.fixture()
def bam(tmp_path):
    path = tmp_path / "sig.bam"
    _build(path)
    return str(path)


def test_ensembl_names_give_the_same_score_as_ucsc_names(bam):
    ucsc = [("chr1", t, "+") for t in TSS]
    ens = [("1", t, "+") for t in TSS]

    a = tss_enrichment(bam, ucsc)
    b = tss_enrichment(bam, ens)

    assert a["tss_enrichment_status"] == "measured"
    assert b["tss_enrichment_status"] == "measured"
    # The number is what rules out a reconciliation that mapped to the wrong
    # contig: that would also avoid `no_tss`.
    assert b["tss_enrichment_score"] == a["tss_enrichment_score"]


def test_the_rename_is_reported_not_silent(bam):
    ucsc = [("chr1", t, "+") for t in TSS]
    ens = [("1", t, "+") for t in TSS]

    assert tss_enrichment(bam, ucsc)["chromosomes_renamed"] == {}
    assert tss_enrichment(bam, ens)["chromosomes_renamed"] == {"1": "chr1"}


def test_a_truly_absent_contig_is_still_reported(bam):
    other = [("chrZ", t, "+") for t in TSS]

    res = tss_enrichment(bam, other)

    assert res["tss_enrichment_status"] == "no_tss"
    assert "chrZ" in res["error"]
    assert res["chromosomes_unmatched"] == ["chrZ"]


def test_mitochondrial_contigs_stay_excluded(bam):
    """Mitochondrial contigs are dropped on purpose, not reconciled.

    `_SKIP_PREFIXES` lists ``MT`` and ``chrM``: a TSS profile built on
    mitochondrial reads measures contamination. The `chr1`-vs-`1` fix must not
    quietly resurrect them, so they still end in `no_tss`.
    """
    for spelling in ("MT", "chrM"):
        res = tss_enrichment(bam, [(spelling, t, "+") for t in TSS])
        assert res["tss_enrichment_status"] == "no_tss", spelling


def test_a_plain_rename_still_works_alongside_the_mito_skip(bam):
    """The skip must not hide a broken reconciler.

    ``1`` is a normal contig: it is renamed, measured, and reported. Only the
    mitochondrial spellings are excluded.
    """
    res = tss_enrichment(bam, [("1", t, "+") for t in TSS])

    assert res["tss_enrichment_status"] == "measured"
    assert res["chromosomes_renamed"] == {"1": "chr1"}
