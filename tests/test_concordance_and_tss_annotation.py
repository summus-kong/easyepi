"""TSS from annotation, the four-state status, and replicate concordance.

The status transitions are the point of the first half: a sparse library must be
reported as measured-and-failed, not as never-measured, or a downstream verdict will
say "not evaluated" for a sample that has an answer.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.concordance import pairwise_concordance, replicate_concordance
from easyepi.tss_enrichment import read_tss_table

pysam = pytest.importorskip("pysam")
from easyepi.tss_enrichment import tss_enrichment  # noqa: E402

HEADER = {"HD": {"VN": "1.6", "SO": "coordinate"},
          "SQ": [{"SN": "chr1", "LN": 1_000_000}]}


def _write_bam(path, starts, reverse=False):
    """Write one 1 bp alignment per start and return the path.

    Returning it matters: every caller hands the result straight to
    ``tss_enrichment``, which opens it with pysam. Leaving the value implicit
    passed pysam ``None``, and the resulting TypeError surfaced three frames
    away as a bare ``{"tss_enrichment_status": "error"}`` -- a failure that
    named nothing that was actually wrong.
    """
    with pysam.AlignmentFile(str(path), "wb", header=HEADER) as out:
        for i, s in enumerate(starts):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A"
            a.flag = 16 if reverse else 0
            a.reference_id = 0
            a.reference_start = int(s)
            a.mapping_quality = 60
            a.cigarstring = "1M"
            a.query_qualities = pysam.qualitystring_to_array("I")
            out.write(a)
    return path


# ------------------------------------------------------------ TSS from annotation
def _write_gtf(path):
    path.write_text(
        'chr1\tsrc\tgene\t1001\t9000\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding";\n'
        'chr1\tsrc\ttranscript\t1001\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; '
        'transcript_type "protein_coding";\n'
        'chr1\tsrc\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tsrc\ttranscript\t7001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; '
        'transcript_type "nonsense_mediated_decay";\n'
        'chr2\tsrc\tgene\t5001\t8000\t.\t-\t.\tgene_id "G2"; gene_type "protein_coding";\n'
        'chr2\tsrc\ttranscript\t5001\t8000\t.\t-\t.\tgene_id "G2"; transcript_id "T3"; '
        'transcript_type "protein_coding";\n',
        encoding="utf-8")


def test_gtf_is_detected_and_read_at_transcript_level(tmp_path):
    gtf = tmp_path / "a.gtf"
    _write_gtf(gtf)
    rows = read_tss_table(str(gtf))
    # Two transcripts on chr1 plus one on chr2.
    assert len(rows) == 3, rows
    # plus strand: TSS is the 0-based start; minus strand: the end
    assert ("chr1", 1000, "+") in rows
    assert ("chr1", 7000, "+") in rows
    assert ("chr2", 8000, "-") in rows


def test_gene_level_gives_one_tss_per_gene(tmp_path):
    gtf = tmp_path / "a.gtf"
    _write_gtf(gtf)
    rows = read_tss_table(str(gtf), level="gene")
    assert len(rows) == 2, rows
    assert ("chr1", 1000, "+") in rows


def test_transcript_type_filter_narrows_the_set(tmp_path):
    gtf = tmp_path / "a.gtf"
    _write_gtf(gtf)
    rows = read_tss_table(str(gtf))
    coding = [r for r in rows]
    assert len(coding) == 3
    from easyepi.annotation_tss import tss_from_annotation

    only = tss_from_annotation(str(gtf), include_types=("protein_coding",))
    assert len(only) == 2, only


def test_a_tsv_is_still_read_as_a_tsv(tmp_path):
    tsv = tmp_path / "tss.tsv"
    tsv.write_text("chrom\ttss\tstrand\nchr1\t500000\t+\n", encoding="utf-8")
    assert read_tss_table(str(tsv)) == [("chr1", 500_000, "+")]


# ------------------------------------------------------------- four-state status
TSS = [("chr1", 500_000, "+")]


def test_normal_library_reports_measured(tmp_path):
    starts = list(range(499_960, 500_041)) * 20 + list(range(498_000, 502_000, 200))
    bam = _write_bam(tmp_path / "ok.bam", starts)
    res = tss_enrichment(bam, TSS, flank=2000, center=50, flank_bg=100)
    assert res["tss_enrichment_status"] == "measured"
    assert res["tss_enrichment_score"] > 5
    assert res["reads_truncated"] is False


def test_zero_background_is_measured_and_failed_not_missing(tmp_path):
    """Signal but no flanks: a sparse library, not an absent measurement."""
    starts = list(range(499_960, 500_041)) * 20      # pile only, zero background
    bam = _write_bam(tmp_path / "sparse.bam", starts)
    res = tss_enrichment(bam, TSS, flank=2000, center=50, flank_bg=100)
    assert res["tss_enrichment_score"] is None
    assert res["tss_enrichment_status"] == "background_zero"
    assert res["signal"] > 0 and res["background"] == 0
    assert "sparse" in res["error"]


def test_an_empty_library_reports_no_reads(tmp_path):
    bam = _write_bam(tmp_path / "empty.bam", [])
    res = tss_enrichment(bam, TSS, flank=2000)
    assert res["tss_enrichment_status"] == "no_reads"


def test_an_unreadable_bam_keeps_the_full_result_shape(tmp_path):
    """A failed read must still return ``reads_truncated``.

    Every return path of `tss_enrichment` promises the same keys. The error path
    once omitted ``reads_truncated``, so the only moment a consumer touched that
    key on a failed sample it raised KeyError -- turning a per-sample data fault
    into a crash of the whole summary.
    """
    res = tss_enrichment(str(tmp_path / "absent.bam"), TSS, flank=2000)
    assert res["tss_enrichment_status"] == "error"
    assert res["tss_enrichment_score"] is None
    assert res["reads_truncated"] is False
    assert res["error"]


def test_coverage_capped_scan_is_flagged(tmp_path):
    """A truncated scan is not interchangeable with a complete one."""
    starts = list(range(499_000, 501_000))
    bam = _write_bam(tmp_path / "big.bam", starts)
    res = tss_enrichment(bam, TSS, flank=2000, max_reads=100)
    assert res["reads_truncated"] is True
    assert res["reads_used"] == 100


# ------------------------------------------------------------------ concordance
def _matrix(seed: int = 3):
    """Six samples in two groups of three, each group sharing its own profile.

    Every sample carries a genome-wide component plus a *group-specific* one, so
    replicates of a group resemble each other more than they resemble the other
    group -- which is the situation the metric exists to detect.

    An earlier version gave each group a constant offset instead
    (``base + offset + noise``). A correlation is mean-invariant, so that offset
    cannot separate the groups at all: the honest margin was then pure noise
    (mean +0.00003, sd 0.00058 over 400 seeds) and ``grouping_supported`` came
    out True on only 48% of them. The two tests below happened to sit on seed 3,
    where the coin landed the right way, so they were passing on luck.
    """
    rng = np.random.default_rng(seed)
    common = rng.normal(size=(200, 1))
    cols = []
    labels = []
    for g in ("A", "B"):
        profile = rng.normal(size=(200, 1))       # shared by the group only
        for r in range(3):
            cols.append(common[:, 0] + profile[:, 0]
                        + rng.normal(scale=0.1, size=200))
            labels.append(f"{g}{r}")
    return np.column_stack(cols), labels


def test_replicates_correlate_more_within_than_between():
    matrix, ids = _matrix()
    groups = {s: s[0] for s in ids}
    res = replicate_concordance(ids, matrix, groups=groups)
    assert res["summary"]["grouping_supported"] is True
    assert res["summary"]["within_group_mean_correlation"] > \
        res["summary"]["between_group_mean_correlation"]


def test_a_swapped_label_lowers_the_margin():
    """Relabelling one replicate must weaken the within-minus-between margin.

    The margin, not the boolean, is what a single swap moves. With three
    replicates per group, moving one sample moves three pairs into ``between``
    and three cross pairs into ``within``; the twelve remaining pairs still
    reflect the true design, so ``grouping_supported`` stays True and *should*
    -- one mislabelled sample among six does not invalidate the grouping, and
    the module deliberately reports a group-level conclusion rather than a
    per-sample verdict.

    Asserting the boolean flipped here is what an earlier version did. It could
    not hold on any fixture with real group structure (measured: True on
    100% of 300 seeds across six designs, before and after the swap), so it was
    passing only because the old fixture had no structure to begin with.
    """
    matrix, ids = _matrix()
    groups = {s: s[0] for s in ids}
    honest = replicate_concordance(ids, matrix, groups=groups)["summary"]

    swapped = dict(groups)
    swapped["A0"] = "B"
    after = replicate_concordance(ids, matrix, groups=swapped)["summary"]

    def margin(summary):
        return (summary["within_group_mean_correlation"]
                - summary["between_group_mean_correlation"])

    # Measured drop is 0.35 on average (min 0.27 over 200 seeds); require a
    # clear fraction of the honest margin so the test states a real effect.
    assert margin(after) < margin(honest) - 0.2 * margin(honest), (honest, after)


def test_region_jaccard_is_reported_when_sets_are_given():
    matrix, ids = _matrix()
    sets = {s: {"r1", "r2", "r3"} for s in ids}
    sets["B0"] = {"r9"}
    res = replicate_concordance(ids, matrix, groups={s: s[0] for s in ids},
                                region_sets=sets)
    w = [p for p in res["pairs"] if p["comparison"] == "within_group"]
    assert all(p["jaccard"] is not None for p in w)
    assert res["summary"]["within_group_mean_jaccard"] is not None


def test_without_regions_jaccard_is_none_not_zero():
    matrix, ids = _matrix()
    pairs = pairwise_concordance(ids, matrix)
    assert all(p["jaccard"] is None for p in pairs)


def test_without_groups_pairs_are_unlabelled():
    matrix, ids = _matrix()
    res = replicate_concordance(ids, matrix)
    assert all(p["comparison"] == "unlabelled" for p in res["pairs"])
    assert "no groups supplied" in res["summary"]["note"]
    assert res["summary"]["grouping_supported"] is None


def test_one_sample_per_group_says_so():
    matrix, ids = _matrix()
    groups = {ids[0]: "A", ids[3]: "B"}
    res = replicate_concordance(ids[:1] + ids[3:4], matrix[:, [0, 3]], groups=groups)
    assert "at least two samples per group" in res["summary"]["note"]


def test_a_single_sample_is_not_an_error():
    matrix, ids = _matrix()
    res = replicate_concordance(ids[:1], matrix[:, :1])
    assert res["pairs"] == []
    assert res["summary"]["within_group_mean_correlation"] is None
