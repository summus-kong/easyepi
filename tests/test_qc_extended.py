"""Unit tests for the extended QC metrics and the threshold verdicts.

Covers the edge cases that matter for batch QC: peaks with no reads, an empty or
missing blacklist, a BAM with no unmapped records, and threshold payloads whose
metrics are absent rather than zero.

Synthetic inputs are built with pysam/pyBigWig into pytest's tmp_path, so the
tests are self-contained and leave nothing behind.
"""

import json

import numpy as np
import pytest

pytest.importorskip("pysam")
pytest.importorskip("pyBigWig")

import pysam  # noqa: E402
import pyBigWig  # noqa: E402

from easyepi.qc_extended import (  # noqa: E402
    blacklist_fraction,
    compute_extended_metrics,
    peak_read_distribution,
    read_flag_metrics,
    tss_enrichment,
)
from easyepi.qc_thresholds import (  # noqa: E402
    evaluate_thresholds,
    annotate_rows_with_thresholds,
    summarize_verdicts,
)

CHROM = "chr1"
LENGTH = 200_000


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _write_bam(path, records, chrom=CHROM, length=LENGTH, mito=None):
    """Write a coordinate-sorted BAM.

    ``records`` items are ``(start, isize, flag_extra, ref_id)`` where ``ref_id``
    0 is the main contig and 1 the optional mitochondrial contig.
    """
    sq = [{"SN": chrom, "LN": length}]
    if mito:
        sq.append({"SN": mito, "LN": 16_569})
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    ordered = sorted(records, key=lambda r: (r[3], r[0]))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (start, isize, extra_flag, ref_id) in enumerate(ordered):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.query_sequence = "A" * 75
            a.query_qualities = pysam.qualitystring_to_array("I" * 75)
            a.flag = (99 if i % 2 == 0 else 147) | extra_flag
            a.reference_id = ref_id
            a.reference_start = start
            a.mapping_quality = 60
            a.cigarstring = "75M"
            a.template_length = isize if i % 2 == 0 else -isize
            out.write(a)
    pysam.index(str(path))


def _write_bed(path, intervals, chrom=CHROM):
    with open(path, "w", encoding="utf-8") as fh:
        for start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\tregion\n")


def _write_gtf(path, tss_positions, chrom=CHROM):
    with open(path, "w", encoding="utf-8") as fh:
        for i, pos in enumerate(tss_positions):
            fh.write(f'{chrom}\tsim\tgene\t{pos + 1}\t{pos + 500}\t.\t+\t.\tgene_id "G{i}";\n')


def _write_bigwig(path, values_by_base, chrom=CHROM):
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(chrom, len(values_by_base))])
    starts = list(range(len(values_by_base)))
    bw.addEntries([chrom] * len(values_by_base), starts,
                  ends=[s + 1 for s in starts],
                  values=[float(v) for v in values_by_base])
    bw.close()


# ---------------------------------------------------------------------------
# peak_read_distribution
# ---------------------------------------------------------------------------
def test_peak_distribution_counts_reads_and_zero_peaks(tmp_path):
    bam = tmp_path / "reads.bam"
    # peak A gets 4 reads, peak B gets 0, peak C gets 1
    _write_bam(bam, [
        (1000, 150, 0, 0), (1010, 150, 0, 0), (1020, 150, 0, 0), (1030, 150, 0, 0),
        (5100, 150, 0, 0),
        (20_000, 150, 0, 0),  # outside every peak (peak C stays unsupported)
    ])
    bed = tmp_path / "peaks.bed"
    _write_bed(bed, [(950, 1150), (4950, 5150), (8950, 9150)])

    res = peak_read_distribution(str(bam), str(bed))
    assert res["n_peaks"] == 3
    assert res["peaks_with_zero_reads"] == 1
    assert res["fraction_peaks_zero"] == pytest.approx(1 / 3)
    assert res["max_reads_per_peak"] == 4
    assert res["median_reads_per_peak"] == 1.0
    # 5 of 6 reads land in peaks (4 in peak A + 1 in peak B)
    assert res["frip"] == pytest.approx(5 / 6)
    assert len(res["hist_counts"]) == len(res["hist_edges"]) - 1


def test_peak_distribution_all_zero_peaks(tmp_path):
    """Every peak unsupported: still a valid result, not a crash."""
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [(50_000, 150, 0, 0)])
    bed = tmp_path / "peaks.bed"
    _write_bed(bed, [(1000, 1200), (3000, 3200)])
    res = peak_read_distribution(str(bam), str(bed))
    assert res["peaks_with_zero_reads"] == 2
    assert res["fraction_peaks_zero"] == 1.0
    assert res["max_reads_per_peak"] == 0
    assert res["top_decile_read_share"] is None  # no reads in peaks at all


def test_peak_distribution_empty_bed(tmp_path):
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [(1000, 150, 0, 0)])
    bed = tmp_path / "peaks.bed"
    _write_bed(bed, [])
    assert peak_read_distribution(str(bam), str(bed))["n_peaks"] == 0


# ---------------------------------------------------------------------------
# blacklist_fraction
# ---------------------------------------------------------------------------
def test_blacklist_fraction_counts_overlapping_reads(tmp_path):
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [
        (1000, 150, 0, 0), (1010, 150, 0, 0),   # inside blacklist
        (50_000, 150, 0, 0), (60_000, 150, 0, 0),  # outside
    ])
    bed = tmp_path / "blacklist.bed"
    _write_bed(bed, [(950, 1200)])
    res = blacklist_fraction(str(bam), str(bed))
    assert res["reads_inspected"] == 4
    assert res["reads_in_blacklist"] == 2
    assert res["blacklist_fraction"] == pytest.approx(0.5)


def test_blacklist_missing_file_returns_error(tmp_path):
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [(1000, 150, 0, 0)])
    res = blacklist_fraction(str(bam), str(tmp_path / "does_not_exist.bed"))
    assert "error" in res
    assert "blacklist_fraction" not in res


def test_blacklist_no_reads_inside(tmp_path):
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [(50_000, 150, 0, 0)])
    bed = tmp_path / "blacklist.bed"
    _write_bed(bed, [(1000, 1200)])
    res = blacklist_fraction(str(bam), str(bed))
    assert res["reads_in_blacklist"] == 0
    assert res["blacklist_fraction"] == 0.0


# ---------------------------------------------------------------------------
# read_flag_metrics / tss_enrichment
# ---------------------------------------------------------------------------
def test_read_flag_metrics_mito_and_duplicates(tmp_path):
    bam = tmp_path / "reads.bam"
    _write_bam(
        bam,
        [(1000, 150, 0, 0), (1100, 150, 0, 0),
         (1200, 150, 1024, 0),          # duplicate flag
         (500, 150, 0, 1)],
        mito="chrM",
    )
    res = read_flag_metrics(str(bam))
    assert res["reads_inspected"] == 4
    assert res["mitochondrial_reads"] == 1
    assert res["mitochondrial_fraction"] == pytest.approx(0.25)
    assert res["duplicate_rate"] == pytest.approx(0.25)
    assert res["duplicate_source"] == "duplicate_flag"
    assert res["median_read_length"] == 75


def test_read_flag_metrics_without_a_mito_contig_reports_undefined(tmp_path):
    """No matching contig is `None`, not 0.0.

    This fixture has only `chr1`. Reporting 0.0 here would claim the library is free
    of mitochondrial reads when the question was never answered -- the failure mode
    this metric exists to avoid. The status carries the reason instead.
    """
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [(1000, 150, 0, 0)])
    res = read_flag_metrics(str(bam))
    assert res["reads_inspected"] == 1
    assert res["mitochondrial_fraction"] is None
    assert res["mitochondrial_status"] == "no_contig_matched"


def test_tss_enrichment_detects_signal(tmp_path):
    bigwig = tmp_path / "signal.bw"
    n = 20_000
    values = np.zeros(n, dtype=float)
    tss = [5000, 12_000]
    for pos in tss:
        values[pos - 100:pos + 100] = 10.0   # strong TSS signal
    values[values == 0] = 0.5                # weak background
    _write_bigwig(bigwig, values)
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf, tss)
    res = tss_enrichment(str(bigwig), str(gtf))
    assert res["n_tss_used"] == 2
    assert res["tss_enrichment"] > 5.0


def test_tss_enrichment_missing_annotation(tmp_path):
    bigwig = tmp_path / "signal.bw"
    _write_bigwig(bigwig, np.full(1000, 1.0))
    res = tss_enrichment(str(bigwig), str(tmp_path / "missing.gtf"))
    assert "error" in res


def test_bigwig_tss_publishes_its_own_key(tmp_path):
    """The BigWig reader must not publish the shared `tss_enrichment` key.

    It scales differently from the BAM reader (+-100 bp over the 1900-2000 bp ring
    versus +-25 bp over the outer 100 bp), so a shared key made the number depend
    on call order. See `test_a_real_annotation_does_not_let_both_readers_overwrite_each_other`.
    """
    bigwig = tmp_path / "signal.bw"
    values = np.zeros(20_000, dtype=float)
    tss = [5000, 12_000]
    for pos in tss:
        values[pos - 100:pos + 100] = 10.0
    values[values == 0] = 0.5
    _write_bigwig(bigwig, values)
    gtf = tmp_path / "genes.gtf"
    _write_gtf(gtf, tss)

    out = compute_extended_metrics("ATAC-seq", bigwig=str(bigwig),
                                   annotation=str(gtf))
    assert out["tss_enrichment_bigwig"] > 5.0
    assert "tss_enrichment" not in out


def test_one_finding_is_not_scored_twice_when_both_tss_readers_run():
    """Both readers present must not count as two pieces of evidence."""
    both = evaluate_thresholds({"tss_enrichment": 1220.0,
                                "tss_enrichment_bigwig": 42.4})
    checks = {c["metric"]: c for c in both["checks"]}
    assert checks["tss_enrichment"]["status"] == "pass"
    assert checks["tss_enrichment_bigwig"]["status"] == "not_applicable"

    # With only the BigWig reading, the BigWig rule is the one that applies.
    only_bw = evaluate_thresholds({"tss_enrichment_bigwig": 1.2})
    checks = {c["metric"]: c for c in only_bw["checks"]}
    assert checks["tss_enrichment_bigwig"]["status"] == "breached"
    assert checks["tss_enrichment"]["status"] == "not_evaluated"


# ---------------------------------------------------------------------------
# threshold verdicts
# ---------------------------------------------------------------------------
def test_threshold_pass_and_warn():
    # "pass" now requires the core ChIP metrics to be present, not just any metric
    # in range: a sheet missing FRiP/enrichment/duplication is unassessed.
    ok = evaluate_thresholds({"frip": 0.6, "signal_enrichment": 5.0,
                              "duplicate_rate": 0.1, "bam_mapping_rate": 0.99,
                              "tss_enrichment": 8.0, "blacklist_fraction": 0.001})
    assert ok["overall"] == "pass"
    assert ok["warnings"] == []
    assert ok["missing_core_metrics"] == []

    bad = evaluate_thresholds({"frip": 0.05, "blacklist_fraction": 0.05})
    assert bad["overall"] == "warn"
    assert "frip" in bad["warnings"]
    assert "blacklist_fraction" in bad["warnings"]


def test_threshold_missing_metric_is_not_evaluated():
    res = evaluate_thresholds({})
    assert res["overall"] == "not_evaluated"
    statuses = {c["metric"]: c["status"] for c in res["checks"]}
    assert statuses["frip"] == "not_evaluated"


def test_threshold_nested_payload_is_flattened():
    res = evaluate_thresholds({"peak_reads": {"fraction_peaks_zero": 0.8}})
    flagged = {c["metric"] for c in res["checks"] if c["status"] == "breached"}
    assert "fraction_peaks_zero" in flagged


def test_threshold_boundary_is_inclusive():
    # exactly at the threshold must NOT breach (strict comparison)
    assert evaluate_thresholds({"frip": 0.20})["overall"] != "warn"
    assert evaluate_thresholds({"blacklist_fraction": 0.01})["overall"] != "warn"


def test_nan_metric_is_ignored():
    res = evaluate_thresholds({"frip": float("nan")})
    statuses = {c["metric"]: c["status"] for c in res["checks"]}
    assert statuses["frip"] == "not_evaluated"


def test_annotate_and_summarize_rows():
    rows = [
        {"sample_id": "s1", "assay_qc": {"frip": 0.5, "signal_enrichment": 5.0,
                                         "duplicate_rate": 0.1,
                                         "bam_mapping_rate": 0.99,
                                         "tss_enrichment": 6.0}},
        {"sample_id": "s2", "assay_qc": {"frip": 0.02, "signal_enrichment": 5.0,
                                         "duplicate_rate": 0.1,
                                         "bam_mapping_rate": 0.99}},
    ]
    annotate_rows_with_thresholds(rows)
    assert rows[0]["qc_verdict"]["overall"] == "pass"
    assert rows[1]["qc_verdict"]["overall"] == "warn"
    summary = summarize_verdicts(rows)
    assert summary["samples"]["pass"] == 1
    assert summary["samples"]["warn"] == 1
    assert summary["warning_metrics"]["frip"] == 1


def test_custom_rules_override_table():
    # Thresholds are advisory: even a rule marked "fail" only warns.
    rules = {"frip": {"op": "lt", "value": 0.5, "severity": "fail",
                      "message": "custom"}}
    res = evaluate_thresholds({"frip": 0.3}, rules=rules)
    assert res["overall"] == "warn"
    assert res["warnings"] == ["frip"]
    assert res["failures"] == []


# ---------------------------------------------------------------------------
# assay name normalisation
# ---------------------------------------------------------------------------
_ALIAS = {
    "chipseq": "chip-qc", "chip": "chip-qc", "cutandtag": "chip-qc",
    "cutrun": "chip-qc", "cutandrun": "chip-qc",
    "atacseq": "atac-qc", "atac": "atac-qc",
    "histone": "histone-qc", "histonemark": "histone-qc", "histoneqc": "histone-qc",
    "mnaseseq": "mnase-qc", "mnase": "mnase-qc",
    "methylation": "methylation-qc", "methylseq": "methylation-qc",
    "wgbs": "methylation-qc", "rrbs": "methylation-qc",
    "bisulfite": "methylation-qc",
}


def _norm(name):
    return name.lower().replace("-", "").replace("_", "").replace(" ", "").replace("&", "and")


@pytest.mark.parametrize("name", [
    "ChIP-seq", "ChIP", "chip", "CUT&Tag", "CUT&RUN",
    "ATAC-seq", "ATAC", "atac",
])
def test_assay_alias_resolves_to_chip_or_atac(name):
    """Regression: assay="ChIP" used to resolve to nothing and silently drop FRiP."""
    from easyepi.assay_metrics import ASSAY_METRIC_BUILDERS

    resolved = _ALIAS.get(_norm(name), name.lower())
    assert resolved in ASSAY_METRIC_BUILDERS, f"{name!r} -> {resolved!r} missing"
    assert "frip" in ASSAY_METRIC_BUILDERS[resolved]


@pytest.mark.parametrize("name,expected", [
    ("histone", "histone-qc"), ("Histone-mark", "histone-qc"),
    ("MNase-seq", "mnase-qc"), ("MNase", "mnase-qc"),
    ("methylation", "methylation-qc"), ("WGBS", "methylation-qc"),
    ("RRBS", "methylation-qc"),
])
def test_assay_alias_for_other_assays(name, expected):
    from easyepi.assay_metrics import ASSAY_METRIC_BUILDERS

    resolved = _ALIAS.get(_norm(name), name.lower())
    assert resolved == expected
    assert resolved in ASSAY_METRIC_BUILDERS


def test_unknown_assay_returns_empty_not_error():
    from easyepi.assay_metrics import compute_assay_metrics

    assert compute_assay_metrics("RNA-seq") == {}
