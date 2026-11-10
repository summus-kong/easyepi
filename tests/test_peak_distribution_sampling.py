"""Peak-read distribution must stride the whole BAM, not the first chromosome.

Same failure mode as FRiP: a coordinate-sorted BAM truncated at N reads leaves
every peak on later chromosomes at zero, so ``top_decile_read_share`` collapses to
1.0. On a real ENCODE H3K4me3 file that produced a meaningless QC row.
"""

from __future__ import annotations

import pytest

pytest.importorskip("pysam")

import pysam  # noqa: E402

from easyepi.qc_extended import peak_read_distribution  # noqa: E402

CHROMS = [("chr1", 200_000), ("chr2", 200_000), ("chr3", 200_000)]
READ_LEN = 50


def _build(tmp_path, per_peak=20, n_peaks_per_chrom=20):
    """Peaks on every chromosome, all equally supported."""
    peaks = {}
    for ref_id, (chrom, size) in enumerate(CHROMS):
        spans = []
        for i in range(n_peaks_per_chrom):
            s = 1000 + i * 2000
            spans.append((s, s + 500))
        peaks[ref_id] = spans

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": n} for c, n in CHROMS]}
    rows = []
    for ref_id, spans in peaks.items():
        for s, e in spans:
            for k in range(per_peak):
                rows.append((ref_id, s + k * (e - s) // max(1, per_peak)))
    rows.sort()

    bam = tmp_path / "s.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for idx, (ref, start) in enumerate(rows):
            rec = pysam.AlignedSegment(out.header)
            rec.query_name = f"r{idx}"
            rec.query_sequence = "A" * READ_LEN
            rec.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            rec.flag = 0
            rec.reference_id = ref
            rec.reference_start = start
            rec.mapping_quality = 60
            rec.cigarstring = f"{READ_LEN}M"
            out.write(rec)
    pysam.index(str(bam))

    bed = tmp_path / "p.bed"
    with open(bed, "w", encoding="utf-8") as fh:
        n = 0
        for ref_id, (chrom, _size) in enumerate(CHROMS):
            for s, e in peaks[ref_id]:
                fh.write(f"{chrom}\t{s}\t{e}\tpk{n}\n")
                n += 1

    return bam, bed, sum(len(v) for v in peaks.values())


def test_exhaustive_distribution_is_uniform(tmp_path):
    bam, bed, n_peaks = _build(tmp_path)
    res = peak_read_distribution(str(bam), str(bed), max_reads=10**9)

    assert res["n_peaks"] == n_peaks
    assert res["stride"] == 1
    assert res["peaks_with_zero_reads"] == 0
    # every peak carries the same count, so the top decile holds ~10%
    assert res["top_decile_read_share"] == pytest.approx(0.1, abs=0.05)


def test_truncated_scan_still_reaches_later_chromosomes(tmp_path):
    """The bug produced top_decile_read_share == 1.0 by starving late peaks."""
    bam, bed, _n = _build(tmp_path)

    full = peak_read_distribution(str(bam), str(bed), max_reads=10**9)
    trunc = peak_read_distribution(str(bam), str(bed), max_reads=200)

    assert trunc["stride"] > 1
    assert trunc["reads_scanned"] > trunc["reads_inspected"]
    # The fix keeps the concentration realistic instead of collapsing to 1.0.
    assert trunc["top_decile_read_share"] < 0.9
    assert trunc["fraction_peaks_zero"] < 0.5
    assert trunc["top_decile_read_share"] == pytest.approx(
        full["top_decile_read_share"], abs=0.35)


def test_distribution_reports_scan_provenance(tmp_path):
    bam, bed, _n = _build(tmp_path)
    res = peak_read_distribution(str(bam), str(bed), max_reads=300)
    for key in ("reads_scanned", "reads_inspected", "stride",
                "top_decile_read_share", "fraction_peaks_zero"):
        assert key in res
