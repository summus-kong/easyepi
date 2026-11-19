"""The un-indexed CPM denominator must not collapse to the mapped count.

`count_mapped` extrapolates the total read count from a sample when the BAM has
no index. A previous version computed

    total = int(round(seen * fraction)) if mapped == seen else 0
    return {"mapped": mapped, "total": total or mapped, ...}

`mapped == seen` is false for any file with even one unmapped read, so `total`
became 0 and `total or mapped` silently substituted the mapped count. CPM then
divided by the numerator's own magnitude, disabling depth normalisation exactly
where it matters. These tests pin the extrapolation down.
"""
from __future__ import annotations

import pysam
import pytest

from easyepi.bam2bigwig import count_mapped

CHROM = "chr1"
CHROM_LEN = 100_000
READ_LEN = 50


def _write(path, n_mapped: int, n_unmapped: int) -> str:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(n_mapped):
            a = pysam.AlignedSegment()
            a.query_name = f"m{i}"
            a.query_sequence = "A" * READ_LEN
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 100 + (i * 37) % (CHROM_LEN - 500)
            a.mapping_quality = 60
            a.cigartuples = [(0, READ_LEN)]
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            out.write(a)
        for i in range(n_unmapped):
            u = pysam.AlignedSegment()
            u.query_name = f"u{i}"
            u.query_sequence = "A" * READ_LEN
            u.flag = 0x4
            u.reference_id = 0
            u.reference_start = 0
            u.mapping_quality = 0
            u.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            out.write(u)
    return str(path)


def test_denominator_equals_mapped_when_everything_mapped(tmp_path):
    info = count_mapped(_write(tmp_path / "all.bam", 500, 0), chrom=None)
    assert info["mapped"] == 500
    assert info["total"] == 500
    assert info["method"] == "sampled_extrapolated"


def test_denominator_grows_with_the_unmapped_fraction(tmp_path):
    """total must exceed mapped: that is the whole point of a denominator."""
    info = count_mapped(_write(tmp_path / "some.bam", 500, 214), chrom=None)
    assert info["mapped"] == 500
    # 714 reads seen, 500 mapped -> estimated total 714, not 500.
    assert info["total"] == pytest.approx(714, abs=1)
    assert info["total"] / info["mapped"] == pytest.approx(1.428, abs=0.01)


def test_denominator_never_falls_below_the_mapped_count(tmp_path):
    info = count_mapped(_write(tmp_path / "mix.bam", 500, 20), chrom=None)
    assert info["total"] >= info["mapped"]


def test_per_chromosome_request_is_reported_as_unsupported(tmp_path):
    """The sampled fraction describes the whole file, so a per-chrom total is
    not extrapolable; it must be flagged rather than silently guessed."""
    info = count_mapped(_write(tmp_path / "chrom.bam", 500, 20), chrom=CHROM)
    assert info["method"] == "sampled_per_chrom_unsupported"
    assert info["total"] == 0
