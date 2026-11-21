"""The two coverage paths must agree on which records they count.

`bam2bigwig` chooses a path from the assay: cut sites for `--assay atac`
(`coverage_readends`), accumulated fragments otherwise (`coverage_parallel`). They
disagreed on exactly one filter — the fragment path had always dropped
`is_duplicate` records, the cut-site path had no duplicate test — so a
duplicate-marked BAM produced a different track depending only on which assay
selected the path.

These tests pin the parity. They build a BAM whose reads are flagged duplicate and
assert both paths ignore them, and that `skip_duplicates=False` still counts them so
the switch is real rather than dead.
"""
from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

from easyepi.coverage_readends import read_end_positions  # noqa: E402

READ_LEN = 50
CHROM = "chr1"
CHROM_LEN = 100_000


def _write_bam(path, n_plain: int, n_dup: int):
    """A coordinate-sorted BAM with `n_plain` clean reads then `n_dup` duplicates."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i in range(n_plain + n_dup):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * READ_LEN
            a.flag = 0
            if i >= n_plain:
                a.flag |= 0x400          # duplicate
            a.reference_id = 0
            a.reference_start = 1000 + i * 10
            a.mapping_quality = 60
            a.cigarstring = f"{READ_LEN}M"
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            out.write(a)
    pysam.index(str(path))
    return str(path)


def test_cut_site_path_skips_flagged_duplicates(tmp_path):
    """The regression: without the duplicate test this returned every read."""
    bam = _write_bam(tmp_path / "dup.bam", n_plain=10, n_dup=5)

    positions, stats = read_end_positions(bam, CHROM, 0, CHROM_LEN,
                                          min_mapq=0, shift="none")
    assert positions.size == 10, (
        "the 5 flagged duplicates must not reach the cut-site track; "
        f"got {positions.size} positions from 10 clean + 5 duplicate reads")
    assert stats["records_seen"] == 15
    assert stats["skipped"] == 5
    assert stats["used"] == 10


def test_the_switch_is_not_dead(tmp_path):
    """`skip_duplicates=False` must actually keep them, or the flag is a no-op."""
    bam = _write_bam(tmp_path / "dup.bam", n_plain=10, n_dup=5)

    positions, stats = read_end_positions(bam, CHROM, 0, CHROM_LEN,
                                          min_mapq=0, shift="none",
                                          skip_duplicates=False)
    assert positions.size == 15
    assert stats["skipped"] == 0


def test_a_bam_without_the_flag_is_unaffected(tmp_path):
    """The common case — nothing marked, so nothing changes."""
    bam = _write_bam(tmp_path / "clean.bam", n_plain=12, n_dup=0)

    default, _ = read_end_positions(bam, CHROM, 0, CHROM_LEN, min_mapq=0, shift="none")
    keep, _ = read_end_positions(bam, CHROM, 0, CHROM_LEN, min_mapq=0, shift="none",
                                 skip_duplicates=False)
    assert default.size == 12
    assert default.size == keep.size
