"""Single-end and paired-end libraries must yield the same per-fragment metrics.

Both fixtures are built from the SAME fragment list, so a correct implementation
has no licence to differ between them. Each assertion here pins one of the three
places where a library layout could leak into a number:

* coverage bodies (paired-end must span the fragment, not just read1),
* cut-site counting (paired-end must contribute one cut site per fragment),
* fragment-length summary (paired-end is defined; single-end is not applicable).
"""
from __future__ import annotations

import pysam
import pytest

from easyepi.assay_metrics import fragment_length_summary
from easyepi.coverage_parallel import CoverageSpec, _accumulate_range

READ_LEN = 50
FRAGMENTS = [(10_000, 150), (10_500, 200), (11_000, 180), (11_500, 150),
             (20_000, 160), (20_500, 190)]


def _header() -> dict:
    return {"HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 100_000}]}


def _write(tmp_path, paired: bool):
    path = tmp_path / ("pe.bam" if paired else "se.bam")
    with pysam.AlignmentFile(str(path), "wb", header=_header()) as out:
        for start, length in FRAGMENTS:
            end = start + length
            if paired:
                for is_read1, flag_extra, pos, tlen in (
                    (True, 0x40 | 0x20, start, length),
                    (False, 0x80 | 0x10, end - READ_LEN, -length),
                ):
                    a = pysam.AlignedSegment()
                    a.query_name = f"f{start}"
                    a.query_sequence = "A" * READ_LEN
                    a.flag = 0x1 | 0x2 | flag_extra
                    a.reference_id = 0
                    a.reference_start = pos
                    a.mapping_quality = 60
                    a.cigartuples = [(0, READ_LEN)]
                    a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
                    a.next_reference_id = 0
                    a.next_reference_start = end - READ_LEN if is_read1 else start
                    a.template_length = tlen
                    out.write(a)
            else:
                a = pysam.AlignedSegment()
                a.query_name = f"f{start}"
                a.query_sequence = "A" * length
                a.flag = 0
                a.reference_id = 0
                a.reference_start = start
                a.mapping_quality = 60
                a.cigartuples = [(0, length)]
                a.query_qualities = pysam.qualitystring_to_array("I" * length)
                out.write(a)
    pysam.index(str(path))
    return path


def _coverage(path, shift: str):
    spec = CoverageSpec(bam_path=str(path), chrom="chr1", chrom_length=100_000,
                        bin_size=1000, min_mapq=30, shift=shift,
                        max_depth=None, chunk_bp=10_000)
    return _accumulate_range(spec, 0, 50_000)


def test_paired_end_body_matches_single_end_bytes(tmp_path):
    """A pair is one fragment: PE coverage must equal the SE coverage body."""
    import numpy as np

    se = _coverage(_write(tmp_path, False), "none")
    pe = _coverage(_write(tmp_path, True), "none")
    assert np.array_equal(se, pe)


def test_paired_end_contributes_one_cut_site_per_fragment(tmp_path):
    """Counting both mates would double the cut-site signal."""
    se = _coverage(_write(tmp_path, False), "tn5")
    pe = _coverage(_write(tmp_path, True), "tn5")
    assert se.sum() == pe.sum() == len(FRAGMENTS)


def test_fragment_length_is_defined_for_pairs_and_flagged_for_single_end(tmp_path):
    pe = fragment_length_summary(str(_write(tmp_path, True)), max_reads=10_000)
    assert pe["n_fragments"] == len(FRAGMENTS)
    lengths = sorted(length for _, length in FRAGMENTS)
    n = len(lengths)
    expected = (lengths[n // 2 - 1] + lengths[n // 2]) / 2
    assert pe["median_fragment_length"] == pytest.approx(expected)

    se = fragment_length_summary(str(_write(tmp_path, False)), max_reads=10_000)
    assert se["n_fragments"] == 0
    assert se["status"] == "not_applicable"
