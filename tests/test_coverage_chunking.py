"""Chunk boundaries must fall on the bin grid.

`_accumulate_range` bins from each chunk's own start, while `_place` maps a chunk to
`first = start // bin_size`. If a chunk began mid-bin the two disagreed and coverage
was written to the wrong bins. Measured on real chr22 data before the fix:

    bin_size=1000 (divides CHUNK_BP): 0 bins wrong
    bin_size=128  (does not divide) : 53,250 bins wrong, 33.2% of total signal
    bin_size=333  (does not divide) : 32,118 bins wrong, 14.8% of total signal

Total signal was conserved in every case, so the defect silently redistributed
coverage rather than losing it. These tests pin the invariant that makes the two
functions agree, without needing a BAM.
"""

from __future__ import annotations

import pytest

from easyepi.coverage_parallel import CoverageSpec, _chunk_tasks

CHROM_LENGTH = 2_500_000  # deliberately not a multiple of 1 Mb


def spec(bin_size: int) -> CoverageSpec:
    return CoverageSpec(
        bam_path="unused.bam",
        chrom="chr1",
        chrom_length=CHROM_LENGTH,
        bin_size=bin_size,
        min_mapq=30,
        shift=None,
        max_depth=None,
    )


@pytest.mark.parametrize("bin_size", [50, 1000, 128, 333, 7])
def test_every_chunk_start_is_a_multiple_of_bin_size(bin_size):
    """The alignment invariant: chunk starts must sit on the bin grid."""
    tasks = _chunk_tasks(spec(bin_size), 1_000_000)
    for _spec, start, _end in tasks:
        assert start % bin_size == 0, (bin_size, start)


@pytest.mark.parametrize("bin_size", [50, 1000, 128, 333])
def test_chunks_cover_the_chromosome_without_gaps(bin_size):
    """Contiguity must survive rounding the step down to a bin multiple."""
    tasks = _chunk_tasks(spec(bin_size), 1_000_000)
    assert tasks[0][1] == 0
    for prev, nxt in zip(tasks, tasks[1:]):
        assert prev[2] == nxt[1], (prev, nxt)
    assert tasks[-1][2] == CHROM_LENGTH


def test_step_still_advances_when_bin_size_exceeds_chunk():
    """A bin larger than the chunk request must not produce a zero-length step."""
    s = spec(5000)
    tasks = _chunk_tasks(s, 1000)  # bin_size > chunk_bp
    starts = [t[1] for t in tasks]
    assert starts == sorted(set(starts)), "duplicate or unordered chunk starts"
    assert starts[0] == 0
    assert len(tasks) < CHROM_LENGTH, "step collapsed to zero"
    for _spec, start, _end in tasks:
        assert start % 5000 == 0
