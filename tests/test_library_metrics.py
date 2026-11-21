"""Pin the three failure modes `library_summary` has already produced once.

Each test below exists because the implementation got it wrong at least once, and
each wrong answer was a plausible-looking number rather than an error:

1. **chrM sorts last.** The first version counted mitochondrial reads while scanning
   records in file order and stopped at `max_reads`. On a coordinate-sorted BAM chrM
   is at the end, so the scan never reached it and a contaminated library reported
   0.0 -- the most dangerous possible answer, because 0 reads as clean.
2. **Numerator and denominator from different populations.** The fix moved the
   numerator to the index but kept the truncated scan total as the denominator, so
   the fraction came out ten-fold high (9026/200000 vs 9026/95302588).
3. **An unmarked BAM looks duplicate-free.** `get_duplicate_flag` is only set if
   something marked duplicates; without that step every read reads as unique, and a
   bare 0.0 claims a clean library that was never assessed.

The fixtures are tiny on purpose: the bugs are about ordering and population, not
magnitude, so two contigs and a handful of reads reproduce all three.
"""
from __future__ import annotations

import pysam
import pytest

from easyepi.assay_metrics import library_summary

READ_LEN = 50
MITO_PATTERN = r"^(chr)?M(T)?$"

#: chr1 reads come first so that chrM is last in a coordinate-sorted file, which is
#: exactly the layout the original scan could not reach.
N_CHR1 = 40
N_CHRM = 10


def _header() -> dict:
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1_000_000},
               {"SN": "chrM", "LN": 16_569}],
    }


def _write(tmp_path, *, name: str = "lib.bam", duplicates: bool = False,
           paired: bool = False):
    """chr1 reads first, chrM last, so chrM is unreachable by a truncated scan."""
    path = tmp_path / name
    with pysam.AlignmentFile(str(path), "wb", header=_header()) as out:
        for i in range(N_CHR1 + N_CHRM):
            on_mito = i >= N_CHR1
            for mate in ((True, False) if paired else (None,)):
                a = pysam.AlignedSegment()
                a.query_name = f"r{i}"
                a.query_sequence = "A" * READ_LEN
                flag = 0
                if paired:
                    flag |= 0x1 | 0x2 | (0x40 if mate else 0x80)
                    a.next_reference_id = 1 if on_mito else 0
                    a.next_reference_start = 100 + (i * 200) + (0 if mate else 100)
                    a.template_length = 100 if mate else -100
                if duplicates:
                    flag |= 0x400  # duplicate
                a.flag = flag
                a.reference_id = 1 if on_mito else 0
                a.reference_start = 100 + (i * 200) + (0 if (mate in (True, None)) else 100)
                a.mapping_quality = 60
                a.cigartuples = [(0, READ_LEN)]
                a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
                out.write(a)
    pysam.index(str(path))
    return str(path)


def test_mito_is_counted_when_chr_m_sorts_last(tmp_path):
    """A truncated scan must not be able to hide chrM.

    With ``max_reads`` far below the record count, a file-order scan sees only chr1
    and reports 0.0. The index-based count must still find every chrM read.
    """
    path = _write(tmp_path)
    res = library_summary(path, max_reads=5, mito_pattern=MITO_PATTERN)

    assert res["mitochondrial_status"] == "measured"
    # The bug produced 0 here; the truth is the full chrM count.
    assert res["mitochondrial_reads"] == N_CHRM
    assert res["mitochondrial_fraction"] > 0


def test_mito_numerator_and_denominator_share_one_population(tmp_path):
    """The fraction must equal chrM over every mapped record, not over the sample.

    Dividing the whole-file numerator by the truncated denominator inflated the
    result ten-fold on real data. Here the sample is 5 reads while the file has 50,
    so a mixed-population ratio would be 10/5 = 2.0 (impossible, >1) and a shared
    one is 10/50 = 0.2.
    """
    path = _write(tmp_path)
    res = library_summary(path, max_reads=5, mito_pattern=MITO_PATTERN)

    # Numerator and denominator both come from the index, so the answer is
    # chrM over every mapped record -- not over the 5 reads this call inspected.
    # Expecting the sample ratio here would re-encode the mixed-population bug.
    assert res["mitochondrial_fraction"] == pytest.approx(N_CHRM / (N_CHR1 + N_CHRM))
    assert res["n_reads"] == 5
    assert 0.0 < res["mitochondrial_fraction"] < 1.0
    assert res["mitochondrial_denominator"] == "index_mapped_records"
    assert res["mitochondrial_mapped_records"] == N_CHR1 + N_CHRM


def test_unmarked_bam_does_not_claim_a_duplicate_free_library(tmp_path):
    """0.0 with no marking step must be labelled, not presented as a measurement."""
    path = _write(tmp_path)
    res = library_summary(path, max_reads=None)

    assert res["duplicate_fraction"] == 0.0
    assert res["duplicate_status"] == "not_marked"
    assert "duplicate_note" in res


def test_marked_duplicates_are_counted(tmp_path):
    """The other side of the same coin: a marked file must report measured."""
    path = _write(tmp_path, name="dup.bam", duplicates=True)
    res = library_summary(path, max_reads=None)

    assert res["duplicate_status"] == "measured"
    assert res["duplicate_fraction"] == pytest.approx(1.0)
    assert res["n_duplicates"] == N_CHR1 + N_CHRM


def test_complexity_counts_one_template_per_pair(tmp_path):
    """Both mates of a pair share a name, so templates are half the records.

    Counting per alignment would make the ratio layout-dependent and put
    single-end and paired-end libraries on different scales.
    """
    path = _write(tmp_path, name="pe.bam", paired=True)
    res = library_summary(path, max_reads=None)

    total_records = 2 * (N_CHR1 + N_CHRM)
    assert res["n_reads"] == total_records
    assert res["n_templates"] == N_CHR1 + N_CHRM
    assert res["library_complexity"] == pytest.approx(0.5)


def test_mito_pattern_matching_nothing_is_not_zero(tmp_path):
    """A pattern that matches no contig means "unknown", not "no mitochondria"."""
    path = _write(tmp_path)
    res = library_summary(path, max_reads=None, mito_pattern=r"^chrMT_only$")

    assert res["mitochondrial_fraction"] is None
    assert res["mitochondrial_status"] == "no_contig_matched"


def test_no_pattern_is_not_applicable(tmp_path):
    """Absent a pattern the contig cannot be identified, so no number is invented."""
    path = _write(tmp_path)
    res = library_summary(path, max_reads=None)

    assert res["mitochondrial_fraction"] is None
    assert res["mitochondrial_status"] == "not_applicable"


def test_invalid_pattern_is_reported_not_raised(tmp_path):
    """A bad regex must come back as data, matching how unreadable BAMs are handled."""
    path = _write(tmp_path)
    res = library_summary(path, max_reads=None, mito_pattern="([")

    assert "error" in res
