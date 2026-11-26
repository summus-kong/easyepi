"""Pin the cut-site port against TOBIAS, using synthetic reads for the soft-clip case.

Measured facts this test encodes:

* `ENCFF327MYY.bam` contains no soft-clipped reads in chr1:1-4 Mb, so the
  `query_length - query_alignment_end` correction never fires on that data.
  The real-data comparison therefore covers the plain branch only.
* The reference alphabet is A=0, T=1, C=2, G=3 (see `GenomicSequence.from_fasta`),
  the complement table being 0<->1, 2<->3. Getting this wrong still yields
  correct cut sites but permutes every k-mer.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_cutsites import (
    N_CODE, cutsite, encode_sequence, kmer_at, revcomp_codes,
)

FASTA = "/data1/GaoLab/houoweibo/reference/human/genome.fa"
BAM = "/data1/GaoLab/houoweibo/easyepi_work/data_atac/ENCFF327MYY.bam"


def _real():
    return pytest.importorskip("pysam")


def test_alphabet_is_alphabetical_not_acgt():
    """A=0, T=1, C=2, G=3: the permutation the reference uses."""
    assert encode_sequence("ATCG").tolist() == [0, 1, 2, 3]
    assert encode_sequence("N").tolist() == [N_CODE]
    assert encode_sequence("atcg").tolist() == [0, 1, 2, 3]


def test_revcomp_uses_the_reference_complement_table():
    """A<->T and C<->G in that alphabet; the sentinel is preserved."""
    # reversed first, then complemented: [4,3,2,1,0] -> [4,2,3,0,1]
    assert revcomp_codes(np.array([0, 1, 2, 3, 4])).tolist() == [4, 2, 3, 0, 1]
    # revcomp of revcomp is the original
    ids = np.arange(5)
    assert np.array_equal(revcomp_codes(revcomp_codes(ids)), ids)


def test_forward_cutsite_without_softclip():
    """No soft clip: reference_start + 1 + 4."""
    assert cutsite(1000, 1076, 0, 76, 76, is_reverse=False) == 1005


def test_reverse_cutsite_without_softclip():
    """No soft clip: reference_end + 1 - 5."""
    assert cutsite(1000, 1076, 0, 76, 76, is_reverse=True) == 1072


def test_soft_clip_shifts_the_cut_site():
    """A trailing soft clip moves the reverse cut site by its length.

    `positions[1] = reference_end + query_length - query_alignment_end`, so a
    10 bp trailing clip that is not in the alignment adds 10 before the shift.
    """
    no_clip = cutsite(1000, 1076, 0, 76, 76, is_reverse=True)
    with_clip = cutsite(1000, 1076, 0, 76, 86, is_reverse=True)
    assert with_clip - no_clip == 10

    # leading clip affects the forward branch symmetrically
    fwd_no = cutsite(1000, 1076, 0, 76, 76, is_reverse=False)
    fwd_clip = cutsite(1000, 1076, 10, 76, 86, is_reverse=False)
    assert fwd_no - fwd_clip == 10


def test_kmer_matches_reference_on_real_reads():
    """End-to-end on real BAM reads: cut sites and k-mers both agree."""
    pytest.importorskip("pysam")
    from tobias.utils.ngs import OneRead
    from tobias.utils.regions import OneRegion
    from tobias.utils.sequences import GenomicSequence

    pysam = _real()
    bam = pysam.AlignmentFile(BAM, "rb")
    fasta = pysam.FastaFile(FASTA)
    k = 12

    region = OneRegion(["chr1", 1_000_000, 1_020_000])
    seq_obj = GenomicSequence(region).from_fasta(fasta)
    codes = encode_sequence(fasta.fetch("chr1", 1_000_000, 1_020_000))
    rcodes = revcomp_codes(codes)

    checked = 0
    for read in bam.fetch("chr1", 1_000_000, 1_020_000):
        ref = OneRead(read)
        ref.get_cutsite((4, -5))
        mine = cutsite(read.reference_start, read.reference_end,
                       read.query_alignment_start, read.query_alignment_end,
                       read.query_length, read.is_reverse)
        assert mine == ref.cutsite, read.query_name

        ref.get_kmer(seq_obj, k)
        mine_kmer = kmer_at(codes, rcodes, 1_000_000, 1_020_000, mine, k,
                            read.is_reverse)
        assert np.array_equal(np.asarray(ref.kmer), mine_kmer), read.query_name
        checked += 1

    bam.close()
    fasta.close()
    assert checked > 0, "no reads were examined"


def test_kmer_outside_the_usable_range_is_all_sentinel():
    """A cut site too close to the region edge yields the reference's filler."""
    codes = np.zeros(100, dtype=np.int64)
    rcodes = revcomp_codes(codes)
    k = 12
    kmer = kmer_at(codes, rcodes, 0, 100, 5, k, False)
    assert kmer.tolist() == [N_CODE] * (2 * k + 1)
