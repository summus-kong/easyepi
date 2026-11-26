"""1:1 reimplementation of ATACorrect's cut-site and k-mer extraction.

This is the one part of ATACorrect that previously lived only in a driver script and
was never isolated. Whole-track comparison showed the forward strand at corr
0.999992 against 0.989637 for reverse; scoring was then ruled out (both strands agree
with TOBIAS to machine precision), so the difference has to be here.

Reference: `tobias/utils/ngs.pyx`, `OneRead.get_cutsite` and `OneRead.get_kmer`.

Two details drive the strand asymmetry, and neither is obvious:

* ``get_cutsite`` corrects for soft clipping. It works from
  ``reference_start - query_alignment_start`` and
  ``reference_end + query_length - query_alignment_end``, where ``query_length``
  includes clipped bases. The reverse branch uses the second expression, which
  carries the ``query_length - query_alignment_end`` term -- so soft-clipped reads
  move reverse cut sites and leave forward ones alone when the clip is on the other
  end. A naive ``reference_start``/``reference_end`` gives the same answer only for
  reads without soft clips.
* ``get_kmer`` refuses cut sites that sit within ``k_flank`` of the region edge, and
  fills a sentinel of ``4`` (N) instead. It also indexes the reverse k-mer from the
  *end*: ``i = seq_end - cutsite``, into the reverse-complement array.
"""
from __future__ import annotations

import numpy as np

#: sentinel filling a k-mer the reference refuses to build
N_CODE = 4

#: the reference alphabet is alphabetical A/T/C/G, not the usual A/C/G/T
BASE_TO_CODE = {"A": 0, "T": 1, "C": 2, "G": 3}

#: complement table in that alphabet: A<->T, C<->G, N->N
CODE_COMPLEMENT = np.array([1, 0, 3, 2, N_CODE], dtype=np.int64)


def encode_sequence(seq: str) -> np.ndarray:
    """A/T/C/G -> 0/1/2/3 and anything else -> 4, matching ``GenomicSequence``."""
    return np.array([BASE_TO_CODE.get(b, N_CODE) for b in seq.upper()],
                    dtype=np.int64)


def revcomp_codes(codes: np.ndarray) -> np.ndarray:
    """Reverse complement in the reference's alphabet; the sentinel stays 4."""
    return CODE_COMPLEMENT[np.asarray(codes, dtype=np.int64)[::-1]]


def cutsite(reference_start: int, reference_end: int, query_alignment_start: int,
            query_alignment_end: int, query_length: int, is_reverse: bool,
            read_shift=(4, -5)) -> int:
    """Cut site as ``OneRead.get_cutsite`` computes it.

    ``positions[0]`` is the reference start pulled back by the leading soft clip;
    ``positions[1]`` is the reference end pushed forward by the trailing soft clip.
    The ``+1`` converts the BAM's 0-based start to the reference's 1-based cut site.
    """
    pos_shift, neg_shift = read_shift
    positions = (
        reference_start - query_alignment_start,
        reference_end + query_length - query_alignment_end,
    )
    if is_reverse:
        return positions[1] + 1 + neg_shift
    return positions[0] + 1 + pos_shift


def kmer_at(codes: np.ndarray, revcomp: np.ndarray, region_start: int,
            region_end: int, cut: int, k_flank: int, is_reverse: bool) -> np.ndarray:
    """The ``2*k_flank+1`` k-mer around a cut site, or an all-N sentinel.

    ``codes`` is the region's forward sequence and ``revcomp`` its reverse complement,
    both indexed from ``region_start``.
    """
    width = 2 * k_flank + 1

    # the reference's guard: the k-mer must fit fully inside the region
    if not (cut > region_start + k_flank + 1 and cut <= region_end - k_flank):
        return np.full(width, N_CODE, dtype=np.int64)

    if is_reverse:
        i = region_end - cut
        return np.asarray(revcomp[i - k_flank:i + k_flank + 1], dtype=np.int64)

    i = cut - 1 - region_start
    return np.asarray(codes[i - k_flank:i + k_flank + 1], dtype=np.int64)


def read_kmers(reads, codes: np.ndarray, revcomp: np.ndarray, region_start: int,
               region_end: int, k_flank: int, read_shift=(4, -5),
               min_per_position: int | None = None):
    """Every usable k-mer in a region, with the reference's deduplication rules.

    ``reads`` are mappings with the pysam attributes used by
    ``OneRead.get_cutsite``. ``min_per_position`` mirrors the cap the reference
    applies when several reads share a cut site (``no_cut = min(n, 10)``); pass
    ``None`` to return one entry per read.
    """
    per_position = {}
    for read in reads:
        if read.cigartuples and read.cigartuples[0][0] == 4:
            # the reference skips reads that start with a soft clip
            continue
        cut = cutsite(read.reference_start, read.reference_end,
                      read.query_alignment_start, read.query_alignment_end,
                      read.query_length, read.is_reverse, read_shift)
        per_position.setdefault(cut, (read.is_reverse, 0))
        is_rev, n = per_position[cut]
        per_position[cut] = (is_rev, n + 1)

    out = []
    for cut, (is_rev, n) in per_position.items():
        weight = n if min_per_position is None else min(n, min_per_position)
        km = kmer_at(codes, revcomp, region_start, region_end, cut, k_flank, is_rev)
        out.append((km, weight, is_rev))
    return out
