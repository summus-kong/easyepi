"""Read-level coverage used when the signal needs read 5' ends (Tn5 shift).

Why a separate path
-------------------
``count_coverage`` histograms the CIGAR-aligned read *body* into per-base depth.
That is the right answer for ChIP-seq signal, and the wrong one for ATAC-seq and
CUT&Tag, where the measured quantity is the transposase insertion site: the 5' end
of each read, corrected by +4 (forward) / -5 (reverse).

Once depth has been accumulated the read ends are gone, so the shift cannot be
applied afterwards. Rather than silently producing an unshifted track, this module
accumulates read-end positions directly.

Why it is still fast
--------------------
Reads are collected in C-level iteration (``fetch``), but only their end positions
are kept, and the histogram is built with ``np.bincount`` instead of a Python loop
per read. The cost is one pass over the file and one vectorised bincount per
chromosome, which is dominated by I/O rather than by Python.

Fragment handling
-----------------
Paired-end BAMs are the norm for ATAC. Using both mates would count each fragment's
two ends twice per strand, so the inter-mate region would be double counted. Only
one end per fragment is used (the 5' end of read 1), which is the standard
"one cut site per molecule" convention.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def _require_pysam():
    from easyepi.core.io import _require
    return _require("pysam", "pip install 'easyepi[bam]'")


def read_end_positions(bam_path: str | Path, chrom: str, start: int = 0,
                       end: int | None = None, *, min_mapq: int = 30,
                       shift: str = "tn5",
                       one_per_fragment: bool = True) -> tuple[np.ndarray, dict[str, Any]]:
    """5' positions of alignments on one chromosome, after the strand shift.

    Returns ``(positions, stats)`` where ``positions`` are 0-based offsets on
    ``chrom`` (clipped to ``[start, end)``). ``stats`` records how many records were
    seen, used, skipped and clipped, so a suspiciously low count is visible.
    """
    pysam = _require_pysam()
    if shift not in ("none", "tn5"):
        raise ValueError("shift must be 'none' or 'tn5'")

    offset = 4 if shift == "tn5" else 0
    offset_rev = -5 if shift == "tn5" else 0

    buf: list[int] = []
    seen = skipped = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if chrom not in bam.references:
            raise ValueError(
                f"chrom {chrom!r} is not in the BAM header "
                f"(has e.g. {list(bam.references)[:5]})")
        if end is None:
            end = bam.get_reference_length(chrom)
        for rec in bam.fetch(chrom, int(start), int(end)):
            seen += 1
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                skipped += 1
                continue
            if rec.mapping_quality < min_mapq:
                skipped += 1
                continue
            # one cut site per fragment: keep read 1 of a pair
            if one_per_fragment and rec.is_paired and not rec.is_read1:
                skipped += 1
                continue
            if rec.is_reverse:
                pos = int(rec.reference_end) + offset_rev
            else:
                pos = int(rec.reference_start) + offset
            if pos < start or pos >= end:
                skipped += 1
                continue
            buf.append(pos - int(start))

    positions = np.fromiter(buf, dtype=np.int64, count=len(buf))
    stats = {"records_seen": seen, "used": int(positions.size),
             "skipped": skipped, "shift": shift,
             "one_per_fragment": one_per_fragment}
    return positions, stats


def read_end_bins(bam_path: str | Path, chrom: str, length: int, bin_size: int,
                  *, min_mapq: int = 30, shift: str = "tn5",
                  one_per_fragment: bool = True,
                  shift_span: int = 1) -> tuple[np.ndarray, dict[str, Any]]:
    """Bin the shifted 5' positions into a fixed bin grid for one chromosome.

    ``shift_span`` widens each cut site into ``[-span, span]`` bins, which
    approximates the footprint footprint of an insertion rather than a single
    base; ``1`` (default) keeps a one-bin spike, matching how the callers consume
    the track.
    """
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")
    n_bins = max(1, int(np.ceil(length / bin_size)))
    positions, stats = read_end_positions(bam_path, chrom, 0, length,
                                          min_mapq=min_mapq, shift=shift,
                                          one_per_fragment=one_per_fragment)
    if positions.size == 0:
        return np.zeros(n_bins), {**stats, "n_bins": n_bins, "bin_size": bin_size}

    bins = positions // bin_size
    bins = bins[bins < n_bins]
    counts = np.bincount(bins, minlength=n_bins).astype(np.float64)

    if shift_span > 1:
        kernel = np.ones(2 * shift_span - 1, dtype=np.float64)
        counts = np.convolve(counts, kernel, mode="same")

    return counts, {**stats, "n_bins": n_bins, "bin_size": bin_size,
                    "shift_span": int(shift_span)}
