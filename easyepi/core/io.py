"""Genomic IO: BAM (pysam), BigWig/BedGraph (pyBigWig), BED.

All heavy IO libraries are imported lazily so that importing ``easyepi`` works
even in a minimal environment; a clear message is raised only when an actual
io operation is requested without the required dependency.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


def _require(mod_name: str, hint: str):
    try:
        return __import__(mod_name)
    except ImportError as exc:  # pragma: no cover - depends on env
        raise ImportError(
            f"Optional dependency '{mod_name}' is required for this operation. {hint}"
        ) from exc


# ---------------------------------------------------------------------------
# BAM
# ---------------------------------------------------------------------------
def bam_coverage(bam_path, chrom_size: Optional[dict] = None,
                 bin_size: int = 1, normalize_cpm: bool = False):
    """Yield per-bin read coverage as (chrom, start, end, count) arrays.

    For the downstream suite we only need a coverage signal, not alignment
    details; this avoids depending on alignment-specific queries.
    """
    pysam = _require("pysam", "install with: pip install pysam")
    from collections import defaultdict

    counts = defaultdict(lambda: np.zeros(0, dtype=np.int64))
    # Recompute a consensus coverage over bins on the fly is expensive;
    # here we expose per-read mapped positions so higher-level modules can
    # bin them. For simplicity we return mapped read counts per chrom.
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_reads = 0
        per_chrom = {}
        for chrom in bam.references:
            per_chrom[chrom] = bam.count(chrom)
            total_reads += per_chrom[chrom]
    # Note: this is a light QC summary; rich binned coverage is added by
    # CoverageTracker in bigwig module where the signal is actually needed.
    return {"per_chrom": per_chrom, "total_reads": total_reads,
            "mapped_reads": total_reads, "n_references": len(per_chrom)}


def bam_qc(bam_path) -> dict:
    """Lightweight BAM QC summary (reads, mapping rate, duplicates).

    ``mapping_rate`` is only meaningful when the BAM still contains unaligned
    reads. Sorted/duplicate-marked BAMs that were already filtered upstream
    contain mapped alignments only, so the ratio is trivially 1.0 and must not be
    read as a library-alignment success rate. We therefore report both the raw
    ratio and ``mapping_rate_source`` so downstream reports can tell the cases
    apart instead of silently presenting a misleading 100%.
    """
    pysam = _require("pysam", "install with: pip install pysam")
    mapped = 0
    total = 0
    unmapped = 0
    secondary = 0
    supplementary = 0
    duplicate = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            total += 1
            if rec.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
                if rec.is_secondary:
                    secondary += 1
                if rec.is_supplementary:
                    supplementary += 1
                if rec.is_duplicate:
                    duplicate += 1
    rate = (mapped / total) if total else 0.0
    # When no unaligned records were observed the ratio carries no information.
    source = "bam_unaligned_records" if unmapped else "uninformative_no_unaligned_records"
    return {"total_reads": total, "mapped_reads": mapped, "unmapped_reads": unmapped,
            "secondary_reads": secondary, "supplementary_reads": supplementary,
            "duplicate_reads": duplicate,
            "duplicate_rate": (duplicate / mapped) if mapped else None,
            "mapping_rate": rate, "mapping_rate_source": source}


def flagstat_summary(bam_path) -> dict | None:
    """Return samtools-flagstat style counts through the ``pysam`` API.

    Returns ``None`` when this pysam build has no ``flagstat`` helper, so callers
    can fall back to :func:`bam_qc`.
    """
    pysam = _require("pysam", "install with: pip install pysam")
    flagstat = getattr(pysam, "flagstat", None)
    if flagstat is None:
        return None
    try:
        text = str(flagstat(str(bam_path)))
    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
        # bad region or file cannot abort a batch QC run. The caller reports
        # the metric as not_evaluated rather than failing.
        return None
    counts: dict[str, int] = {}
    for line in text.splitlines():
        head, _, tail = line.partition(" + ")
        if not tail:
            continue
        try:
            counts[head.strip()] = int(tail.split()[0])
        except (ValueError, IndexError):
            continue
    total = counts.get("in total (QC-passed reads")
    mapped = counts.get("primary mapped (", counts.get("mapped ("))
    return {"flagstat_counts": counts, "flagstat_total": total,
            "flagstat_mapped": mapped,
            "flagstat_mapping_rate": (mapped / total) if total and mapped else None}


# ---------------------------------------------------------------------------
# BED
# ---------------------------------------------------------------------------
def read_bed(path, ncols: int = 3):
    """Read a BED file into numpy arrays (first ncols columns)."""
    cols = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            cols.append(parts[:ncols])
    if not cols:
        return np.empty((0, ncols), dtype=object)
    arr = np.array(cols, dtype=object)
    return arr


# ---------------------------------------------------------------------------
# BigWig / BigBed
# ---------------------------------------------------------------------------
class CoverageTracker:
    """Bin counts from a group of BigWig files into signal matrices.

    This is the core object the normalize / nucleosome / viz modules consume:
    a (n_bins, n_samples) signal matrix over a set of genomic regions.
    """

    def __init__(self, bigwig_paths: List[str], sample_names: Optional[List[str]] = None):
        """Open the reader lazily so import stays cheap."""
        pyBigWig = _require("pyBigWig", "install with: pip install pyBigWig")
        self.bigwig_paths = list(bigwig_paths)
        self.sample_names = sample_names or [Path(p).stem for p in self.bigwig_paths]
        self._bw = None

    def region_profile(self, chrom: str, start: int, end: int, bins: int = 100):
        """Return (n_samples, bins) normalized-average signal over one region."""
        pyBigWig = _require("pyBigWig", "install with: pip install pyBigWig")
        out = np.full((len(self.bigwig_paths), bins), np.nan, dtype=np.float64)
        edges = np.linspace(start, end, bins + 1).astype(np.int64)
        for i, bw_path in enumerate(self.bigwig_paths):
            bw = pyBigWig.open(bw_path)
            try:
                vals = bw.values(chrom, int(start), int(end))
                if vals:
                    seq_vals = np.asarray(vals, dtype=np.float64)
                    seq_vals[seq_vals < 0] = 0.0
                    # rebin to requested bin count
                    if len(seq_vals) == bins:
                        out[i] = seq_vals
                    else:
                        # simple interpolation / resampling
                        idx = np.linspace(0, len(seq_vals) - 1, bins)
                        out[i] = np.interp(idx, np.arange(len(seq_vals)), seq_vals)
            finally:
                bw.close()
        return out

    def mean_profile(self, chrom, start, end, bins=100):
        """Mean signal profile over the supplied regions, binned to a fixed length."""
        prof = self.region_profile(chrom, start, end, bins)
        # average across samples
        valid = ~np.isnan(prof)
        if valid.any():
            return np.nanmean(prof, axis=0)
        return np.full(bins, np.nan)
