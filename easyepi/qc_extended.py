"""Extended QC metrics: TSS enrichment, duplication, MAPQ, mitochondrial and read length.

Complements :mod:`easyepi.assay_metrics` with the metrics that every
ATAC/ChIP/CUT&Tag QC report is expected to show. All functions are defensive:
an unreadable input yields ``{"error": ...}`` instead of raising, so a single bad
sample never aborts a multi-sample run. Readers are imported lazily.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

# Background for TSS enrichment: flanks 1900-2000 bp from the TSS, matching the
# ENCODE ATAC-seq definition (signal within +-100 bp of the TSS over signal in
# the +-1900..2000 bp flanks).
_TSS_FLANK_INNER = 100
_TSS_FLANK_OUTER = 2000


def _ratio(num, den):
    return float(num) / float(den) if den else None


# ---------------------------------------------------------------------------
# Read-level flags: duplication, MAPQ, mitochondrial fraction, read length
# ---------------------------------------------------------------------------
def read_flag_metrics(bam_path: str, max_reads: int = 500_000,
                      mito_prefixes: tuple[str, ...] = ("chrM", "MT", "chrMT", "M")) -> dict[str, Any]:
    """Duplication rate, MAPQ distribution, mitochondrial fraction, read length.

    ``duplicate_rate`` uses the ``0x400`` duplicate flag when present and falls
    back to a positional exact-duplicate estimate (same reference, start, strand,
    mate start) for BAMs where duplicates were never marked, which is the common
    case for CUT&Tag/ATAC alignments straight out of the aligner.
    ``mitochondrial_fraction`` is a load-bearing QC signal for ATAC/CUT&Tag:
    high values indicate cell lysis or cytoplasmic carry-over.
    """
    import pysam

    total = 0
    flagged_dup = 0
    mito = 0
    lengths: list[int] = []
    mapq_bins = {"mapq_ge_30": 0, "mapq_10_29": 0, "mapq_lt_10": 0}
    seen: set[tuple] = set()
    positional_dup = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                total += 1
                if rec.is_duplicate:
                    flagged_dup += 1
                ref = rec.reference_name or ""
                if ref in mito_prefixes or ref.startswith("chrM"):
                    mito += 1
                key = (rec.reference_id, rec.reference_start, rec.is_reverse,
                       rec.next_reference_id, rec.next_reference_start)
                if key in seen:
                    positional_dup += 1
                else:
                    seen.add(key)
                mq = rec.mapping_quality
                if mq >= 30:
                    mapq_bins["mapq_ge_30"] += 1
                elif mq >= 10:
                    mapq_bins["mapq_10_29"] += 1
                else:
                    mapq_bins["mapq_lt_10"] += 1
                qlen = rec.query_length or 0
                if qlen:
                    lengths.append(qlen)
                if total >= max_reads:
                    break
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}

    if not total:
        return {"reads_inspected": 0}

    lengths.sort()
    n = len(lengths)
    median_read_length = None
    if n:
        median_read_length = (lengths[n // 2] if n % 2
                              else (lengths[n // 2 - 1] + lengths[n // 2]) / 2)
    dup_source = "duplicate_flag" if flagged_dup else "positional_estimate"
    dup_count = flagged_dup if flagged_dup else positional_dup
    return {
        "reads_inspected": total,
        "duplicate_rate": _ratio(dup_count, total),
        "duplicate_source": dup_source,
        "mitochondrial_reads": mito,
        "mitochondrial_fraction": _ratio(mito, total),
        "median_read_length": median_read_length,
        **{k: _ratio(v, total) for k, v in mapq_bins.items()},
    }


# ---------------------------------------------------------------------------
# TSS enrichment (ENCODE ATAC-seq metric)
# ---------------------------------------------------------------------------
def read_tss(annotation: str, max_genes: int = 60_000) -> list[tuple[str, int]]:
    """Extract (chrom, tss) from a GTF/GFF annotation."""
    tss: list[tuple[str, int]] = []
    path = Path(annotation)
    if not path.exists():
        return tss
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2].lower() != "gene":
                    continue
                try:
                    start, end = int(fields[3]) - 1, int(fields[4])
                except ValueError:
                    continue
                strand = fields[6]
                tss.append((fields[0], start if strand != "-" else end - 1))
                if len(tss) >= max_genes:
                    break
    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
        # bad region or file cannot abort a batch QC run. The caller reports
        # the metric as not_evaluated rather than failing.
        return tss
    return tss


def tss_enrichment(bigwig_path: str, annotation: str, max_genes: int = 4000) -> dict[str, Any]:
    """ENCODE-style TSS enrichment score.

    Computes mean signal in +-100 bp around each TSS, divided by the mean signal
    in the +-1900..2000 bp flanking windows. A good ATAC library scores well
    above 2; a failed library sits near 1.
    """
    import pyBigWig

    tss = read_tss(annotation, max_genes=max_genes * 2)
    if not tss:
        return {"error": "no TSS parsed from annotation"}
    if len(tss) > max_genes:
        step = max(1, len(tss) // max_genes)
        tss = tss[::step][:max_genes]

    inner_sum = inner_n = outer_sum = outer_n = 0
    used = 0
    try:
        with pyBigWig.open(bigwig_path) as bw:
            chroms = bw.chroms() or {}
            for chrom, pos in tss:
                length = chroms.get(chrom)
                if length is None:
                    continue
                lo, hi = max(0, pos - _TSS_FLANK_INNER), min(length, pos + _TSS_FLANK_INNER)
                if hi <= lo:
                    continue
                try:
                    vals = bw.values(chrom, lo, hi)
                except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
                    # bad region or file cannot abort a batch QC run. The caller reports
                    # the metric as not_evaluated rather than failing.
                    continue
                vals = [v for v in (vals or []) if v is not None and not math.isnan(v)]
                used += 1
                if vals:
                    inner_sum += sum(vals)
                    inner_n += len(vals)
                for fs, fe in ((max(0, pos - _TSS_FLANK_OUTER), max(0, pos - (_TSS_FLANK_OUTER - _TSS_FLANK_INNER))),
                               (min(length, pos + (_TSS_FLANK_OUTER - _TSS_FLANK_INNER)), min(length, pos + _TSS_FLANK_OUTER))):
                    if fe <= fs:
                        continue
                    try:
                        fvals = bw.values(chrom, fs, fe)
                    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
                        # bad region or file cannot abort a batch QC run. The caller reports
                        # the metric as not_evaluated rather than failing.
                        continue
                    fvals = [v for v in (fvals or []) if v is not None and not math.isnan(v)]
                    if fvals:
                        outer_sum += sum(fvals)
                        outer_n += len(fvals)
    except Exception as exc:
        return {"error": f"bigwig: {type(exc).__name__}: {exc}"}

    inner_mean = _ratio(inner_sum, inner_n)
    outer_mean = _ratio(outer_sum, outer_n)
    score = None
    if inner_mean is not None and outer_mean:
        score = inner_mean / outer_mean
    return {"n_tss_used": used, "tss_mean_signal": inner_mean,
            "flank_mean_signal": outer_mean, "tss_enrichment": score}


# ---------------------------------------------------------------------------
# Peak-level read distribution (replaces the single FRiP summary)
# ---------------------------------------------------------------------------
def _load_intervals(bed_path: str, max_intervals: int = 300_000):
    """Load BED intervals as {chrom: [(start, end), ...]} sorted per chromosome."""
    intervals: dict[str, list[tuple[int, int]]] = {}
    try:
        with Path(bed_path).open(encoding="utf-8", errors="replace") as fh:
            n = 0
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    continue
                if end > start:
                    intervals.setdefault(fields[0], []).append((start, end))
                    n += 1
                    if n >= max_intervals:
                        break
    except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
        # bad region or file cannot abort a batch QC run. The caller reports
        # the metric as not_evaluated rather than failing.
        return {}
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def _find_interval(ivs: list[tuple[int, int]], pos: int) -> int | None:
    """Binary search: index of the interval containing ``pos``, else None."""
    if not ivs:
        return None
    lo, hi = 0, len(ivs) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        s, e = ivs[mid]
        if pos < s:
            hi = mid - 1
        elif pos >= e:
            lo = mid + 1
        else:
            return mid
    return None


def peak_read_distribution(bam_path: str, peak_path: str, max_reads: int = 500_000,
                           max_peaks: int = 50_000, n_bins: int = 20) -> dict[str, Any]:
    """Per-peak read counts: the distribution behind a single FRiP number.

    A lone FRiP value cannot distinguish "every peak is weakly supported" from
    "a few peaks carry all the signal". This reports the per-peak counts, how many
    peaks have zero reads, the concentration of reads in the top decile, and a
    log10 histogram suitable for plotting.
    """
    import pysam

    peaks: list[tuple[str, int, int]] = []
    try:
        with Path(peak_path).open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    s, e = int(fields[1]), int(fields[2])
                except ValueError:
                    continue
                if e > s:
                    peaks.append((fields[0], s, e))
                    if len(peaks) >= max_peaks:
                        break
    except Exception as exc:
        return {"error": f"peak: {type(exc).__name__}: {exc}"}

    if not peaks:
        return {"n_peaks": 0}

    # Group peak indices by chromosome for O(log n) lookup while streaming reads.
    by_chrom: dict[str, list[int]] = {}
    for idx, (chrom, _s, _e) in enumerate(peaks):
        by_chrom.setdefault(chrom, []).append(idx)
    chrom_peaks: dict[str, list[tuple[int, int]]] = {}
    for chrom, idxs in by_chrom.items():
        ordered = sorted(idxs, key=lambda i: peaks[i][1])
        by_chrom[chrom] = ordered
        chrom_peaks[chrom] = [(peaks[i][1], peaks[i][2]) for i in ordered]

    counts = [0] * len(peaks)
    total = 0
    in_peak = 0
    scanned = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Same trap as FRiP: a BAM is coordinate-sorted, so stopping after N
            # reads samples only the first chromosomes and leaves most peaks at
            # zero (observed top_decile_read_share = 1.0 on a real file). Stride
            # across the whole file instead.
            try:
                header_count = bam.mapped
            except Exception:
                header_count = 0
            stride = max(1, int(header_count / max_reads)) if header_count and max_reads else 1
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                scanned += 1
                if stride > 1 and (scanned % stride):
                    continue
                total += 1
                chrom = rec.reference_name or ""
                ivs = chrom_peaks.get(chrom)
                if ivs:
                    hit = _find_interval(ivs, rec.reference_start)
                    if hit is not None:
                        counts[by_chrom[chrom][hit]] += 1
                        in_peak += 1
    except Exception as exc:
        return {"error": f"bam: {type(exc).__name__}: {exc}"}

    import math

    ordered_counts = sorted(counts)
    n = len(ordered_counts)
    zero = sum(1 for c in ordered_counts if c == 0)
    total_peak_reads = sum(ordered_counts)
    top_decile_n = max(1, n // 10)
    top_decile_reads = sum(ordered_counts[-top_decile_n:])
    # Concentration: share of in-peak reads held by the top 10% of peaks.
    concentration = _ratio(top_decile_reads, total_peak_reads)
    median = (ordered_counts[n // 2] if n % 2
              else (ordered_counts[n // 2 - 1] + ordered_counts[n // 2]) / 2) if n else None

    # log10(1+count) histogram for the figure.
    import numpy as np
    logs = np.log10(np.asarray(ordered_counts, dtype=float) + 1.0)
    hist, edges = np.histogram(logs, bins=n_bins)
    return {
        "n_peaks": n,
        "reads_inspected": total,
        "reads_scanned": scanned,
        "stride": stride,
        "reads_in_peaks": in_peak,
        "frip": _ratio(in_peak, total),
        "peaks_with_zero_reads": zero,
        "fraction_peaks_zero": _ratio(zero, n),
        "median_reads_per_peak": median,
        "mean_reads_per_peak": _ratio(total_peak_reads, n),
        "max_reads_per_peak": ordered_counts[-1] if n else None,
        "top_decile_read_share": concentration,
        "hist_counts": [int(x) for x in hist],
        "hist_edges": [float(x) for x in edges],
    }


# ---------------------------------------------------------------------------
# Blacklist fraction
# ---------------------------------------------------------------------------
def blacklist_fraction(bam_path: str, blacklist_path: str,
                       max_reads: int = 500_000) -> dict[str, Any]:
    """Fraction of mapped reads overlapping ENCODE-style blacklist regions.

    High values flag artefact-prone libraries (hyper-chimeric or high-signal
    regions) that should be excluded before peak calling.
    """
    import pysam

    regions = _load_intervals(blacklist_path)
    if not regions:
        return {"error": "blacklist empty or unreadable"}
    total = 0
    in_bl = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                total += 1
                ivs = regions.get(rec.reference_name or "")
                if ivs and _find_interval(ivs, rec.reference_start) is not None:
                    in_bl += 1
                if total >= max_reads:
                    break
    except Exception as exc:
        return {"error": f"bam: {type(exc).__name__}: {exc}"}
    return {"reads_inspected": total, "reads_in_blacklist": in_bl,
            "blacklist_fraction": _ratio(in_bl, total)}


def compute_extended_metrics(assay: str, *, bam: str | None = None,
                             bigwig: str | None = None, annotation: str | None = None,
                             peak: str | None = None, blacklist: str | None = None,
                             max_reads: int = 200_000) -> dict[str, Any]:
    """Assemble the extra QC metrics that apply to ``assay``."""
    out: dict[str, Any] = {}
    if bam:
        flags = read_flag_metrics(bam, max_reads=max_reads)
        if "error" not in flags:
            out["read_flags"] = flags
            for key in ("duplicate_rate", "duplicate_source", "mitochondrial_fraction",
                        "median_read_length", "mapq_ge_30", "mapq_lt_10"):
                if key in flags and flags[key] is not None:
                    out[key] = flags[key]
    if bigwig and annotation:
        tss = tss_enrichment(bigwig, annotation)
        if "error" not in tss:
            out["tss"] = tss
            if tss.get("tss_enrichment") is not None:
                out["tss_enrichment"] = tss["tss_enrichment"]
    if bam and peak:
        prd = peak_read_distribution(bam, peak, max_reads=max_reads)
        if "error" not in prd:
            out["peak_reads"] = prd
            for key in ("peaks_with_zero_reads", "fraction_peaks_zero",
                        "median_reads_per_peak", "max_reads_per_peak",
                        "top_decile_read_share"):
                if prd.get(key) is not None:
                    out[key] = prd[key]
    if bam and blacklist:
        bl = blacklist_fraction(bam, blacklist, max_reads=max_reads)
        if "error" not in bl:
            out["blacklist"] = bl
            if bl.get("blacklist_fraction") is not None:
                out["blacklist_fraction"] = bl["blacklist_fraction"]
    return out
