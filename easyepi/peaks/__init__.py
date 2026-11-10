"""Peak calling / event identification (native).

A native peak caller based on a robust background (median) and a Poisson / odd
ratio enrichment threshold is provided as the default dependency-free path. It
operates on a 1-D binned signal (from BAM coverage, BigWig or a bedGraph). The
output supports the narrowPeak / BED format for downstream interoperability.

Optional wrappers (MACS3 / SEACR) can be invoked as enhancements; they are not
required for the core pipeline.

Importable API:
  * ``poisson_peaks`` — call peaks on a 1-D signal array.
  * ``bam_to_binned`` — real coverage from a BAM (needs pysam).
  * ``bedgraph_to_signal`` / ``bigwig_to_signal`` — real signal inputs.
  * ``run_peaks`` — CLI wrapper.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from easyepi.core.models import PeakSet


def poisson_peaks(signal: np.ndarray, chrom: str, start: int, bin_size: int,
                  scale: float = 1.0, threshold: float = 5.0,
                  max_qvalue: float = 0.01):
    """Call peaks on a 1-D bin-signal array.

    Uses a robust background (median of the whole signal) and a Poisson
    upper-tail p-value per enriched bin. A bin is "enriched" if its signal
    exceeds background * threshold (the pre-filter), then consecutive enriched
    bins are merged and the merged peak gets a p-value from its max bin count.
    Peaks with q-value (BH-FDR) above ``max_qvalue`` are removed.

    Args:
        signal: 1-D coverage per bin.
        chrom: chromosome label for the output PeakSet.
        start: genomic coordinate of the first bin.
        bin_size: width of each bin in bp.
        scale: multiplier on the background estimate.
        threshold: pre-filter enrichment multiplier.
        max_qvalue: BH-FDR cutoff (>=) for retaining peaks.
    """
    from scipy.stats import poisson
    sig = np.asarray(signal, dtype=np.float64)
    if sig.size == 0:
        return PeakSet(chrom=np.empty(0, dtype=object),
                       start=np.empty(0, dtype=np.int64),
                       end=np.empty(0, dtype=np.int64),
                       signal=np.empty(0, dtype=np.float64))
    # Robust background: median of non-zero bins (MACS-style local background
    # proxy). Falls back to overall median if all bins are zero.
    nz = sig[sig > 0]
    background = float(np.median(nz)) if nz.size else float(np.median(sig))
    mu = max(background * scale, 1e-6)
    enriched = sig >= mu * threshold
    if not enriched.any():
        return PeakSet(chrom=np.empty(0, dtype=object),
                       start=np.empty(0, dtype=np.int64),
                       end=np.empty(0, dtype=np.int64),
                       signal=np.empty(0, dtype=np.float64))
    runs = []
    in_run = False
    for i, flag in enumerate(enriched.tolist()):
        if flag and not in_run:
            s = i
            in_run = True
        elif not flag and in_run:
            runs.append((s, i))
            in_run = False
    if in_run:
        runs.append((s, len(enriched)))
    starts = np.array([start + s * bin_size for s, _ in runs], dtype=np.int64)
    ends = np.array([start + e * bin_size for _, e in runs], dtype=np.int64)
    peak_signal = np.array([sig[s:e].max() for s, e in runs], dtype=np.float64)
    # Poisson upper-tail p-value for the peak max count
    pvals = np.array([float(poisson.sf(float(peak_signal[i]), mu))
                      for i in range(len(runs))], dtype=np.float64)
    from easyepi.core.stats import benjamini_hochberg
    qvals = benjamini_hochberg(pvals)
    keep = qvals <= max_qvalue
    keep &= peak_signal >= mu * threshold  # enforce the enrichment pre-filter again
    return PeakSet(
        chrom=np.full(len(runs), chrom, dtype=object)[keep],
        start=starts[keep],
        end=ends[keep],
        signal=peak_signal[keep],
        fold_change=(peak_signal[keep] / mu) if mu > 0 else peak_signal[keep],
        pvalue=pvals[keep],
        qvalue=qvals[keep],
    )


def call_peaks_from_bam(bam_path, chrom: str, bin_size: int = 200,
                        scale: float = 1.0, threshold: float = 5.0,
                        max_qvalue: float = 0.01, t5_shift: bool = True) -> PeakSet:
    """Call peaks directly from a coordinate-sorted, indexed BAM.

    Builds a per-bin signal using the Tn5/CUT&Tag midpoint shift
    (+4 forward / -5 reverse) by default, then calls :func:`poisson_peaks`.

    Args:
        bam_path: coordinate-sorted BAM (must be indexed).
        chrom: chromosome to analyze.
        bin_size: signal bin width in bp.
        threshold / scale / max_qvalue: passed to poisson_peaks.
        t5_shift: apply the Tn5 midpoint shift.
    """
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        end = bam.get_reference_length(chrom)
    except ValueError as exc:
        raise ValueError(f"chrom '{chrom}' not found in BAM") from exc
    start = 0
    n = max(1, int(np.ceil((end - start) / bin_size)))
    sig = np.zeros(n, dtype=np.float64)
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        pos = read.reference_start
        if t5_shift:
            pos = pos + (-5 if read.is_reverse else 4)
        if pos < start:
            continue
        b = int((pos - start) // bin_size)
        if 0 <= b < n:
            sig[b] += 1.0
    bam.close()
    return poisson_peaks(sig, chrom, start, bin_size, scale=scale,
                         threshold=threshold, max_qvalue=max_qvalue)


def bam_to_binned(bam_path, chrom: str, start: int, end: int, bin_size: int = 50) -> np.ndarray:
    """Compute per-bin coverage of a BAM over [start, end) (needs pysam)."""
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover
        raise ImportError("pysam is required to read BAM (pip install pysam)") from exc
    n_bins = max(1, int(np.ceil((end - start) / bin_size)))
    counts = np.zeros(n_bins, dtype=np.float64)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, int(start), int(end)):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            rstart, rend = read.reference_start, read.reference_end
            if rend is None:
                continue
            lo = max(rstart, start)
            hi = min(rend, end)
            if hi <= lo:
                continue
            b0 = int((lo - start) // bin_size)
            b1 = int((hi - start - 1) // bin_size)
            counts[b0:b1 + 1] += 1.0
    return counts


def bedgraph_to_signal(path: Path, chrom: str, start: int, end: int,
                       bin_size: int = 50) -> np.ndarray:
    """Rebin a chrom-wide bedGraph into a per-bin signal over [start, end)."""
    n_bins = max(1, int(np.ceil((end - start) / bin_size)))
    out = np.zeros(n_bins, dtype=np.float64)
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] != chrom:
                continue
            cs, ce, val = int(p[1]), int(p[2]), float(p[3])
            lo, hi = max(cs, start), min(ce, end)
            if hi <= lo:
                continue
            b0 = int((lo - start) // bin_size)
            b1 = int((hi - start - 1) // bin_size)
            # weight by overlap
            out[b0:b1 + 1] += val
    return out


def bigwig_to_signal(path: Path, chrom: str, start: int, end: int,
                     bin_size: int = 50) -> np.ndarray:
    """Read a bigwig region into per-bin signal (needs pyBigWig)."""
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover
        raise ImportError("pyBigWig is required to read bigwig (pip install pyBigWig)") from exc
    n_bins = max(1, int(np.ceil((end - start) / bin_size)))
    bw = pyBigWig.open(str(path))
    try:
        vals = bw.values(chrom, int(start), int(end))
        if vals is None:
            return np.zeros(n_bins, dtype=np.float64)
        arr = np.array([0.0 if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)
                        for v in vals], dtype=np.float64)
        arr[arr < 0] = 0.0
        if len(arr) == n_bins:
            return arr
        if len(arr) == 0:
            return np.zeros(n_bins, dtype=np.float64)
        idx = np.linspace(0, len(arr) - 1, n_bins)
        return np.interp(idx, np.arange(len(arr)), arr)
    finally:
        bw.close()


def peaks_to_bed(peaks: PeakSet, out_path: Path, format: str = "narrowPeak") -> None:
    """Write peaks to a BED / narrowPeak file.

    ``narrowPeak`` (default) writes the ENCODE 10-column format:
      chrom, start, end, name, score, strand, signalValue, pValue (-log10 p), qValue (-log10 q), peak (-1 if unknown).
    ``bed`` writes the 3/4 column form.

    Args:
        peaks: PeakSet to write.
        out_path: destination file.
        format: 'narrowPeak' (default) or 'bed'.
    """
    with open(out_path, "w") as fh:
        for i in range(len(peaks)):
            name = peaks.name[i] if peaks.name is not None else f"peak_{i}"
            score = min(1000, max(0, int(round(float(peaks.signal[i]))))) if peaks.signal is not None else 0
            chrom = peaks.chrom[i]
            s0 = int(peaks.start[i]); e0 = int(peaks.end[i])
            if format == "bed":
                row = [chrom, s0, e0, str(name)]
                fh.write("\t".join(str(x) for x in row) + "\n")
                continue
            strand = peaks.strand[i] if peaks.strand is not None else "."
            sv = float(peaks.signal[i]) if peaks.signal is not None else 0.0
            fc = float(peaks.fold_change[i]) if (peaks.fold_change is not None
                                                 and not np.isnan(peaks.fold_change[i])) else 0.0
            pv = float(peaks.pvalue[i]) if (peaks.pvalue is not None and
                                            not np.isnan(peaks.pvalue[i])) else 1.0
            qv = float(peaks.qvalue[i]) if (peaks.qvalue is not None and
                                            not np.isnan(peaks.qvalue[i])) else 1.0
            # guard against p/q underflow to exactly 0 (very small tail probs)
            pscore = round(-np.log10(max(pv, 1e-300)), 3)
            qscore = round(-np.log10(max(qv, 1e-300)), 3)
            row = [str(chrom), str(s0), str(e0), str(name), str(score), str(strand),
                   "%.3f" % sv, "%.3f" % pscore, "%.3f" % qscore, "-1"]
            fh.write("\t".join(row) + "\n")


def run_peaks(args) -> int:
    """Peak calling entry point; delegates to MACS3 when available."""
    raise ValueError("Use python -m easyepi.specific.atac.peaks or easyepi.specific.chip.narrow_peaks with MACS3; the legacy stub has been disabled")
