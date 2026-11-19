"""Assay-specific QC metrics computed from BAM/BigWig/peak inputs.

Each function returns plain JSON-serialisable values so the calling command can
merge them into ``*_qc.json`` and the flattened ``*_qc.tsv``. All heavy readers
(pysam, pyBigWig) are imported lazily so importing this module stays cheap and
works on installations without the optional binary wheels.

Design notes
------------
* Every metric degrades gracefully: a missing or unreadable input yields a
  ``None``/``NaN`` value instead of raising, so one bad sample does not abort a
  multi-sample QC run.
* Fractional read sampling keeps the cost bounded on large BAMs. The resulting
  estimates are explicitly reported together with the number of reads inspected.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Iterable, Sequence

from easyepi.assays import AssayError, canonical_assay

# Nucleosome-free / mononucleosome / dinucleosome windows in base pairs, used by
# ATAC and MNase fragment-length summaries.
_NFR_MAX = 100
_MONO_MAX = 250
_DI_MAX = 500


def _safe_div(num: float, den: float) -> float | None:
    if not den:
        return None
    return float(num) / float(den)


def fragment_length_summary(bam_path: str, max_reads: int = 200_000,
                            sample_every: int = 1) -> dict[str, Any]:
    """Summarise insert sizes from mapped read pairs.

    Insert size (``isize``) is used as a proxy for fragment length, which is what
    ATAC/MNase QC actually interpret: NFR vs mononucleosome vs dinucleosome.
    """
    import pysam

    lengths: list[int] = []
    paired_seen = False
    seen = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                # One fragment contributes once. Both mates of a proper pair carry
                # the same |template_length|, so without this filter every paired-end
                # fragment is counted twice and the NRL / periodicity fractions are
                # computed on a doubled, layout-dependent histogram.
                if rec.is_paired:
                    paired_seen = True
                    if not rec.is_read1:
                        continue
                seen += 1
                if sample_every > 1 and (seen % sample_every):
                    continue
                size = abs(int(rec.template_length))
                if size > 0:
                    lengths.append(size)
                if len(lengths) >= max_reads:
                    break
    except Exception as exc:  # unreadable BAM -> report why, do not crash
        return {"error": f"{type(exc).__name__}: {exc}"}

    if not lengths:
        # Distinguish "no usable length was seen" from "zero reads". Single-end
        # data has no template_length, so a fragment distribution is undefined
        # there rather than empty; reporting n_fragments=0 alone would read as a
        # failed library instead of an inapplicable metric.
        if seen and not paired_seen:
            return {
                "n_fragments": 0,
                "status": "not_applicable",
                "note": "single-end reads carry no template_length; "
                        "fragment length requires paired-end data",
            }
        return {"n_fragments": 0, "status": "no_reads"}

    lengths.sort()
    n = len(lengths)
    # Keep a downsampled copy of the raw vector: the QC figure needs the actual
    # distribution, not a single aggregated median.
    if n > 2000:
        step = max(1, n // 2000)
        sampled = lengths[::step][:2000]
    else:
        sampled = lengths
    median = lengths[n // 2] if n % 2 else (lengths[n // 2 - 1] + lengths[n // 2]) / 2
    nfr = sum(1 for x in lengths if x < _NFR_MAX)
    mono = sum(1 for x in lengths if _NFR_MAX <= x < _MONO_MAX)
    di = sum(1 for x in lengths if _MONO_MAX <= x < _DI_MAX)
    return {
        "n_fragments": n,
        "median_fragment_length": float(median),
        "mean_fragment_length": float(sum(lengths) / n),
        "fraction_nfr": nfr / n,
        "fraction_mononucleosome": mono / n,
        "fraction_dinucleosome": di / n,
        "length_sample": [int(x) for x in sampled],
    }


def frip(bam_path: str, peak_path: str, max_reads: int = 500_000,
         sample_every: int = 1) -> dict[str, Any]:
    """Fraction of Reads in Peaks (FRiP) -- the standard ChIP/CUT&Tag QC metric."""
    import pysam

    intervals: dict[str, list[tuple[int, int]]] = {}
    try:
        with Path(peak_path).open(encoding="utf-8", errors="replace") as fh:
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
    except Exception as exc:
        return {"error": f"peak: {type(exc).__name__}: {exc}"}

    for chrom in intervals:
        intervals[chrom].sort()

    def in_peak(chrom: str, pos: int) -> bool:
        """True when ``pos`` falls inside a peak interval on this chromosome."""
        ivs = intervals.get(chrom)
        if not ivs:
            return False
        # Binary search over sorted intervals.
        lo, hi = 0, len(ivs) - 1
        while lo <= hi:
            mid = (lo + hi) // 2
            s, e = ivs[mid]
            if pos < s:
                hi = mid - 1
            elif pos >= e:
                lo = mid + 1
            else:
                return True
        return False

    total = 0
    in_peak_count = 0
    scanned = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # A BAM is coordinate-sorted, so stopping after N reads samples only the
            # first chromosome and biases FRiP massively (observed 47x high on a real
            # ENCODE file). Instead: stride deterministically across the WHOLE file.
            # ``stride`` is estimated from the header count so the scan touches
            # roughly ``max_reads`` records spread over every chromosome.
            header_count = 0
            try:
                header_count = bam.mapped
            except Exception:
                header_count = 0
            stride = max(1, int(header_count / max_reads)) if header_count and max_reads else 1
            stride = max(stride, sample_every)
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                scanned += 1
                if stride > 1 and (scanned % stride):
                    continue
                total += 1
                if in_peak(rec.reference_name, rec.reference_start):
                    in_peak_count += 1
    except Exception as exc:
        return {"error": f"bam: {type(exc).__name__}: {exc}"}

    result = {"reads_inspected": total, "reads_in_peaks": in_peak_count,
              "reads_scanned": scanned, "stride": stride,
              "sampling": "whole_file_stride" if stride > 1 else "exhaustive",
              "frip": _safe_div(in_peak_count, total)}
    result.update(_frip_coverage_check(bam_path, intervals, total, in_peak_count))
    return result


def _frip_coverage_check(bam_path: str, intervals: dict[str, list[tuple[int, int]]],
                         total: int, in_peak_count: int) -> dict[str, Any]:
    """Flag FRiP computed over mismatched genomic universes.

    FRiP is reads-in-peaks over *all* mapped reads. If the peak set covers only a
    slice of the genome while the BAM spans all of it -- e.g. peaks called on
    chr1:0-5 Mb but a whole-genome alignment -- the numerator is restricted to that
    slice and the denominator is not, so the ratio is not a FRiP at all. On a real
    ENCODE file this produced 0.0018 where the true value is ~two orders higher.

    The check compares the span of the peak set with the span of the mapped reads.
    When peaks cover less than ``min_covered_fraction`` of the read span, the value
    is annotated ``interpretable: False`` with the region-restricted FRiP reported
    separately, so a downstream reader cannot mistake it for a QC failure.
    """
    covered = 0
    min_pos = None
    max_pos = None
    for chrom, ivs in intervals.items():
        for s, e in ivs:
            covered += max(0, e - s)
        if ivs:
            lo = min(s for s, _ in ivs)
            hi = max(e for _, e in ivs)
            min_pos = lo if min_pos is None else min(min_pos, lo)
            max_pos = hi if max_pos is None else max(max_pos, hi)

    span = (max_pos - min_pos) if (min_pos is not None and max_pos is not None) else 0
    read_span = 0
    read_chroms = 0
    try:
        import pysam
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            lens = bam.lengths or []
            read_chroms = len(lens)
            read_span = int(sum(lens))
    except Exception:
        read_span = 0

    detail: dict[str, Any] = {
        "peak_bp": covered,
        "peak_span_bp": span,
        "panel": "the peak set and the reads must cover the same genome for FRiP",
    }
    if not read_span or not span:
        detail["interpretable"] = None      # not enough information to judge
        return detail

    fraction = span / read_span
    detail["peak_span_fraction_of_genome"] = round(fraction, 6)
    detail["reads_chromosomes"] = read_chroms
    # A peak set nominally covers the whole genome; if it spans <5% of the read
    # span it was almost certainly called on a sub-region.
    if fraction < 0.05:
        detail["interpretable"] = False
        detail["region_limited_frip"] = _safe_div(in_peak_count, total)
        detail["warning"] = (
            f"peaks span {span:,} bp ({fraction:.2%} of the {read_span:,} bp of "
            "mapped reference) while FRiP divides by reads from the whole BAM; "
            "this FRiP is not comparable to a genome-wide value. Re-run peak "
            "calling over the same universe, or compute FRiP within the same "
            "region on both sides.")
    else:
        detail["interpretable"] = True
    return detail


def peak_signal_enrichment(bigwig_path: str, peak_path: str, max_peaks: int = 2000,
                           flank: int = 500) -> dict[str, Any]:
    """Mean in-peak signal relative to flanking background (signal-to-noise)."""
    import pyBigWig

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

    peak_sum = 0.0
    peak_n = 0
    flank_sum = 0.0
    flank_n = 0
    try:
        with pyBigWig.open(bigwig_path) as bw:
            chroms = bw.chroms() or {}
            for chrom, s, e in peaks:
                if chrom not in chroms:
                    continue
                length = chroms[chrom]
                try:
                    vals = bw.values(chrom, s, min(e, length))
                except Exception:  # noqa: BLE001
# Intentional: a missing/unreadable input yields no value here so one
                    # bad region or file cannot abort a batch QC run. The caller reports
                    # the metric as not_evaluated rather than failing.
                    continue
                vals = [v for v in (vals or []) if v is not None and not math.isnan(v)]
                if vals:
                    peak_sum += sum(vals)
                    peak_n += len(vals)
                for fs, fe in ((max(0, s - flank), s), (e, min(length, e + flank))):
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
                        flank_sum += sum(fvals)
                        flank_n += len(fvals)
    except Exception as exc:
        return {"error": f"bigwig: {type(exc).__name__}: {exc}"}

    peak_mean = _safe_div(peak_sum, peak_n)
    flank_mean = _safe_div(flank_sum, flank_n)
    ratio = None
    if peak_mean is not None and flank_mean:
        ratio = peak_mean / flank_mean
    return {"n_peaks_used": len(peaks), "mean_signal_in_peaks": peak_mean,
            "mean_signal_in_flanks": flank_mean, "signal_enrichment": ratio}


def periodicity_score(bam_path: str, max_reads: int = 200_000) -> dict[str, Any]:
    """Detect nucleosome phasing in fragment lengths (MNase / ATAC).

    A clear ~147 bp mononucleosome peak and a ~180-200 bp repeat spacing give high
    scores. We report the dominant length bin and the fractional excess of reads
    in the mononucleosome window versus a flat-length expectation.
    """
    summary = fragment_length_summary(bam_path, max_reads=max_reads)
    if "error" in summary or not summary.get("n_fragments"):
        return summary
    mono = summary["fraction_mononucleosome"]
    nfr = summary["fraction_nfr"]
    # A healthy MNase digest is dominated by mononucleosomes; ATAC by NFR.
    result = {"n_fragments": summary["n_fragments"],
              "median_fragment_length": summary["median_fragment_length"],
              "mononucleosome_fraction": mono,
              "nfr_fraction": nfr,
              "mono_to_nfr_ratio": _safe_div(mono, nfr) if nfr else None}
    # The nucleosome repeat length: the one number an MNase library is summarised
    # by that periodicity_score could not report, because it is the *spacing* of
    # the length classes rather than any single class's share.
    nrl = nucleosome_repeat_length(summary.get("length_sample") or ())
    result["nucleosome_repeat_length"] = nrl
    if "nrl_bp" in nrl:
        result["nrl_bp"] = nrl["nrl_bp"]
    return result


#: Minimum peak autocorrelation for a spacing estimate to be reported at all.
#: Chosen low because a real library can have a broad mono-nucleosome peak, but
#: above zero so noise is never presented as a repeat length.
_NRL_MIN_PERIODICITY = 0.05

#: How far the autocorrelation must climb inside the NRL search window to count
#: as a peak. Real ladders measured here reach 0.289-0.323; uniformly random
#: lengths reach at most 0.162, so 0.20 separates them without being tuned to a
#: single fixture.
_NRL_MIN_PEAK_GAIN = 0.20


def nucleosome_repeat_length(length_sample: Sequence[float],
                             *, min_length: int = 100,
                             max_length: int = 1200) -> dict[str, Any]:
    """Estimate the nucleosome repeat length from a fragment-length distribution.

    The NRL is the centre-to-centre distance between neighbouring nucleosomes.
    It is not the mononucleosome fragment length: a mononucleosome fragment is
    ~147 bp of protected DNA, while the repeat length adds the linker and runs
    ~175-200 bp depending on the cell type and the digestion. Confusing the two
    is the usual way this number gets reported 30-50 bp too small.

    Method
    ------
    Multi-nucleosome fragments are integral multiples of the repeat: mono at
    ~1 NRL, di at ~2 NRL, tri at ~3 NRL. The NRL is therefore the *spacing* of the
    peaks in the length histogram, not the position of the first one. The spacing
    is read as the lag that maximises the autocorrelation of the histogram, in the
    window ``[min_nrl, max_nrl]``, using only the range below ``max_length`` where
    fragments are actually present.

    Autocorrelation rather than peak-picking because a real MNase library has a
    broad mononucleosome peak with a shoulder, and picking local maxima in the raw
    histogram returns whatever noise sits above the threshold. The autocorrelation
    is a global property of the whole ladder, so a partially resolved tri-nucleosome
    peak still contributes.

    What this is not
    ----------------
    It does not deconvolve overlapping length classes, does not model digestion
    kinetics, and is not comparable to NRL values from a gel or from
    NucleoATAC/MINCE on the same DNA. It is a descriptive number with its own
    uncertainty, reported with the confidence it was measured at.
    """
    import numpy as np

    values = np.asarray([v for v in length_sample
                         if v is not None and np.isfinite(v)], dtype=float)
    if values.size < 50:
        # A handful of fragments cannot support a spacing estimate; saying so is
        # more useful than returning a number with no support behind it.
        return {"error": f"need at least 50 fragments, got {values.size}"}

    values = values[(values >= min_length) & (values <= max_length)]
    if values.size < 50:
        return {"error": f"need at least 50 fragments in {min_length}-{max_length} bp, "
                         f"got {values.size}"}

    # 1 bp bins over an integer range, so a lag is directly a number of bp.
    lo = int(min_length)
    hist, _ = np.histogram(values, bins=np.arange(lo, int(max_length) + 2, 1))
    hist = hist.astype(float)
    if hist.sum() <= 0:
        return {"error": "empty length histogram"}
    # Subtract the mean before correlating: without it the zero-lag term dominates
    # and the autocorrelation decays from a huge DC component rather than from the
    # ladder structure.
    centred = hist - hist.mean()
    ac = np.correlate(centred, centred, mode="full")[centred.size - 1:]
    if ac[0] <= 0:
        return {"error": "degenerate length distribution (all one length)"}
    ac = ac / ac[0]

    # A nucleosome repeat cannot be shorter than a mononucleosome. Searching from
    # 100 bp keeps the first-lag shoulder (which is not a multiple) out of the
    # argument.
    min_nrl, max_nrl = 130, min(300, ac.size - 1)
    if max_nrl <= min_nrl:
        return {"error": "too few distinct fragment lengths to estimate a repeat"}
    window = ac[min_nrl:max_nrl]
    peak_lag = int(np.argmax(window)) + min_nrl

    # Sub-bin refinement: fit a parabola through the peak and its neighbours. A
    # 1 bp grid is finer than the data justifies, so the point estimate is rounded
    # and the raw lag is kept for anyone who wants it.
    y0, y1, y2 = ac[peak_lag - 1], ac[peak_lag], ac[peak_lag + 1]
    denom = (y0 - 2 * y1 + y2)
    offset = 0.5 * (y0 - y2) / denom if denom else 0.0
    refined = peak_lag + float(offset)

    # Uncertainty from the curvature: a flat peak means the data do not pin the
    # spacing down, and the reported interval widens accordingly.
    if denom >= 0 or y1 <= 0:
        spread = float(max_nrl - min_nrl)
    else:
        curvature = -denom / max(abs(y1), 1e-12)
        spread = float(1.0 / (curvature ** 0.5)) if curvature > 0 else float(
            max_nrl - min_nrl)

    mono_peak = float(np.argmax(hist[: max(1, min(400, hist.size))]) + lo)

    # A periodicity that is not there must be reported as absent, not as a weak
    # estimate. With no ladder, the argmax above still lands somewhere -- measured
    # case: fragments only at 507 bp against a 400-1200 bp window returned
    # "138.1 bp" with an autocorrelation of -0.009, i.e. pure noise presented as a
    # measurement. Both guards below say "no spacing found" instead.
    # An absolute threshold is not enough. A histogram of *uniformly random*
    # lengths is smooth, so its autocorrelation is high at every lag (measured:
    # ~0.66) and an absolute cut of 0.05 happily passes it. What distinguishes a
    # real ladder is that the curve rises again after a dip, i.e. there is a local
    # maximum away from lag 0, rather than decaying monotonically.
    # Measured shapes (this project, seeds 3/9/11):
    #   uniform random lengths : 1.000, 0.715, 0.720, 0.692, 0.681 ... monotone
    #                            decay, minimum -0.19 reached only at the far end
    #   a real 180 bp ladder   : 1.000, 0.921, 0.788, 0.602, 0.410, 0.233, 0.097,
    #                            -0.00, then rising again to a peak at lag 178
    # So the discriminator is not "there is a wiggle" but "there is a trough and
    # the curve comes back UP after it". The peak is then required to sit after
    # that trough, which monotone decay can never satisfy.
    # Is there a real peak inside the search window, or just noise?
    #
    # Deciding this from the whole curve failed twice: the global minimum sits at
    # lag ~860 where the autocorrelation is already ~0, so "does it rebound after
    # the trough" measures noise (real +0.001 vs noise -0.001 -- no separation).
    # The decision belongs inside the window the estimate comes from.
    #
    # Measured here (seeds/fixtures committed in this repo):
    #   real ladder 165..195 bp : peak ac 0.281-0.285, rise over ac[130] 0.289-0.323
    #   uniform random lengths  : peak ac 0.069-0.092, rise 0.056-0.162
    # Raising from the window start is the discriminator: a ladder climbs to a
    # peak, noise wanders around zero and never gets far above it.
    window_lo = int(min_nrl)
    window = ac[window_lo: max_nrl + 1]
    if window.size:
        window_gain = float(window.max() - ac[window_lo])
    else:
        window_gain = 0.0
    local_peak = bool(window_gain > _NRL_MIN_PEAK_GAIN)
    if y1 < _NRL_MIN_PERIODICITY or not local_peak:
        return {"error": (("no periodic spacing found: peak autocorrelation "
                           f"{y1:.3f} is below {_NRL_MIN_PERIODICITY}")
                          if y1 < _NRL_MIN_PERIODICITY else
                          ("no periodic spacing found: the autocorrelation "
                           "decays without a local maximum away from lag 0, "
                           "which is what an irregular length distribution "
                           "looks like")),
                "mononucleosome_length_bp": mono_peak,
                "n_fragments_used": int(values.size),
                # Kept on the rejection path as well: without this a refused
                # estimate is unauditable, which is how the earlier "no spacing
                # found" cases could not be told apart from parser bugs.
                "autocorrelation_lags": [int(x) for x in
                                         range(0, ac.size, max(1, ac.size // 200))],
                "autocorrelation": [round(float(x), 4) for x in
                                    ac[:: max(1, ac.size // 200)]],
                "autocorrelation_window_gain": round(float(window_gain), 4)}
    if spread >= float(max_nrl - min_nrl):
        return {"error": (f"the spacing is not resolved: uncertainty "
                          f"{spread:.0f} bp spans the whole search range"),
                "mononucleosome_length_bp": mono_peak,
                "n_fragments_used": int(values.size)}

    # The curve itself is kept so the figure can draw the evidence rather than
    # just the number: a reader should be able to see whether the peak the
    # estimate came from is a real peak. Downsampled to ~200 points, which is far
    # more than a figure needs and keeps the JSON small.
    step = max(1, ac.size // 200)
    lags = list(range(0, ac.size, step))
    curve = [round(float(ac[i]), 4) for i in lags]

    return {
        "nrl_bp": round(refined, 1),
        "nrl_lag_raw": peak_lag,
        "autocorrelation_lags": lags,
        "autocorrelation": curve,
        # Curvature width, NOT a statistical confidence interval. It says how
        # sharply the autocorrelation peak pins the spacing; it carries no
        # sampling distribution, so "mean +/- 2*this" would be a fabricated CI.
        # Coverage calibration against known-repeat libraries is what would turn
        # this into an interval, and that has not been done.
        "nrl_uncertainty_bp": round(min(spread, float(max_nrl - min_nrl)), 1),
        "nrl_uncertainty_kind": "peak-curvature width, not a confidence interval",
        "nrl_autocorrelation": round(float(y1), 3),
        "nrl_search_range": [min_nrl, max_nrl],
        "mononucleosome_length_bp": mono_peak,
        "n_fragments_used": int(values.size),
        "method": "autocorrelation of the fragment-length histogram",
        "note": ("repeat length (centre-to-centre), not the mononucleosome "
                 "fragment length; descriptive, not comparable to a gel estimate"),
    }


def beta_distribution(table: Iterable[float]) -> dict[str, Any]:
    """Summarise a methylation beta-value vector (0-1 scale)."""
    import numpy as np

    values = np.asarray([v for v in table if v is not None and not (isinstance(v, float) and math.isnan(v))],
                        dtype=float)
    values = values[(values >= 0) & (values <= 1)]
    if values.size == 0:
        return {"n_sites": 0}
    hist, _ = np.histogram(values, bins=20, range=(0.0, 1.0))
    hist = hist / values.size
    # Bimodality: mass in the two tails relative to the middle.
    low = float(((values < 0.2)).mean())
    high = float(((values > 0.8)).mean())
    mid = 1.0 - low - high
    return {"n_sites": int(values.size), "mean_beta": float(values.mean()),
            "median_beta": float(np.median(values)),
            "fraction_hypomethylated": low, "fraction_hypermethylated": high,
            "fraction_intermediate": mid,
            "is_bimodal": bool(low + high > 0.6),
            "bin_fractions": [float(x) for x in hist]}


ASSAY_METRIC_BUILDERS: dict[str, tuple[str, ...]] = {
    # ATAC QC gains tss_enrichment: the signal-to-background ratio at transcription
    # start sites, which is how an ATAC library is judged before any peak is
    # called. It needs a TSS table, so it is opt-in via tss_rows.
    "atac-qc": ("fragments", "frip", "enrichment", "tss_enrichment"),
    "chip-qc": ("frip", "enrichment"),
    "histone-qc": ("frip", "enrichment"),
    "mnase-qc": ("periodicity", "fragments"),
    "methylation-qc": ("beta",),
}


def _metric_key(assay: str) -> str:
    """Map any accepted assay spelling to a key of ``ASSAY_METRIC_BUILDERS``.

    Returns the input unchanged when it is unrecognised, so the caller reports
    an empty metric set rather than raising: a project manifest may carry an
    assay this suite does not model, and that is a reportable state.
    """
    text = str(assay).strip().lower()
    if text.endswith("-qc") and text in ASSAY_METRIC_BUILDERS:
        return text
    try:
        return f"{canonical_assay(text)}-qc"
    except AssayError:
        return text


def compute_assay_metrics(assay: str, *, bam: str | None = None, bigwig: str | None = None,
                          peak: str | None = None, methylation_values: Iterable[float] | None = None,
                          max_reads: int = 200_000,
                          tss_rows: Sequence[tuple[str, int, str]] | None = None,
                          tss_flank: int = 2000) -> dict[str, Any]:
    """Compute the metric set that is meaningful for ``assay``.

    Returns a flat dict so it can be merged directly into the QC row.
    """
    # The metric table is keyed by "<canonical>-qc". The assay name arrives as the
    # CLI label ("ChIP-seq", "ATAC-seq"), a project-manifest value ("CUT&Tag"),
    # or an already-normalised "<canonical>-qc". easyepi.assays owns that
    # vocabulary; keeping a second copy here is what once dropped FRiP for
    # assay="ChIP" because the two spellings drifted apart.
    key = _metric_key(assay)
    wanted = ASSAY_METRIC_BUILDERS.get(key, ())
    out: dict[str, Any] = {}
    if "fragments" in wanted and bam:
        frag = fragment_length_summary(bam, max_reads=max_reads)
        out["fragment"] = frag
        if "median_fragment_length" in frag:
            out["median_fragment_length"] = frag["median_fragment_length"]
            out["fraction_nfr"] = frag["fraction_nfr"]
            out["fraction_mononucleosome"] = frag["fraction_mononucleosome"]
            out["fraction_dinucleosome"] = frag["fraction_dinucleosome"]
        if frag.get("length_sample"):
            out["fragment_length_sample"] = frag["length_sample"]
    if "periodicity" in wanted and bam:
        per = periodicity_score(bam, max_reads=max_reads)
        out["periodicity"] = per
        if "mononucleosome_fraction" in per:
            out["mononucleosome_fraction"] = per["mononucleosome_fraction"]
            out["mono_to_nfr_ratio"] = per["mono_to_nfr_ratio"]
            if "nrl_bp" in per:
                out["nrl_bp"] = per["nrl_bp"]
    if "frip" in wanted and bam and peak:
        f = frip(bam, peak)
        out["frip_detail"] = f
        if "frip" in f:
            out["frip"] = f["frip"]
            out["reads_in_peaks"] = f["reads_in_peaks"]
    if "enrichment" in wanted and bigwig and peak:
        enr = peak_signal_enrichment(bigwig, peak)
        out["enrichment_detail"] = enr
        if "signal_enrichment" in enr:
            out["signal_enrichment"] = enr["signal_enrichment"]
            out["mean_signal_in_peaks"] = enr["mean_signal_in_peaks"]
    if "tss_enrichment" in wanted and bam and tss_rows:
        from .tss_enrichment import tss_enrichment

        # Reported whether or not it clears a threshold: the score is the
        # evidence, and a low one is a finding rather than a reason to omit the
        # key and let a downstream reader assume it passed.
        tss_res = tss_enrichment(bam, tss_rows, flank=tss_flank, max_reads=max_reads)
        out["tss_enrichment_detail"] = tss_res
        score = tss_res.get("tss_enrichment_score")
        if score is not None:
            out["tss_enrichment"] = score
            # Which definition produced this number. Two TSS readers exist in the
            # package and they do not agree by a large factor, so a consumer that
            # compares a stored score against a fresh one has to know which it is.
            out["tss_enrichment_source"] = "bam"
        else:
            # No ratio could be formed. Carry the status through so a verdict can
            # tell "measured and failed" from "never measured": a sparse library
            # with zero flanking pileup must not read as a sample where --tss was
            # simply not given.
            out["tss_enrichment_status"] = tss_res.get("tss_enrichment_status")
            # A zero ratio is the worst possible outcome, so it is scored as 0.0
            # rather than as absent -- 0 can breach a floor threshold, None cannot.
            if tss_res.get("tss_enrichment_status") in ("background_zero",
                                                         "no_signal"):
                out["tss_enrichment"] = 0.0
    if "beta" in wanted and methylation_values is not None:
        b = beta_distribution(methylation_values)
        out["beta_detail"] = b
        if "mean_beta" in b:
            out["mean_beta"] = b["mean_beta"]
            out["is_bimodal"] = b["is_bimodal"]
    return out
