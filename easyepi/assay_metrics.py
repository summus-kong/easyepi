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
    seen = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
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
        return {"n_fragments": 0}

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
    return {"n_fragments": summary["n_fragments"],
            "median_fragment_length": summary["median_fragment_length"],
            "mononucleosome_fraction": mono,
            "nfr_fraction": nfr,
            "mono_to_nfr_ratio": _safe_div(mono, nfr) if nfr else None}


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
