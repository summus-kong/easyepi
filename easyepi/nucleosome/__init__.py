"""Nucleosome occupancy / positioning analysis (native, dependency-free).

Provides occupancy profiling, spacing analysis and a Fourier-based periodicity
estimate for MNase-seq / ATAC-seq fragment data. It is the differentiating
module of the suite: it couples nucleosome phasing with accessibility (ATAC)
and TF footprints.

Importable API:
  * ``occupancy_profile`` / ``periodicity`` — pure numpy/scipy.
  * ``nucleosome_scores`` / ``nfr_score`` — occupancy-based scoring.
  * ``bigwig_signal`` — read a bigwig into a per-bin array (needs pyBigWig).
  * ``run_nucleosome`` — CLI wrapper (reads a matrix/bigwig, writes a report).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


def occupancy_profile(signal: np.ndarray, bin_size: int = 1) -> np.ndarray:
    """Return a normalized occupancy profile (z-scored) of a 1-D signal."""
    sig = np.asarray(signal, dtype=np.float64)
    if sig.size == 0:
        return sig
    mu, sd = float(sig.mean()), float(sig.std())
    if sd <= 0:
        return np.zeros_like(sig)
    return (sig - mu) / sd


def periodicity(signal: np.ndarray, fs_hint: float = 0.1) -> dict:
    """Estimate the dominant spatial periodicity of an occupancy signal.

    Uses the power spectral density (Welch) and returns the dominant period in
    signal bins. For MNase nucleosome arrays the ~10.5 bp phasing manifests as a
    periodic component; with bin_size given in bp this maps to the spacing.
    """
    from scipy.signal import welch
    sig = np.asarray(signal, dtype=np.float64)
    if sig.size < 8:
        return {"dominant_period_bins": None, "n_bins": int(sig.size)}
    nperseg = min(sig.size, 256)
    f, psd = welch(sig, fs=fs_hint, nperseg=nperseg)
    idx = np.argmax(psd[1:]) + 1
    period = 1.0 / f[idx] if f[idx] > 0 else None
    return {"dominant_period_bins": float(period) if period else None,
            "n_bins": int(sig.size)}


def nucleosome_repeat(signal: np.ndarray, bin_size: int = 1,
                     lo_bp: float = 120.0, hi_bp: float = 400.0) -> dict:
    """Nucleosome repeat length from a signal, with an explicit verdict.

    ``periodicity`` returns whatever the PSD peaks at. That is only meaningful
    when the input is phased: a whole-chromosome occupancy track carries a 1/f
    trend and heterogeneous spacing, so the peak lands on the lowest bins. On
    GSE247785 (5 Mb of real MNase data) the raw PSD reported 1280 bp for every
    window size from 20 kb to 5 Mb, and a detrended autocorrelation found no
    ~185 bp peak at all.

    So the number is returned together with a plain-language verdict, and a
    caller must not quote ``period_bp`` without checking ``reliable``.
    """
    sig = np.asarray(signal, dtype=np.float64)
    bins = periodicity(sig, fs_hint=1.0) if sig.size >= 8 else {
        "dominant_period_bins": None, "n_bins": int(sig.size)}
    bins_per = bins.get("dominant_period_bins")
    period_bp = round(float(bins_per) * int(bin_size), 1) if bins_per else None
    reliable = bool(period_bp is not None and lo_bp <= period_bp <= hi_bp)
    out = {"dominant_period_bins": bins_per,
           "dominant_period_bp": period_bp,
           "search_range_bp": [lo_bp, hi_bp],
           "reliable": reliable,
           "method": "welch PSD argmax over the whole supplied signal"}
    if not reliable:
        out["note"] = ("no repeat peak inside {0:g}-{1:g} bp; measure the repeat "
                       "on aligned regions or a metagene profile, not on an "
                       "arbitrary interval".format(lo_bp, hi_bp))
    return out


def nucleosome_scores(signal: np.ndarray, bin_size: int = 1,
                      nfr_quantile: float = 0.2) -> dict:
    """Compute per-bin nucleosome occupancy summary scores.

    Returns the z-scored occupancy, the mean occupancy of above-mean bins, and
    the nucleosome-free fraction.

    The NFR fraction is defined by a *quantile*, not by a fixed z threshold.
    MNase coverage is strongly right-skewed (on GSE247785 chr1:0-5 Mb, 38-82% of
    10 bp bins are exactly zero and the tail reaches 17x), so after z-scoring
    the minimum is only about -0.9 and a ``z < -1`` rule reports 0% nucleosome
    free regions on every real input. The lowest ``nfr_quantile`` of bins is the
    population that is depleted relative to the rest of the same sample.
    """
    sig = np.asarray(signal, dtype=np.float64)
    if sig.size == 0:
        return {"occupancy": sig, "nfr_score": np.nan, "mean_occupancy": np.nan}
    occ = occupancy_profile(sig)
    mean_occ = float(occ[occ > 0].mean()) if (occ > 0).any() else float(occ.mean())
    if not (0.0 < nfr_quantile < 1.0):
        raise ValueError("nfr_quantile must be in (0, 1)")
    cutoff = float(np.quantile(sig, nfr_quantile))
    nfr = float((sig <= cutoff).mean())
    return {"occupancy": occ, "nfr_score": nfr,
            "nfr_cutoff": cutoff, "nfr_quantile": float(nfr_quantile),
            "mean_occupancy": mean_occ, "n_bins": int(occ.size)}


def bigwig_signal(path: Path, chrom: str, start: int, end: int,
                  bins: int = 2000) -> np.ndarray:
    """Read a bigwig region into a per-bin signal array (needs pyBigWig)."""
    try:
        import pyBigWig
    except ImportError as exc:  # pragma: no cover
        raise ImportError("pyBigWig is required to read bigwig (pip install pyBigWig)") from exc
    bw = pyBigWig.open(str(path))
    try:
        vals = bw.values(chrom, int(start), int(end))
        if vals is None:
            return np.zeros(bins, dtype=np.float64)
        arr = np.array([0.0 if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)
                        for v in vals], dtype=np.float64)
        arr[arr < 0] = 0.0
        if len(arr) == bins:
            return arr
        if len(arr) == 0:
            return np.zeros(bins, dtype=np.float64)
        idx = np.linspace(0, len(arr) - 1, bins)
        return np.interp(idx, np.arange(len(arr)), arr)
    finally:
        bw.close()


def run_nucleosome(args):
    """Occupancy and NFR summary for a bigWig region.

    This is the signal-track counterpart to ``nucleosome-positioning``, which
    needs paired-end BAM and reports dyad calls. Here the input is an existing
    occupancy track, so the output is descriptive: a z-scored profile, the
    nucleosome-free fraction, and a repeat estimate that states whether it is
    trustworthy.
    """
    import json

    bigwig = getattr(args, "bigwig", None)
    if not bigwig:
        raise ValueError("nucleosome needs --bigwig (an occupancy track); "
                         "for dyad calls from paired-end BAM use nucleosome-positioning")
    chrom, start, end = args.chrom, int(args.start), int(args.end)
    if end <= start:
        raise ValueError("--end must be greater than --start")

    bin_size = int(getattr(args, "bin_size", 10) or 10)
    n_bins = max(1, (end - start) // bin_size)
    sig = bigwig_signal(Path(bigwig), chrom, start, end, bins=n_bins)
    if sig.size == 0:
        raise ValueError("bigWig returned no values for {0}:{1}-{2}".format(chrom, start, end))

    scores = nucleosome_scores(sig, bin_size=bin_size)
    repeat = nucleosome_repeat(sig, bin_size=bin_size)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    report = {
        "bigwig": str(bigwig),
        "region": "{0}:{1}-{2}".format(chrom, start, end),
        "bin_size": bin_size,
        "n_bins": int(sig.size),
        "mean_occupancy_z": round(float(np.asarray(scores["occupancy"]).mean()), 6),
        "nfr_score": round(float(scores["nfr_score"]), 6),
        "nfr_definition": "{0:g} lowest-quantile bins of the raw signal".format(
            scores["nfr_quantile"]),
        "nucleosome_repeat": repeat,
        "method": "z-scored occupancy + nucleosome-free fraction by quantile",
    }
    (out / "nucleosome_summary.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print("[easyepi] nucleosome: {0} bins in {1}; NFR fraction {2:.3f}".format(
        sig.size, report["region"], report["nfr_score"]))
    if repeat["reliable"]:
        print("[easyepi] repeat {0} bp".format(repeat["dominant_period_bp"]))
    else:
        print("[easyepi] repeat NOT reported: {0}".format(repeat["note"]))
    return {
        "n_bins": int(sig.size),
        "nfr_score": round(float(scores["nfr_score"]), 6),
        "repeat_bp": repeat["dominant_period_bp"],
        "repeat_reliable": repeat["reliable"],
    }
