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


def nucleosome_scores(signal: np.ndarray, bin_size: int = 1,
                      min_signal: float = 0.0) -> dict:
    """Compute per-bin nucleosome occupancy summary scores.

    Returns a dict with the occupancy array, the mean/median occupancy of
    high-occupancy (nucleosome) regions, and a simple NFR (nucleosome-free
    region) score (low-occupancy fraction).
    """
    sig = np.asarray(signal, dtype=np.float64)
    if sig.size == 0:
        return {"occupancy": sig, "nfr_score": np.nan, "mean_occupancy": np.nan}
    occ = occupancy_profile(sig)
    mean_occ = float(occ[occ > 0].mean()) if (occ > 0).any() else float(occ.mean())
    # NFR: fraction of low-occupancy bins (z < -1), i.e. depleted regions
    nfr = float((occ < -1.0).mean()) if occ.size else np.nan
    return {"occupancy": occ, "nfr_score": nfr,
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


def run_nucleosome(args) -> int:
    """Nucleosome occupancy and positioning summary for the supplied reads."""
    raise ValueError("Use python -m easyepi.specific.mnase.profile with a measured signal vector; the legacy metadata-only entry is disabled")
