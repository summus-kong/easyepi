"""Signal / count normalization: library size, spike-in, RPGC, CPM.

Implemented with numpy so it needs no alignment library; the raw per-region
counts are passed in or read from a precomputed matrix.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

from easyepi.core.models import CountMatrix, NormalizationFactors


def cpm_normalize(matrix: np.ndarray, lib_sizes: np.ndarray) -> np.ndarray:
    """Counts-per-million normalization. lib_sizes in reads; divide by total*1e6."""
    lib = _positive_vector(matrix, lib_sizes, "library sizes")
    return matrix / lib[None, :] * 1_000_000.0


def spikein_normalize(matrix: np.ndarray, spikein_factors: np.ndarray) -> np.ndarray:
    """Normalize by user-provided spike-in scale factors (multiply)."""
    f = _positive_vector(matrix, spikein_factors, "spike-in factors")
    return matrix * f[None, :]


def spikein_factors_from_counts(spike_counts: np.ndarray) -> np.ndarray:
    """Compute spike-in scale factors from per-sample spike-in read counts.

    Higher spike-in count -> lower scale (less genomic DNA per spike read):
    factor = median(spike_counts) / spike_counts (so the median sample keeps ~1x,
    and samples with fewer spike reads are scaled up). Returns a (n_samples,) array.
    """
    sp = np.asarray(spike_counts, dtype=np.float64)
    if sp.ndim != 1 or sp.size == 0 or not np.isfinite(sp).all() or (sp <= 0).any():
        raise ValueError("spike-in counts must be positive finite measurements")
    return np.median(sp) / sp


def apply_spike_norm(matrix: np.ndarray, mode: str = "none",
                     spikein_factors: Optional[np.ndarray] = None,
                     spike_counts: Optional[np.ndarray] = None) -> np.ndarray:
    """Normalize a count matrix by spike-in, per the requested mode.

    Modes:
      * ``"none"`` (default): return the matrix unchanged.
      * ``"manual"``: use provided ``spikein_factors`` (per sample) -> scale matrix.
      * ``"auto"``: compute factors from per-sample ``spike_counts`` and normalize.
    Returns the (possibly) normalized matrix (regions x samples).
    """
    m = np.asarray(matrix, dtype=np.float64)
    if mode == "none":
        return m
    if mode == "manual":
        if spikein_factors is None:
            raise ValueError("spike-norm=manual requires --spike-factors")
        return spikein_normalize(m, spikein_factors)
    if mode == "auto":
        if spike_counts is None:
            raise ValueError("spike-norm=auto requires spike-in counts (e.g. from BAM spike reads)")
        return spikein_normalize(m, spikein_factors_from_counts(spike_counts))
    raise ValueError(f"unknown spike-norm mode '{mode}' (none/manual/auto)")


def rpgc_normalize(matrix: np.ndarray, effective_genome_size: int,
                   total_reads: np.ndarray, fragment_length: float | None = None) -> np.ndarray:
    """Scale per-base coverage to 1x genome coverage using sequenced fragment bases.

    The matrix must contain coverage values, not arbitrary region counts.
    total_reads refers to the number of contributing fragments.
    """
    reads = _positive_vector(matrix, total_reads, "fragment counts")
    if fragment_length is None or fragment_length <= 0 or effective_genome_size <= 0:
        raise ValueError("RPGC requires positive effective genome size and mean fragment length")
    return np.asarray(matrix, dtype=float) * (effective_genome_size / (reads * fragment_length))[None,:]


def _positive_vector(matrix, vector, label):
    values = np.asarray(matrix,dtype=float)
    factors = np.asarray(vector,dtype=float)
    if values.ndim != 2 or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("matrix must be finite, nonnegative and two-dimensional")
    if factors.ndim != 1 or len(factors) != values.shape[1] or not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError(f"{label} must be positive, finite and match sample columns")
    return factors


def compute_normalization_factors(samples, lib_sizes=None, spikeins=None) -> NormalizationFactors:
    """Build a NormalizationFactors object from per-sample library sizes and
    optional spike-in factors.

    Args:
        samples: list of sample ids.
        lib_sizes: optional dict sample->total reads.
        spikeins: optional dict sample->spike-in scale factor.
    Returns:
        NormalizationFactors.
    """
    lib_sizes = lib_sizes or {}
    spikeins = spikeins or {}
    factor = NormalizationFactors()
    for s in samples:
        if s not in lib_sizes:
            raise ValueError(f"Missing measured library size for {s}")
        lib = float(lib_sizes[s])
        if not np.isfinite(lib) or lib <= 0:
            raise ValueError("Library sizes must be positive and finite")
        factor.lib_sizes[s] = lib
        if s in spikeins:
            if not np.isfinite(float(spikeins[s])) or float(spikeins[s]) <= 0:
                raise ValueError("Spike-in multipliers must be positive and finite")
            factor.spikein[s] = float(spikeins[s])
            factor.scale_factors[s] = 1_000_000.0 / lib * float(spikeins[s])
        else:
            factor.scale_factors[s] = 1_000_000.0 / lib
    return factor


def run_normalize(args) -> int:
    """Normalise a counts/signal matrix with explicit factors, recording provenance."""
    import pandas as pd
    from easyepi.differential import read_count_matrix
    regions, matrix, samples = read_count_matrix(args.matrix)
    factors = pd.read_csv(args.factors, sep="\t")
    if not {"sample_id","value"}.issubset(factors.columns) or factors.sample_id.duplicated().any():
        raise ValueError("Factor TSV requires unique sample_id and value")
    if set(samples) != set(factors.sample_id):
        raise ValueError("Factor/matrix sample mismatch")
    vector = factors.set_index("sample_id").loc[samples,"value"].to_numpy(float)
    if args.method == "cpm":
        normalized = cpm_normalize(matrix,vector)
    elif args.method == "spikein":
        normalized = spikein_normalize(matrix,vector)
    else:
        normalized = rpgc_normalize(matrix,args.effective_genome_size,vector,args.fragment_length)
    out = Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(normalized,index=pd.Index(regions,name="region"),columns=samples).to_csv(out,sep="\t")
    print(f"[easyepi] normalized matrix written to {out}")
    return 0
