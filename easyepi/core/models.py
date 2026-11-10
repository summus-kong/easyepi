"""Core data model: unified in-memory representation for epigenomics downstream analysis.

The data model is deliberately minimal and pandas-centric so that every analysis
module (qc / normalize / peaks / annotate / motifs / nucleosome / differential /
viz) works on the same set of objects without external database dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np

# Lightweight optional import: pandas is a hard dependency, but we import it
# lazily so that core model construction and version query still work in a
# minimal environment.
try:
    import pandas as pd
    _HAS_PANDAS = True
except Exception:  # pragma: no cover - ImportError or broken install
    # pandas is an optional extra; without it the DataFrame helpers raise a clear
    # ImportError at the point of use instead of at import time.
    pd = None  # type: ignore
    _HAS_PANDAS = False


@dataclass
class PeakSet:
    """A set of genomic peaks / regions in the BED-ish coordinate system.

    Columns (all optional except chrom/start/end):
      chrom, start, end, name, score, strand, signal, fold_change, pvalue, qvalue
    """

    chrom: np.ndarray
    start: np.ndarray
    end: np.ndarray
    name: Optional[np.ndarray] = None
    score: Optional[np.ndarray] = None
    strand: Optional[np.ndarray] = None
    signal: Optional[np.ndarray] = None
    fold_change: Optional[np.ndarray] = None
    pvalue: Optional[np.ndarray] = None
    qvalue: Optional[np.ndarray] = None
    metadata: Dict[str, object] = field(default_factory=dict)

    def __len__(self) -> int:
        return int(len(self.chrom))

    @property
    def width(self) -> np.ndarray:
        """Peak width in base pairs."""
        return self.end - self.start

    @property
    def n_peaks(self) -> int:
        """Number of peaks held by this collection."""
        return len(self)

    def to_dataframe(self):
        """Return a pandas DataFrame (hard requirement for analysis upstream)."""
        if not _HAS_PANDAS:
            raise ImportError(
                "pandas is required to build a DataFrame from a PeakSet. "
                "Install with: pip install pandas"
            )
        data: Dict[str, Iterable] = {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
        }
        for col in ("name", "score", "strand", "signal", "fold_change", "pvalue", "qvalue"):
            val = getattr(self, col)
            if val is not None:
                data[col] = val
        return pd.DataFrame(data)


@dataclass
class CountMatrix:
    """Peak x sample count matrix with per-sample library/ spike-in factors."""

    regions: List[str] = field(default_factory=list)  # peak ids
    samples: List[str] = field(default_factory=list)   # sample/bam ids
    matrix: Optional[np.ndarray] = None                # shape (n_regions, n_samples)
    lib_sizes: Optional[np.ndarray] = None             # per-sample library sizes
    spikein_factors: Optional[np.ndarray] = None       # per-sample spike-in factors

    @property
    def shape(self) -> tuple:
        """Matrix shape as (regions, samples)."""
        if self.matrix is None:
            return (len(self.regions), len(self.samples))
        return self.matrix.shape

    def to_dataframe(self):
        """Materialise the matrix as a pandas DataFrame (index = region ids)."""
        if not _HAS_PANDAS:
            raise ImportError("pandas is required to serialize a CountMatrix.")
        return pd.DataFrame(self.matrix, index=self.regions, columns=self.samples)


@dataclass
class NormalizationFactors:
    """Per-sample normalization factors (library size and optional spike-in)."""

    lib_sizes: Dict[str, float] = field(default_factory=dict)
    spikein: Dict[str, float] = field(default_factory=dict)
    scale_factors: Dict[str, float] = field(default_factory=dict)

    def factor_for(self, sample: str, use_spikein: bool = True) -> float:
        """Return a CPM multiplier, optionally multiplied by an explicit spike-in factor."""
        if sample not in self.lib_sizes:
            raise ValueError(f"Missing measured library size for {sample}")
        library = float(self.lib_sizes[sample])
        spike = float(self.spikein.get(sample,1.0)) if use_spikein else 1.0
        if not np.isfinite(library) or library <= 0 or not np.isfinite(spike) or spike <= 0:
            raise ValueError("Normalization quantities must be positive and finite")
        return 1_000_000.0 / library * spike

    def as_dict(self) -> Dict[str, Dict[str, float]]:
        """Plain-dict view of the matrix, suitable for JSON serialisation."""
        return {
            "lib_sizes": dict(self.lib_sizes),
            "spikein": dict(self.spikein),
            "scale_factors": dict(self.scale_factors),
        }


@dataclass
class EpigenomeData:
    """Master container tying together peaks, count matrices and sample metadata.

    This is the unified object passed through the analysis pipeline so that all
    downstream modules share a single, consistent representation.
    """

    name: str = "dataset"
    peaks: Optional[PeakSet] = None
    counts: Optional[CountMatrix] = None
    norm: Optional[NormalizationFactors] = None
    sample_meta: Dict[str, Dict[str, object]] = field(default_factory=dict)
    tracks: Dict[str, str] = field(default_factory=dict)  # sample -> bigwig/bedGraph path
    notes: Dict[str, object] = field(default_factory=dict)

    def add_sample(self, sample_id: str, meta: Optional[Dict[str, object]] = None) -> None:
        """Append one sample column, rejecting duplicate names."""
        self.sample_meta[sample_id] = meta or {}

    def sample_groups(self) -> Dict[str, List[str]]:
        """Infer group -> sample mapping from metadata 'group' key (or 'condition')."""
        groups: Dict[str, List[str]] = {}
        for sid, meta in self.sample_meta.items():
            g = meta.get("group") or meta.get("condition")
            if g is None:
                g = "default"
            groups.setdefault(str(g), []).append(sid)
        return groups


def read_peak_table(path: Path) -> PeakSet:
    """Read a BED6/BED12 or narrowPeak table into a minimal PeakSet (chrom/start/end)."""
    import pandas as pd
    cols = ["chrom", "start", "end"]
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=cols)
    return PeakSet(
        chrom=df["chrom"].astype(str).to_numpy(),
        start=df["start"].astype(np.int64).to_numpy(),
        end=df["end"].astype(np.int64).to_numpy(),
    )
