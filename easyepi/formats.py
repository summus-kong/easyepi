"""One definition of what file extensions mean.

Why this module exists
---------------------
``BAM_EXTS``, ``BIGWIG_EXTS``, ``PEAK_EXTS`` and ``METHYL_EXTS`` had been declared in
both ``qc_inputs`` and ``quickstart``, and the two copies had already diverged:
``quickstart`` accepted ``.cov`` as a methylation table while ``qc_inputs`` did not, so
the same folder produced different sample lists depending on which command scanned it.

Anything that classifies a file by its name goes through :func:`classify` or the
constants below, so a format added once is recognised everywhere.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

#: Track kind -> the extensions that identify it. Longer extensions must be checked
#: first (``.bed.gz`` before ``.gz``), which :func:`classify` guarantees.
KIND_EXTENSIONS: dict[str, tuple[str, ...]] = {
    "bam": (".bam",),
    "bigwig": (".bw", ".bigWig", ".bigwig"),
    "peak": (".narrowPeak", ".broadPeak", ".bed.gz", ".bed"),
    "bedgraph": (".bedGraph.gz", ".bedGraph", ".bdg"),
    "methylation": (".bedMethyl", ".cov", ".cov.gz"),
    "annotation": (".gtf.gz", ".gtf", ".gff3.gz", ".gff3", ".gff.gz", ".gff"),
    "fasta": (".fa.gz", ".fasta.gz", ".fa", ".fasta", ".fna"),
    "motif": (".pfm", ".jaspar", ".meme", ".motif"),
    "matrix": (".tsv", ".csv"),
}

#: Index and sidecar suffixes that are never an analysis input.
IGNORED_SUFFIXES: tuple[str, ...] = (".bai", ".csi", ".tbi", ".fai", ".json")

# Flat names kept for the modules that already import them.
BAM_EXTS = KIND_EXTENSIONS["bam"]
BIGWIG_EXTS = KIND_EXTENSIONS["bigwig"]
PEAK_EXTS = KIND_EXTENSIONS["peak"]
METHYL_EXTS = KIND_EXTENSIONS["methylation"]
BEDGRAPH_EXTS = KIND_EXTENSIONS["bedgraph"]

#: Longest extension per kind, so a stem can be recovered by stripping it.
_LONGEST = {kind: max(exts, key=len) for kind, exts in KIND_EXTENSIONS.items()}


def classify(name: str | Path) -> str | None:
    """Return the track kind of a filename, or ``None``.

    Longest-suffix-wins: ``a.bed.gz`` is a peak file, not an unknown ``.gz``.
    Index files return ``None`` regardless of their name.
    """
    text = str(name)
    if any(text.endswith(s) for s in IGNORED_SUFFIXES):
        return None
    best_kind, best_len = None, -1
    for kind, exts in KIND_EXTENSIONS.items():
        for ext in exts:
            if text.endswith(ext) and len(ext) > best_len:
                best_kind, best_len = kind, len(ext)
    return best_kind


def stem_of(name: str | Path) -> str:
    """Filename without its known extension, for use as a sample id.

    ``a.bed.gz`` becomes ``a``; an unknown suffix falls back to ``Path.stem``.
    """
    text = Path(name).name
    kind = classify(text)
    if kind:
        ext = _LONGEST[kind]
        if text.endswith(ext):
            return text[: -len(ext)]
    for exts in KIND_EXTENSIONS.values():
        for ext in sorted(exts, key=len, reverse=True):
            if text.endswith(ext):
                return text[: -len(ext)]
    return Path(text).stem


def extensions_for(kinds: Iterable[str]) -> tuple[str, ...]:
    """Every extension for the given kinds, longest first."""
    out: list[str] = []
    for kind in kinds:
        out.extend(KIND_EXTENSIONS.get(kind, ()))
    return tuple(sorted(set(out), key=len, reverse=True))
