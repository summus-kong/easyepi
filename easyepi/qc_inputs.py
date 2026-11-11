"""Resolve QC inputs from explicit flags, comma lists and a scanned folder.

Positional pairing of ``--bam`` and ``--bigwig`` used to rely on list order, which
silently mismatches as soon as one sample lacks a track or the flags are written in
a different order. Inputs are therefore resolved to a ``{sample_id: paths}`` map
first, and every track is attached by sample name.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from easyepi.formats import (BAM_EXTS, BIGWIG_EXTS, METHYL_EXTS, PEAK_EXTS,
                            stem_of)


def split_values(items: Sequence[str] | None) -> list[str]:
    """Flatten repeated and comma-separated flag values."""
    out: list[str] = []
    for raw in items or []:
        for piece in str(raw).split(","):
            value = piece.strip()
            if value:
                out.append(value)
    return out


def _stem(path: str | Path, exts: Sequence[str]) -> str:
    """Sample id for a file: the name with a known track extension removed.

    ``a.bed.gz`` must reduce to ``a``, not ``a.bed``, so longer extensions are
    stripped first.
    """
    name = Path(path).name
    for ext in sorted(exts, key=len, reverse=True):
        if name.endswith(ext):
            return name[: -len(ext)]
    return Path(path).stem


def scan_folder(directory: str | Path) -> dict[str, dict[str, str]]:
    """Group track files in one folder into ``{sample: {kind: path}}``.

    Files are matched by extension only; nothing is inferred from a filename
    prefix. A sample that ends up with two files of the same kind is an error,
    because picking one arbitrarily would analyse the wrong track.
    """
    base = Path(directory)
    if not base.is_dir():
        raise NotADirectoryError(f"--input-dir is not a directory: {directory}")

    samples: dict[str, dict[str, str]] = {}
    for path in sorted(base.iterdir()):
        if not path.is_file():
            continue
        for kind, exts in (("bam", BAM_EXTS), ("bigwig", BIGWIG_EXTS),
                           ("peak", PEAK_EXTS), ("methylation", METHYL_EXTS)):
            if not any(path.name.endswith(ext) for ext in exts):
                continue
            # a .bai index is not a sample
            if path.name.endswith((".bai", ".csi")):
                break
            sample = _stem(path, exts)
            entry = samples.setdefault(sample, {})
            if kind in entry:
                raise ValueError(
                    f"sample '{sample}' has two {kind} files: "
                    f"{Path(entry[kind]).name} and {path.name}")
            entry[kind] = str(path)
            break
    if not samples:
        raise ValueError(
            f"no BAM/BigWig/peak files found in {directory}")
    return samples


def resolve_inputs(args) -> dict[str, dict[str, str]]:
    """Merge flag inputs and a scanned folder into ``{sample: {kind: path}}``.

    Explicit flags win over a folder scan for the same sample and kind, so a user
    can scan a directory and then override one file.
    """
    resolved: dict[str, dict[str, str]] = {}

    # ``--project`` is the documented way to run QC over a whole study. It was
    # never read here, so ``qc --project`` reported "no inputs". The project is
    # resolved first and explicit flags may still override individual tracks.
    project_path = getattr(args, "project", None)
    if project_path:
        from easyepi.project import load_project

        config = load_project(project_path)
        if not config.samples:
            raise ValueError(f"{project_path} defines no samples")
        for sample in config.samples:
            entry = resolved.setdefault(sample.sample_id, {})
            for kind, value in (("bam", sample.bam), ("bigwig", sample.bigwig),
                                ("peak", sample.peak),
                                ("methylation", sample.methylation)):
                if value:
                    entry[kind] = str(value)

    for kind, attr, exts in (("bam", "bam", BAM_EXTS),
                             ("bigwig", "bigwig", BIGWIG_EXTS),
                             ("peak", "peak", PEAK_EXTS)):
        for path in split_values(getattr(args, attr, None)):
            if not Path(path).exists():
                raise FileNotFoundError(f"--{attr} file not found: {path}")
            sample = _stem(path, exts)
            entry = resolved.setdefault(sample, {})
            if kind in entry and entry[kind] != path:
                raise ValueError(
                    f"sample '{sample}' given two {kind} files: "
                    f"{entry[kind]} and {path}")
            entry[kind] = path

    directory = getattr(args, "input_dir", None)
    if directory:
        for sample, tracks in scan_folder(directory).items():
            entry = resolved.setdefault(sample, {})
            for kind, path in tracks.items():
                entry.setdefault(kind, path)

    if not resolved:
        raise ValueError(
            "no inputs: pass --bam/--bigwig/--peak or --input-dir")
    return resolved
