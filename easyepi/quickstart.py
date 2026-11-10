"""One-command and guided entry points.

Why this module exists
----------------------
The individual subcommands are correct but demanding: ``footprint`` needs six
arguments, ``nucleosome-positioning`` five, and every path has to be spelled out.
That is fine for a scripted pipeline and hostile for a first run. These entry
points sit *on top of* the same implementations and change nothing about them:

* ``easyepi init``  — write a project file from a folder of BAM/BigWig/peak files
  by detecting what is there, so the user edits one file instead of building a
  command line.
* ``easyepi quick`` — run the standard chain (QC -> peaks -> annotation -> motifs
  -> report) from a project file or a folder, skipping steps whose inputs are
  absent, and print one summary.

Neither hides a decision: every step that is skipped or defaulted is named in the
report. ``quick`` never invents values; when it cannot proceed it says which input
is missing.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

ASSAY_HINTS = (
    # (regex on filename, assay label); first match wins, longest/most specific first
    (r"(?i)(^|[^a-z])(atac)([^a-z]|$)", "ATAC-seq"),
    (r"(?i)(^|[^a-z])cut(&|and)?tag([^a-z]|$)", "CUT&Tag"),
    (r"(?i)(^|[^a-z])cut(&|and)?run([^a-z]|$)", "CUT&RUN"),
    (r"(?i)(^|[^a-z])mnase([^a-z]|$)", "MNase-seq"),
    (r"(?i)(^|[^a-z])(h3k27me3|h3k9me3|h3k36me3)([^a-z]|$)", "histone"),
    (r"(?i)(^|[^a-z])(h3k4me1|h3k4me3|h3k27ac|h3k9ac|h2az)([^a-z]|$)", "histone"),
    (r"(?i)(^|[^a-z])(meth|wgbs|wgbbs|bisulfite|cov)([^a-z]|$)", "methylation"),
    (r"(?i)(^|[^a-z])(chip|tf)([^a-z]|$)", "ChIP-seq"),
)

BAM_EXTS = (".bam",)
BIGWIG_EXTS = (".bw", ".bigWig", ".bigwig")
PEAK_EXTS = (".narrowPeak", ".broadPeak", ".bed", ".bed.gz")
METHYL_EXTS = (".cov", ".bedMethyl")

GENOME_ALIASES = {
    "hg38": "hg38", "grch38": "hg38", "hg19": "hg19", "grch37": "hg19",
    "mm10": "mm10", "grcm38": "mm10", "mm39": "mm39", "grcm39": "mm39",
}


def guess_assay(name: str) -> str:
    """Infer an assay label from a sample/file name.

    Returns the first matching hint, else ``"unknown"``. The caller must surface
    the guess so a wrong inference is visible rather than silently propagated into
    QC thresholds, which differ per assay.
    """
    for pattern, assay in ASSAY_HINTS:
        if re.search(pattern, name):
            return assay
    return "unknown"


def guess_genome(text: str) -> str | None:
    """Resolve a genome build from free text (a path, a filename, or a value)."""
    lowered = str(text).lower()
    for alias, build in GENOME_ALIASES.items():
        if alias in lowered:
            return build
    return None


def _strip_ext(name: str, exts: Sequence[str]) -> str:
    for ext in sorted(exts, key=len, reverse=True):
        if name.endswith(ext):
            return name[: -len(ext)]
    return Path(name).stem


def scan_samples(directory: str | Path) -> dict[str, dict[str, str]]:
    """Group files in one folder into ``{sample: {kind: path}}``.

    Extension-only matching; a sample with two files of the same kind is an error,
    because choosing one arbitrarily would analyse the wrong track.
    """
    base = Path(directory)
    if not base.is_dir():
        raise NotADirectoryError(f"not a directory: {directory}")
    out: dict[str, dict[str, str]] = {}
    for path in sorted(base.iterdir()):
        if not path.is_file() or path.name.endswith((".bai", ".csi", ".tbi")):
            continue
        for kind, exts in (("bam", BAM_EXTS), ("bigwig", BIGWIG_EXTS),
                           ("peak", PEAK_EXTS), ("methylation", METHYL_EXTS)):
            if not any(path.name.endswith(e) for e in exts):
                continue
            sample = _strip_ext(path.name, exts)
            entry = out.setdefault(sample, {})
            if kind in entry:
                raise ValueError(
                    f"sample '{sample}' has two {kind} files: "
                    f"{Path(entry[kind]).name} and {path.name}")
            entry[kind] = str(path)
            break
    if not out:
        raise ValueError(f"no BAM/BigWig/peak files found in {directory}")
    return out


def build_project(
    directory: str | Path,
    *,
    project_id: str | None = None,
    genome_build: str | None = None,
    annotation: str | None = None,
    motif_database: str | None = None,
    assays: Mapping[str, str] | None = None,
) -> tuple[dict[str, Any], list[str]]:
    """Build a project dict from a folder, plus the notes a user must see.

    Returns ``(project, notes)``. ``notes`` records every guess (assay, genome) so
    the caller can print them; nothing inferred is presented as certain.
    """
    samples = scan_samples(directory)
    notes: list[str] = []

    detected_genome = genome_build or guess_genome(str(directory))
    if detected_genome is None and annotation:
        detected_genome = guess_genome(str(annotation))
    if detected_genome:
        notes.append(f"genome build inferred: {detected_genome}")
    else:
        notes.append("genome build not detected; set genome_build before peak calling")

    rows = []
    for sample, tracks in samples.items():
        assay = (assays or {}).get(sample) or guess_assay(sample)
        if assay == "unknown":
            notes.append(f"{sample}: assay not detected; set 'assay' for this sample")
        rows.append({"sample_id": sample, "assay": assay, **tracks})

    project: dict[str, Any] = {
        "project": {
            "project_id": project_id or Path(directory).resolve().name,
            "root": str(Path(directory).resolve()),
        },
        "samples": rows,
    }
    if detected_genome:
        project["project"]["genome_build"] = detected_genome
    refs: dict[str, str] = {}
    if annotation:
        refs["annotation"] = str(annotation)
    if motif_database:
        refs["motif_database"] = str(motif_database)
    if refs:
        project["references"] = refs
    return project, notes


def _toml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def project_to_toml(project: Mapping[str, Any]) -> str:
    """Serialise a project dict to TOML without a third-party writer.

    Written by hand because ``tomli-w`` is not a dependency and the schema is
    small: one ``[project]``/``[references]`` table plus repeated ``[[samples]]``.
    """
    lines: list[str] = ["# Generated by `easyepi init`; edit paths and assay labels."]
    project_table = dict(project.get("project") or {})
    if project_table:
        lines.append("")
        lines.append("[project]")
        for key, value in project_table.items():
            lines.append(f"{key} = {_toml_scalar(value)}")
    refs = dict(project.get("references") or {})
    if refs:
        lines.append("")
        lines.append("[references]")
        for key, value in refs.items():
            lines.append(f"{key} = {_toml_scalar(value)}")
    defaults = dict(project.get("defaults") or {})
    if defaults:
        lines.append("")
        lines.append("[defaults]")
        for key, value in defaults.items():
            lines.append(f"{key} = {_toml_scalar(value)}")
    for sample in project.get("samples") or []:
        lines.append("")
        lines.append("[[samples]]")
        for key, value in sample.items():
            if value is None:
                continue
            lines.append(f"{key} = {_toml_scalar(value)}")
    return "\n".join(lines) + "\n"
