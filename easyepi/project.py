"""Project and sample input contracts. Explicit CLI reference overrides take precedence over project values; resolved inputs are recorded for each run."""

from __future__ import annotations

import json
import os
import re
import tomllib
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping


class ProjectError(ValueError):
    """Invalid project configuration."""


_PATH_KEYS = {
    "bam", "control_bam", "spikein_bam", "bigwig", "peak", "peaks",
    "peak_file", "genome", "genome_fasta", "fasta", "annotation", "gtf",
    "gff", "motif_database", "motif_db", "background_bed", "background_fasta",
    "methylation", "methylation_table", "bedmethyl", "coverage",
}


def _is_remote(value: str) -> bool:
    return bool(re.match(r"^(?:https?|s3|gs)://", value))


def _resolve(value: Any, base: Path) -> Any:
    if not isinstance(value, str) or not value or _is_remote(value):
        return value
    p = Path(os.path.expandvars(os.path.expanduser(value)))
    return str(p if p.is_absolute() else (base / p).resolve())


@dataclass
class SampleRecord:
    """One sample with optional assay-specific inputs."""

    sample_id: str
    assay: str = "unknown"
    bam: str | None = None
    control_bam: str | None = None
    spikein_bam: str | None = None
    bigwig: str | None = None
    peak: str | None = None
    methylation: str | None = None
    group: str | None = None
    time: str | float | None = None
    batch: str | None = None
    subject_id: str | None = None
    genome_build: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any], base: Path) -> "SampleRecord":
        """Build a sample record from a TOML table or manifest row."""
        raw = dict(value)
        sample_id = str(raw.pop("sample_id", raw.pop("name", ""))).strip()
        if not sample_id or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample_id) or sample_id in {".",".."}:
            raise ProjectError("sample_id must start with a letter/digit and contain only ASCII letters, digits, underscore, period or hyphen")
        for key in list(raw):
            if key in _PATH_KEYS and raw[key] is not None:
                raw[key] = _resolve(raw[key], base)
        peak = raw.pop("peak", raw.pop("peaks", raw.pop("peak_file", None)))
        methylation = raw.pop("methylation", raw.pop("methylation_table", raw.pop("bedmethyl", raw.pop("coverage", None))))
        bigwig = raw.pop("bigwig", raw.pop("track", None))
        known = {
            "assay", "bam", "control_bam", "spikein_bam", "group", "time",
            "batch", "subject_id", "genome_build",
        }
        kwargs = {k: raw.pop(k, None) for k in known}
        kwargs.update({"sample_id": sample_id, "peak": peak, "bigwig": bigwig, "methylation": methylation})
        return cls(**kwargs, metadata=raw)


@dataclass
class ProjectConfig:
    """Resolved project configuration."""

    project_id: str = "easyepi-project"
    root: str = "."
    genome_build: str | None = None
    genome: str | None = None
    annotation: str | None = None
    motif_database: str | None = None
    # Metric -> rule overrides for QC warn thresholds; populated from the
    # [qc_thresholds] TOML table or a JSON file referenced by qc_thresholds_file.
    qc_thresholds: dict[str, Any] = field(default_factory=dict)
    samples: list[SampleRecord] = field(default_factory=list)
    defaults: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    source: str | None = None

    def to_dict(self) -> dict[str, Any]:
        """Plain-dict view of the resolved project, used for the run record."""
        return asdict(self)

    def write_resolved(self, out_dir: str | Path) -> Path:
        """Persist the resolved configuration as ``resolved_config.json``."""
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "resolved_config.json"
        path.write_text(json.dumps(self.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8")
        return path


def _load_tsv_manifest(path: Path) -> list[dict[str, Any]]:
    import csv
    with path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader((line for line in fh if line.strip() and not line.startswith("#")), delimiter="\t"))
    if not rows:
        raise ProjectError(f"manifest is empty: {path}")
    return [dict(row) for row in rows]


def load_project(path: str | Path) -> ProjectConfig:
    """Read a TOML project or a tab-separated sample manifest."""
    source = Path(path).expanduser().resolve()
    if source.suffix.lower() in {".tsv", ".txt"}:
        samples = [SampleRecord.from_mapping(x, source.parent) for x in _load_tsv_manifest(source)]
        cfg = ProjectConfig(root=str(source.parent), samples=samples, source=str(source))
        validate_project(cfg, check_files=False)
        return cfg
    if not source.exists():
        raise FileNotFoundError(f"project file not found: {source}")
    with source.open("rb") as fh:
        raw = tomllib.load(fh)
    project = dict(raw.get("project", {}))
    refs = dict(raw.get("references", {}))
    defaults = dict(raw.get("defaults", {}))
    base = source.parent
    samples_raw = raw.get("samples", [])
    if isinstance(samples_raw, Mapping):
        samples_raw = [dict(v, sample_id=k) if isinstance(v, Mapping) else {"sample_id": k} for k, v in samples_raw.items()]
    if not isinstance(samples_raw, list):
        raise ProjectError("[samples] must be an array of tables ([[samples]])")
    samples: list[SampleRecord] = []
    for item in samples_raw:
        if not isinstance(item, Mapping):
            raise ProjectError("each [[samples]] entry must be a table")
        merged = dict(defaults)
        merged.update(item)
        samples.append(SampleRecord.from_mapping(merged, base))
    def ref(*keys: str) -> str | None:
        """Resolve a reference path from the project, falling back to ``[references]``."""
        for key in keys:
            if refs.get(key) is not None:
                return _resolve(refs[key], base)
            if project.get(key) is not None:
                return _resolve(project[key], base)
        return None
    cfg = ProjectConfig(
        project_id=str(project.get("project_id", project.get("name", source.stem))),
        root=str(_resolve(project.get("root", base), base)),
        genome_build=project.get("genome_build", refs.get("genome_build")),
        genome=ref("genome", "genome_fasta", "fasta"),
        annotation=ref("annotation", "gtf", "gff"),
        motif_database=ref("motif_database", "motif_db"),
        qc_thresholds=dict(raw.get("qc_thresholds") or {}),
        samples=samples,
        defaults=defaults,
        metadata={k: v for k, v in raw.items() if k not in {"project", "references", "defaults", "samples"}},
        source=str(source),
    )
    validate_project(cfg, check_files=False)
    return cfg


def validate_project(config: ProjectConfig, *, check_files: bool = True) -> list[str]:
    """Validate sample uniqueness, assays and input paths; return nonfatal warnings."""
    if not config.samples:
        raise ProjectError("project contains no samples")
    ids = [x.sample_id for x in config.samples]
    if len(ids) != len(set(ids)):
        raise ProjectError("sample_id values must be unique")
    warnings: list[str] = []
    for sample in config.samples:
        if not sample.assay or sample.assay.lower() == "unknown":
            warnings.append(f"{sample.sample_id}: assay is not specified")
        for label, value in (("bam", sample.bam), ("control_bam", sample.control_bam), ("spikein_bam", sample.spikein_bam), ("bigwig", sample.bigwig), ("peak", sample.peak), ("methylation", sample.methylation)):
            if value and check_files and _is_remote(value):
                raise ProjectError(f"{label}: download remote files locally before analysis")
            if value and check_files and not Path(value).is_file():
                raise FileNotFoundError(f"{sample.sample_id}: {label} not found: {value}")
    for label, value in (("genome", config.genome), ("annotation", config.annotation), ("motif_database", config.motif_database)):
        if value and check_files and _is_remote(value):
            raise ProjectError(f"{label}: download remote references locally before analysis")
        if value and check_files and not Path(value).is_file():
            raise FileNotFoundError(f"{label} not found: {value}")
    return warnings


def apply_overrides(config: ProjectConfig, overrides: Mapping[str, Any]) -> ProjectConfig:
    """Apply non-None CLI reference overrides to a project."""
    data = config.to_dict()
    for key, value in overrides.items():
        if value is not None and key in {"project_id", "root", "genome_build", "genome", "annotation", "motif_database"}:
            data[key] = value
    return ProjectConfig(
        **{k: data[k] for k in ("project_id", "root", "genome_build", "genome", "annotation", "motif_database", "qc_thresholds", "defaults", "metadata", "source")},
        samples=[SampleRecord(**x) for x in data["samples"]],
    )


def project_inputs(config: ProjectConfig) -> list[dict[str, str]]:
    """Flatten sample inputs for the shared inputs.tsv record."""
    rows: list[dict[str, str]] = []
    for sample in config.samples:
        for key in ("bam", "control_bam", "spikein_bam", "bigwig", "peak", "methylation"):
            value = getattr(sample, key)
            if value:
                rows.append({"sample_id": sample.sample_id, "kind": key, "path": str(value)})
    for kind, value in (("genome", config.genome), ("annotation", config.annotation), ("motif_database", config.motif_database)):
        if value:
            rows.append({"sample_id": "__project__", "kind": kind, "path": str(value)})
    return rows


def write_inputs_tsv(config: ProjectConfig, out_dir: str | Path) -> Path:
    """Flatten every sample input into ``inputs.tsv`` for the run record."""
    import csv
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "inputs.tsv"
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "kind", "path"], delimiter="\t")
        writer.writeheader()
        writer.writerows(project_inputs(config))
    return path

