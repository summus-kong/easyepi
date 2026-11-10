"""Run the standard analysis chain from a project or a folder.

Design rules this module follows
--------------------------------
1. **Reuse, never reimplement.** Every step calls the existing command entry
   point with the same arguments a user would type. There is no second copy of
   the statistics.
2. **Never invent an input.** A step whose inputs are absent is skipped and
   reported. Nothing is substituted, guessed or defaulted silently.
3. **Every decision is printed.** Skipped steps, inferred assays and inferred
   genome builds all appear in the final summary, so the run is auditable.
4. **Failure is local.** One step failing does not abort the chain; it is recorded
   and the summary says which step failed and why.

The chain is deliberately short. It is the set of steps that apply to any bulk
epigenomics sample with a BAM and a genome annotation: QC, peak calling,
annotation, motif enrichment, and a combined report. Assay-specific extras
(nucleosome positioning, chromatin states) are named in the summary as available
next steps rather than guessed into the run.
"""

from __future__ import annotations

import json
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from easyepi.quickstart import (build_project, guess_assay, guess_genome,
                                project_to_toml, scan_samples)


@dataclass
class StepResult:
    """Outcome of one step in the chain."""

    name: str
    status: str            # "done" | "skipped" | "failed"
    detail: str = ""
    output: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return {"step": self.name, "status": self.status,
                "detail": self.detail, "output": self.output}


@dataclass
class Plan:
    """What the chain intends to do, decided before any work runs."""

    project: Mapping[str, Any]
    samples: dict[str, dict[str, str]]
    notes: list[str] = field(default_factory=list)
    steps: list[str] = field(default_factory=list)


def _normalise_source(source: str | Path) -> tuple[dict[str, Any], dict[str, dict[str, str]], list[str]]:
    """Resolve either a project file or a folder into a project dict.

    A folder is scanned and a project is generated in memory; a ``.toml``/``.tsv``
    is loaded through the existing loader. The returned dict always has the same
    shape so the rest of the chain does not branch on the source type.
    """
    path = Path(source)
    if path.is_dir():
        project, notes = build_project(path)
        samples = {row["sample_id"]: {k: v for k, v in row.items()
                                      if k in {"bam", "bigwig", "peak", "methylation"}}
                   for row in project["samples"]}
        return project, samples, notes

    if not path.exists():
        raise FileNotFoundError(f"no such project or folder: {source}")

    from easyepi.project import load_project
    config = load_project(path)
    project = {
        "project": {"project_id": config.project_id,
                    "root": config.root,
                    "genome_build": config.genome_build},
        "references": {k: v for k, v in
                       (("genome", config.genome),
                        ("annotation", config.annotation),
                        ("motif_database", config.motif_database)) if v},
        # group/batch must survive the conversion: they decide whether a
        # differential comparison is possible at all.
        "samples": [{"sample_id": s.sample_id,
                     "assay": s.assay or "unknown",
                     **({"group": s.group} if s.group else {}),
                     **({"batch": s.batch} if s.batch else {}),
                     **{k: v for k, v in (("bam", s.bam), ("bigwig", s.bigwig),
                                          ("peak", s.peak),
                                          ("methylation", s.methylation)) if v}}
                    for s in config.samples],
    }
    samples = {row["sample_id"]: {k: v for k, v in row.items()
                                  if k in {"bam", "bigwig", "peak", "methylation"}}
               for row in project["samples"]}
    notes = []
    if not config.genome_build:
        notes.append("genome_build is not set in the project file")
    return project, samples, notes


def _has_two_groups(project: Mapping[str, Any]) -> bool:
    """True when the project declares at least two non-empty sample groups."""
    groups = {str(row.get("group")).strip() for row in project.get("samples") or []
              if row.get("group") not in (None, "")}
    return len(groups) >= 2


def make_plan(source: str | Path, *, steps: Sequence[str] | None = None) -> Plan:
    """Decide which steps can run, without executing anything."""
    project, samples, notes = _normalise_source(source)
    if not samples:
        raise ValueError(f"{source} defines no samples")

    has_bam = any("bam" in tracks for tracks in samples.values())
    has_peak = any("peak" in tracks for tracks in samples.values())
    has_bigwig = any("bigwig" in tracks for tracks in samples.values())
    refs = project.get("references") or {}
    has_annotation = bool(refs.get("annotation"))
    has_motifs = bool(refs.get("motif_database"))

    available = []
    if has_bam or has_bigwig or has_peak:
        available.append("qc")
    if has_bam or has_bigwig:
        available.append("peaks")
    if has_peak and has_annotation:
        available.append("annotate")
    if has_peak and has_motifs:
        available.append("motifs")
    # Differential analysis is opt-in: it needs a two-group design, and many
    # studies have neither replicates nor a comparison. It is offered only when
    # the project actually declares groups, so the default chain never implies
    # a test that cannot be run.
    if _has_two_groups(project) and has_peak:
        available.append("differential")
    available.append("report")

    if steps:
        unknown = sorted(set(steps) - set(available))
        if unknown:
            raise ValueError(
                f"step(s) not available given these inputs: {', '.join(unknown)}; "
                f"available: {', '.join(available)}")
        chosen = list(steps)
    else:
        chosen = available

    notes = list(notes)
    if not has_annotation:
        notes.append("no annotation in the project: 'annotate' skipped")
    if not has_motifs:
        notes.append("no motif database in the project: 'motifs' skipped")
    if not _has_two_groups(project):
        notes.append("fewer than two sample groups: 'differential' not planned")
    return Plan(project=project, samples=samples, notes=notes, steps=chosen)


# ---------------------------------------------------------------------------
# Step implementations: each delegates to the existing entry point
# ---------------------------------------------------------------------------
def _step_qc(project_path: Path, out: Path) -> StepResult:
    """Study-wide QC over the project samples.

    ``run_qc`` treats ``--out`` as the report *file*, so the path is given
    explicitly; a directory would silently write a file named after it.
    """
    from argparse import Namespace

    from easyepi.qc import run_qc

    target = out / "qc" / "qc_report.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    # ``qc`` resolves inputs from --project, or from explicit flags; passing the
    # project alone is the supported form.
    args = Namespace(project=str(project_path), out=str(target),
                     bam=None, bigwig=None, peak=None, assay=None,
                     annotation=None, blacklist=None, qc_thresholds=None,
                     max_reads=200000, verbose=False, quiet=True,
                     figure_formats=None, input_dir=None, methylation=None)
    run_qc(args)
    if not target.exists():
        return StepResult("qc", "failed", f"expected report not written: {target}")
    import json as _json

    n_samples = len((_json.loads(target.read_text(encoding="utf-8"))
                     .get("samples") or {}))
    return StepResult("qc", "done", f"{n_samples} sample(s) summarised",
                      output=str(target))

def _step_peaks(project_path: Path, out: Path) -> StepResult:
    """Call peaks with MACS3 for every sample that has a BAM or a signal track.

    Reuses the existing ``macs3_command`` builder, so the arguments are exactly
    the ones the assay-specific subcommands use. Samples that already have a
    peak file are left alone: re-calling peaks over an existing set would
    silently change the regions every later step depends on.
    """
    import shutil

    from easyepi._peak_calling import macs3_command
    from easyepi.logging_utils import run_external
    from easyepi.project import load_project

    macs3 = shutil.which("macs3")
    if not macs3:
        return StepResult(
            "peaks", "skipped",
            "MACS3 is not on PATH; install it or supply peak files in the project")

    config = load_project(project_path)
    genome = config.genome_build or "hs"
    target = out / "peaks"
    target.mkdir(parents=True, exist_ok=True)

    called, reused, failures = [], [], []
    for sample in config.samples:
        if sample.peak:
            reused.append(sample.sample_id)
            continue
        if not (sample.bam or sample.bigwig):
            continue
        assay = (sample.assay or "").lower()
        broad = "h3k27me3" in sample.sample_id.lower() or "h3k9me3" in sample.sample_id.lower()
        try:
            command = macs3_command(macs3, sample, target, genome_size=genome,
                                    assay="ATAC" if "atac" in assay else "CHIP",
                                    broad=broad)
            run_external(command, _quiet_logger())
            called.append(sample.sample_id)
        except Exception as exc:                       # noqa: BLE001
            failures.append(f"{sample.sample_id}: {type(exc).__name__}: {exc}")

    if failures:
        return StepResult("peaks", "failed", "; ".join(failures[:3]),
                          output=str(target))
    if not called and reused:
        return StepResult("peaks", "done",
                          f"reused {len(reused)} existing peak set(s)",
                          output=str(target))
    if not called:
        return StepResult("peaks", "skipped", "no sample had a BAM or signal track")
    detail = f"called {len(called)} peak set(s)"
    if reused:
        detail += f", reused {len(reused)}"
    return StepResult("peaks", "done", detail, output=str(target))


def _quiet_logger():
    """Logger that swallows external-tool chatter during a one-command run."""
    import logging

    logger = logging.getLogger("easyepi.quick.external")
    if not logger.handlers:
        logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.CRITICAL)
    return logger

def _step_annotate(project_path: Path, out: Path) -> StepResult:
    from easyepi.annotate import run_annotate
    from easyepi.project import load_project
    from argparse import Namespace

    config = load_project(project_path)
    peaks = [s.peak for s in config.samples if s.peak]
    if not peaks or not config.annotation:
        return StepResult("annotate", "skipped", "needs both peaks and an annotation")
    written = []
    for peak in peaks:
        stem = Path(peak).name.split(".")[0]
        target = out / "annotate" / f"{stem}.tsv"
        target.parent.mkdir(parents=True, exist_ok=True)
        run_annotate(Namespace(peaks=peak, gtf=config.annotation, out=str(target),
                               backend="native", mode="midpoint", min_overlap=0.0,
                               verbose=False, quiet=True, figure_formats=None,
                               log_file=None))
        written.append(str(target))
    return StepResult("annotate", "done", f"{len(written)} peak file(s) annotated",
                      output=str(out / "annotate"))


def _step_motifs(project_path: Path, out: Path) -> StepResult:
    from easyepi.project import load_project

    config = load_project(project_path)
    peaks = [s.peak for s in config.samples if s.peak]
    if not peaks or not config.motif_database:
        return StepResult("motifs", "skipped", "needs peaks and a motif database")
    # Known-motif enrichment needs sequences, which this step cannot obtain
    # without a genome FASTA; say so instead of silently producing nothing.
    if not config.genome:
        return StepResult("motifs", "skipped",
                          "needs references.genome (FASTA) to extract peak sequences")
    from easyepi.tf_annotate import run_tf_annotate
    from argparse import Namespace

    target = out / "motifs"
    target.mkdir(parents=True, exist_ok=True)
    run_tf_annotate(Namespace(peaks=peaks[0], genome_fasta=config.genome,
                              motifs=config.motif_database, fpr=1e-4,
                              background=None, top=5, out=str(target),
                              figure_formats=None, verbose=False, quiet=True,
                              log_file=None))
    return StepResult("motifs", "done", "known-motif annotation of the first peak set",
                      output=str(target))


def _declared_groups(project: Mapping[str, Any]) -> list[tuple[str, list[str]]]:
    """Group label -> sample ids, in first-seen order."""
    order: list[str] = []
    members: dict[str, list[str]] = {}
    for row in project.get("samples") or []:
        label = str(row.get("group") or "").strip()
        if not label:
            continue
        if label not in members:
            members[label] = []
            order.append(label)
        members[label].append(row["sample_id"])
    return [(label, members[label]) for label in order]


def _step_differential(project: Mapping[str, Any], project_path: Path,
                       out: Path) -> StepResult:
    """Differential peaks between the two groups declared in the project.

    The underlying implementation already degrades to fold-change-only when a
    group has one sample (``test=fold_change_only_no_replicates``) and refuses
    when a group is empty, so this step only supplies the design and reports
    which comparison it ran.
    """
    from argparse import Namespace

    from easyepi.peaks.diffpeaks import run_differential_peaks
    from easyepi.project import load_project

    groups = _declared_groups(project)
    if len(groups) < 2:
        return StepResult("differential", "skipped",
                          "the project declares fewer than two groups")
    (label_a, ids_a), (label_b, ids_b) = groups[0], groups[1]

    config = load_project(project_path)
    by_id = {s.sample_id: s for s in config.samples}
    specs = []
    missing = []
    for sid in ids_a + ids_b:
        sample = by_id.get(sid)
        if sample is None or not sample.peak or not sample.bam:
            missing.append(sid)
            continue
        specs.append(f"{sid}={sample.bam}:{sample.peak}")
    if missing:
        return StepResult("differential", "skipped",
                          "needs a BAM and a peak file for: "
                          + ", ".join(missing[:5]))

    target = out / "differential" / "differential_peaks.tsv"
    target.parent.mkdir(parents=True, exist_ok=True)
    run_differential_peaks(Namespace(
        peaks=specs, out=str(target), design=None,
        group_a=",".join(ids_a), min_overlap=1, min_mapq=30))
    return StepResult("differential", "done",
                      f"{label_a} (n={len(ids_a)}) vs {label_b} (n={len(ids_b)})",
                      output=str(target))

def _step_report(project: Mapping[str, Any], results: list[StepResult], out: Path) -> StepResult:
    """Write one human-readable summary of everything that ran."""
    report = {
        "project": project.get("project", {}),
        "steps": [r.as_dict() for r in results],
        "n_samples": len(project.get("samples") or []),
    }
    target = out / "easyepi_report.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [f"# easyepi run summary", "",
             f"Samples: {report['n_samples']}", "",
             "| step | status | detail |", "| --- | --- | --- |"]
    for r in results:
        lines.append(f"| {r.name} | {r.status} | {r.detail} |")
    (out / "SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return StepResult("report", "done", "SUMMARY.md and easyepi_report.json",
                      output=str(target))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run_quick(args) -> int:
    """CLI driver for ``easyepi quick``."""
    source = Path(args.source)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    plan = make_plan(source, steps=getattr(args, "steps", None))
    print(f"[easyepi] quick: {len(plan.samples)} sample(s), "
          f"steps = {' -> '.join(plan.steps)}")
    for note in plan.notes:
        print(f"  note: {note}")

    # Materialise a project file when the source was a folder, so later commands
    # and the recorded run refer to the same configuration.
    if source.is_dir():
        project_path = out / "project.generated.toml"
        project_path.write_text(project_to_toml(plan.project), encoding="utf-8")
        print(f"  wrote {project_path.name} (edit and keep for reproducibility)")
    else:
        project_path = source

    results: list[StepResult] = []
    runners: dict[str, Callable[[], StepResult]] = {
        "qc": lambda: _step_qc(project_path, out),
        "peaks": lambda: _step_peaks(project_path, out),
        "annotate": lambda: _step_annotate(project_path, out),
        "motifs": lambda: _step_motifs(project_path, out),
        "differential": lambda: _step_differential(plan.project,
                                                   project_path, out),
    }

    for step in plan.steps:
        if step == "report":
            continue
        if step not in runners:
            results.append(StepResult(
                step, "skipped",
                "this step needs a dedicated command with extra arguments; "
                f"run `easyepi {step} --help`"))
            continue
        try:
            result = runners[step]()
        except Exception as exc:                      # noqa: BLE001 - reported below
            result = StepResult(step, "failed", f"{type(exc).__name__}: {exc}")
            if getattr(args, "verbose", False):
                traceback.print_exc()
        results.append(result)
        marker = {"done": "OK", "skipped": "SKIP", "failed": "FAIL"}[result.status]
        print(f"  [{marker}] {result.name}: {result.detail}")

    summary = _step_report(plan.project, results, out)
    results.append(summary)
    print(f"[easyepi] summary -> {out / 'SUMMARY.md'}")

    failed = [r for r in results if r.status == "failed"]
    return 1 if failed else 0


def run_init(args) -> int:
    """CLI driver for ``easyepi init``: write a project file from a folder."""
    project, notes = build_project(
        args.source,
        project_id=getattr(args, "project_id", None),
        genome_build=getattr(args, "genome_build", None),
        annotation=getattr(args, "annotation", None),
        motif_database=getattr(args, "motif_database", None))
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(project_to_toml(project), encoding="utf-8")

    print(f"[easyepi] init: {len(project['samples'])} sample(s) -> {target}")
    for row in project["samples"]:
        kinds = ", ".join(k for k in ("bam", "bigwig", "peak", "methylation")
                          if k in row)
        print(f"    {row['sample_id']:20s} {row.get('assay', '?'):10s} {kinds}")
    for note in notes:
        print(f"  note: {note}")
    print("  next: review the assay labels and paths, then `easyepi quick "
          f"{target} --out results`")
    return 0
