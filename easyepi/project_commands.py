"""Independent project validation and QC-only reporting."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .logging_utils import RunRecorder, setup_logging
from .project import (ProjectConfig, ProjectError, SampleRecord, apply_overrides,
                      load_project, validate_project, write_inputs_tsv)
from ._plots import selected_formats


def _config_from_args(args: Any) -> ProjectConfig:
    if getattr(args, "project", None):
        config = load_project(args.project)
    else:
        samples = []
        bams = getattr(args, "bam", None) or []
        bigwigs = getattr(args, "bigwig", None) or []
        peaks = getattr(args, "peak", None) or []
        if not bams and not bigwigs and not peaks:
            raise ProjectError("provide --project or at least one of --bam/--bigwig/--peak")
        n = max(len(bams), len(bigwigs), len(peaks))
        if any(len(values) not in (0,n) for values in (bams,bigwigs,peaks)):
            raise ProjectError("Repeated BAM/BigWig/peak inputs must have matching lengths; use a manifest for partial inputs")
        for i in range(n):
            samples.append(SampleRecord(
                sample_id=f"sample_{i + 1}",
                assay=getattr(args, "assay", None) or "unknown",
                bam=bams[i] if i < len(bams) else None,
                bigwig=bigwigs[i] if i < len(bigwigs) else None,
                peak=peaks[i] if i < len(peaks) else None,
            ))
        config = ProjectConfig(
            project_id="command-line-project",
            genome_build=getattr(args, "genome_build", None),
            genome=getattr(args, "genome", None),
            annotation=getattr(args, "annotation", None),
            samples=samples,
        )
    overrides = {
        "genome": getattr(args, "genome", None),
        "annotation": getattr(args, "annotation", None),
        "genome_build": getattr(args, "genome_build", None),
    }
    return apply_overrides(config, overrides)


def validate_project_command(args: Any) -> int:
    """Validate a project file and write the resolved configuration."""
    out = Path(args.out)
    logger = setup_logging("PROJECT", out, verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    with RunRecorder("PROJECT", out, logger, figure_formats=selected_formats(args)) as run:
        config = _config_from_args(args)
        warnings = validate_project(config, check_files=not args.no_check_files)
        config.defaults.update({"threads": args.threads, "seed": args.seed})
        resolved = config.write_resolved(out)
        write_inputs_tsv(config, out)
        run.set_metrics(samples=len(config.samples), warnings=len(warnings), checked_files=not args.no_check_files)
        if args.print_config:
            print(resolved.read_text(encoding="utf-8"))
        logger.info("Samples detected: %d", len(config.samples))
        for warning in warnings:
            logger.warning(warning)
        logger.info("Wrote resolved_config.json and inputs.tsv")
    return 0


def make_qc_report_command(args: Any) -> int:
    """Summarize input QC only; no differential, motif or RNA conclusions."""
    from .qc import bam_qc_summary
    from .qc_payload import build_sample_payload
    from .qc_thresholds import load_rules, summarize_verdicts

    out = Path(args.out)
    logger = setup_logging("QC-REPORT", out, verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    with RunRecorder("QC-REPORT", out, logger, figure_formats=selected_formats(args)) as run:
        config = _config_from_args(args)
        validate_project(config, check_files=not args.no_check_files)
        config.defaults.update({"threads": args.threads, "seed": args.seed})
        resolved = config.write_resolved(out)
        write_inputs_tsv(config, out)
        if args.print_config:
            print(resolved.read_text(encoding="utf-8"))
        if args.dry_run:
            run.set_metrics(samples=len(config.samples), dry_run=True)
            logger.info("Dry run: input validation completed")
            return 0
        result: dict[str, Any] = {"module": "qc-report", "samples": {}}
        for sample in config.samples:
            logger.info("Processing sample %s", sample.sample_id)
            metrics: dict[str, Any] = {"sample_id": sample.sample_id, "assay": sample.assay}
            if sample.bam:
                metrics["bam"] = bam_qc_summary(sample.bam)
            if sample.bigwig:
                p = Path(sample.bigwig)
                metrics["bigwig"] = {"path": sample.bigwig, "exists": p.exists(), "size_bytes": p.stat().st_size if p.exists() else None}
            if sample.peak:
                p = Path(sample.peak)
                metrics["peak"] = {"path": sample.peak, "exists": p.exists(), "size_bytes": p.stat().st_size if p.exists() else None}
            # Same builder the assay commands and `qc` use, so this batch-level
            # verdict cannot disagree with the per-assay report for a sample.
            _rules = load_rules(getattr(args, "qc_thresholds", None),
                                getattr(config, "qc_thresholds", None))
            # Command-line inputs override the manifest, matching the assay
            # commands. Without this the batch report scored fewer metrics than
            # `chip-qc` did for the same sample and reported not_evaluated.
            _cli_peaks = getattr(args, "peak", None) or []
            _cli_bigwigs = getattr(args, "bigwig", None) or []
            _cli_bams = getattr(args, "bam", None) or []
            if _cli_peaks:
                sample.peak = str(_cli_peaks[0])
            if _cli_bigwigs:
                sample.bigwig = str(_cli_bigwigs[0])
            if _cli_bams:
                sample.bam = str(_cli_bams[0])
            _built = build_sample_payload(
                sample.assay,
                bam=sample.bam,
                bigwig=sample.bigwig,
                peak=sample.peak,
                methylation=sample.methylation,
                annotation=getattr(args, "annotation", None) or config.annotation,
                blacklist=getattr(args, "blacklist", None) or getattr(config, "blacklist", None),
                bam_qc=metrics.get("bam"),
                max_reads=getattr(args, "max_reads", 200_000),
                rules=_rules,
            )
            if _built["assay_qc"]:
                metrics["assay_qc"] = _built["assay_qc"]
            if _built["qc_verdict"] is not None:
                metrics["qc_verdict"] = _built["qc_verdict"]
            result["samples"][sample.sample_id] = metrics
        # Batch-level roll-up across all samples, so a multi-batch run reads its
        # overall health from one place instead of summing per-sample warnings.
        result["verdict_summary"] = summarize_verdicts(result["samples"].values())
        report = out / "qc_report.json"
        report.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
        run.set_metrics(samples=len(config.samples), report=str(report))
        logger.info("Wrote qc_report.json")
    return 0

