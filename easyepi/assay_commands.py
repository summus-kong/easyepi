"""Independent assay input QC using the shared project and logging contract. These summaries do not replace dedicated analysis backends."""

from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path
from typing import Any

from .logging_utils import RunRecorder, setup_logging
from .project import project_inputs, validate_project, write_inputs_tsv
from .project_commands import _config_from_args
from .qc import bam_qc_summary
from ._plots import selected_formats
from .assay_metrics import compute_assay_metrics
from .qc_payload import build_sample_payload
from .qc_thresholds import annotate_rows_with_thresholds, resolve_rules, summarize_verdicts


def _peak_metrics(path: str) -> dict[str, Any]:
    n = 0
    widths: list[int] = []
    malformed = 0
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                start, end = int(parts[1]), int(parts[2])
                if end < start:
                    raise ValueError
                n += 1
                widths.append(end - start)
            except (IndexError, ValueError):
                malformed += 1
    return {"path": path, "n_intervals": n, "median_width": statistics.median(widths) if widths else None,
            "malformed_rows": malformed}


def _bigwig_metrics(path: str) -> dict[str, Any]:
    import pyBigWig
    with pyBigWig.open(path) as bw:
        chroms = bw.chroms() or {}
    return {"path":path,"n_chromosomes":len(chroms),"genome_bases":int(sum(chroms.values()))}


def _generic_input_metrics(sample: Any) -> dict[str, Any] | None:
    value = sample.methylation or sample.metadata.get("input") or sample.metadata.get("methylation") or sample.metadata.get("coverage")
    if not value:
        return None
    p = Path(str(value))
    return {"path": str(value), "exists": p.exists(), "size_bytes": p.stat().st_size if p.exists() else None}


def run_assay_qc(args: Any, assay: str) -> int:
    """Run input QC for one assay and write JSON, TSV and optional figures."""
    # The alias commands (`atac-qc`, ...) pass a canonical label, but the value is
    # normalised anyway so `make_assay_runner("cut&tag")` and any future caller
    # resolve through one vocabulary rather than through the spelling they used.
    # `assay` is normalised only to the canonical label for logic (metric table and
    # thresholds). The output filename keeps the caller-supplied label because
    # `atac-qc` has always written `atac_seq_qc.tsv`, and `merge-qc-reports`
    # documents globbing `*_seq_qc.json`; silently renaming it would break both.
    from .qc_assay import canonical_assay
    canonical = canonical_assay(assay)
    out = Path(args.out)
    label = assay.upper().replace("-SEQ", "") + "-QC"
    logger = setup_logging(label, out, verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    with RunRecorder(label, out, logger, figure_formats=selected_formats(args)) as run:
        config = _config_from_args(args)
        # Rules must be known before building each payload: the payload carries
        # context (e.g. partial peak coverage) that the verdict needs, and that
        # context is lost if the verdict is recomputed later from assay_qc alone.
        _marker = getattr(args, "marker", None) or getattr(config, "marker", None)
        _rules = resolve_rules(_marker,
                               getattr(args, "qc_thresholds", None),
                               getattr(config, "qc_thresholds", None))
        for sample_index, sample in enumerate(config.samples):
            if sample.assay == "unknown":
                sample.assay = canonical
            # Command-line inputs override the manifest so a single sample can be
            # re-QC'd against a freshly called peak file without editing the
            # manifest. Without this, --peak/--bigwig/--bam were silently ignored.
            if getattr(args, "peak", None):
                sample.peak = str(args.peak[sample_index] if sample_index < len(args.peak) else args.peak[0])
            if getattr(args, "bigwig", None):
                sample.bigwig = str(args.bigwig[sample_index] if sample_index < len(args.bigwig) else args.bigwig[0])
            if getattr(args, "bam", None):
                sample.bam = str(args.bam[sample_index] if sample_index < len(args.bam) else args.bam[0])
        validate_project(config, check_files=not args.no_check_files)
        config.defaults.update({"threads": args.threads, "seed": args.seed})
        config.write_resolved(out)
        write_inputs_tsv(config, out)
        if args.print_config:
            print((out / "resolved_config.json").read_text(encoding="utf-8"))
        if args.dry_run:
            run.set_metrics(assay=canonical, samples=len(config.samples), dry_run=True)
            logger.info("Dry run: input validation completed")
            return 0
        rows: list[dict[str, Any]] = []
        for sample in config.samples:
            logger.info("Processing sample %s", sample.sample_id)
            row: dict[str, Any] = {"sample_id": sample.sample_id, "assay": sample.assay,
                                   "group": sample.group, "time": sample.time,
                                   "batch": sample.batch, "subject_id": sample.subject_id}
            if sample.bam:
                row["bam_qc"] = bam_qc_summary(sample.bam)
            if sample.bigwig:
                row["bigwig_qc"] = _bigwig_metrics(sample.bigwig)
            if sample.peak and Path(sample.peak).exists():
                row["peak_qc"] = _peak_metrics(sample.peak)
            elif sample.peak:
                row["peak_qc"] = {"path": sample.peak, "exists": False}
            generic = _generic_input_metrics(sample)
            if generic:
                row["input_qc"] = generic
            if canonical.startswith("methyl") and sample.methylation and Path(sample.methylation).exists():
                from .methylation import read_methylation_table, summarize_methylation
                row["methylation_qc"] = summarize_methylation(read_methylation_table(sample.methylation))
            # Assay-specific + extended metrics come from the shared builder so this
            # command, the generic `qc` command and make-qc-report cannot drift apart.
            _annotation = getattr(args, "annotation", None) or config.annotation
            _blacklist = getattr(args, "blacklist", None) or getattr(config, "blacklist", None)
            _built = build_sample_payload(
                canonical,
                bam=sample.bam,
                bigwig=sample.bigwig,
                peak=sample.peak,
                methylation=sample.methylation,
                annotation=_annotation,
                blacklist=_blacklist,
                bam_qc=row.get("bam_qc"),
                max_reads=getattr(args, "max_reads", 200000),
                rules=_rules,
            )
            if _built["assay_qc"]:
                row["assay_qc"] = _built["assay_qc"]
            if _built["qc_verdict"] is not None:
                row["qc_verdict"] = _built["qc_verdict"]
            rows.append(row)
        # Verdicts were produced per sample above with the resolved rule table;
        # this call is a safety net for rows that somehow lack one.
        rows = annotate_rows_with_thresholds(rows)
        verdict_summary = summarize_verdicts(rows)
        json_path = out / f"{assay.lower().replace('-', '_')}_qc.json"
        json_path.write_text(json.dumps({"assay": assay,
                                         "verdict_summary": verdict_summary,
                                         "samples": rows},
                                        indent=2, ensure_ascii=False), encoding="utf-8")
        tsv_path = out / f"{assay.lower().replace('-', '_')}_qc.tsv"
        with tsv_path.open("w", encoding="utf-8", newline="") as fh:
            base_fields = ["sample_id", "assay", "group", "time", "batch", "subject_id",
                           "bam_total_reads", "bam_mapping_rate", "peak_n_intervals", "peak_median_width"]
            # Discover assay-specific scalar columns so FRiP / periodicity / beta
            # metrics are visible in the TSV, not only in the JSON payload.
            extra_fields: list[str] = []
            for _row in rows:
                for _k, _v in (_row.get("assay_qc") or {}).items():
                    if isinstance(_v, (int, float, str, bool)) and _k not in extra_fields:
                        extra_fields.append(_k)
            fields = base_fields + sorted(extra_fields)
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                flat = dict(row)
                b = row.get("bam_qc") or {}
                p = row.get("peak_qc") or {}
                flat.update({"bam_total_reads": b.get("total_reads"), "bam_mapping_rate": b.get("mapping_rate"),
                             "peak_n_intervals": p.get("n_intervals"), "peak_median_width": p.get("median_width")})
                for _k, _v in (row.get("assay_qc") or {}).items():
                    if isinstance(_v, (int, float, str, bool)):
                        flat[_k] = _v
                writer.writerow(flat)
        run.set_metrics(assay=canonical, samples=len(rows), outputs=[str(tsv_path), str(json_path)])
        logger.info("Wrote %s and %s", tsv_path.name, json_path.name)
    return 0


def make_assay_runner(assay: str):
    """Bind an assay label to its QC command for argparse dispatch."""
    def runner(args: Any) -> int:
        """argparse entry point bound to one assay."""
        return run_assay_qc(args, assay)
    runner.__name__ = f"run_{assay.lower().replace('-', '_')}_qc"
    return runner

