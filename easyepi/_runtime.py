"""Shared project parsing and run initialization for independent scripts."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from easyepi.logging_utils import RunRecorder, setup_logging
from easyepi.project import ProjectConfig, load_project, validate_project, write_inputs_tsv
from easyepi._plots import add_figure_args, selected_formats


def add_project_args(parser: argparse.ArgumentParser, *, require_project: bool = True) -> None:
    """Add the shared project/output/threads arguments to a parser."""
    parser.add_argument("--project", required=require_project, help="Project TOML or TSV sample manifest")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--no-check-files", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)


def prepare_run(args: Any, module: str) -> tuple[ProjectConfig, RunRecorder]:
    """Resolve the project, write provenance and return an open RunRecorder."""
    out = Path(args.out)
    logger = setup_logging(module, out, verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    config = load_project(args.project)
    validate_project(config, check_files=not args.no_check_files)
    if args.threads < 1:
        raise ValueError("threads must be positive")
    config.defaults.update({"threads": args.threads, "seed": args.seed})
    config.write_resolved(out)
    write_inputs_tsv(config, out)
    if args.dry_run:
        with RunRecorder(module, out, logger, render_figures=False) as run:
            run.set_metrics(dry_run=True, samples=len(config.samples))
        raise SystemExit(0)
    return config, RunRecorder(module, out, logger, figure_formats=selected_formats(args))

