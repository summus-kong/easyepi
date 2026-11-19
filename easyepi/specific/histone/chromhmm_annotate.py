"""Summarize ChromHMM state coverage and optional user-reviewed labels."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from easyepi.logging_utils import RunRecorder, setup_logging
from easyepi._plots import add_figure_args, selected_formats


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Summarize ChromHMM state intervals; optional user-supplied biological labels")
    parser.add_argument("--states", required=True, help="ChromHMM state BED/TSV")
    parser.add_argument("--labels", help="TSV columns: state, label; manually reviewed state interpretation")
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    out = Path(args.out); work_dir = out.parent if out.suffix.lower() == ".tsv" else out
    if out.suffix.lower() != ".tsv":
        out = work_dir / "states.tsv"
    logger = setup_logging("CHROMHMM-ANNOTATE", work_dir, verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    with RunRecorder("CHROMHMM-ANNOTATE", work_dir, logger, figure_formats=selected_formats(args)) as run:
        frame = pd.read_csv(args.states, sep="\t", comment="#", header=None)
        if frame.shape[1] < 4:
            raise ValueError("State BED requires chrom, start, end, state")
        frame = frame.iloc[:,:4]
        frame.columns = ["chrom", "start", "end", "state"]
        frame["start"] = pd.to_numeric(frame["start"],errors="raise")
        frame["end"] = pd.to_numeric(frame["end"],errors="raise")
        if (frame.start < 0).any() or (frame.end <= frame.start).any():
            raise ValueError("Invalid state interval")
        frame["length_bp"] = frame.end-frame.start
        if args.labels:
            labels = pd.read_csv(args.labels,sep="\t",dtype={"state":str})
            frame = frame.merge(labels[["state","label"]],on="state",how="left",validate="many_to_one")
        frame.groupby("state")["length_bp"].agg(["count","sum"]).to_csv(work_dir/"state_coverage.tsv",sep="\t")
        out.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out, sep="\t", index=False)
        run.set_metrics(rows=len(frame), output=str(out)); logger.info("Wrote %s", out.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
