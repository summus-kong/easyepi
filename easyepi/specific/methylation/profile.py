"""Assay-specific signal or methylation profile summaries."""
from __future__ import annotations

import argparse
from pathlib import Path

from easyepi.methylation import profile_methylation, read_methylation_table
from easyepi._runtime import add_project_args, prepare_run


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi methylation profile")
    add_project_args(parser)
    parser.add_argument("--min-coverage", type=float, default=5)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "METHYLATION-PROFILE")
    with recorder as run:
        out = Path(args.out)
        n = 0
        for sample in config.samples:
            if not sample.methylation:
                continue
            profile_methylation(read_methylation_table(sample.methylation), out / f"{sample.sample_id}.methylation.tsv", args.min_coverage)
            n += 1
        run.set_metrics(samples=n, min_coverage=args.min_coverage)
        run.logger.info("Wrote %d methylation profiles", n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
