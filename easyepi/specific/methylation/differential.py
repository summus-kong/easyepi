"""Exploratory replicated two-group comparison."""
from __future__ import annotations

import argparse
from pathlib import Path

from easyepi.methylation import differential_methylation, read_methylation_table
from easyepi._runtime import add_project_args, prepare_run


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi methylation differential analysis")
    add_project_args(parser)
    parser.add_argument("--min-coverage", type=float, default=5)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "METHYLATION-DIFFERENTIAL")
    with recorder as run:
        samples = [s for s in config.samples if s.methylation]
        if len(samples) < 2 or any(not s.group for s in samples):
            raise ValueError("methylation differential requires at least two samples with group")
        frames = [read_methylation_table(s.methylation) for s in samples]
        result = differential_methylation(frames, [s.sample_id for s in samples], [s.group for s in samples], args.min_coverage)
        out = Path(args.out)
        result.to_csv(out / "methylation_differential.tsv", sep="\t", index=False)
        run.set_metrics(samples=len(samples), rows=len(result), output=str(out / "methylation_differential.tsv"))
        run.logger.info("Wrote methylation_differential.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

