"""Construct consensus intervals across samples."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from easyepi._runtime import add_project_args, prepare_run


def _read(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip() and not line.startswith(("#", "track")):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    rows.append((p[0], int(p[1]), int(p[2])))
    return rows


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi ATAC consensus peaks")
    add_project_args(parser)
    parser.add_argument("--min-samples", type=int, default=1)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "ATAC-CONSENSUS")
    with recorder as run:
        if not 1 <= args.min_samples <= len(config.samples):
            raise ValueError("min-samples must be between one and the number of samples")
        events = defaultdict(list)
        for sample_index, sample in enumerate(config.samples):
            if not sample.peak:
                raise ValueError(f"{sample.sample_id} requires a peak BED")
            for chrom, start, end in _read(sample.peak):
                if start < 0 or end <= start:
                    raise ValueError("Invalid BED coordinates")
                events[chrom].append((start, 1, sample_index))
                events[chrom].append((end, -1, sample_index))
        intervals = []
        for chrom, points in sorted(events.items()):
            counts = defaultdict(int)
            previous = None
            for position, delta, sample_index in sorted(points):
                support = sum(value > 0 for value in counts.values())
                if previous is not None and position > previous and support >= args.min_samples:
                    if intervals and intervals[-1][0] == chrom and intervals[-1][2] == previous:
                        c, start, _, old_support = intervals[-1]
                        intervals[-1] = (c, start, position, min(old_support, support))
                    else:
                        intervals.append((chrom, previous, position, support))
                counts[sample_index] += delta
                previous = position
        out = Path(args.out) / "atac_consensus.bed"
        with out.open("w", encoding="utf-8") as fh:
            for i, (chrom, start, end, n) in enumerate(intervals, 1):
                fh.write(f"{chrom}\t{start}\t{end}\tconsensus_{i}\t{n}\n")
        run.set_metrics(intervals=len(intervals), min_samples=args.min_samples)
        run.logger.info("Wrote %s", out.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
