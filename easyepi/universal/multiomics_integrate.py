"""Integrate normalized epigenomic tracks and RNA result tables by gene."""
from __future__ import annotations

import argparse
from pathlib import Path

from easyepi.integration import run_integration
from easyepi._entry import run_entry
from easyepi._plots import add_figure_args


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi epigenome-RNA integration")
    parser.add_argument("--gene-bed", required=True)
    parser.add_argument("--track", action="append", required=True, metavar="NAME=BIGWIG")
    parser.add_argument("--out", required=True)
    parser.add_argument("--clusters", type=int, default=4)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--peaks", default=None)
    parser.add_argument("--rna-results", default=None)
    parser.add_argument("--rna-gene-column", default="gene_id")
    parser.add_argument("--max-link-distance", type=int, default=100000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    out = Path(args.out)
    with run_entry('RNA-INTEGRATION', args) as run:
        tracks = {}
        for item in args.track:
            if "=" not in item:
                raise ValueError("--track must be NAME=BIGWIG")
            name, path = item.split("=", 1)
            if not name or not path or name in tracks:
                raise ValueError("track names must be unique and non-empty")
            tracks[name] = path
        logger.info("Integrating %d epigenome tracks", len(tracks))
        result = run_integration(args.gene_bed, tracks, out, clusters=args.clusters, seed=args.seed,
                                 peaks=args.peaks, rna_results=args.rna_results,
                                 rna_gene_column=args.rna_gene_column,
                                 max_link_distance=args.max_link_distance)
        run.set_metrics(genes_input=result["genes_input"], genes_clustered=result["genes_clustered"], tracks=len(tracks),
                        rna_matched=result.get("rna_genes_matched", 0), peaks_linked=result.get("peaks_linked", 0))
        logger.info("Wrote gene clusters, optional RNA merge and peak-gene links")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

