"""Normalize a measured matrix with explicit factors and run provenance."""
import argparse

from easyepi._entry import redirect_out_to_file, run_entry
from easyepi._plots import add_figure_args
from easyepi.normalize import run_normalize


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for this entry point."""
    parser = argparse.ArgumentParser(description="Normalize an observed features-by-samples matrix")
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--factors", required=True, help="TSV: sample_id,value in matrix sample order by ID")
    parser.add_argument("--method", choices=["cpm", "spikein", "rpgc"], default="cpm")
    parser.add_argument("--effective-genome-size", type=int, default=0)
    parser.add_argument("--fragment-length", type=float)
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    return parser


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    args = build_parser().parse_args(argv)
    redirect_out_to_file(args, "normalized_matrix.tsv")
    with run_entry("NORMALIZE", args, method=args.method, output=args.out):
        run_normalize(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
