"""Prepare RNA result identifiers with explicit logging and no implicit model fitting."""
import argparse
from pathlib import Path
from easyepi.integration import read_rna_results
from easyepi._entry import run_entry
from easyepi._plots import add_figure_args


def main(argv=None):
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Normalize Ensembl-style RNA gene ID versions")
    parser.add_argument("--input", required=True)
    parser.add_argument("--gene-column", default="gene_id")
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    output = Path(args.out)
    work_dir = output.parent if output.suffix == ".tsv" else output
    if output.suffix != ".tsv":
        output = output / "rna_prepared.tsv"
    with run_entry('RNA-PREPARE', args) as run:
        frame = read_rna_results(args.input, args.gene_column)
        frame.to_csv(output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
