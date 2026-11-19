"""Estimate TF activity from measured normalized features using decoupler ULM."""
import argparse
from pathlib import Path
import pandas as pd
from easyepi.activity import estimate_tf_activity
from easyepi._entry import run_entry
from easyepi._plots import add_figure_args


def main(argv=None):
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="decoupler ULM activity; requires measured normalized accessibility or expression")
    parser.add_argument("--matrix",required=True,help="TSV: feature ID then sample columns, normalized beforehand")
    parser.add_argument("--network",required=True,help="TSV: source (TF), target (feature ID), weight")
    parser.add_argument("--min-targets",type=int,default=5)
    parser.add_argument("--out",required=True)
    parser.add_argument("--verbose",action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("--quiet",action="store_true")
    add_figure_args(parser)
    args=parser.parse_args(argv)
    out=Path(args.out)
    work_dir=out.parent if out.suffix == ".tsv" else out
    if out.suffix != ".tsv":
        out=out/"tf_activity.tsv"
    with run_entry('TF-ACTIVITY', args) as run:
        matrix=pd.read_csv(args.matrix,sep="\t",index_col=0)
        network=pd.read_csv(args.network,sep="\t")
        result=estimate_tf_activity(matrix,network,args.min_targets)
        result.to_csv(out,sep="\t",index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
