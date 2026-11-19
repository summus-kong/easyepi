"""Compare externally supported TF activity estimates between two replicated groups."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu
from easyepi.core.stats import benjamini_hochberg
from easyepi._entry import run_entry
from easyepi._plots import add_figure_args


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi TF activity differential")
    parser.add_argument("--activity", required=True)
    parser.add_argument("--design", required=True, help="sample_id/group TSV")
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    out = Path(args.out); work_dir = out.parent if out.suffix.lower() == ".tsv" else out
    if out.suffix.lower() != ".tsv":
        out = work_dir / "tf_activity_differential.tsv"
    with run_entry('TF-ACTIVITY-DIFFERENTIAL', args) as run:
        act, design = pd.read_csv(args.activity, sep="\t"), pd.read_csv(args.design, sep="\t")
        if not {"sample_id","tf","activity"}.issubset(act.columns):
            raise ValueError("Activity input requires sample_id, tf, activity")
        if act.duplicated(["sample_id","tf"]).any():
            raise ValueError("Duplicate sample/TF estimates")
        if not {"sample_id","group"}.issubset(design.columns) or design.sample_id.duplicated().any():
            raise ValueError("Design requires unique sample_id and group")
        if set(act.sample_id) != set(design.sample_id) or design.group.isna().any():
            raise ValueError("Activity/design sample mismatch or missing group")
        frame = act.merge(design[["sample_id", "group"]], on="sample_id", validate="many_to_one")
        levels = list(dict.fromkeys(frame["group"]))
        if len(levels) != 2: raise ValueError("TF activity differential requires exactly two groups")
        rows = []
        for tf, sub in frame.groupby("tf"):
            a, b = sub.loc[sub.group == levels[0], "activity"].dropna(), sub.loc[sub.group == levels[1], "activity"].dropna()
            if min(len(a),len(b)) < 2:
                raise ValueError("At least two biological replicates per group and TF are required")
            p = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
            rows.append({"tf": tf, "group1": levels[0], "group2": levels[1], "delta": float(b.mean() - a.mean()) if len(a) and len(b) else None, "pvalue": p})
        if not rows:
            raise ValueError("No TF activity estimates")
        result = pd.DataFrame(rows); result["padj"] = benjamini_hochberg(result["pvalue"].to_numpy(float))
        out.parent.mkdir(parents=True, exist_ok=True); result.to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
