"""Exploratory replicated two-group comparison."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from easyepi.differential_matrix import build_count_matrix, differential_matrix
from easyepi._runtime import add_project_args, prepare_run


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Exploratory two-group comparison of normalized BigWig signal; not a count model")
    add_project_args(parser)
    parser.add_argument("--method", choices=["mannwhitney", "ttest"], default="mannwhitney")
    parser.add_argument("--max-regions", type=int, default=None)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "ATAC-DIFFERENTIAL")
    with recorder as run:
        samples = [s for s in config.samples if s.bigwig and s.peak]
        if len(samples) < 2 or any(not s.group for s in samples):
            raise ValueError("ATAC differential requires at least two samples with bigwig, peak and group")
        groups = list(dict.fromkeys(s.group for s in samples))
        if len(groups) != 2:
            raise ValueError("ATAC differential requires exactly two groups")
        data = build_count_matrix([s.bigwig for s in samples], [s.peak for s in samples],
                                  [s.sample_id for s in samples], max_regions=args.max_regions, seed=args.seed)
        result = differential_matrix(data["matrix"], np.array([s.group == groups[0] for s in samples]), args.method)
        rows = []
        for region, lfc, p, q, sig in zip(data["regions"], result["log2fc"], result["pvalue"], result["padj"], result["significant"]):
            rows.append({"chrom": region[0], "start": region[1], "end": region[2], "log2fc": lfc, "pvalue": p, "padj": q, "significant": bool(sig)})
        out = Path(args.out)
        pd.DataFrame(rows).to_csv(out / "atac_differential.tsv", sep="\t", index=False)
        comparison={"group_a":groups[0],"group_b":groups[1],"log2fc_direction":"group_a / group_b",
                    "units":"integrated normalized BigWig signal", "method":args.method,
                    "max_regions":args.max_regions,"seed":args.seed}
        (out/"comparison.json").write_text(json.dumps(comparison,indent=2),encoding="utf-8")
        run.set_metrics(samples=len(samples), regions=len(rows), **comparison)
        run.logger.info("Wrote atac_differential.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

