"""Independent assay input QC."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from easyepi.methylation import read_methylation_table, summarize_methylation
from easyepi._runtime import add_project_args, prepare_run


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi methylation table QC")
    add_project_args(parser)
    parser.add_argument("--min-coverage", type=float, default=0)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "METHYLATION-QC")
    with recorder as run:
        rows = []
        for sample in config.samples:
            if not sample.methylation:
                run.logger.warning("%s has no methylation table; skipped", sample.sample_id)
                continue
            run.logger.info("Processing sample %s", sample.sample_id)
            metrics = summarize_methylation(read_methylation_table(sample.methylation), args.min_coverage)
            metrics.update({"sample_id": sample.sample_id, "group": sample.group, "time": sample.time, "input": sample.methylation})
            rows.append(metrics)
        out = Path(args.out)
        (out / "methylation_qc.json").write_text(json.dumps({"samples": rows}, indent=2, ensure_ascii=False), encoding="utf-8")
        import pandas as pd
        pd.DataFrame(rows).to_csv(out / "methylation_qc.tsv", sep="\t", index=False)
        run.set_metrics(samples=len(rows), min_coverage=args.min_coverage)
        run.logger.info("Wrote methylation_qc.json and methylation_qc.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

