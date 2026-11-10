"""Independent BAM input QC."""
from __future__ import annotations
import argparse, json
from pathlib import Path
from easyepi.qc import bam_qc_summary
from easyepi.logging_utils import RunRecorder, setup_logging
from easyepi._plots import add_figure_args, selected_formats


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi BAM QC")
    parser.add_argument("--bam", action="append", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    out = Path(args.out)
    work_dir = out.parent if out.suffix.lower() == ".json" else out
    logger = setup_logging("BAM-QC", work_dir, verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    with RunRecorder("BAM-QC", work_dir, logger, figure_formats=selected_formats(args)) as run:
        result = {Path(p).stem: bam_qc_summary(p) for p in args.bam}
        target = out if out.suffix.lower() == ".json" else out / "bam_qc.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps({"samples": result}, indent=2), encoding="utf-8")
        import pandas as pd
        rows = [{"sample_id": name, **metrics} for name, metrics in result.items()]
        pd.DataFrame(rows).to_csv(target.parent / "bam_qc.tsv", sep="\t", index=False)
        run.set_metrics(samples=len(rows), output=str(target))
        logger.info("Wrote %s", target.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
