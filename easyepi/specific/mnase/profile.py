"""Assay-specific signal or methylation profile summaries."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from easyepi.nucleosome import nucleosome_scores, periodicity
from easyepi.logging_utils import RunRecorder, setup_logging
from easyepi._plots import add_figure_args, selected_formats


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi MNase profile")
    parser.add_argument("--signal", required=True, help="Single-column tab-separated numeric signal")
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    out = Path(args.out); work_dir = out.parent if out.suffix.lower() == ".json" else out
    if out.suffix.lower() != ".json":
        out = work_dir / "mnase_profile.json"
    logger = setup_logging("MNASE-PROFILE", work_dir, verbose=args.verbose, quiet=args.quiet, log_file=args.log_file)
    with RunRecorder("MNASE-PROFILE", work_dir, logger, figure_formats=selected_formats(args)) as run:
        values = np.loadtxt(args.signal, delimiter="\t", ndmin=1)
        if values.size == 0 or not np.isfinite(values).all():
            raise ValueError("signal must contain finite numeric values")
        result = nucleosome_scores(values)
        result["periodicity"] = periodicity(values)
        # occupancy is an ndarray; JSON has no array type, so export it explicitly
        result["occupancy"] = np.asarray(result["occupancy"], dtype=float).tolist()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        run.set_metrics(points=int(values.size), output=str(out)); logger.info("Wrote %s", out.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
