"""Known-motif enrichment using Biopython PWM thresholds, Fisher exact tests and BH correction."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact
from easyepi.core.stats import benjamini_hochberg
from easyepi.motifs import read_pwm, scan_sequences
from easyepi.universal.motif_scan import _fasta
from easyepi._entry import run_entry
from easyepi._plots import add_figure_args


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi motif enrichment")
    parser.add_argument("--target", required=True)
    parser.add_argument("--background", required=True)
    parser.add_argument("--pwm", action="append", required=True, help="PFM; repeatable")
    parser.add_argument("--fpr", type=float, default=1e-4, help="PWM score false positive rate under a uniform background")
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    out = Path(args.out)
    work_dir = out.parent if out.suffix.lower() == ".json" else out
    if out.suffix.lower() != ".json":
        out = work_dir / "motif_enrichment.json"
    with run_entry('MOTIF-ENRICHMENT', args) as run:
        target, background = _fasta(args.target), _fasta(args.background)
        if not target or not background or Path(args.target).resolve() == Path(args.background).resolve():
            raise ValueError("Provide nonempty, distinct target and matched background FASTA files")
        if not 0 < args.fpr < 1:
            raise ValueError("fpr must lie between zero and one")
        from Bio.motifs.matrix import PositionSpecificScoringMatrix
        rows = []
        for pwm_path in args.pwm:
            pwm = read_pwm(Path(pwm_path))
            model = PositionSpecificScoringMatrix("ACGT", {base: pwm[i].tolist() for i,base in enumerate("ACGT")})
            threshold = model.distribution(background={base:0.25 for base in "ACGT"},precision=1000).threshold_fpr(args.fpr)
            ts = scan_sequences(pwm, target)["scores"] >= threshold
            bs = scan_sequences(pwm, background)["scores"] >= threshold
            table = [[int(ts.sum()), int((~ts).sum())], [int(bs.sum()), int((~bs).sum())]]
            odds, p = fisher_exact(table, alternative="greater")
            rows.append({"motif": Path(pwm_path).stem, "target_hits": int(ts.sum()), "target_n": len(ts), "background_hits": int(bs.sum()), "background_n": len(bs), "score_threshold": float(threshold), "fpr_per_strand":args.fpr, "odds_ratio": float(odds), "pvalue": float(p)})
        q = benjamini_hochberg(np.array([r["pvalue"] for r in rows]))
        for row, value in zip(rows, q): row["padj"] = float(value)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps({"motifs": rows}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
