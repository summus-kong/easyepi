"""Scan known JASPAR/PFM motifs and return best sequence scores."""
from __future__ import annotations
import argparse, json
from pathlib import Path
from easyepi.motifs import read_pwm, scan_sequences
from easyepi._entry import run_entry
from easyepi._plots import add_figure_args


def _fasta(path):
    seqs, current = [], []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current)); current = []
            else:
                current.append(line.strip())
        if current: seqs.append("".join(current))
    return seqs


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi native PWM motif scan")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--pwm", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    out = Path(args.out)
    work_dir = out.parent if out.suffix.lower() == ".json" else out
    if out.suffix.lower() != ".json":
        out = work_dir / "motif_scan.json"
    with run_entry('MOTIF-SCAN', args) as run:
        result = scan_sequences(read_pwm(Path(args.pwm)), _fasta(args.fasta))
        serial = {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in result.items()}
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(serial, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
