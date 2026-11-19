"""Batch BigWig -> bedGraph -> narrowPeak over a directory.

Why `bdgpeakcall`
-----------------
A BigWig is already a depth track, so MACS3's `bdgpeakcall` is the matching
subcommand: it thresholds signal directly. `callpeak` expects reads and would
apply fragment modelling that is meaningless here.

Cutoff choice
-------------
The cutoff is an absolute signal value, not a q-value. On these imputed tracks
the chr1 signal quantiles were p50=0.23, p95=3.1, p99=5.86, p99.9=30.95. A
cutoff of 5.0 sits near p99 and gave ~3k peaks on chr1 for H3K27ac, which is the
right order of magnitude for an active mark. It is recorded per run so the
threshold is never implicit.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pyBigWig

MACS3 = "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi/bin/macs3"


def to_bedgraph(bw_path: Path, out: Path, chrom: str, bin_size: int) -> int:
    """Write a fixed-width bedGraph for one chromosome; returns bin count."""
    bw = pyBigWig.open(str(bw_path))
    try:
        length = bw.chroms().get(chrom)
        if not length:
            raise KeyError(f"{chrom} absent from {bw_path.name}")
        values = bw.stats(chrom, 0, length, type="mean", nBins=length // bin_size)
        with out.open("w") as fh:
            for i, v in enumerate(values):
                if v is None or v <= 0:
                    continue
                fh.write(f"{chrom}\t{i * bin_size}\t{(i + 1) * bin_size}\t{v:.6f}\n")
        return len(values)
    finally:
        bw.close()


def call_peaks(bedgraph: Path, out: Path, cutoff: float,
               max_gap: int, min_len: int) -> Path:
    """Run `bdgpeakcall`; returns the produced peak file.

    `-o` names the output file itself. Unlike `callpeak`, `bdgpeakcall` does not
    append a suffix, so the caller's path is used verbatim in both the command
    and the existence check.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        MACS3, "bdgpeakcall",
        "-i", str(bedgraph),
        "-o", str(out),
        "-c", str(cutoff),
        "-g", str(max_gap),
        "-l", str(min_len),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"bdgpeakcall failed ({proc.returncode})\n{proc.stderr[-1500:]}"
        )
    if not out.exists():
        raise RuntimeError(f"bdgpeakcall produced no peaks file at {out}")
    return out


def count_peaks(path: Path) -> int:
    with path.open() as fh:
        return sum(1 for line in fh if line.strip() and not line.startswith(("#", "track")))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--bin-size", type=int, default=100)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--max-gap", type=int, default=200)
    ap.add_argument("--min-len", type=int, default=200)
    args = ap.parse_args()

    indir, outdir = Path(args.indir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bgdir = outdir / "bedgraph"
    bgdir.mkdir(exist_ok=True)

    bws = sorted(indir.glob("*.bw"))
    if not bws:
        print(f"no .bw files in {indir}", file=sys.stderr)
        return 1

    failures = 0
    for bw in bws:
        name = bw.stem
        bg = bgdir / f"{name}.bedGraph"
        try:
            n = to_bedgraph(bw, bg, args.chrom, args.bin_size)
            peaks = call_peaks(bg, outdir / f"{name}.narrowPeak",
                               args.cutoff, args.max_gap, args.min_len)
            print(f"{name:28s} bins={n:7d}  peaks={count_peaks(peaks):6d}")
        except Exception as exc:  # noqa: BLE001 - reported, not swallowed
            failures += 1
            print(f"{name:28s} FAILED: {exc}", file=sys.stderr)

    print(f"\n{len(bws) - failures}/{len(bws)} succeeded"
          f"  (chrom={args.chrom} bin={args.bin_size} cutoff={args.cutoff})")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
