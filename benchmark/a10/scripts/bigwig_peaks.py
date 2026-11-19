"""BigWig -> bedGraph -> peak calls, using MACS3's bedGraph input path.

Why bedGraph
------------
MACS3 cannot read BigWig. Its `--format BEDGRAPH` accepts a 4-column
chrom/start/end/value file, which is exactly what a BigWig stores internally, so
the conversion is lossless and avoids inventing a signal model.

Why a fixed bin
----------------
BigWig intervals are already discrete, but their widths vary wildly (1 bp to
megabases). MACS3's local-lambda window is defined in bp, so feeding raw
variable-width intervals would make the effective number of points per window
depend on how the file happened to be written. Rebinning to a fixed width makes
the poisson test see a uniform grid, which is what its model assumes.

Axis convention
---------------
bedGraph and narrowPeak are both 0-based half-open, and pyBigWig also reports
0-based half-open, so no coordinate conversion happens here. This is deliberate:
an off-by-one introduced at the I/O boundary would silently shift every peak.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pyBigWig

# Fixed bin for the poisson grid. 100 bp is the scale at which H3K27ac and
# H3K4me3 peaks are resolved without exploding the point count: chr1 at 100 bp
# is ~2.5 M points, which MACS3 handles in seconds.
BIN_BP = 100


def read_bigwig(path: Path, chroms: List[str], bin_bp: int = BIN_BP) -> List[Tuple[str, int, int, float]]:
    """Fixed-width bedGraph rows, skipping bins with no coverage."""
    bw = pyBigWig.open(str(path))
    if bw is None:
        raise RuntimeError(f"cannot open {path}")
    sizes = bw.chroms()
    rows: List[Tuple[str, int, int, float]] = []
    try:
        for chrom in chroms:
            if chrom not in sizes:
                continue
            length = sizes[chrom]
            # `stats` is the vectorised accessor; one call per chromosome beats
            # per-bin `values` by orders of magnitude.
            vals = bw.stats(chrom, 0, length, type="mean", nBins=(length // bin_bp) + 1)
            for i, v in enumerate(vals):
                if v is None:
                    continue
                start = i * bin_bp
                rows.append((chrom, start, start + bin_bp, float(v)))
    finally:
        bw.close()
    return rows


def write_bedgraph(rows: List[Tuple[str, int, int, float]], out: Path) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        for chrom, start, end, value in rows:
            if value > 0:
                fh.write(f"{chrom}\t{start}\t{end}\t{value:.6f}\n")
    return sum(1 for r in rows if r[3] > 0)


def call_peaks(bedgraph: Path, out: Path, cutoff: float,
               max_gap: int, min_len: int) -> Path:
    """Run MACS3's `bdgpeakcall` on a bedGraph; raises MACS3's own message.

    `bdgpeakcall` is the macs3 subcommand for a depth track: it thresholds the
    signal rather than stacking reads, so `-c` takes a cutoff in signal units and
    `-g` the maximum gap between merged regions. A q-value has no meaning here.
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    exe = Path(sys.executable).with_name("macs3")
    if not exe.exists():
        raise RuntimeError(
            f"macs3 console script not found at {exe}; "
            "install with `micromamba install -c bioconda macs3`"
        )
    cmd = [
        str(exe), "bdgpeakcall",
        "-i", str(bedgraph),
        "-c", str(cutoff),
        "-g", str(max_gap),
        "-l", str(min_len),
        "-o", str(out),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not out.exists():
        raise RuntimeError(
            f"macs3 bdgpeakcall failed (exit {proc.returncode})\n"
            f"cmd: {' '.join(cmd)}\n{proc.stderr[-1500:]}"
        )
    return out


def count_peaks(path: Path) -> int:
    with path.open(encoding="utf-8") as fh:
        return sum(1 for line in fh if line.strip() and not line.startswith(("#", "track")))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigwig", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cutoff", type=float, default=1.0,
                    help="bdgpeakcall signal cutoff, in track units")
    ap.add_argument("--max-gap", type=int, default=200,
                    help="merge regions separated by at most this many bp")
    ap.add_argument("--min-len", type=int, default=200)
    ap.add_argument("--chroms", default="chr1")
    args = ap.parse_args()

    bw_path = Path(args.bigwig)
    out_dir = Path(args.out)
    chroms = [c for c in args.chroms.split(",") if c]
    name = bw_path.stem

    rows = read_bigwig(bw_path, chroms)
    bg = out_dir / f"{name}.bedGraph"
    kept = write_bedgraph(rows, bg)
    print(f"{name}: {len(rows)} bins read, {kept} with signal -> {bg.name}")

    if kept == 0:
        print("  no signal in the requested region; skipping peak calling")
        return 0

    peaks = call_peaks(bg, out_dir / f"{name}.narrowPeak",
                       args.cutoff, args.max_gap, args.min_len)
    print(f"  peaks: {count_peaks(peaks):,} -> {peaks.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
