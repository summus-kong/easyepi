"""Benchmark the chunked/parallel depth path against the previous pileup path.

The comparison is the point: a speed claim without a measured baseline is not a
claim. Both implementations are run on the same chromosome of the same real BAM and
the resulting arrays are checked for equality, so "faster" cannot silently mean
"different".

Run:
    python scripts/bench_coverage.py --bam <bam> --chrom chr1 --length 20000000
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def reference_pileup(bam_path: str, chrom: str, length: int, bin_size: int,
                     min_mapq: int) -> np.ndarray:
    """The previous implementation: per-column pileup into Python lists."""
    import pysam

    n_bins = max(1, int(np.ceil(length / bin_size)))
    out = np.zeros(n_bins, dtype=np.float64)
    positions: list[int] = []
    depths: list[int] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for col in bam.pileup(chrom, 0, length, truncate=True, stepper="all",
                              min_base_quality=0, min_mapping_quality=min_mapq,
                              max_depth=1_000_000):
            if col.nsegments == 0:
                continue
            positions.append(int(col.reference_pos))
            depths.append(int(col.nsegments))
    if not positions:
        return out
    pos = np.asarray(positions, dtype=np.int64)
    dep = np.asarray(depths, dtype=np.float64)
    bins = pos // bin_size
    keep = (bins >= 0) & (bins < n_bins)
    np.add.at(out, bins[keep], dep[keep])
    return out


def main() -> int:
    from easyepi.coverage_parallel import (CoverageSpec, available_workers,
                                           coverage_chromosome)
    import pysam

    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--length", type=int, default=20_000_000)
    ap.add_argument("--bin-size", type=int, default=25)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--chunk-mb", type=float, default=1.0)
    ap.add_argument("--skip-reference", action="store_true",
                    help="skip the slow baseline (use once it is known)")
    args = ap.parse_args()

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        full_length = bam.get_reference_length(args.chrom)
    length = min(args.length, full_length)

    spec = CoverageSpec(bam_path=args.bam, chrom=args.chrom,
                        chrom_length=length, bin_size=args.bin_size,
                        min_mapq=args.min_mapq, shift="none", max_depth=None)
    chunk_bp = int(args.chunk_mb * 1_000_000)

    results = {}
    if not args.skip_reference:
        t0 = time.time()
        ref = reference_pileup(args.bam, args.chrom, length, args.bin_size,
                               args.min_mapq)
        results["baseline_pileup"] = (time.time() - t0, ref)

    for workers in (1, available_workers(args.workers)):
        t0 = time.time()
        got = coverage_chromosome(spec, workers=workers, chunk_bp=chunk_bp)
        elapsed = time.time() - t0
        results[f"chunked_w{workers}"] = (elapsed, got)

    print(f"region      : {args.chrom}:0-{length:,}  bin={args.bin_size} bp")
    print(f"cores visible: {available_workers(args.workers)}")
    print()
    baseline = results.get("baseline_pileup")
    for name, (elapsed, values) in results.items():
        agree = ""
        if baseline is not None and name != "baseline_pileup":
            same = np.allclose(values, baseline[1], rtol=0, atol=1e-9)
            agree = "  identical" if same else "  DIFFERS"
        speed = ""
        if baseline is not None and name != "baseline_pileup":
            speed = f"  speedup x{baseline[0] / elapsed:.1f}"
        print(f"{name:18s} {elapsed:7.2f}s{speed}{agree}")

    if baseline is not None:
        for name, (elapsed, values) in results.items():
            if name == "baseline_pileup":
                continue
            if not np.allclose(values, baseline[1], rtol=0, atol=1e-9):
                worst = float(np.max(np.abs(values - baseline[1])))
                print(f"\nMISMATCH in {name}: max abs difference {worst:g}")
                return 1
    print("\nall computed tracks agree")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
