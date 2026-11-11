#!/usr/bin/env python3
"""Compare easyepi's built-in peak caller against MACS3 on the same BAM.

Why this is the right comparison to make
---------------------------------------
The MACS3 wrapper was already shown to be argument-transparent (byte-identical peak
files). That says nothing about ``call_peaks_from_bam``, which is a separate
implementation: a Poisson test on binned depth. Reviewers ask whether the pure-Python
path finds the same regions as the tool everyone uses, and until now there was no
number for that.

What is measured, and why these statistics
------------------------------------------
Peak callers are not expected to be identical. MACS3 models local background with a
dynamic lambda and merges enriched bins; the built-in caller uses a global Poisson
null. So exact coordinate equality is the wrong test, and reporting it would
misrepresent both tools.

The comparison therefore reports:

1. **Reciprocal overlap** -- the fraction of each tool's peaks that touch a peak of
   the other. This is the number a biologist means by "do they agree".
2. **Base-pair Jaccard** -- width-robust set similarity, so a caller that emits many
   small peaks inside one large MACS3 peak is not penalised for fragmenting.
3. **Coverage agreement** -- correlation of per-bin signal inside the union of both
   peak sets, which checks that both are picking out the same enriched regions rather
   than the same coordinates by accident.

Gates are set at levels that would catch a broken implementation (a sign error, an
off-by-one in binning, a wrong null) rather than at levels that demand the two models
coincide.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

WORK = Path("/data1/GaoLab/houoweibo/easyepi_work")
OUT = WORK / "validation_peak_concordance"
# macs3 lives in the analysis environment, not the R validation environment
MACS = "/data1/GaoLab/houoweibo/soft/micromamba/root/envs/easyepi/bin/macs3"


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue
        try:
            rows.append((f[0], int(f[1]), int(f[2])))
        except ValueError:
            continue
    return sorted(rows)


def merge_and_index(rows: list[tuple[str, int, int]]):
    """Merge overlapping intervals per chromosome and return sorted arrays."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e in rows:
        if e > s:
            by_chrom.setdefault(chrom, []).append((s, e))
    merged = {}
    for chrom, spans in by_chrom.items():
        spans.sort()
        out = [list(spans[0])]
        for s, e in spans[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[chrom] = np.asarray(out, dtype=np.int64)
    return merged


def reciprocal_overlap(a: list[tuple[str, int, int]],
                       b: list[tuple[str, int, int]], min_bp: int = 1) -> dict:
    """Fraction of each set that overlaps the other."""
    idx_b = merge_and_index(b)
    idx_a = merge_and_index(a)

    def covered(rows, index):
        hit = 0
        for chrom, s, e in rows:
            spans = index.get(chrom)
            if spans is None or not len(spans):
                continue
            hi = int(np.searchsorted(spans[:, 0], e, side="left"))
            if hi == 0:
                continue
            shared = np.minimum(spans[:hi, 1], e) - np.maximum(spans[:hi, 0], s)
            if np.any(shared >= min_bp):
                hit += 1
        return hit / len(rows) if rows else float("nan")

    return {"a_covered_by_b": round(covered(a, idx_b), 4),
            "b_covered_by_a": round(covered(b, idx_a), 4)}


def bp_jaccard(a: list[tuple[str, int, int]],
               b: list[tuple[str, int, int]]) -> float:
    ma, mb = merge_and_index(a), merge_and_index(b)
    shared = 0
    for chrom, sa in ma.items():
        sb = mb.get(chrom)
        if sb is None:
            continue
        i = j = 0
        while i < len(sa) and j < len(sb):
            lo = max(sa[i, 0], sb[j, 0])
            hi = min(sa[i, 1], sb[j, 1])
            if hi > lo:
                shared += int(hi - lo)
            if sa[i, 1] <= sb[j, 1]:
                i += 1
            else:
                j += 1
    bp_a = int(sum(int(e - s) for arr in ma.values() for s, e in arr))
    bp_b = int(sum(int(e - s) for arr in mb.values() for s, e in arr))
    denom = bp_a + bp_b - shared
    return float(shared / denom) if denom else float("nan")


def bins_covered(rows, chrom: str, bin_size: int, n_bins: int) -> np.ndarray:
    """Binary vector: is bin ``i`` inside any peak of ``rows``?

    Projecting both peak sets onto the same grid is what makes them comparable. One
    value per peak cannot be: the two callers emit different numbers of peaks at
    different widths (measured median 1000 bp vs 191 bp here), so an index-wise
    comparison of per-peak values measures the pairing, not the agreement.
    """
    vec = np.zeros(n_bins, dtype=np.int8)
    for c, s, e in rows:
        if c != chrom:
            continue
        lo = max(0, s // bin_size)
        hi = min(n_bins, max(s // bin_size + 1, int(np.ceil(e / bin_size))))
        if hi > lo:
            vec[lo:hi] = 1
    return vec


def depth_on_grid(bam: Path, chrom: str, chrom_length: int, bin_size: int) -> np.ndarray:
    """Depth per bin over the whole chromosome, computed once."""
    from easyepi.coverage_parallel import CoverageSpec, coverage_chromosome

    spec = CoverageSpec(bam_path=str(bam), chrom=chrom, chrom_length=chrom_length,
                        bin_size=bin_size, min_mapq=30, shift="none", max_depth=None)
    return coverage_chromosome(spec, workers=1)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", default=str(WORK / "data_atac/ENCFF327MYY.bam"))
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--bin-size", type=int, default=1000)
    ap.add_argument("--threshold", type=float, default=5.0)
    ap.add_argument("--max-qvalue", type=float, default=0.01)
    ap.add_argument("--min-overlap", type=float, default=0.5,
                    help="gate: reciprocal overlap fraction both ways")
    args = ap.parse_args()

    from easyepi.peaks import call_peaks_from_bam, peaks_to_bed

    OUT.mkdir(parents=True, exist_ok=True)
    bam = Path(args.bam)
    if not bam.exists():
        raise SystemExit(f"input BAM not found: {bam}")

    # ---- built-in caller --------------------------------------------------
    mine_bed = OUT / "easyepi.narrowPeak"
    peak_set = call_peaks_from_bam(str(bam), args.chrom, bin_size=args.bin_size,
                                   threshold=args.threshold,
                                   max_qvalue=args.max_qvalue, t5_shift=True)
    peaks_to_bed(peak_set, mine_bed, format="narrowPeak")
    mine = read_bed(mine_bed)
    print(f"easyepi : {len(mine)} peaks on {args.chrom}")

    # ---- MACS3 ------------------------------------------------------------
    macs_dir = OUT / "macs3"
    macs_dir.mkdir(exist_ok=True)
    cmd = [MACS, "callpeak", "-t", str(bam), "-f", "BAMPE", "-g", "hs",
           "-n", "macs", "--outdir", str(macs_dir), "-q", "0.01"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"MACS3 failed:\n{proc.stderr[-2000:]}")
    produced = sorted(macs_dir.glob("*_peaks.narrowPeak"))
    if not produced:
        raise SystemExit("MACS3 produced no narrowPeak file")
    theirs = [r for r in read_bed(produced[0]) if r[0] == args.chrom]
    print(f"MACS3   : {len(theirs)} peaks on {args.chrom}")

    # ---- overlap statistics ----------------------------------------------
    if not mine or not theirs:
        report = {"error": "one caller produced no peaks on this chromosome",
                  "easyepi_peaks": len(mine), "macs3_peaks": len(theirs)}
        (OUT / "concordance.json").write_text(json.dumps(report, indent=2))
        print("\nFAIL: a caller produced no peaks; nothing to compare")
        return 1

    overlap = reciprocal_overlap(mine, theirs)
    jaccard = bp_jaccard(mine, theirs)

    # Common grid over the chromosome, so both sets describe the same coordinates.
    import pysam

    with pysam.AlignmentFile(str(bam), "rb") as handle:
        chrom_length = int(handle.get_reference_length(args.chrom))
    n_bins = max(1, chrom_length // args.bin_size)
    mask_mine = bins_covered(mine, args.chrom, args.bin_size, n_bins)
    mask_theirs = bins_covered(theirs, args.chrom, args.bin_size, n_bins)

    both = int(((mask_mine == 1) & (mask_theirs == 1)).sum())
    either = int(((mask_mine == 1) | (mask_theirs == 1)).sum())
    bin_jaccard = both / either if either else float("nan")

    depth = depth_on_grid(bam, args.chrom, chrom_length, args.bin_size)
    depth = depth[:n_bins]
    # Correlate depth in bins selected by each caller against the shared selection:
    # this asks whether both prefer the same bins, on a common axis.
    sel_mine = np.flatnonzero(mask_mine)
    sel_theirs = np.flatnonzero(mask_theirs)
    from scipy.stats import spearmanr

    shared_sel = np.intersect1d(sel_mine, sel_theirs)
    rho = None
    if shared_sel.size > 5:
        # within the shared bins, do the two masks still rank them consistently?
        # both are 1 there, so instead correlate depth with the *agreement* pattern
        # across the union: 2 = both call it, 1 = one calls it.
        union_sel = np.union1d(sel_mine, sel_theirs)
        agree = (mask_mine[union_sel] + mask_theirs[union_sel]).astype(float)
        rho = float(spearmanr(agree, depth[union_sel]).statistic)

    gates = {
        "easyepi_covered_by_macs3": overlap["a_covered_by_b"] >= args.min_overlap,
        "macs3_covered_by_easyepi": overlap["b_covered_by_a"] >= args.min_overlap,
        "bp_jaccard_above_chance": jaccard >= 0.10,
        # The shared-bin Jaccard is the interpretable agreement measure. A caller
        # using a global Poisson null will not reproduce MACS3's boundaries, so the
        # gate asks that the majority of occupied bins are common.
        "bin_jaccard_above_0_5": (not np.isnan(bin_jaccard)) and bin_jaccard >= 0.5,
    }

    report = {
        "inputs": {"bam": str(bam), "chrom": args.chrom,
                   "bin_size": args.bin_size, "threshold": args.threshold,
                   "max_qvalue": args.max_qvalue},
        "easyepi_peaks": len(mine), "macs3_peaks": len(theirs),
        "reciprocal_overlap": overlap, "bp_jaccard": round(jaccard, 4),
        "bin_jaccard": round(bin_jaccard, 4),
        "bins_shared": both, "bins_union": either,
        "signal_spearman_within_union": None if rho is None or np.isnan(rho) else round(rho, 4),
        "median_width_bp": {
            "easyepi": float(np.median([e - s for _, s, e in mine])),
            "macs3": float(np.median([e - s for _, s, e in theirs])),
        },
        "method_note": "the two callers use different null models (global Poisson vs "
                       "MACS3 dynamic lambda), so coordinate equality is not expected; "
                       "the gates check that they select the same enriched regions",
        "gates": {k: bool(v) for k, v in gates.items()},
    }
    (OUT / "concordance.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"\nreciprocal overlap : easyepi in MACS3 {overlap['a_covered_by_b']:.3f}, "
          f"MACS3 in easyepi {overlap['b_covered_by_a']:.3f}")
    print(f"base-pair Jaccard  : {jaccard:.4f}")
    print(f"shared-bin Jaccard : {bin_jaccard:.4f} "
          f"({both} common bins of {either})")
    print(f"median peak width  : easyepi {np.median([e - s for _, s, e in mine]):.0f} bp, "
          f"MACS3 {np.median([e - s for _, s, e in theirs]):.0f} bp")
    print("\ngates:")
    for name, ok in gates.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    print(f"\nreport: {OUT / 'concordance.json'}")
    return 0 if all(gates.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
