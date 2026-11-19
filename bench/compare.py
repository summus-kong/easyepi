#!/usr/bin/env python3
"""Compare ChIPseeker and easyepi peak annotations, peak by peak.

Both runs saw the same peak file in the same order and both used a TxDb built from
GENCODE v50, so index i in each output is the same peak. Anything that differs is
therefore a difference in classification, not in inputs.

Reports the five quantities agreed before the run:
  1. per-peak exact agreement over the six shared bins
  2. coarse agreement (Promoter / Genic / Intergenic)
  3. the full confusion matrix
  4. per-bin frequency difference
  5. how the disagreements are distributed by TSS distance, which is how a systematic
     cause is told apart from noise
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

SHARED = ["Promoter", "5'UTR", "3'UTR", "Exon", "Intron", "Intergenic"]


def coarse(label: str) -> str:
    if label == "Promoter":
        return "Promoter"
    if label == "Intergenic":
        return "Intergenic"
    return "Genic"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chipseeker", required=True)
    ap.add_argument("--easyepi", required=True)
    args = ap.parse_args()

    cs = pd.read_csv(args.chipseeker, sep="\t")
    ee = pd.read_csv(args.easyepi, sep="\t")

    if len(cs) != len(ee):
        print(f"row count differs: chipseeker {len(cs)} vs easyepi {len(ee)}")
        return 2

    # Both are in input order; assert the coordinates agree so a silent reordering
    # cannot masquerade as agreement.
    #
    # ChIPseeker returns GenomicRanges coordinates: 1-based, closed interval. easyepi
    # writes BED convention: 0-based, half-open. The same peak therefore reads
    # (start+1, end) on the R side, so normalise before comparing.
    same_coords = ((cs["seqnames"].to_numpy() == ee["chrom"].to_numpy()) &
                   (cs["start"].to_numpy() - 1 == ee["start"].to_numpy()) &
                   (cs["end"].to_numpy() == ee["end"].to_numpy()))
    if not same_coords.all():
        bad = int((~same_coords).sum())
        print(f"coordinate mismatch on {bad} rows; refusing to compare")
        # Show a couple so a systematic offset is visible rather than mysterious.
        for i in (~same_coords).nonzero()[0][:3]:
            print(f"   row {i}: chipseeker {cs['seqnames'][i]}:{cs['start'][i]}-{cs['end'][i]}"
                  f"  vs easyepi {ee['chrom'][i]}:{ee['start'][i]}-{ee['end'][i]}")
        return 2

    cs_bin = cs["bin_mapped"].to_numpy()
    ee_bin = ee["annotation"].to_numpy()
    unknown = sorted(set(cs_bin) - set(SHARED))
    if unknown:
        print("unexpected chipseeker bins after mapping:", unknown)
        return 2

    n = len(cs_bin)
    exact = (cs_bin == ee_bin).sum()
    cs_coarse = [coarse(x) for x in cs_bin]
    ee_coarse = [coarse(x) for x in ee_bin]
    coarse_agree = sum(a == b for a, b in zip(cs_coarse, ee_coarse))

    print(f"peaks compared: {n}")
    print()
    print(f"1. per-peak exact agreement : {exact}/{n} = {100*exact/n:.2f}%")
    print(f"2. coarse agreement        : {coarse_agree}/{n} = "
          f"{100*coarse_agree/n:.2f}%")
    print()

    print("3. confusion matrix (rows = ChIPseeker, cols = easyepi)")
    conf = pd.crosstab(pd.Series(cs_bin, name="chipseeker"),
                       pd.Series(ee_bin, name="easyepi"))
    conf = conf.reindex(index=SHARED, columns=SHARED, fill_value=0)
    print(conf.to_string())
    print()

    print("4. per-bin frequency difference (percentage points)")
    cs_pct = conf.sum(axis=1) / n * 100
    ee_pct = conf.sum(axis=0) / n * 100
    for b in SHARED:
        print(f"   {b:11s} chipseeker {cs_pct[b]:6.2f}%  easyepi {ee_pct[b]:6.2f}%"
              f"  delta {ee_pct[b]-cs_pct[b]:+6.2f}")
    max_delta = max(abs(ee_pct[b] - cs_pct[b]) for b in SHARED)
    print(f"   max |delta| = {max_delta:.2f} pp")
    print()

    # 5. where do the disagreements sit?
    disagree = cs_bin != ee_bin
    print("5. disagreements by TSS distance (ChIPseeker distanceToTSS)")
    dist = cs["distanceToTSS"].abs()
    bins = [0, 500, 1000, 2000, 3000, 10000, float("inf")]
    labels = ["<=500", "500-1k", "1k-2k", "2k-3k", "3k-10k", ">10k"]
    cut = pd.cut(dist, bins=bins, labels=labels)
    total = cut.value_counts().reindex(labels)
    bad = cut[disagree].value_counts().reindex(labels).fillna(0).astype(int)
    for lab in labels:
        t, b = int(total.get(lab, 0)), int(bad.get(lab, 0))
        pct = (100 * b / t) if t else 0.0
        print(f"   {lab:8s} {b:6d}/{t:6d} disagree ({pct:5.1f}%)")

    print()
    print("top confusion pairs (chipseeker -> easyepi):")
    pairs = pd.Series(list(zip(cs_bin[disagree], ee_bin[disagree]))).value_counts()
    for (a, b), c in pairs.head(8).items():
        print(f"   {a:11s} -> {b:11s} {c:6d}  ({100*c/n:5.2f}% of all peaks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
