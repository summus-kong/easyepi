"""Where does the remaining 17.6% disagreement actually live?

Variant B (whole-peak overlap + transcript TSS) reaches 82.42%. This script takes
that variant and dissects the failures, so the next fix targets evidence rather than
a guess.
"""
from __future__ import annotations

import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, "/data1/GaoLab/houoweibo/easyepi_work/easyepiv1")

from easyepi.annotate import (_make_tx_feature_sets, _read_gtf_txs,
                              classify_peaks_to_genomic_bins, read_intervals)

GTF = "/data1/GaoLab/houoweibo/reference/human/gencode.v50.primary_assembly.annotation.gtf.gz"
PEAKS = ("/data1/GaoLab/houoweibo/easyepi_work/encode_h3k4me3/peaks_full/"
         "H3K4me3/H3K4me3/H3K4me3_peaks.narrowPeak")
CS = "/tmp/port/cs.tsv"


def main():
    cs = pd.read_csv(CS, sep="\t")

    # Pull ChIPseeker's own distanceToTSS and the rich annotation string, which the
    # baseline file already carries -- no need to re-run R.
    pc, ps, pe, _ = read_intervals(PEAKS)
    txs = _read_gtf_txs(GTF)
    fs = _make_tx_feature_sets(txs)
    labels, _ = classify_peaks_to_genomic_bins(pc, ps, pe, fs, mode="overlap")
    mine = np.array(labels, dtype=object)
    theirs = cs["bin_mapped"].to_numpy()

    n = len(mine)
    bad = mine != theirs
    print(f"peaks {n}, disagree {bad.sum()} ({100*bad.sum()/n:.2f}%)")
    print()

    print("=== top confusion pairs ===")
    for (t, m), c in Counter(zip(theirs[bad], mine[bad])).most_common(10):
        print(f"  chipseeker {t:11s} -> mine {m:11s} {c:6d} ({100*c/n:5.2f}%)")
    print()

    # Does ChIPseeker's own distance explain the Promoter disagreements?
    dist = cs["distanceToTSS"].to_numpy()
    in_win = (dist >= -3000) & (dist <= 3000)
    print("=== promoter rule check ===")
    print(f"  ChIPseeker calls Promoter on {int((theirs=='Promoter').sum())} peaks")
    print(f"  peaks with |distanceToTSS| <= 3000: {int(in_win.sum())}")
    print(f"  match: {int(((theirs=='Promoter') == in_win).sum())}/{n}")
    disagree_promo = (theirs == "Promoter") & ~in_win
    miss_promo = (theirs != "Promoter") & in_win
    print(f"  called Promoter but |dist|>3000: {int(disagree_promo.sum())}")
    print(f"  |dist|<=3000 but not Promoter  : {int(miss_promo.sum())}")
    print()

    # For peaks ChIPseeker calls Promoter that we do not: how far is our nearest TSS?
    mask = (theirs == "Promoter") & (mine != "Promoter")
    print("=== peaks ChIPseeker calls Promoter, we do not ===")
    print(f"  count: {int(mask.sum())}")
    sub = np.abs(dist[mask])
    if len(sub):
        print(f"  their |distanceToTSS|: median {np.median(sub):.0f}"
              f"  p90 {np.percentile(sub,90):.0f}")
    print(f"  what we called them: {Counter(mine[mask]).most_common(5)}")
    print()

    # Reverse: we call promoter, they do not.
    mask2 = (mine == "Promoter") & (theirs != "Promoter")
    print("=== peaks we call Promoter, ChIPseeker does not ===")
    print(f"  count: {int(mask2.sum())}")
    if mask2.sum():
        print(f"  their |distanceToTSS|: median {np.median(np.abs(dist[mask2])):.0f}")
        print(f"  they called them: {Counter(theirs[mask2]).most_common(5)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
