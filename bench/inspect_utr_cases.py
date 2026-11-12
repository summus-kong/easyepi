"""The 105 remaining disagreements are UTR-vs-Exon. Find out why.

Reading GENCODE's `UTR` rows instead of deriving them changes nothing: for a sampled
transcript the two interval sets are identical. So the difference is not which UTR
intervals exist, it is how they are split into 5' and 3', or which one wins when both
overlap a peak.

This inspects the actual failing peaks against both implementations' feature sets.
"""
from __future__ import annotations

import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, "/data1/GaoLab/houoweibo/easyepi_work/easyepiv1")

import easyepi.annotate as A

GTF = "/data1/GaoLab/houoweibo/reference/human/gencode.v50.primary_assembly.annotation.gtf.gz"
PEAKS = ("/data1/GaoLab/houoweibo/easyepi_work/encode_h3k4me3/peaks_full/"
         "H3K4me3/H3K4me3/H3K4me3_peaks.narrowPeak")
CS = "/tmp/port/cs.tsv"


def main():
    cs = pd.read_csv(CS, sep="\t")
    pc, ps, pe, _ = A.read_intervals(PEAKS)
    txs = A._read_gtf_txs(GTF)
    fs = A._make_tx_feature_sets(txs)
    labels, _ = A.classify_peaks_to_genomic_bins(pc, ps, pe, fs, mode="overlap",
                                                 promoter_window=3000)
    mine = np.array(labels, dtype=object)
    theirs = cs.bin_mapped.to_numpy()

    bad = np.nonzero(mine != theirs)[0]
    print(f"total disagreements: {len(bad)}")
    print("pairs:", Counter(zip(theirs[bad], mine[bad])).most_common())
    print()

    # Inspect a few: does the peak overlap BOTH a UTR and an Exon interval?
    uc, us, ue, ust = fs["5'UTR"]
    tc, ts, te, tst = fs["3'UTR"]
    ec, es, ee, est = fs["Exon"]

    def overlaps(chrom, s, e, fc, fst_, fen):
        sel = np.nonzero((fc == chrom) & (fen > s) & (fst_ < e))[0]
        return [(int(fst_[k]), int(fen[k])) for k in sel]

    for i in bad[:8]:
        c, s, e = pc[i], int(ps[i]), int(pe[i])
        print(f"--- {c}:{s}-{e}")
        print(f"    chipseeker={theirs[i]!r}  mine={mine[i]!r}")
        print(f"    overlaps 5'UTR: {len(overlaps(c,s,e,uc,us,ue))}"
              f"  3'UTR: {len(overlaps(c,s,e,tc,ts,te))}"
              f"  Exon: {len(overlaps(c,s,e,ec,es,ee))}")
        print(f"    tss_distance(cs)={cs.distanceToTSS[i]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
