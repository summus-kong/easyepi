"""Trace how a peak that overlaps an Intron ends up labelled Intergenic.

The previous diagnostic showed all 1,786 peaks do overlap an intron interval. If the
label comes out Intergenic, either the intron pass never saw them or something wrote
over the label afterwards. This walks one case through the actual code path.
"""
from __future__ import annotations

import sys

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
    fs = A._make_tx_feature_sets(A._read_gtf_txs(GTF))

    labels, _ = A.classify_peaks_to_genomic_bins(pc, ps, pe, fs, mode="overlap",
                                                 promoter_window=3000)
    mine = np.array(labels, dtype=object)
    theirs = cs.bin_mapped.to_numpy()

    idx = np.nonzero((theirs == "Intron") & (mine == "Intergenic"))[0]
    print(f"cases: {len(idx)}")
    print("feature set sizes:", {k: len(v[1]) for k, v in fs.items()})

    # Pick one case and ask the classifier directly with the intron set alone.
    for i in idx[:3]:
        c = pc[i]
        s, e = int(ps[i]), int(pe[i])
        print(f"\n--- peak {i}: {c}:{s}-{e}  our label={mine[i]!r}")

        # count intron overlaps manually for this chromosome
        fc, fst, fen, fstr = fs["Intron"]
        sel = np.nonzero((fc == c) & (fen > s) & (fst < e))[0]
        print(f"    intron intervals overlapping: {len(sel)}")
        if len(sel):
            print(f"    first: {int(fst[sel[0]])}-{int(fen[sel[0]])} strand={fstr[sel[0]]}")

        # what does the intron pass alone say?
        only = {"Intron": fs["Intron"]}
        lab2, _ = A.classify_peaks_to_genomic_bins(
            np.array([c], dtype=object), np.array([s]), np.array([e]),
            only, mode="overlap", promoter_window=None)
        print(f"    intron-only classification: {lab2[0]!r}")

        # does it also overlap an exon?
        ec, est, een, estr = fs["Exon"]
        esel = np.nonzero((ec == c) & (een > s) & (est < e))[0]
        print(f"    exon intervals overlapping: {len(esel)}")
        if len(esel):
            print(f"    first exon: {int(est[esel[0]])}-{int(een[esel[0]])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
