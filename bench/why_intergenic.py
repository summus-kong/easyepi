"""Why are 1,786 ChIPseeker Intron peaks called Intergenic here?

Three candidate explanations, tested rather than assumed:
  A. the peak overlaps nothing at all in our feature set (a coverage hole)
  B. the peak overlaps an exon of a *different* transcript, so our priority order
     assigns Exon and then... no, we would say Exon. If we say Intergenic, we matched
     no exon and no intron.
  C. our intron intervals are computed from exon gaps within one transcript, so a
     peak inside an intron of transcript A but inside an *exon* of transcript B is
     Exon for us, not Intergenic -- again not this case.

Since we report Intergenic, the peak matched nothing: neither exon nor intron nor
promoter. So the question is which feature ChIPseeker has that we do not.
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


def coverage_by_chrom(chrom, start, end):
    """Sorted interval list per chromosome for a feature set."""
    out = {}
    for c, s, e in zip(chrom, start, end):
        out.setdefault(c, []).append((int(s), int(e)))
    for c in out:
        out[c].sort()
    return out


def hits(index, chrom, s, e):
    import bisect
    items = index.get(chrom)
    if not items:
        return 0
    hi = bisect.bisect_left(items, (e,))
    n = 0
    for k in range(hi - 1, -1, -1):
        a, b = items[k]
        if b <= s:
            continue
        n += 1
    return n


def main():
    cs = pd.read_csv(CS, sep="\t")
    pc, ps, pe, _ = read_intervals(PEAKS)
    txs = _read_gtf_txs(GTF)
    fs = _make_tx_feature_sets(txs)
    labels, _ = classify_peaks_to_genomic_bins(pc, ps, pe, fs, mode="overlap",
                                               promoter_window=3000)
    mine = np.array(labels, dtype=object)
    theirs = cs.bin_mapped.to_numpy()

    mask = (theirs == "Intron") & (mine == "Intergenic")
    idx = np.nonzero(mask)[0]
    print(f"peaks ChIPseeker calls Intron, we call Intergenic: {len(idx)}")

    idx_exon = coverage_by_chrom(*fs["Exon"][:3])
    idx_intron = coverage_by_chrom(*fs["Intron"][:3])
    idx_utr5 = coverage_by_chrom(*fs["5'UTR"][:3])
    idx_utr3 = coverage_by_chrom(*fs["3'UTR"][:3])

    tally = Counter()
    for i in idx:
        c, s, e = pc[i], int(ps[i]), int(pe[i])
        tally["any_exon"] += hits(idx_exon, c, s, e) > 0
        tally["any_intron"] += hits(idx_intron, c, s, e) > 0
        tally["any_utr5"] += hits(idx_utr5, c, s, e) > 0
        tally["any_utr3"] += hits(idx_utr3, c, s, e) > 0
    print("  overlapping at least one ...", dict(tally))

    # What transcripts do they fall inside? Use the gene-level span of our records.
    spans = {}
    for t in txs:
        gs = t["gene_start"] if t["gene_start"] is not None else t["start"]
        ge = t["gene_end"] if t["gene_end"] is not None else t["end"]
        spans.setdefault(t["chrom"], []).append((gs, ge))
    for c in spans:
        spans[c].sort()
    inside_gene = 0
    for i in idx:
        if hits(spans, pc[i], int(ps[i]), int(pe[i])) > 0:
            inside_gene += 1
    print(f"  inside some transcript's genomic span: {inside_gene} "
          f"({100*inside_gene/len(idx):.1f}%)")

    # Are these peaks perhaps overlapping exons that our parser missed because the
    # transcript has a single exon and we therefore created no intron?
    single = sum(1 for t in txs if len(t["exons"]) == 1)
    print(f"\nsingle-exon transcripts: {single} of {len(txs)}")
    print(f"intron intervals we build: {len(fs['Intron'][1])}")
    print(f"exon intervals we build  : {len(fs['Exon'][1])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
