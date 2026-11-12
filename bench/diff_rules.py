"""Differential test: change ONE rule at a time and measure against ChIPseeker.

The point is to attribute the remaining disagreement to specific rules rather than
guess. Each variant reuses the same peaks, the same GTF and the same ChIPseeker
output, so the only thing that moves is the rule under test.

Variants:
  A. midpoint-in-interval, transcript TSS        (what 0.4.0 ships)
  B. whole-peak overlap, transcript TSS          (ChIPseeker's overlap rule)
  C. midpoint-in-interval, gene TSS              (what 0.3.0 shipped)
  D. whole-peak overlap, gene TSS                (closest to ChIPseeker)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/data1/GaoLab/houoweibo/easyepi_work/easyepiv1")

from easyepi.annotate import (_make_tx_feature_sets, _read_gtf_txs,
                              classify_peaks_to_genomic_bins, read_intervals)

GTF = "/data1/GaoLab/houoweibo/reference/human/gencode.v50.primary_assembly.annotation.gtf.gz"
PEAKS = ("/data1/GaoLab/houoweibo/easyepi_work/encode_h3k4me3/peaks_full/"
         "H3K4me3/H3K4me3/H3K4me3_peaks.narrowPeak")
CS = Path("/tmp/port/cs.tsv")

SHARED = ["Promoter", "5'UTR", "3'UTR", "Exon", "Intron", "Intergenic"]


def normalise(label: str) -> str:
    """Map a rich label to one of the six shared bins."""
    if isinstance(label, str) and label.startswith("Promoter"):
        return "Promoter"
    for name, out in (("5' UTR", "5'UTR"), ("3' UTR", "3'UTR"), ("Exon", "Exon"),
                      ("Intron", "Intron"), ("Downstream", "Intergenic"),
                      ("Intergenic", "Intergenic")):
        if isinstance(label, str) and label.startswith(name):
            return out
    return "Intergenic"


def score(name, mine, theirs):
    a = np.array([normalise(x) for x in mine], dtype=object)
    b = np.array([normalise(x) for x in theirs], dtype=object)
    n = len(a)
    exact = (a == b).sum()

    def coarse(x):
        return "Promoter" if x == "Promoter" else ("Intergenic" if x == "Intergenic" else "Genic")

    coarse_agree = sum(coarse(x) == coarse(y) for x, y in zip(a, b))
    conf = pd.crosstab(pd.Series(b, name="cs"), pd.Series(a, name="mine"))
    conf = conf.reindex(index=SHARED, columns=SHARED, fill_value=0)
    cs_pct = conf.sum(axis=1) / n * 100
    mine_pct = conf.sum(axis=0) / n * 100
    max_delta = max(abs(mine_pct[k] - cs_pct[k]) for k in SHARED)
    print(f"{name:34s} exact {100*exact/n:6.2f}%  coarse {100*coarse_agree/n:6.2f}%"
          f"  maxΔbin {max_delta:6.2f}pp")
    return exact / n, coarse_agree / n, max_delta


def main():
    if not CS.exists():
        print("missing", CS, "- run the ChIPseeker side first")
        return 2
    cs = pd.read_csv(CS, sep="\t")
    pc, ps, pe, _ = read_intervals(PEAKS)
    if len(cs) != len(ps):
        print(f"row count differs: cs {len(cs)} vs peaks {len(ps)}")
        return 2

    txs = _read_gtf_txs(GTF)
    print(f"transcript records: {len(txs)}")

    # Collapse to one record per gene, taking the gene-level span, for variants C/D.
    genes = {}
    for t in txs:
        gid = t["gene_id"]
        if gid not in genes:
            genes[gid] = {
                "gene_id": gid, "transcript_id": gid,
                "chrom": t["chrom"],
                "start": t["gene_start"] if t["gene_start"] is not None else t["start"],
                "end": t["gene_end"] if t["gene_end"] is not None else t["end"],
                "strand": t["strand"],
                "exons": list(t["exons"]), "cds": list(t["cds"]),
            }
        else:
            g = genes[gid]
            g["exons"].extend(t["exons"])
            g["cds"].extend(t["cds"])
    gene_recs = []
    for g in genes.values():
        strand = g["strand"]
        exons = g["exons"] or [(g["start"], g["end"])]
        g["tss"] = min(s for s, _ in exons) if strand != "-" else max(e for _, e in exons)
        g["tes"] = max(e for _, e in exons) if strand != "-" else min(s for s, _ in exons)
        gene_recs.append(g)
    print(f"gene records      : {len(gene_recs)}")

    fs_tx = _make_tx_feature_sets(txs)
    fs_gene = _make_tx_feature_sets(gene_recs)
    print()

    results = {}
    for name, recs, fs, mode in (
        ("A. midpoint + transcript TSS", txs, fs_tx, "midpoint"),
        ("B. whole-peak + transcript TSS", txs, fs_tx, "overlap"),
        ("C. midpoint + gene TSS", gene_recs, fs_gene, "midpoint"),
        ("D. whole-peak + gene TSS", gene_recs, fs_gene, "overlap"),
    ):
        t0 = time.time()
        labels, _ = classify_peaks_to_genomic_bins(pc, ps, pe, fs, mode=mode)
        results[name] = score(name, labels, cs.bin_mapped.values)
        print(f"    ({time.time()-t0:.1f}s)")

    print()
    print("ChIPseeker reference bins:", cs.bin_mapped.value_counts().to_dict())
    return 0


if __name__ == "__main__":
    sys.exit(main())
