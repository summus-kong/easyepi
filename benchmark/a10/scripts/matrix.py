"""A10 comparison matrix: label and distance agreement across gene-model variants.

The first A10 pass reported ~91% label agreement and blamed the remainder on
"ChIPseeker picks a different nearest TSS". That diagnosis was incomplete: the two
annotators were reading *different gene models*, so no rule-level conclusion could be
drawn. This driver measures every combination so the effect of each input is visible
on its own:

  EasyEpi GTF   x  ChIPseeker gene model

Variants
--------
  ee_txdb_full   EasyEpi  --gtf txdb_full.gtf          (GTF whose exon rows carry gene
                                                        ids instead of transcript ids)
  ee_fixed       EasyEpi  --gtf txdb_fixed.gtf         (same GTF, gene_id backfilled)
  cs_ucsc        ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene  (the shipped package)
  cs_gtf_full    ChIPseeker, TxDb built from txdb_full.gtf
  cs_gtf_fixed   ChIPseeker, TxDb built from txdb_fixed.gtf

Only pairs that share a gene model are a fair test of the annotation rules. The
mismatched pairs are kept because they quantify how much of the original 91% was
gene-model disagreement rather than rule disagreement.

Label spellings are folded (ChIPseeker "5' UTR"/"Distal Intergenic" vs EasyEpi
"5'UTR"/"Intergenic"); "Downstream" has no EasyEpi counterpart and is counted
separately rather than folded.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10")
OUT = ROOT / "results" / "matrix.json"

ALIASES = {"5' UTR": "5'UTR", "3' UTR": "3'UTR", "Distal Intergenic": "Intergenic"}

EE_VARIANTS = {
    "ee_txdb_full": ROOT / "heart_annot",
    "ee_fixed": ROOT / "heart_annot_fixed",
}
CS_VARIANTS = {
    "cs_ucsc": ROOT / "heart_annot",
    "cs_gtf_full": ROOT / "heart_annot_txdb",
    "cs_gtf_fixed": ROOT / "heart_annot_txdb2",
}


def norm(label: str) -> str:
    base = label.split(" (")[0].strip()
    return ALIASES.get(base, base)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def compare(ee_path: Path, cs_path: Path) -> dict:
    ee, cs = rows(ee_path), rows(cs_path)
    n = min(len(ee), len(cs))
    ee, cs = ee[:n], cs[:n]
    a = [norm(r["annotation"]) for r in ee]
    b = [norm(r["annotation"]) for r in cs]
    agree = sum(1 for x, y in zip(a, b) if x == y)
    dd = [abs(float(ee[i]["tss_distance"]) - float(cs[i]["distanceToTSS"]))
          for i in range(n)]
    conf = Counter((x, y) for x, y in zip(a, b) if x != y)
    return {
        "n": n,
        "label_agree": agree,
        "label_frac": round(agree / n, 6) if n else None,
        "dist_exact": sum(1 for v in dd if v == 0),
        "dist_within_1": sum(1 for v in dd if v <= 1),
        "dist_within_10": sum(1 for v in dd if v <= 10),
        "top_disagreements": [
            {"easyepi": x, "chipseeker": y, "n": c} for (x, y), c in conf.most_common(5)
        ],
    }


def main() -> int:
    datasets = sorted(p.name for p in (ROOT / "heart_annot").iterdir() if p.is_dir())
    matrix: dict[str, dict] = {}

    for ds in datasets:
        for etag, eroot in EE_VARIANTS.items():
            ee_path = eroot / ds / "ee" / "annotated.tsv"
            if not ee_path.exists():
                continue
            for ctag, croot in CS_VARIANTS.items():
                cs_path = croot / ds / "cs.tsv"
                if not cs_path.exists():
                    continue
                matrix.setdefault(f"{etag}|{ctag}", {})[ds] = compare(ee_path, cs_path)

    print(f"{'pair':40s} {'datasets':>8s} {'label agree':>18s} {'dist<=1bp':>18s}")
    for pair, per in sorted(matrix.items()):
        n = sum(v["n"] for v in per.values())
        la = sum(v["label_agree"] for v in per.values())
        d1 = sum(v["dist_within_1"] for v in per.values())
        flag = "  <- matched gene model" if pair in (
            "ee_fixed|cs_gtf_fixed",) else ""
        print(f"{pair:40s} {len(per):8d} {la:8,d}/{n:<8,d} {100*la/n:6.2f}% "
              f"{d1:8,d}/{n:<8,d} {100*d1/n:6.2f}%{flag}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"matrix": matrix, "aliases": ALIASES}, indent=2),
                   encoding="utf-8")
    print(f"\nwritten -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
