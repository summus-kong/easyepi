"""A10 on independent data: per-peak annotation agreement, EasyEpi vs ChIPseeker.

Input layout (produced by the A10 heart pipeline):
  heart_annot/<dataset>/ee/annotated.tsv   chrom start end width annotation tss_distance
  heart_annot/<dataset>/cs.tsv             seqnames start end annotation distanceToTSS

Both annotators read the same BED and resolve features from the same gene model,
so a label difference can only come from the annotation rules. Rows are in peak
order on both sides; the only coordinate offset is the BED 0-based -> GRanges
1-based start, which is checked rather than assumed (see check_alignment).

Label vocabulary: ChIPseeker splits UTRs ("5' UTR") and intergenic
("Distal Intergenic"), and can emit "Downstream"; EasyEpi writes "5'UTR",
"3'UTR", "Intergenic" and has no Downstream bin. Aliases are folded so the
comparison is on biology, not spelling. "Downstream" has no EasyEpi counterpart
and is reported separately instead of being folded into Intergenic.

Reported per dataset, no thresholded "pass":
  * exact label match after alias folding
  * distance agreement (exact bp, and within 1/10/100 bp)
  * per-label recall and the top confusion pairs
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10")
ANNOT = ROOT / "heart_annot"
OUT = ROOT / "results" / "heart_agreement.json"

# ChIPseeker spelling -> EasyEpi spelling. "Downstream" is deliberately absent.
ALIASES = {
    "5' UTR": "5'UTR",
    "3' UTR": "3'UTR",
    "Distal Intergenic": "Intergenic",
}


def norm(label: str) -> str:
    """Fold a ChIPseeker annotation to its bin name, keeping sub-bin detail out."""
    base = label.split(" (")[0].strip()
    return ALIASES.get(base, base)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_alignment(ee: list[dict], cs: list[dict]) -> dict:
    """Verify row order matches and the only coordinate offset is start + 1."""
    problems = []
    start_offsets = Counter()
    for i, (a, b) in enumerate(zip(ee, cs)):
        if a["chrom"] != b["seqnames"]:
            problems.append(f"row {i}: chrom {a['chrom']} != {b['seqnames']}")
        if int(a["end"]) != int(b["end"]):
            problems.append(f"row {i}: end {a['end']} != {b['end']}")
        start_offsets[int(b["start"]) - int(a["start"])] += 1
        if len(problems) >= 5:
            break
    return {"row_count_match": len(ee) == len(cs),
            "start_offset_histogram": dict(start_offsets),
            "problems": problems}


def compare(name: str) -> dict:
    d = ANNOT / name
    ee_rows = read_rows(d / "ee" / "annotated.tsv")
    cs_rows = read_rows(d / "cs.tsv")
    align = check_alignment(ee_rows, cs_rows)
    n = min(len(ee_rows), len(cs_rows))
    ee_rows, cs_rows = ee_rows[:n], cs_rows[:n]

    ee_lab = [norm(r["annotation"]) for r in ee_rows]
    cs_lab = [norm(r["annotation"]) for r in cs_rows]

    agree = sum(1 for a, b in zip(ee_lab, cs_lab) if a == b)
    conf = Counter((a, b) for a, b in zip(ee_lab, cs_lab) if a != b)

    per_label = {}
    for want in sorted(set(ee_lab) | set(cs_lab)):
        tot = sum(1 for x in ee_lab if x == want)
        hit = sum(1 for x, y in zip(ee_lab, cs_lab) if x == want and y == want)
        per_label[want] = {
            "n_easyepi": tot,
            "n_chipseeker": sum(1 for y in cs_lab if y == want),
            "matched": hit,
            "recall": round(hit / tot, 4) if tot else None,
        }

    d_ee = [float(r["tss_distance"]) for r in ee_rows]
    d_cs = [float(r["distanceToTSS"]) for r in cs_rows]
    dd = [abs(a - b) for a, b in zip(d_ee, d_cs)]
    dist = {
        "exact": sum(1 for x in dd if x == 0),
        "within_1bp": sum(1 for x in dd if x <= 1),
        "within_10bp": sum(1 for x in dd if x <= 10),
        "within_100bp": sum(1 for x in dd if x <= 100),
        "max_abs_diff": max(dd) if dd else None,
        "median_abs_diff": sorted(dd)[len(dd) // 2] if dd else None,
    }

    return {
        "n_peaks": n,
        "alignment": align,
        "exact_label_match": agree,
        "exact_label_frac": round(agree / n, 6) if n else None,
        "per_label": per_label,
        "distance": dist,
        "chipseeker_only_downstream": sum(1 for r in cs_rows
                                          if r["annotation"].startswith("Downstream")),
        "top_disagreements": [
            {"easyepi": a, "chipseeker": b, "n": c} for (a, b), c in conf.most_common(8)
        ],
    }


def main() -> int:
    if not ANNOT.is_dir():
        print(f"missing {ANNOT}", file=sys.stderr)
        return 1
    names = sorted(p.name for p in ANNOT.iterdir() if p.is_dir())
    results, skipped = {}, {}

    for name in names:
        d = ANNOT / name
        missing = [str(f) for f in (d / "ee" / "annotated.tsv", d / "cs.tsv")
                   if not f.exists()]
        if missing:
            skipped[name] = {"missing": missing}
            print(f"{name}: SKIPPED (missing {', '.join(missing)})")
            continue
        r = compare(name)
        results[name] = r
        n = r["n_peaks"]
        print(f"{name}: label {r['exact_label_match']:,}/{n:,} "
              f"({100 * r['exact_label_frac']:.4f}%)  "
              f"dist exact {r['distance']['exact']:,} "
              f"<=10bp {r['distance']['within_10bp']:,}  "
              f"maxdiff {r['distance']['max_abs_diff']:g}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    payload = {"datasets": results, "skipped": skipped,
               "aliases": ALIASES,
               "note": "ChIPseeker 1.46.1, level=transcript, tssRegion +/-3kb, "
                       "TxDb.Hsapiens.UCSC.hg38.knownGene; EasyEpi from the same TxDb-derived GTF"}
    OUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"\n{len(results)}/{len(names)} datasets compared -> {OUT}")
    return 0 if results else 1


if __name__ == "__main__":
    raise SystemExit(main())
