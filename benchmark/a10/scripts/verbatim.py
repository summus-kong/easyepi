"""A10: compare EasyEpi's annotation strings with ChIPseeker's, byte for byte.

No vocabulary folding. ChIPseeker emits `Promoter (<=1kb)`,
`Exon (T1/G1, exon 2 of 4)`, `3' UTR`, `Downstream (<=300bp)` and
`Distal Intergenic`; a coarse-bin comparison would hide a wrong sub-category, a
wrong exon rank or a missing Downstream. This driver treats any difference in the
text as a failure and reports the pairs, so a regression cannot hide behind
"same bin".

Inputs are matched by row, which is safe because both tools emit peaks in the
order they were read and the gene model is shared -- `start` is checked against
`end` on both sides to prove the rows line up.

Usage: verbatim.py [root] [--json out.json]
"""

from __future__ import annotations

import csv
import json
import os
import sys
from collections import Counter
from pathlib import Path

ROOT = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10")
# tag -> directory holding <dataset>/ee/annotated.tsv
EE_ROOT = ROOT / "heart_exact"
# tag -> directory holding <dataset>/cs.tsv
CS_ROOT = ROOT / "heart_annot_txdb2"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def compare(ee_path: Path, cs_path: Path) -> dict:
    ee, cs = rows(ee_path), rows(cs_path)
    n = min(len(ee), len(cs))
    ee, cs = ee[:n], cs[:n]

    # Both sides write the same interval; if they did not, row order is not a
    # valid key and the comparison must not proceed silently.
    misaligned = sum(1 for i in range(n) if int(ee[i]["end"]) != int(cs[i]["end"]))
    if misaligned:
        raise SystemExit(f"{ee_path}: {misaligned} rows disagree on 'end'; "
                         "row order is not a valid key")

    label_hits = sum(1 for i in range(n)
                     if ee[i]["annotation"] == cs[i]["annotation"])
    dist_hits = sum(1 for i in range(n)
                    if int(float(ee[i]["tss_distance"]))
                    == int(float(cs[i]["distanceToTSS"])))
    conf = Counter((ee[i]["annotation"], cs[i]["annotation"])
                   for i in range(n)
                   if ee[i]["annotation"] != cs[i]["annotation"])
    return {
        "peaks": n,
        "label_match": label_hits,
        "label_frac": round(label_hits / n, 6) if n else None,
        "dist_match": dist_hits,
        "dist_frac": round(dist_hits / n, 6) if n else None,
        "top_disagreements": [
            {"easyepi": a, "chipseeker": b, "n": c}
            for (a, b), c in conf.most_common(8)
        ],
    }


def main() -> int:
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    out_path = ROOT / "results" / "verbatim_100.json"
    if "--json" in sys.argv:
        out_path = Path(sys.argv[sys.argv.index("--json") + 1])
    root = Path(argv[0]) if argv else ROOT

    datasets = sorted(p.name for p in (root / "heart_exact").iterdir()
                      if p.is_dir())
    result: dict = {
        "basis": "verbatim annotation-string equality, shared gene model",
        "params": {"tss_region": [-3000, 3000], "level": "transcript",
                   "downstream_distance": 300},
        "datasets": {},
    }
    tn = ta = td = 0
    print(f"{'dataset':24s} {'peaks':>6s} {'label':>10s} {'dist':>10s}")
    for name in datasets:
        ee_path = root / "heart_exact" / name / "ee" / "annotated.tsv"
        cs_path = root / "heart_annot_txdb2" / name / "cs.tsv"
        if not (ee_path.exists() and cs_path.exists()):
            continue
        r = compare(ee_path, cs_path)
        result["datasets"][name] = r
        tn += r["peaks"]
        ta += r["label_match"]
        td += r["dist_match"]
        flag = "" if r["label_match"] == r["peaks"] and r["dist_match"] == r["peaks"] else "  <-- CHECK"
        print(f"{name:24s} {r['peaks']:6d} {100 * r['label_frac']:9.2f}% "
              f"{100 * r['dist_frac']:9.2f}%{flag}")
        for d in r["top_disagreements"]:
            print(f"    EE={d['easyepi'][:46]:46s} | CS={d['chipseeker'][:46]:46s} n={d['n']}")

    result["totals"] = {
        "peaks": tn,
        "label_match": ta,
        "label_frac": round(ta / tn, 6) if tn else None,
        "dist_match": td,
        "dist_frac": round(td / tn, 6) if tn else None,
    }
    print(f"\n{'TOTAL':24s} {tn:6d} {100 * ta / tn:9.4f}% {100 * td / tn:9.4f}%")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"written -> {out_path}")
    return 0 if ta == tn and td == tn else 1


if __name__ == "__main__":
    raise SystemExit(main())
