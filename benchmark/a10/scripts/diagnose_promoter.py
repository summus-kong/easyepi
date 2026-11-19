"""A10 diagnostic: is the EasyEpi "Promoter" call on disputed peaks supported by the gene model?

Motivation
----------
The agreement run shows a strictly one-directional disagreement: EasyEpi labels a
peak Promoter where ChIPseeker labels it Intron/Exon/UTR/Intergenic, and never the
reverse. Two explanations fit that asymmetry:

  (a) EasyEpi invents promoters  -> its nearest-TSS search is wrong, or
  (b) ChIPseeker misses TSSs     -> its nearest-TSS search is wrong.

The two are distinguishable without re-running either annotator: recompute the true
strand-aware nearest TSS for every disputed peak straight from the GTF and compare
to the distance each tool reported. If a transcript TSS really lies within 3 kb,
`tssRegion` should have produced a Promoter label, so ChIPseeker omitted a candidate.

This script recomputes nothing that either tool reported; it reads the GTF itself and
is therefore independent evidence rather than a restatement of EasyEpi's output.
"""

from __future__ import annotations

import bisect
import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10")
ANNOT = ROOT / "heart_annot"
GTF = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_annotation/annot/txdb_full.gtf")
OUT = ROOT / "results" / "heart_promoter_diagnosis.json"
TSS_REGION = 3000

ALIASES = {"5' UTR": "5'UTR", "3' UTR": "3'UTR", "Distal Intergenic": "Intergenic"}


def norm(label: str) -> str:
    base = label.split(" (")[0].strip()
    return ALIASES.get(base, base)


def load_tss() -> dict[str, list[int]]:
    """Strand-aware TSS of every transcript in the GTF, per chromosome (sorted)."""
    tss: dict[str, list[int]] = {}
    with GTF.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            start, end, strand = int(f[3]), int(f[4]), f[6]
            pos = start if strand == "+" else end
            tss.setdefault(f[0], []).append(pos)
    for chrom in tss:
        tss[chrom] = sorted(set(tss[chrom]))
    return tss


def nearest_tss(points: list[int], start: int, end: int) -> int | None:
    """Smallest |TSS - peak| over the interval [start, end]; None if chromosome empty."""
    if not points:
        return None
    lo = bisect.bisect_left(points, start)
    cands = points[max(0, lo - 1): lo + 1]
    hi = bisect.bisect_right(points, end)
    cands += points[max(0, hi - 1): hi + 1]
    return min((abs(p - start) for p in cands), default=None) if cands else None


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    if not GTF.exists():
        print(f"missing {GTF}", file=sys.stderr)
        return 1
    tss = load_tss()
    print(f"GTF transcripts: {sum(len(v) for v in tss.values()):,} distinct TSS")

    verdicts = Counter()
    detail = []
    off_by_one = Counter()

    for name in sorted(p.name for p in ANNOT.iterdir() if p.is_dir()):
        d = ANNOT / name
        if not (d / "cs.tsv").exists():
            continue
        ee_rows, cs_rows = rows(d / "ee" / "annotated.tsv"), rows(d / "cs.tsv")
        for a, b in zip(ee_rows, cs_rows):
            la, lb = norm(a["annotation"]), norm(b["annotation"])
            if la == lb:
                # Record the coordinate convention from peaks where both agree.
                off_by_one[int(float(b["distanceToTSS"])) - int(float(a["tss_distance"]))] += 1
                continue
            if la != "Promoter":
                verdicts["other_direction"] += 1
                continue
            s, e = int(a["start"]), int(a["end"])
            # -1: EasyEpi's start is 1-based on the GRanges row; BED columns are 0-based.
            truth = nearest_tss(tss.get(a["chrom"], []), s - 1, e - 1)
            ee_d = abs(float(a["tss_distance"]))
            cs_d = abs(float(b["distanceToTSS"]))
            if truth is None:
                verdicts["no_tss_on_chrom"] += 1
                continue
            within = truth <= TSS_REGION
            verdicts["gtf_tss_within_3kb" if within else "gtf_tss_beyond_3kb"] += 1
            if len(detail) < 12:
                detail.append({
                    "dataset": name,
                    "peak": f"{a['chrom']}:{s}-{e}",
                    "easyepi_label": la, "easyepi_dist": ee_d,
                    "chipseeker_label": b["annotation"], "chipseeker_dist": cs_d,
                    "gtf_nearest_tss_dist": truth,
                    "easyepi_matches_gtf": ee_d == truth or ee_d == truth + 1,
                    "chipseeker_matches_gtf": cs_d == truth or cs_d == truth + 1,
                })

    total = sum(verdicts.values())
    print(f"\ndisputed Promoter calls: {total:,}")
    for k, v in verdicts.most_common():
        print(f"  {k:24s} {v:7,d}  {100 * v / total:6.2f}%")
    print("\ndistance delta (CS - EE) on peaks where labels AGREE:")
    for k, v in off_by_one.most_common(5):
        print(f"  {k:+4d} bp : {v:7,d}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"verdicts": dict(verdicts), "examples": detail}, indent=2),
                   encoding="utf-8")
    print("\nexamples:")
    for e in detail[:8]:
        print(f"  {e['peak']:22s} EE={e['easyepi_dist']:8.0f} CS={e['chipseeker_dist']:9.0f} "
              f"GTF={e['gtf_nearest_tss_dist']:6d}  {e['chipseeker_label'].split(' (')[0]}")
    print(f"\nwritten -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
