"""Diagnose the residual label mismatches against ChIPseeker.

The distance function now matches R on all 1,254 peaks, so every remaining
difference must come from child-feature matching. For each mismatched peak this
prints the peak span, the R label, and what the exon/intron tables actually
contain there, which distinguishes "table missing data" from "lookup wrong".
"""

from __future__ import annotations

import bisect
import csv
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from chipseeker_port import annotate, load_features  # noqa: E402

WORK = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10")


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm(text: str) -> str:
    text = text.strip()
    for prefix, value in (
        ("Promoter", "Promoter"), ("5' UTR", "5'UTR"), ("5'UTR", "5'UTR"),
        ("3' UTR", "3'UTR"), ("3'UTR", "3'UTR"), ("Exon", "Exon"),
        ("Intron", "Intron"), ("Intergenic", "Intergenic"),
        ("Distal Intergenic", "Intergenic"), ("Downstream", "Downstream"),
    ):
        if text.startswith(prefix):
            return value
    return text


def main() -> int:
    feats = load_features(WORK / "feats_v2")
    peaks = []
    with (WORK / "peaks/H3K4me3/H3K4me3/H3K4me3_peaks.narrowPeak").open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                peaks.append((f[0], int(f[1]) + 1, int(f[2])))

    labels, _, _near = annotate(feats, [p[0] for p in peaks],
                         [p[1] for p in peaks], [p[2] for p in peaks])

    # `cs_transcript.tsv` rows are in peak order, so index alignment is valid.
    cs_rows = rows(WORK / "annot/H3K4me3/cs_transcript.tsv")
    cs_labels = [norm(r["annotation"]) for r in cs_rows]

    bad = [(i, peaks[i], norm(labels[i]), cs_labels[i])
           for i in range(len(peaks)) if norm(labels[i]) != cs_labels[i]]
    print("mismatches:", len(bad))
    print("pairs:", Counter((a, b) for _, _, a, b in bad).most_common(6))
    print()

    for i, (chrom, s, e), mine, theirs in bad[:10]:
        print(f"{chrom}:{s}-{e}  py={mine}  R={theirs}")
        for kind in ("exon", "intron"):
            table = feats.child.get(kind, {}).get(chrom, [])
            if not table:
                print(f"   {kind}: no rows for {chrom}")
                continue
            starts = [r[0] for r in table]
            hi = bisect.bisect_right(starts, e)
            hits = [r for r in table[:hi] if r[1] >= s]
            print(f"   {kind}: {len(table)} rows on {chrom}, {len(hits)} overlap")
            if hits:
                print(f"      first hit: {hits[0]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
