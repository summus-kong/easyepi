"""Check the Python port against ChIPseeker's own output, peak by peak."""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/data1/GaoLab/houoweibo/easyepi_work/easyepiv1")

from easyepi.annotate.chipseeker_port import annotate, load_features  # noqa: E402


def rows(path: Path):
    with Path(path).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm(text: str) -> str:
    text = text.strip()
    for prefix, value in (
        ("Promoter", "Promoter"),
        ("5' UTR", "5'UTR"),
        ("5'UTR", "5'UTR"),
        ("3' UTR", "3'UTR"),
        ("3'UTR", "3'UTR"),
        ("Exon", "Exon"),
        ("Intron", "Intron"),
        ("Intergenic", "Intergenic"),
        ("Distal Intergenic", "Intergenic"),
        ("Downstream", "Downstream"),
    ):
        if text.startswith(prefix):
            return value
    return text


def main() -> int:
    work = Path("/data1/GaoLab/houoweibo/easyepi_work/bench_a10")
    feat_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else work / "feats_transcript"
    label = sys.argv[2] if len(sys.argv) > 2 else "transcript"
    name = sys.argv[3] if len(sys.argv) > 3 else "H3K4me3"

    feats = load_features(feat_dir)
    print(f"loaded {len(feats.strand)} features, "
          f"{sum(len(v) for v in feats.tss.values())} TSS rows")

    peaks = []
    with (work / f"peaks/{name}/{name}/{name}_peaks.narrowPeak").open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                # BED 0-based half-open -> GRanges 1-based inclusive
                peaks.append((f[0], int(f[1]) + 1, int(f[2])))

    ref = rows(work / f"annot/{name}/dist_{label}.tsv")
    ref_by_key = {r["key"]: r for r in ref}

    labels, dists, _near = annotate(
        feats,
        [p[0] for p in peaks],
        [p[1] for p in peaks],
        [p[2] for p in peaks],
    )

    d_ok = d_bad = d_missing = 0
    examples = []
    for (chrom, s, e), mine in zip(peaks, dists):
        r = ref_by_key.get(f"{chrom}:{s}:{e}")
        if r is None:
            d_missing += 1
            continue
        theirs = float(r["distance"])
        if mine is not None and abs(mine - theirs) < 1e-6:
            d_ok += 1
        else:
            d_bad += 1
            if len(examples) < 6:
                examples.append((f"{chrom}:{s}:{e}", mine, theirs,
                                 r["tx_strand"], r["tx_start"], r["tx_end"]))

    print(f"\ndistance: exact={d_ok}/{len(peaks)} diff={d_bad} missing={d_missing}")
    for key, mine, theirs, st, ts, te in examples:
        print(f"   {key}  py={mine}  R={theirs}  strand={st} tx=[{ts},{te}]")

    # Label comparison uses ChIPseeker's full annotatePeak output, not just dist.
    cs_rows = rows(work / f"annot/{name}/cs_{label}.tsv")
    cs_labels = [norm(r["annotation"]) for r in cs_rows]
    mine_labels = [norm(x) for x in labels]
    same = sum(1 for a, b in zip(mine_labels, cs_labels) if a == b)
    print(f"labels  : exact={same}/{len(cs_labels)} ({100*same/len(cs_labels):.2f}%)")
    print(f"   py:", dict(Counter(mine_labels).most_common()))
    print(f"   R :", dict(Counter(cs_labels).most_common()))
    print("   top:", Counter(
        (a, b) for a, b in zip(mine_labels, cs_labels) if a != b
    ).most_common(5))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
