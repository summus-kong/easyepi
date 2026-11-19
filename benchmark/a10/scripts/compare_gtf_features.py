"""Do gtf_features.py and export_features_from_gtf.R produce the same tables?

Why this test exists: `chipseeker_port` reproduces annotatePeak exactly *given*
the feature tables. On a GTF the tables have to be built by `gtf_features.py`,
and if that builder disagrees with the TxDb accessors then the port's 100%
agreement is about the algorithm only, not about what a GTF user actually gets.

Both sides read the same GTF: R builds a TxDb with makeTxDbFromGFF and then
applies ChIPseeker's accessors; Python parses the GTF directly. The comparison
is on intervals, not on ids, because the two sides name transcripts differently
(tx_name vs a positional index) and the ids carry no scientific content.

Usage: compare_gtf_features.py <feature_dir_py> <feature_dir_R>
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def spans_by_kind(feature_tsv: Path) -> dict[str, Counter]:
    """Multiset of (chrom, start, end) per kind.

    A multiset rather than a set: the same interval can legitimately be emitted
    twice (a shared exon under `by = "gene"` merging, or two transcripts with
    identical introns), and collapsing that would hide a real discrepancy.
    """
    out: dict[str, Counter] = {}
    for row in read_tsv(feature_tsv):
        key = (row["chrom"], int(row["start"]), int(row["end"]))
        out.setdefault(row["kind"], Counter())[key] += 1
    return out


def tss_multiset(tss_tsv: Path) -> Counter:
    return Counter(
        (r["chrom"], int(r["pos"])) for r in read_tsv(tss_tsv)
    )


def describe(name: str, py: Counter, r: Counter) -> int:
    """Report the comparison for one table; return the number of mismatches."""
    only_py = py - r
    only_r = r - py
    total = sum(py.values()) or 1
    bad = sum(only_py.values()) + sum(only_r.values())
    print(f"\n{name}: python={sum(py.values())} R={sum(r.values())} "
          f"exact={bad == 0} ({100 * (1 - bad / total):.4f}%)")
    if only_py:
        print(f"   only in python ({sum(only_py.values())}):",
              list(only_py.items())[:5])
    if only_r:
        print(f"   only in R ({sum(only_r.values())}):",
              list(only_r.items())[:5])
    return bad


def main() -> int:
    py_dir, r_dir = Path(sys.argv[1]), Path(sys.argv[2])

    print("=== TSS ===")
    bad = describe("tss", tss_multiset(py_dir / "tss.tsv"),
                   tss_multiset(r_dir / "tss.tsv"))

    py_kind = spans_by_kind(py_dir / "feature.tsv")
    r_kind = spans_by_kind(r_dir / "feature.tsv")
    for kind in sorted(set(py_kind) | set(r_kind)):
        bad += describe(kind, py_kind.get(kind, Counter()),
                        r_kind.get(kind, Counter()))

    print(f"\nVERDICT: {'IDENTICAL' if bad == 0 else f'{bad} interval mismatches'}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
