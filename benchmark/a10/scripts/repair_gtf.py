"""Repair a GTF whose feature rows are missing `gene_id`, so a TxDb can be built.

Why this is needed for A10
--------------------------
`txdb_full.gtf` (the file EasyEpi annotates with) is unusable as a gene model:
its `transcript` rows carry ENST ids while its `exon`/`CDS`/`UTR` rows carry the
UCSC numeric *gene* ids and no attribute links the two namespaces. Consequences
measured on the A10 heart data:

  * a TxDb built from it reports 824,088 transcripts for 412,044 transcript rows,
    because every numeric id also becomes a transcript;
  * EasyEpi, which builds introns from the gaps between one transcript's exons,
    then forms 112 kb "introns" spanning whole genes and swallows real exons;
  * that produced 1,257 peak-level `EasyEpi=Intron / ChIPseeker=Exon` disputes.

`txdb_hg38.gffread.gtf` does link the two namespaces through `transcript_id`
(412,044 of 412,044 exon ids resolve to a transcript row), but 196,282 of its
5'UTR/3'UTR/exon rows omit `gene_id`, which makes `makeTxDbFromGFF` abort with
"'genes$gene_id' must be a character vector with no NAs".

The repair is deliberately minimal and auditable: fill a missing `gene_id` from
the `transcript_id` -> `gene_id` map learned from the transcript rows of the same
file. No coordinate is touched and no row is added or dropped.

Usage: repair_gtf.py <in.gtf> <out.gtf>
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

GENE_ID_QUOTED = re.compile(r'gene_id\s+"([^"]*)"')
TX_ID_QUOTED = re.compile(r'transcript_id\s+"([^"]*)"')
GENE_ID_BARE = re.compile(r"gene_id=([^;]*)")
TX_ID_BARE = re.compile(r"transcript_id=([^;]*)")


def get(attr: str, quoted: re.Pattern, bare: re.Pattern) -> str | None:
    m = quoted.search(attr)
    if m:
        return m.group(1)
    m = bare.search(attr)
    return m.group(1).strip() if m else None


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])

    # Pass 1: transcript_id -> gene_id, from rows that carry both.
    t2g: dict[str, str] = {}
    with src.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attr = parts[8]
            tx = get(attr, TX_ID_QUOTED, TX_ID_BARE)
            gene = get(attr, GENE_ID_QUOTED, GENE_ID_BARE)
            if tx and gene:
                t2g[tx] = gene
    print(f"learned {len(t2g):,} transcript->gene mappings")

    # Pass 2: copy through, filling only the rows that are missing gene_id.
    n_fill = n_unfillable = n_rows = 0
    with src.open(encoding="utf-8", errors="replace") as fin, \
            dst.open("w", encoding="utf-8", newline="") as fout:
        for line in fin:
            if line.startswith("#") or "\t" not in line:
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 9:
                n_rows += 1
                attr = parts[8]
                if not get(attr, GENE_ID_QUOTED, GENE_ID_BARE):
                    tx = get(attr, TX_ID_QUOTED, TX_ID_BARE)
                    gene = t2g.get(tx or "")
                    if gene is None:
                        n_unfillable += 1
                    else:
                        sep = "" if attr.rstrip().endswith(";") else ";"
                        parts[8] = f'{attr}{sep} gene_id "{gene}";'
                        line = "\t".join(parts) + "\n"
                        n_fill += 1
            fout.write(line)

    print(f"rows: {n_rows:,}  gene_id filled: {n_fill:,}  still missing: {n_unfillable:,}")
    print(f"written -> {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
