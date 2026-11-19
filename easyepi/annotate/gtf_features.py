"""Build the ChIPseeker feature tables from a GTF, without R.

Why this exists
---------------
`chipseeker_port.py` reproduces `annotatePeak` exactly, but it consumes the
tables `export_features.R` writes, i.e. the TxDb objects themselves. Users have
a GTF, not a TxDb, so this module reproduces the same TxDb accessors from GTF
rows:

    tss.tsv      <- getGene(txdb, by = "transcript")
    exons        <- exonsBy(txdb)                  # note: by = "gene"
    introns      <- intronsByTranscript(txdb)
    5'/3' UTR    <- fiveUTRsByTranscript / threeUTRsByTranscript

The `by = "gene"` default for `exonsBy` is not a detail: ChIPseeker calls it
without `by`, so exons are grouped per gene. Using `by = "tx"` instead leaves
~28 large peaks unannotated on a 1,254-peak test set, which is why the
port honours the default here.

Coordinate convention: GTF is 1-based inclusive, the same as GRanges, so no
shift is applied. BED input is converted by the caller.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .chipseeker_port import Features


class _Tx:
    """One transcript row, in the order the GTF presented it."""

    __slots__ = ("chrom", "start", "end", "strand", "gene", "exons", "cds",
                 "utr5", "utr3", "stops")

    def __init__(self, chrom: str, start: int, end: int, strand: str,
                 gene: str) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.gene = gene
        self.exons: List[Tuple[int, int]] = []
        self.cds: List[Tuple[int, int]] = []
        self.utr5: List[Tuple[int, int]] = []
        self.utr3: List[Tuple[int, int]] = []
        self.stops: List[Tuple[int, int]] = []


def _attrs(text: str) -> Dict[str, str]:
    """Parse GTF attributes; GFF3-style `key=value` is accepted too."""
    out: Dict[str, str] = {}
    if "=" in text and ";" not in text.split("=")[0]:
        for part in text.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                out[k.strip()] = v.strip().strip('"')
        return out
    for part in text.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        k, v = part.split(" ", 1)
        out[k.strip()] = v.strip().strip('"')
    return out


def _merge(spans: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Collapse overlapping/adjacent intervals, as a TxDb group would."""
    if not spans:
        return []
    ordered = sorted(spans)
    out = [list(ordered[0])]
    for s, e in ordered[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def _introns(exons: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Gaps between a transcript's exons, matching intronsByTranscript.

    TxDb orders the exons but does not collapse overlapping ones, so the gaps
    are taken between consecutive distinct spans.
    """
    spans = sorted(set(exons))
    return [(spans[i][1] + 1, spans[i + 1][0] - 1)
            for i in range(len(spans) - 1)
            if spans[i + 1][0] - 1 >= spans[i][1] + 1]


def _id_sort_key(tx_id: str) -> Tuple[int, int, str]:
    """Sort transcript ids the way TxDb numbers them.

    TxDb assigns ``tx_id`` in ascending id order, so a purely numeric id sorts
    numerically ("512" before "9245"); ids with a prefix (``TxID:1234``) sort
    after the bare numbers and then numerically, and anything else falls back to
    lexical order. Sorting these lexically would put "9245" before "512" and
    change which transcript wins a multi-transcript overlap.
    """
    if tx_id.isdigit():
        return (0, int(tx_id), "")
    head, _, tail = tx_id.rpartition(":")
    if tail.isdigit():
        return (1, int(tail), head)
    return (2, 0, tx_id)


def _index_child(rows: Sequence[Tuple]) -> Tuple[List[Tuple], List[Tuple[int, int]]]:
    """Keep ``rows`` in array order and build the ``(start, index)`` lookup view.

    ChIPseeker's ``getFirstHitIndex`` keeps the overlapping feature with the
    lowest *subject index*, and the subject index is the array position produced
    by ``unlist(...)``. The row list therefore must not be sorted; the sort lives
    only in this auxiliary view, which is used to find candidates quickly.
    """
    return list(rows), sorted((r[0], j) for j, r in enumerate(rows))


def load_from_gtf(gtf_path: Path) -> Features:
    """Reproduce annotatePeak's TxDb view from a GTF file."""
    txs: Dict[str, _Tx] = {}
    order: List[str] = []

    with Path(gtf_path).open(encoding="utf-8") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start, end, _s, strand = parts[:7]
            attrs = _attrs(parts[8])
            tx_id = attrs.get("transcript_id")
            if not tx_id:
                continue
            s, e = int(start), int(end)
            rec = txs.get(tx_id)
            if rec is None:
                rec = txs[tx_id] = _Tx(chrom, s, e, strand,
                                       attrs.get("gene_id", ""))
                order.append(tx_id)
            else:
                # A transcript row may come after its children in some GTFs.
                rec.start = min(rec.start, s)
                rec.end = max(rec.end, e)
            if feature == "exon":
                rec.exons.append((s, e))
            elif feature == "CDS":
                rec.cds.append((s, e))
            elif feature in ("five_prime_utr", "5UTR"):
                rec.utr5.append((s, e))
            elif feature in ("three_prime_utr", "3UTR"):
                rec.utr3.append((s, e))
            elif feature == "stop_codon":
                rec.stops.append((s, e))

    feats = Features()
    # TxDb numbers transcripts by ascending id, and ChIPseeker's
    # `GRegion <- unlist(intronsByTranscript(txdb))` therefore orders features in
    # that numeric order. `getFirstHitIndex` then keeps the hit with the lowest
    # subject index, i.e. the smallest transcript id -- measured on
    # chr1:174127900, where tx 13017 sits at subject 62198 and tx 30843 at
    # 149486, so R reports 13017. GTF file order is not that order, so it is
    # sorted here rather than trusted.
    order.sort(key=_id_sort_key)
    feats.gene_of = {tx: txs[tx].gene for tx in order if txs[tx].gene}

    for tx_id in order:
        rec = txs[tx_id]
        # resize(width = 1) is strand aware: the left edge on +, right on -.
        tss = rec.start if rec.strand != "-" else rec.end
        idx = len(feats.strand)
        feats.strand.append(rec.strand if rec.strand in ("+", "-") else "+")
        feats.tss.setdefault(rec.chrom, []).append((int(tss), idx))
        feats.span.append((int(rec.start), int(rec.end)))
        feats.chrom.append(rec.chrom)

        # ChIPseeker reports "Intron (<tx_name>/<gene_id>, intron k of N)".
        # TXID2TXEG pastes tx_name and gene_id, and when a transcript has no
        # gene_id it first substitutes tx_name -- which yields the *duplicated*
        # "TxID:9245/TxID:9245", because the "/NA" strip only fires on a literal
        # NA that the substitution has already replaced. Reproducing that
        # duplication is the difference between matching and not.
        gene_tag = rec.gene if rec.gene else tx_id
        label = f"{tx_id}/{gene_tag}"

        introns = _introns(rec.exons)
        n_int = len(introns)
        for i, (s, e) in enumerate(introns):
            # intron_rank is 1-based in transcript direction: on the minus
            # strand R does rev(seq_len(N)), so the first row gets N.
            rank = n_int - i if rec.strand == "-" else i + 1
            per = feats.child.setdefault("intron", {}).setdefault(rec.chrom, [])
            # Final field is insertion order, mirroring unlist(intronsByTranscript()):
            # it decides which transcript wins a multi-transcript overlap.
            per.append((s, e, idx, label, rank, n_int, len(per)))
            if e - s > feats.child_max_len.get("intron", 0):
                feats.child_max_len["intron"] = e - s

        # R stores a transcript's exons the way exonsBy() returns them: on the
        # minus strand that is DESCENDING genomic order, with exon_rank = 1..N
        # aligned to that array. Verified on tx 18518 (n = 15, strand -):
        # starts are 1348896, 1342689, 1342363, ... and exon_rank is 1..15.
        # Keeping the array in that same order matters twice over: the reported
        # rank, and which row getFirstHitIndex keeps when several exons overlap
        # one peak (it keeps the lowest array index). An ascending array with
        # N-i ranking gives the right rank but picks the wrong row.
        n_ex = len(rec.exons)
        exons = sorted(rec.exons, reverse=(rec.strand == "-"))
        for i, (s, e) in enumerate(exons):
            per = feats.child.setdefault("exon", {}).setdefault(rec.chrom, [])
            per.append((s, e, idx, label, i + 1, n_ex, len(per)))
            if e - s > feats.child_max_len.get("exon", 0):
                feats.child_max_len["exon"] = e - s

        for kind, spans in (
            ("utr5", rec.utr5 or _derived_utr(rec, "5")),
            ("utr3", _trim_stop_codon(rec.utr3 or _derived_utr(rec, "3"),
                                      rec.stops, rec.strand)),
        ):
            per = feats.child.setdefault(kind, {})
            for s, e in spans:
                rows = per.setdefault(rec.chrom, [])
                rows.append((s, e, idx, label, 0, 0, len(rows)))
                if e - s > feats.child_max_len.get(kind, 0):
                    feats.child_max_len[kind] = e - s

    for chrom, rows in feats.tss.items():
        rows.sort()
    # Row lists keep the assembly (array) order, which is ChIPseeker's subject
    # index; sorting them would silently reassign getFirstHitIndex's choice.
    for per in feats.child.values():
        for chrom, rows in per.items():
            per[chrom] = _index_child(rows)

    # genes(TxDb) is range(transcriptsBy(txdb, by = "gene")). A transcript with
    # no gene_id is its own group and is named by its transcript id, which is
    # why genes() names are a GeneID/TxID mixture (35,356 + 85,737 = 121,093).
    # Dropping the gene_id-less transcripts loses 85,737 of them and shifts
    # every gene index with it.
    _build_gene_view(feats, [
        (txs[tx_id].gene or tx_id, txs[tx_id].chrom,
         int(txs[tx_id].start), int(txs[tx_id].end), txs[tx_id].strand)
        for tx_id in order
    ])
    return feats


def _build_gene_view(feats: Features, records: Sequence[Tuple[str, str, int, int, str]]) -> None:
    """Populate R's ``genes(TxDb)`` view, used by ``level = "gene"``.

    ``records`` is ``(gene_id, chrom, start, end, strand)`` per transcript, in
    ``unlist(transcriptsBy())`` order.

    ``genes()`` is ``range(transcriptsBy(txdb, by = "gene"))`` with the groups
    that span more than one seqlevel or more than one strand dropped. Measured on
    the GTF-derived TxDb: 123,262 groups -> 121,093 genes, 2,169 dropped for
    spanning seqlevels and 238 for spanning strands.

    The gene span is the transcript min-start/max-end, exactly as ``range()``
    builds it; the gene's strand-side edge (start on "+", end on "-") is the TSS
    ``follow()`` sees in gene mode.
    """
    per_gene: Dict[str, Dict[str, object]] = {}
    for gene, chrom, s, e, strand in records:
        g = per_gene.get(gene)
        if g is None:
            per_gene[gene] = {"chrom": chrom, "strand": strand,
                              "start": s, "end": e, "multi": False}
            continue
        if g["chrom"] != chrom or g["strand"] != strand:
            g["multi"] = True
            continue
        g["start"] = min(int(g["start"]), s)
        g["end"] = max(int(g["end"]), e)

    for gene, g in per_gene.items():
        if g["multi"]:
            continue
        chrom = g["chrom"]
        s, e = int(g["start"]), int(g["end"])
        strand = g["strand"]
        tss = s if strand == "+" else e
        feats.tss_gene.setdefault(chrom, []).append((tss, gene))
        feats.span_gene.setdefault(chrom, []).append((s, e, strand, gene))

    for chrom, rows in feats.tss_gene.items():
        rows.sort()
    for chrom, rows in feats.span_gene.items():
        feats.span_gene[chrom] = sorted(rows, key=lambda r: r[1])


def write_tables(gtf_path: Path, outdir: Path) -> Dict[str, int]:
    """Write tx.tsv/tss.tsv/feature.tsv, mirroring ``export_features.R``.

    Exists so the GTF path can be compared row by row against the TxDb path.
    Without this, "the two agree" would rest on the peaks happening to be
    annotated identically, which is a weaker statement than the tables matching.
    """
    import csv

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    feats = load_from_gtf(gtf_path)

    # Recover the per-transcript rows in file order, which is the order the
    # feature indices in `feats` were assigned.
    txs: Dict[int, List] = {}
    for chrom, rows in feats.tss.items():
        for pos, idx in rows:
            txs[idx] = [chrom, pos]

    with (outdir / "tss.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "chrom", "pos", "strand"])
        for idx in sorted(txs):
            chrom, pos = txs[idx]
            w.writerow([f"T{idx}", chrom, pos, feats.strand[idx]])

    n_feat = 0
    with (outdir / "feature.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "kind", "chrom", "start", "end", "gene", "rank", "total"])
        for kind, per in sorted(feats.child.items()):
            for chrom, entry in sorted(per.items()):
                # child[chrom] is (rows, order); rows keep array order.
                rows = entry[0]
                for s, e, owner, gene, rank, total, _order in rows:
                    w.writerow([f"T{owner}", kind, chrom, s, e, gene, rank, total])
                    n_feat += 1

    # tx.tsv carries the full span: the Downstream pass follows the unresized
    # feature, so a TSS-only table cannot reproduce it.
    with (outdir / "tx.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "chrom", "start", "end", "strand"])
        for idx in sorted(txs):
            chrom, _pos = txs[idx]
            s, e = feats.span[idx]
            w.writerow([f"T{idx}", chrom, s, e, feats.strand[idx]])

    return {"tss": len(txs), "feature": n_feat}


def _trim_stop_codon(
    spans: Sequence[Tuple[int, int]],
    stops: Sequence[Tuple[int, int]],
    strand: str,
) -> List[Tuple[int, int]]:
    """Remove the 3' UTR bases that overlap the stop codon.

    TxDb's threeUTRsByTranscript excludes the stop codon, but GENCODE's UTR
    row starts at the first base of it: a UTR of 70006-71585 with a stop codon
    at 70006-70008 comes back as 70009-71585. Without this trim the interval is
    three bases too wide on every such transcript, which is enough to move a
    boundary peak from 3'UTR to Exon.
    """
    if not stops:
        return list(spans)
    out: List[Tuple[int, int]] = []
    for s, e in spans:
        pieces = [(s, e)]
        for ss, se in stops:
            nxt: List[Tuple[int, int]] = []
            for ps, pe in pieces:
                if se < ps or ss > pe:
                    nxt.append((ps, pe))
                    continue
                # Keep the part of the UTR beyond the stop codon: on the plus
                # strand that is the right-hand remainder, on minus the left.
                if strand != "-":
                    if se < pe:
                        nxt.append((se + 1, pe))
                else:
                    if ss > ps:
                        nxt.append((ps, ss - 1))
            pieces = nxt
        out.extend(p for p in pieces if p[1] >= p[0])
    return out


def _derived_utr(rec: _Tx, end: str) -> List[Tuple[int, int]]:
    """UTR as exon minus CDS, when the GTF has no explicit UTR rows."""
    if not rec.cds or not rec.exons:
        return []
    cds_lo = min(s for s, _ in rec.cds)
    cds_hi = max(e for _, e in rec.cds)
    out: List[Tuple[int, int]] = []
    for s, e in rec.exons:
        if end == "5":
            hi = min(e, cds_lo - 1) if rec.strand != "-" else e
            lo = s if rec.strand != "-" else max(s, cds_hi + 1)
        else:
            lo = max(s, cds_hi + 1) if rec.strand != "-" else s
            hi = e if rec.strand != "-" else min(e, cds_lo - 1)
        if hi >= lo:
            out.append((lo, hi))
    return out
