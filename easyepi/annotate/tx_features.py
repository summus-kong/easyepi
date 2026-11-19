"""A transcription table in the shape ChIPseeker's TxDb presents.

Why this module exists
----------------------
ChIPseeker does not work from a flat GTF. It queries a ``TxDb`` for five
separate feature sets:

    intronsByTranscript(txdb)
    exonsBy(txdb)
    threeUTRsByTranscript(txdb)
    fiveUTRsByTranscript(txdb)
    transcripts(txdb)            (TSS = strand-aware start)

Each set is a GRangesList keyed by transcript, and a peak is assigned by
``findOverlaps`` against it. Two properties of that arrangement matter and are
lost in a flat, per-bin coordinate list:

* a feature belongs to a *transcript*, so the promoter window and the
  ``distanceToTSS`` both refer to a specific transcript start, not to the
  nearest start among all isoforms collapsed together;
* UTRs are explicit intervals supplied by the annotation, not something
  inferred by subtracting CDS from exons.

This reader therefore returns, per transcript, the five interval lists plus the
strand-aware TSS, so the classifier can reproduce ChIPseeker's assignment order
instead of approximating it.

Input formats
-------------
GTF (``gene_id``/``transcript_id`` attributes) and GFF3 (``ID``/``Parent``, with
comma-separated multi-parent children) are both accepted. A GFF3 child naming
several parents contributes to each of them: on a real TxDb dump 149,499 of
412,044 transcripts lost exons when only the first parent was honoured, which
silently emptied the Intron and UTR bins downstream.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

Interval = Tuple[int, int]


class TranscriptRecord:
    """One transcript, carrying the feature lists ChIPseeker queries."""

    __slots__ = ("chrom", "start", "end", "strand", "tss", "gene_id",
                 "transcript_id", "exons", "cds", "utr5", "utr3", "introns")

    def __init__(self, chrom: str, start: int, end: int, strand: str,
                 transcript_id: str, gene_id: str):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand or "+"
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        # TSS is the 5' end in the direction of transcription. Storing it here
        # (rather than recomputing from `start`) is what makes the promoter
        # window strand-aware.
        self.tss = self.start if self.strand != "-" else self.end
        self.exons: List[Interval] = []
        self.cds: List[Interval] = []
        self.utr5: List[Interval] = []
        self.utr3: List[Interval] = []
        self.introns: List[Interval] = []

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (f"<{self.transcript_id} {self.chrom}:{self.start}-{self.end}"
                f"{self.strand} exons={len(self.exons)}>")


# ---------------------------------------------------------------------------
# attribute parsing
# ---------------------------------------------------------------------------

def _split_gtf_attrs(field: str) -> Dict[str, str]:
    """Parse ``key "value";`` attributes.

    GTF carries no ``=``; the value is quoted. rtracklayer's ``export(...,
    "gtf")`` writes GFF2-style ``ID "TxID:1";`` / ``Parent "TxID:1";`` under a
    ``##gff-version 2`` header, so the same reader has to understand ``ID`` and
    ``Parent`` as well as ``gene_id``/``transcript_id``.
    """
    out: Dict[str, str] = {}
    for chunk in field.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if " " in chunk:
            key, _, val = chunk.partition(" ")
            out[key.strip()] = val.strip().strip('"')
    return out


def _split_gff3_attrs(field: str) -> Dict[str, str]:
    """Parse ``key=value`` attributes; values may be comma-separated lists."""
    out: Dict[str, str] = {}
    for chunk in field.split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        key, _, val = chunk.partition("=")
        out[key.strip()] = val.strip()
    return out


def detect_format(path) -> str:
    """Return ``"gff3"`` or ``"gtf"`` for an annotation file.

    The GFF3/GTF distinction is decided on attribute syntax rather than the file
    suffix, because exports are routinely mislabelled: a file named ``.gtf`` may
    contain ``key=value`` attributes and vice versa.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##gff-version"):
                return "gff3"
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attr = parts[8]
            # `key=value` (GFF3) never contains the `key "value"` shape that GTF
            # uses for every attribute, so the presence of `=` before the first
            # quote is decisive.
            return "gff3" if "=" in attr.split('"')[0] else "gtf"
    return "gtf"


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------

def read_transcripts(path) -> List[TranscriptRecord]:
    """Read a GTF or GFF3 into per-transcript records.

    Returns one record per transcript, with ``exons``/``cds``/``utr5``/``utr3``
    populated and ``introns`` derived as the gaps between consecutive exons.
    """
    fmt = detect_format(path)
    parse = _split_gff3_attrs if fmt == "gff3" else _split_gtf_attrs
    return _read_gff3(path, parse) if fmt == "gff3" else _read_gtf(path, parse)


def _new_record(chrom: str, start: int, end: int, strand: str,
                tx_id: str, gene_id: str,
                by_tx: Dict[str, TranscriptRecord]) -> TranscriptRecord:
    rec = by_tx.get(tx_id)
    if rec is None:
        rec = TranscriptRecord(chrom, start, end, strand, tx_id, gene_id)
        by_tx[tx_id] = rec
    else:
        # The record may have been created from a child line before its own
        # transcript line was seen, so widen rather than trust the first value.
        rec.start = min(rec.start, int(start))
        rec.end = max(rec.end, int(end))
        rec.tss = rec.start if rec.strand != "-" else rec.end
        if gene_id and not rec.gene_id:
            rec.gene_id = gene_id
    return rec


def _derive_utrs(rec: TranscriptRecord) -> None:
    """Fill ``utr5``/``utr3`` when the annotation carried no explicit UTR rows.

    gffread and most GTF exporters emit only ``exon`` and ``CDS`` (a real
    gffread export of a TxDb contained 2,305,263 exons and 983,861 CDS rows and
    zero ``UTR`` rows). ChIPseeker queries ``fiveUTRsByTranscript`` on the TxDb,
    which returns exonic bases outside the CDS blocks, so the same intervals have
    to be reconstructed here or the UTR bins come out empty and every promoter
    peak falls through to Intron.

    The subtraction is per exon against the merged CDS blocks, not against the
    global CDS span: an internal exon between the first and last CDS block has no
    coding bases yet is neither 5' nor 3' of the global span.
    """
    if rec.utr5 or rec.utr3 or not rec.exons or not rec.cds:
        return
    blocks = _merge_intervals(rec.cds)
    remaining = []
    for s, e in sorted(rec.exons):
        remaining.extend(_subtract_interval((s, e), blocks))
    if not remaining:
        return
    lo = min(b for _, b in rec.cds)
    hi = max(a for a, _ in rec.cds)
    five, three = [], []
    for s, e in remaining:
        if e <= lo:
            five.append((s, e))
        elif s >= hi:
            three.append((s, e))
    if rec.strand == "-":
        five, three = three, five
    rec.utr5, rec.utr3 = sorted(five), sorted(three)


def _merge_intervals(ivs):
    out = []
    for s, e in sorted((int(a), int(b)) for a, b in ivs):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _subtract_interval(interval, blocks):
    s, e = interval
    pieces = [(s, e)]
    for bs, be in blocks:
        nxt = []
        for ps, pe in pieces:
            if be <= ps or bs >= pe:
                nxt.append((ps, pe))
                continue
            if ps < bs:
                nxt.append((ps, bs))
            if be < pe:
                nxt.append((be, pe))
        pieces = nxt
        if not pieces:
            break
    return pieces


def _finish(by_tx: Dict[str, TranscriptRecord]) -> List[TranscriptRecord]:
    """Derive introns and sort, once all children are attached."""
    for rec in by_tx.values():
        _derive_utrs(rec)
        if len(rec.exons) > 1:
            ex = sorted(rec.exons)
            rec.introns = [(ex[i][1], ex[i + 1][0])
                           for i in range(len(ex) - 1)
                           if ex[i + 1][0] > ex[i][1]]
    return [r for r in by_tx.values() if r.exons or r.start <= r.end]


def _read_gtf(path, parse) -> List[TranscriptRecord]:
    by_tx: Dict[str, TranscriptRecord] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start, end, _s, strand = parts[:7]
            attrs = parse(parts[8])
            # GTF names the transcript explicitly. Falling back to `Parent` here
            # is wrong: a `gene` row carries `ID "1";` and no transcript id, so
            # the numeric gene id would become a transcript and every gene would
            # appear twice in the record set (824,088 records for 412,044
            # transcripts on a real TxDb export).
            tx_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id") or ""
            if not tx_id:
                continue
            s, e = int(start) - 1, int(end)      # GTF is 1-based closed
            feature = feature.lower()
            if feature in ("transcript", "mrna", "gene"):
                _new_record(chrom, s, e, strand, tx_id, gene_id, by_tx)
                continue
            rec = _new_record(chrom, s, e, strand, tx_id, gene_id, by_tx)
            if feature == "exon":
                rec.exons.append((s, e))
            elif feature == "cds":
                rec.cds.append((s, e))
            elif feature == "five_prime_utr":
                rec.utr5.append((s, e))
            elif feature == "three_prime_utr":
                rec.utr3.append((s, e))
            elif feature == "utr":
                # GTF has a single `UTR` feature; the side is decided by CDS.
                rec.utr5.append((s, e)) if _is_five_prime(rec, s, e) else rec.utr3.append((s, e))
    return _finish(by_tx)


def _is_five_prime(rec: TranscriptRecord, s: int, e: int) -> bool:
    if not rec.cds:
        return False
    lo = min(c[0] for c in rec.cds)
    hi = max(c[1] for c in rec.cds)
    return e <= lo if rec.strand != "-" else s >= hi


def _read_gff3(path, parse) -> List[TranscriptRecord]:
    by_tx: Dict[str, TranscriptRecord] = {}
    children: List[tuple] = []          # (feature, chrom, s, e, strand, parents)
    gene_of_tx: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start, end, _s, strand = parts[:7]
            attrs = parse(parts[8])
            s, e = int(start) - 1, int(end)      # GFF3 is 1-based closed
            feature = feature.lower()
            if feature in ("transcript", "mrna"):
                tx_id = attrs.get("ID")
                if not tx_id:
                    continue
                gene = attrs.get("Parent", "")
                gene_of_tx[tx_id] = gene
                _new_record(chrom, s, e, strand, tx_id, gene, by_tx)
                continue
            parent = attrs.get("Parent", "")
            if parent:
                # A shared exon names several parents; every one of them owns it.
                children.append((feature, chrom, s, e, strand,
                                 [p for p in parent.split(",") if p]))
    for feature, chrom, s, e, strand, parents in children:
        for p in parents:
            rec = _new_record(chrom, s, e, strand, p,
                              gene_of_tx.get(p, ""), by_tx)
            if feature == "exon":
                rec.exons.append((s, e))
            elif feature == "cds":
                rec.cds.append((s, e))
            elif feature in ("five_prime_utr", "utr"):
                if feature == "utr" and not _is_five_prime(rec, s, e):
                    rec.utr3.append((s, e))
                else:
                    rec.utr5.append((s, e))
            elif feature == "three_prime_utr":
                rec.utr3.append((s, e))
    return _finish(by_tx)
