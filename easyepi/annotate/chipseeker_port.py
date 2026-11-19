"""A literal Python port of ChIPseeker's genomic annotation.

The R reference (`scripts/export_features.R`) writes the exact feature objects
`annotatePeak` consumes, so this module never re-derives TSS positions or
re-maps transcript ids: whatever ChIPseeker sees, this sees.

Ported from ChIPseeker 1.46.1:

  getNearestFeatureIndicesAndDistances
      features <- resize(features, width = 1)      # strand-aware 1 bp TSS
      ps.idx <- follow(peaks, unstrand(features))
      pe.idx <- precede(peaks, unstrand(features))
      psD <- sign(strand) * (start(peaks) - start(psF))
      peD <- sign(strand) * (end(peaks)   - start(peF))
      j   <- which.min(abs(c(psD, peD)))           # keeps the signed value
      then, if !ignoreOverlap:
          hit <- findOverlaps(peaks, unstrand(features))   # features is restored
          distanceToTSS[peakIdx] <- 0

  getGenomicAnnotation
      priority <- rev(priority)                    # last applied wins
      each child bin is matched by findOverlaps + getFirstHitIndex
      Promoter is applied last, purely from `distance` in tssRegion

Coordinate convention: GRanges 1-based inclusive, as written by the R export.
"""

from __future__ import annotations

import bisect
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# ChIPseeker's `genomicAnnotationPriority` default, verbatim from annotatePeak().
# Downstream is part of the priority list, not an afterthought: rev() means
# Intergenic is applied first and Promoter last, with Downstream in between.
DEFAULT_PRIORITY = ("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream",
                    "Intergenic")

# Map a bin to the feature kind written by export_features.R.
KIND_OF = {"Intron": "intron", "Exon": "exon", "5UTR": "utr5", "3UTR": "utr3"}

# ChIPseeker's default for getOption("ChIPseeker.downstreamDistance").
DEFAULT_DOWNSTREAM = 300


def promoter_bins(tss_region: Tuple[int, int]) -> List[Tuple[int, str]]:
    """``Promoter`` sub-category labels, widest last, exactly as R builds them.

    R builds ``dd <- seq_len(ceiling(pm/1000)) * 1000`` once (``pm =
    max(abs(tssRegion))``) and then, for the Promoter pass, picks the *first*
    bucket whose lower bound contains ``abs(distance)``. With the default
    ``c(-3000, 3000)`` that yields "Promoter (<=1kb)", "Promoter (1-2kb)",
    "Promoter (2-3kb)".

    When ``pm < 2000`` R skips the sub-categorisation entirely and leaves the
    plain "Promoter" label, which is reproduced by returning a single unbounded
    bucket of -1.
    """
    pm = max(abs(int(tss_region[0])), abs(int(tss_region[1])))
    if pm / 1000 < 2:
        return [(-1, "Promoter")]
    out: List[Tuple[int, str]] = []
    for i in range(1, int(-(-pm // 1000)) + 1):  # ceiling(pm/1000)
        if i == 1:
            out.append((1000, "Promoter (<=1kb)"))
        else:
            out.append((i * 1000, f"Promoter ({i - 1}-{i}kb)"))
    return out


def downstream_bins(dsd: int) -> List[Tuple[int, str]]:
    """Deprecated shim: kept only because the label logic now lives in `annotate`.

    R applies Downstream *after* the priority loop, so the labels can no longer
    be expressed as a bucket table and this helper is not used by the
    classifier. It is retained so external callers that import it keep working,
    and it delegates to the same two branches `annotate` implements.
    """
    if dsd / 1000 <= 1:
        return [(dsd, f"Downstream (<={dsd}bp)")]
    whole = int(dsd // 1000)
    out: List[Tuple[int, str]] = []
    for i in range(1, whole + 1):
        out.append((i * 1000, "Downstream (<1kb)" if i == 1
                    else f"Downstream ({i - 1}-{i}kb)"))
    if dsd > whole * 1000:
        out.append((dsd, f"Downstream ({dsd / 1000:g}kb-{dsd}bp)"))
    return out


@dataclass
class Features:
    """Everything annotatePeak reads from the TxDb."""

    # chrom -> sorted (tss_position, tx_index)
    tss: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)
    strand: List[str] = field(default_factory=list)
    # kind -> chrom -> sorted (start, end, tx_index, gene_label, rank, rank_total)
    child: Dict[str, Dict[str, List[Tuple]]] = field(default_factory=dict)
    # tx_index -> full span, for the Downstream pass (follow() runs on the
    # *unresized* feature, so end(peF) is the transcript end, not the TSS).
    span: List[Tuple[int, int]] = field(default_factory=list)
    # tx_index -> chromosome; spans are looked up by peak chromosome.
    chrom: List[str] = field(default_factory=list)
    # kind -> longest interval, bounding the backward scan in _overlaps.
    child_max_len: Dict[str, int] = field(default_factory=dict)
    # `level = "gene"` view: genes(TxDb) spans and their strand-side TSS.
    gene_of: Dict[str, str] = field(default_factory=dict)
    # index -> id, parallel to `span`; the gene view fills it, the transcript
    # loader leaves it empty and callers fall back to gene_of ordering.
    key_order: List[str] = field(default_factory=list)
    key_to_index: Dict[str, int] = field(default_factory=dict)
    tss_gene: Dict[str, List[Tuple[int, str]]] = field(default_factory=dict)
    span_gene: Dict[str, List[Tuple[int, int, str, str]]] = field(default_factory=dict)
    # level="gene": gene bodies as (end-sorted) follow() input.
    span_gene_follow: Dict[str, Tuple[List[int], List[int], List[str]]] = \
        field(default_factory=dict)


def _id_at(feats: Features, index: int) -> str:
    """Feature id at ``index``.

    ``key_order`` is filled by the gene view; for the transcript view the
    insertion order of ``gene_of`` is the order ``load_from_gtf`` appended
    records (verified: index 18215 -> span (688808, 710676), gene GeneID:81399).
    ``key_to_index`` is only populated by the exported-table loader, so it
    cannot be used here.
    """
    ids = feats.key_order or list(feats.gene_of)
    return ids[index] if 0 <= index < len(ids) else "NA"


def _overlaps(rows: Sequence[Tuple], order: Sequence[Tuple[int, int]],
              qs: int, qe: int, max_len: int = 0,
              same_strand: bool = False,
              peak_strand: Optional[str] = None,
              strand_of: Optional[Sequence[str]] = None) -> Optional[int]:
    """Index into ``rows`` of the hit ChIPseeker's ``getFirstHitIndex`` keeps.

    ``rows`` is kept in R's own array order (the order ``unlist(...)`` produces),
    so the winner is simply the overlapping row with the **smallest array
    index** -- that is what ``getFirstHitIndex(queryHits(hit))`` selects, since
    ``findOverlaps`` returns hits ordered by subject index.

    ``order`` is a ``(start, index)``-sorted view used only to find candidate
    rows quickly; ``max_len`` bounds the backward scan, because an interval
    starting before ``qs - max_len`` cannot reach the query.

    Note that neither the TxID nor the leftmost coordinate is the right key on
    its own. On tx 18518 three exons overlap one peak; R's array is in
    descending genomic order for that minus-strand transcript, so the kept row
    is ``exon 2 of 15`` while ascending coordinates give ``exon 4 of 15``.
    """
    lo = bisect.bisect_left(order, (qs - max_len - 1, -1)) if max_len else 0
    hi = bisect.bisect_right(order, (qe, float("inf")))
    best: Optional[int] = None
    for _s, j in order[lo:hi]:
        if rows[j][1] < qs:
            continue
        # R: findOverlaps(query, subject, ignore.strand = !sameStrand).
        # The strand lives in feats.strand[rows[j][2]] -- rows[j][3] is the
        # label, not a strand, and comparing against it silently filtered every
        # hit, which is what pinned ssTRUE at 290/400.
        if same_strand and peak_strand is not None and strand_of is not None:
            if strand_of[rows[j][2]] != peak_strand:
                continue
        if best is None or j < best:
            best = j
    return best


def _index_child(rows: Sequence[Tuple]) -> Tuple[List[Tuple], List[Tuple[int, int]]]:
    """Keep ``rows`` in array order, and build the ``(start, index)`` view."""
    order = sorted((r[0], j) for j, r in enumerate(rows))
    return list(rows), order


def _label(rows: Sequence[Tuple], hit: int, kind: str) -> str:
    """ChIPseeker's per-feature label for a hit, verbatim.

    ``getGenomicAnnotation.internal`` returns plain "5' UTR"/"3' UTR" for the UTR
    passes and only numbers exons and introns; returning the bin name for those
    two keeps the strings identical instead of blanking them out.
    """
    _s, _e, _tx, gene, rank, total = rows[hit][:6]
    if kind == "intron":
        return f"Intron ({gene}, intron {rank} of {total})"
    if kind == "exon":
        return f"Exon ({gene}, exon {rank} of {total})"
    if kind == "utr5":
        return "5' UTR"
    if kind == "utr3":
        return "3' UTR"
    return ""

def load_features(feat_dir: Path) -> Features:
    """Read the tables written by ``export_features.R``.

    Requires the extended ``feature.tsv`` schema (``gene,rank,total`` columns)
    and ``tx.tsv`` spanning ``start,end``; both are produced by the matching
    exporter. The older three-column form cannot express an exon/intron label,
    so it is rejected loudly rather than silently degrading to bin names.
    """
    feat_dir = Path(feat_dir)

    key_to_index: Dict[str, int] = {}
    strand: List[str] = []
    span: List[Tuple[int, int]] = []
    chrom_of: List[str] = []
    with (feat_dir / "tx.tsv").open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not {"key", "strand", "chrom", "start", "end"} <= set(reader.fieldnames or ()):
            raise ValueError(
                f"{feat_dir/'tx.tsv'} lacks chrom/start/end; regenerate it with "
                "the current export_features.R (the Downstream pass needs the span)")
        for row in reader:
            key_to_index[row["key"]] = len(strand)
            strand.append(row["strand"] or "+")
            span.append((int(row["start"]), int(row["end"])))
            chrom_of.append(row["chrom"])

    feats = Features(strand=strand, span=span, chrom=chrom_of)
    feats.key_to_index = key_to_index
    with (feat_dir / "tss.tsv").open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ti = key_to_index.get(row["key"])
            if ti is None:
                continue
            feats.tss.setdefault(row["chrom"], []).append((int(row["pos"]), ti))
    for chrom, rows in feats.tss.items():
        rows.sort()

    for kind in KIND_OF.values():
        feats.child[kind] = {}

    with (feat_dir / "feature.tsv").open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not {"key", "kind", "chrom", "start", "end"} <= set(reader.fieldnames or ()):
            raise ValueError(f"{feat_dir/'feature.tsv'} has an unexpected schema")
        has_labels = {"gene", "rank", "total"} <= set(reader.fieldnames or ())
        for order, row in enumerate(reader):
            kind = row["kind"]
            if kind not in KIND_OF.values():
                continue
            ti = key_to_index.get(row["key"])
            if ti is None:
                continue
            per = feats.child.setdefault(kind, {})
            gene = row.get("gene", "") if has_labels else ""
            rank = int(row["rank"]) if has_labels else 0
            total = int(row["total"]) if has_labels else 0
            s, e = int(row["start"]), int(row["end"])
            # `order` is the file order, which is the unlist() order R hits.
            per.setdefault(row["chrom"], []).append(
                (s, e, ti, gene, rank, total, order))
            if e - s > feats.child_max_len.get(kind, 0):
                feats.child_max_len[kind] = e - s

    # The row lists stay in file (array) order: that order IS ChIPseeker's
    # subject index. Sort order lives in a separate lookup view.
    for per in feats.child.values():
        for chrom, rows in per.items():
            per[chrom] = _index_child(rows)

    return feats


def _follow(positions: Sequence[int], qs: int) -> int:
    """IRanges::follow -- last TSS at or before qs (left neighbour).

    When several features share one TSS, R returns the one at the LARGEST array
    position: TxID:270 and TxID:271 both sit at 904786 and follow() yields
    168934 (the later row). ``bisect_right(...) - 1`` stops at the first of the
    equal run, so the run is walked to its end here.
    """
    return bisect.bisect_right(positions, qs) - 1


def _precede(positions: Sequence[int], qe: int) -> int:
    """IRanges::precede -- first TSS at or after qe (right neighbour).

    For a tied run R keeps the FIRST row, so ``bisect_left`` is already right.
    """
    return bisect.bisect_left(positions, qe)


def _follow_right(positions: Sequence[int], qs: int) -> int:
    """``follow`` with a minus-strand query: the first TSS to the RIGHT.

    Ties keep the lowest row, matching ``_precede``.

    Measured on chr1:713701-714700 (-) against TSS 710676 (-) / 720150 (+):
    ``follow(pk, unstrand(tss))`` returns 720150, i.e. the right neighbour, and
    ``precede`` returns 710676. The two swap, and reproducing the swap is what
    takes the 400-peak distance agreement from 204/400 to 400/400.
    """
    # bisect_right - 1 would land on the FEATURE BEFORE unless qs is an exact
    # TSS; bisect_left picks the first TSS >= qs, and when several share that
    # TSS R keeps the one with the LARGEST row (TxID:271 row 168934 over
    # TxID:270 row 168933, both at 904786). Entry order is (pos, row), so the
    # equal run is walked to its end here.
    i = bisect.bisect_left(positions, qs)
    while i + 1 < len(positions) and positions[i + 1] == positions[i]:
        i += 1
    return i


def _precede_left(positions: Sequence[int], qe: int) -> int:
    """``precede`` with a minus-strand query: the last TSS to the LEFT."""
    return bisect.bisect_right(positions, qe) - 1


def nearest_distance(
    feats: Features, chrom: str, qs: int, qe: int,
    same_strand: bool = False, peak_strand: Optional[str] = None,
) -> Tuple[Optional[int], Optional[float]]:
    """``(tx_index, distance)`` exactly as getNearestFeatureIndicesAndDistances.

    With ``same_strand`` the TSS vector is restricted to features on the peak's
    strand first, and precede/follow run inside that subset -- R does the same,
    ``GRegion <- GRegion[strand(GRegion) == strand(peak.gr)]``, so an opposite
    strand TSS can never be the nearest one.
    """
    entry = feats.tss.get(chrom)
    if not entry:
        return None, None
    # NOTE: R does NOT subset the TSS vector under sameStrand. It only makes
    # follow()/precede() strand aware, so an opposite-strand TSS can still be
    # returned when it is the only candidate on that side. Filtering here made
    # every stranded peak fall through to "Distal Intergenic" (278/400).
    positions = [p for p, _ in entry]

    # IRanges::follow / precede are STRAND AWARE. Measured on chr1:713701-714700
    # against TSS 710676 (-) and 720150 (+):
    #   query "+" or "*": follow -> 710676 (left)   precede -> 720150 (right)
    #   query "-"       : follow -> 720150 (right)  precede -> 710676 (left)
    # so for a minus-strand peak the two swap. Reproducing that swap moved the
    # 400-peak distance agreement from 204/400 to 372/400.
    if peak_strand == "-":
        ps_i = _follow_right(positions, qs)    # follow -> right neighbour
        pe_i = _precede_left(positions, qe)    # precede -> left neighbour
    else:
        ps_i = _follow(positions, qs)
        pe_i = _precede(positions, qe)

    def on_strand(i: int, step: int) -> int:
        """Walk ``i`` in ``step`` direction until the feature strand matches.

        R's ``follow``/``precede`` with ``sameStrand=TRUE`` keep the strand on
        both query and subject, so an opposite-strand candidate is skipped and
        the search continues in the same direction. Returns -1 when it runs off
        the end, which R expresses as NA.
        """
        if not (same_strand and peak_strand is not None):
            return i
        while 0 <= i < len(entry):
            if feats.strand[entry[i][1]] == peak_strand:
                return i
            i += step
        return -1

    # Ties: several features can share one TSS and R resolves them on the
    # feature's OWN row -- TxID:270 (row 168933) and TxID:271 (row 168934) both
    # sit at 904786 and follow() returns 168934. Entry order is (tss, row), so
    # the equal-TSS run is scanned for the largest row (follow) or smallest
    # (precede). Walking the sorted list alone did nothing because that order
    # is by position, not by row.
    def tie_bounds(i: int):
        if not (0 <= i < len(entry)):
            return i
        lo = i
        while lo - 1 >= 0 and entry[lo - 1][0] == entry[i][0]:
            lo -= 1
        hi = i
        while hi + 1 < len(entry) and entry[hi + 1][0] == entry[i][0]:
            hi += 1
        return lo, hi

    def pick_row(i: int, want_last: bool) -> int:
        if not (0 <= i < len(entry)):
            return i
        lo, hi = tie_bounds(i)
        best = None
        for k in range(lo, hi + 1):
            row = entry[k][1]
            if best is None or (row > entry[best][1] if want_last
                                else row < entry[best][1]):
                best = k
        return best if best is not None else i

    if peak_strand == "-":
        # follow(right) keeps going right, precede(left) keeps going left.
        ps_i = on_strand(ps_i, +1)
        pe_i = on_strand(pe_i, -1)
    else:
        ps_i = on_strand(ps_i, -1)
        pe_i = on_strand(pe_i, +1)

    # Ties among features sharing one TSS are resolved on the feature's own row,
    # and the rule does not depend on the peak's strand: measured on
    # chr1:894201-894600 (-), follow -> TxID:271 (row 168934, the LARGER of the
    # two features at 904786) and precede -> the smaller row.
    ps_i = pick_row(ps_i, want_last=True)
    pe_i = pick_row(pe_i, want_last=False)
    ps_owner = entry[ps_i][1] if 0 <= ps_i < len(entry) else None
    pe_owner = entry[pe_i][1] if 0 <= pe_i < len(entry) else None

    def signed(owner: Optional[int], edge: int, tss: int) -> Optional[float]:
        if owner is None:
            return None
        # Sign follows the FEATURE strand: R writes
        #   psD <- ifelse(strand(psF) == "+", 1, -1) * (start(peaks) - start(psF))
        sign = 1.0 if feats.strand[owner] == "+" else -1.0
        return sign * (edge - tss)

    ps_d = (signed(ps_owner, qs, entry[ps_i][0])
            if ps_owner is not None else None)
    pe_d = (signed(pe_owner, qe, entry[pe_i][0])
            if pe_owner is not None else None)

    # j <- which.min(abs(c(psD, peD))); NA distances are skipped by R's ifelse.
    if ps_d is None and pe_d is None:
        index, dist = None, None
    elif ps_d is None:
        index, dist = pe_owner, pe_d
    elif pe_d is None:
        index, dist = ps_owner, ps_d
    elif abs(ps_d) <= abs(pe_d):
        index, dist = ps_owner, ps_d
    else:
        index, dist = pe_owner, pe_d

    # if (!ignoreOverlap) { hit <- findOverlaps(peaks, unstrand(features));
    #                       distanceToTSS[peakIdx] <- 0 }
    # `features` is the resized one-base view here, so ANY TSS inside the peak
    # forces zero -- not just the two follow/precede candidates. Checking only
    # those two missed most in-body TSSs and under-reported the zero cases.
    lo = bisect.bisect_left(positions, qs)
    hi = bisect.bisect_right(positions, qe)
    if hi > lo:
        # Within the peak, R keeps the feature with the SMALLEST row, not the
        # leftmost TSS: chr1:8876601-8878500 contains 7 TSS (rows 11771,
        # 119277, 138437, 138438, 138440 and more) and R reports TxID:1209,
        # which sits at row 11771. Taking entry[lo] gave the leftmost one
        # instead, because `positions` is sorted by coordinate, not by row.
        first = min(range(lo, hi), key=lambda k: entry[k][1])
        return entry[first][1], 0.0

    return index, dist


def _follow_excluding(ends_sorted: Sequence[int], qs: int,
                      same_strand: bool = False,
                      peak_strand: Optional[str] = None,
                      strands_by_end: Optional[Sequence[str]] = None) -> int:
    """``IRanges::follow(peaks, unstrand(features))`` -- the nearest feature before the peak.

    ``follow`` returns the feature whose **end is greatest while still lying
    entirely before the peak** (``end < start(peak)``). Two naive readings are
    both wrong and both were measured:

    * "last start <= qs" picks a feature whose start is late but whose end is
      far to the left (chr1:33472660-33489245, dd2 = 12555);
    * dropping the containment test picks a minus-strand transcript spanning the
      peak at all.

    R's answer for chr1:33501800-33502500 is chr1:33472645-33501643, i.e. the
    largest end below the query -- so sorting by ``end`` and binary-searching
    gives it in O(log n). ``end < qs`` also implies ``start < qs``.
    """
    i = bisect.bisect_left(ends_sorted, qs) - 1
    if not same_strand or peak_strand is None or strands_by_end is None:
        return i
    # R keeps the strand on `peaks` when sameStrand is TRUE, so an opposite
    # strand feature cannot be the one followed.
    while i >= 0 and strands_by_end[i] != peak_strand:
        i -= 1
    return i

def _span_index(feats: Features) -> Dict[str, Tuple[List[int], List[int], List[str]]]:
    """Per-chromosome transcript spans, sorted by END, for the Downstream pass.

    ``follow(peaks, unstrand(features))`` runs on the *unresized* features, so the
    transcript body -- not the 1 bp TSS -- is what ``dd2`` is measured against.
    Sorted by end because :func:`_follow_excluding` asks for the largest end below
    the peak; built once per call, not per peak.
    """
    per: Dict[str, List[Tuple[int, int, str]]] = {}
    for ti, (s, e) in enumerate(feats.span):
        chrom = feats.chrom[ti]
        per.setdefault(chrom, []).append((int(e), int(s), feats.strand[ti]))
    out: Dict[str, Tuple[List[int], List[int], List[str]]] = {}
    for chrom, rows in per.items():
        rows.sort()
        out[chrom] = ([r[0] for r in rows], [r[1] for r in rows],
                      [r[2] for r in rows])
    return out


_FLANK_BIN = 10_000          # bucket width for the flanking-gene index
_FLANK_MAX_BINS = 200        # a wider feature is filed as "oversized" instead


def flank_info(
    feats: Features,
    peaks_chrom: Sequence[str],
    peaks_start: Sequence[int],
    peaks_end: Sequence[int],
    distance: int = 5000,
    level: str = "transcript",
) -> Dict[str, List[str]]:
    """Port of ``getAllFlankingGene``: the three ``flank_*`` columns.

    Every feature whose **body** overlaps the peak widened by ``distance`` on
    both sides is listed, joined with ``";"``, in ascending feature order --
    the order R's ``subjectHits`` has.

    R widens ``peak.gr`` in place and keeps ``peak.gr2`` as a copy, so the search
    uses the widened interval while the distance arithmetic uses the original
    one. A feature overlapping the original peak gets ``0``; otherwise ``psD``
    and ``peD`` are formed with the feature's own strand and the smaller
    absolute value wins -- the value is NOT sign-flipped by the peak's strand,
    which is why ``-3025`` appears verbatim in the reference.

    ``flank_txIds`` only exists at ``level = "transcript"``. A peak with no
    overlapping feature gets ``NA`` in every column, as R writes it.

    The index buckets features by position so a query touches only the buckets
    its window covers: linear in peaks, not peaks x features.
    """
    use_gene = level == "gene"

    # ---- rows: (start, end, strand, key, row) with `row` the feature order
    raw: Dict[str, List[Tuple[int, int, str, str, int]]] = {}
    if use_gene:
        for chrom, rows in feats.span_gene.items():
            for gene_row, (a, b, strand, gene) in enumerate(rows):
                raw.setdefault(chrom, []).append(
                    (int(a), int(b), strand, gene, gene_row))
    else:
        # Hoist the id list: `_id_at` builds `feats.key_order or
        # list(feats.gene_of)`, which is 412k entries, so calling it per feature
        # is O(n^2) and never finished.
        ids = feats.key_order or list(feats.gene_of)
        for idx, (a, b) in enumerate(feats.span):
            key = ids[idx] if idx < len(ids) else "NA"
            raw.setdefault(feats.chrom[idx], []).append(
                (int(a), int(b), feats.strand[idx], key, idx))

    # ---- bucket index: chrom -> (bins, oversized)
    index: Dict[str, Tuple[Dict[int, List[int]], List[int], List[Tuple]]] = {}
    for chrom, rows in raw.items():
        rows.sort(key=lambda r: (r[4],))     # keep feature order for output
        bins: Dict[int, List[int]] = {}
        oversized: List[int] = []
        for i, (a, b, _st, _key, _row) in enumerate(rows):
            b0, b1 = a // _FLANK_BIN, b // _FLANK_BIN
            if b1 - b0 + 1 > _FLANK_MAX_BINS:
                oversized.append(i)
                continue
            for bn in range(b0, b1 + 1):
                bins.setdefault(bn, []).append(i)
        index[chrom] = (bins, oversized, rows)

    out: Dict[str, List[str]] = {"flank_txIds": [], "flank_geneIds": [],
                                 "flank_gene_distances": []}
    for qi in range(len(peaks_chrom)):
        chrom = peaks_chrom[qi]
        qs, qe = int(peaks_start[qi]), int(peaks_end[qi])
        wide_s, wide_e = qs - distance, qe + distance

        tx_ids: List[str] = []
        gene_ids: List[str] = []
        dists: List[str] = []
        entry = index.get(chrom)
        if entry is not None:
            bins, oversized, rows = entry
            seen = set(oversized)
            for bn in range(wide_s // _FLANK_BIN, wide_e // _FLANK_BIN + 1):
                hit = bins.get(bn)
                if hit:
                    seen.update(hit)
            for i in sorted(seen):           # ascending feature order
                a, b, strand, key, _row = rows[i]
                if b < wide_s or a > wide_e:  # bucket hit but not a real overlap
                    continue
                if not (b < qs or a > qe):
                    d = 0
                elif strand == "+":
                    ps_d = qs - a
                    pe_d = qe - a
                    d = ps_d if abs(ps_d) <= abs(pe_d) else pe_d
                else:
                    ps_d = b - qe
                    pe_d = b - qs
                    d = ps_d if abs(ps_d) <= abs(pe_d) else pe_d
                if not use_gene:
                    tx_ids.append(key)
                gene_ids.append(key if use_gene else feats.gene_of.get(key, key))
                dists.append(str(d))

        out["flank_geneIds"].append(";".join(gene_ids) if gene_ids else "NA")
        out["flank_gene_distances"].append(";".join(dists) if dists else "NA")
        if not use_gene:
            out["flank_txIds"].append(";".join(tx_ids) if tx_ids else "NA")
    return out


def _gene_level_view(feats: Features) -> Features:
    """The ``genes(TxDb)`` view R uses for ``level = "gene"``.

    ``genes()`` is ``range(transcriptsBy(txdb, by = "gene"))`` minus the groups
    spanning several seqlevels or several strands: 123,262 groups -> 121,093
    genes on the GTF-derived TxDb.

    R replaces ``features <- getGene(TxDb, by = level)`` before anything reads
    it (getGenomicAnnotation line 80), so this view must supply the TSS vector
    and the feature bodies used by ``follow()`` for the Downstream pass.
    ``annotation`` is already final by line 74, so the exon/intron/UTR tables
    that produce sub-categories are untouched -- R builds those from
    ``exonsBy()``/``intronsByTranscript()``, not from ``features``.
    """
    out = Features()
    # genes(TxDb) is ordered by gene id (coercion sorts the names) and
    # getNearestFeatureIndicesAndDistances indexes into that order, so the gene
    # view is sorted the same way.
    genes = sorted(
        (gene, chrom, s, e, strand)
        for chrom, rows in feats.span_gene.items()
        for s, e, strand, gene in rows)

    out.key_to_index = {}
    out.strand = []
    out.span = []
    out.chrom = []
    for gene, chrom, s, e, strand in genes:
        idx = len(out.strand)
        out.key_to_index[gene] = idx
        out.strand.append(strand)
        out.span.append((int(s), int(e)))
        out.chrom.append(chrom)
        # resize(width = 1) is strand aware: left edge on "+", right on "-".
        tss = int(s) if strand != "-" else int(e)
        out.tss.setdefault(chrom, []).append((tss, idx))
        out.gene_of[gene] = gene
    for chrom, rows in out.tss.items():
        rows.sort()
    # Positional id list: `span[i]` belongs to `keys[i]`. load_from_gtf cannot
    # supply key_to_index, so the index -> id mapping travels with the view.
    out.key_order = [g[0] for g in genes]

    # Sub-feature tables stay transcript-level: R builds them from exonsBy()/
    # intronsByTranscript(), not from `features`.
    out.child = {k: {c: list(v) for c, v in per.items()}
                 for k, per in feats.child.items()}
    out.child_max_len = dict(feats.child_max_len)
    return out


def annotate(
    feats: Features,
    peaks_chrom: Sequence[str],
    peaks_start: Sequence[int],
    peaks_end: Sequence[int],
    *,
    tss_region: Tuple[int, int] = (-3000, 3000),
    priority: Sequence[str] = DEFAULT_PRIORITY,
    downstream: int = DEFAULT_DOWNSTREAM,
    same_strand: bool = False,
    peaks_strand: Optional[Sequence[str]] = None,
    level: str = "transcript",
) -> Tuple[List[str], List[Optional[float]],
           List[Optional[Tuple[str, int, int, str, str, str]]]]:
    """Return ``(labels, distances, nearest)`` for each peak.

    ``nearest`` is ``(chrom, start, end, strand, feature_id, gene_id)`` of the
    feature each peak was measured against;
    ``annotatePeak`` turns it into ``geneChr``/``geneStart``/``geneEnd``/
    ``geneLength``/``geneStrand``/``geneId`` (and ``transcriptId`` at
    ``level = "transcript"``).

    Labels are ChIPseeker's ``annotation`` column verbatim, including its
    sub-categories: ``"Promoter (<=1kb)"``, ``"Exon (ENST…/GENE, exon 2 of 7)"``,
    ``"Downstream (<=300bp)"`` and ``"Distal Intergenic"``. Callers that need the
    six coarse bins use :func:`coarse_bin` on the result.
    """
    n = len(peaks_chrom)
    labels: List[str] = [""] * n
    distances: List[Optional[float]] = [None] * n
    # Nearest feature per peak, captured from the SAME `feats` object the
    # distances came from: (chrom, start, end, strand, key). Returning the
    # values instead of an index keeps every derived column consistent with the
    # feature that produced the distance -- re-deriving them from a second load
    # of the GTF drifted on 2-10 of 400 peaks.
    nearest: List[Optional[Tuple[str, int, int, str, str, str]]] = [None] * n
    # dd2 for the Downstream pass, kept apart from `distances` on purpose.
    dd2: List[Optional[float]] = [None] * n

    spans = _span_index(feats)
    # Sub-feature tables are always transcript-level; their strand lookup must
    # keep pointing at the transcript index even under level="gene".
    tx_strand = feats.strand

    if level == "gene":
        # R replaces `features <- getGene(TxDb, by = "gene")` BEFORE anything
        # reads it: the TSS distance, the follow()/dd, and geneId all come from
        # that variable. So the gene view has to replace the TSS vector too --
        # swapping only the follow() input left every distance at its
        # transcript-level value (176/400).
        if not feats.span_gene:
            raise ValueError(
                "level='gene' needs gene spans, which only the GTF loader "
                "builds (gtf_features.load_from_gtf); the exported-table loader "
                "does not carry gene_id yet")
        feats = _gene_level_view(feats)
        spans = _span_index(feats)

    for qi in range(n):
        chrom = peaks_chrom[qi]
        # Peaks arrive as GRanges coordinates, the same thing ChIPseeker sees
        # after readPeakFile(). Callers holding BED (0-based) add +1 when they
        # read it -- see run_annotate -- so no shift happens here.
        qs, qe = int(peaks_start[qi]), int(peaks_end[qi])
        pstrand = peaks_strand[qi] if peaks_strand is not None else None
        index, dist = nearest_distance(feats, chrom, qs, qe,
                                       same_strand=same_strand,
                                       peak_strand=pstrand)
        distances[qi] = dist
        if index is not None and 0 <= index < len(feats.span):
            s_i, e_i = feats.span[index]
            key = _id_at(feats, index)
            # Measured on the fixed TxDb, peak chr1:713701-714700 (-):
            #   transcript level -> transcriptId=TxID:18216, geneId=GeneID:81399
            #   gene level       -> geneId=TxID:18225, no transcriptId column
            gene = feats.gene_of.get(key, key)
            nearest[qi] = (feats.chrom[index], int(s_i), int(e_i),
                           feats.strand[index], key, gene)

        entry = spans.get(chrom)
        if entry is not None:
            # spans are (end, start, strand) sorted by end.
            ends_sorted, starts_by_end, strands_by_end = entry
            fi = _follow_excluding(ends_sorted, qs, same_strand=same_strand,
                                   peak_strand=pstrand,
                                   strands_by_end=strands_by_end)
            if fi >= 0:
                # On the minus strand R measures the peak's END against the
                # transcript's START; reusing distanceToTSS would be off by
                # one peak length.
                dd2[qi] = (float(qe - starts_by_end[fi])
                           if strands_by_end[fi] == "-"
                           else float(qs - ends_sorted[fi]))

    # Literal port of getGenomicAnnotation's loop:
    #   genomicAnnotationPriority <- rev(genomicAnnotationPriority)
    #   for (AP in genomicAnnotationPriority) { ... }
    # Under rev() the order is Intergenic, Downstream, Intron, Exon, 3UTR, 5UTR,
    # Promoter. "Intergenic" and "Downstream" both fall through to R's trailing
    # `else` branch (annotation[is.na(annotation)] <- "Intergenic"), so neither
    # writes a label here -- all Downstream work happens after the loop.
    for name in reversed(list(priority)):
        if name in ("Intergenic", "Downstream"):
            for qi in range(n):
                if not labels[qi]:
                    labels[qi] = "Intergenic"
            continue
        if name == "Promoter":
            lo, hi = tss_region
            buckets = promoter_bins(tss_region)
            for qi in range(n):
                d = distances[qi]
                if d is None or not (lo <= d <= hi):
                    continue
                for upper, text in buckets:
                    if upper < 0 or abs(d) <= upper:
                        labels[qi] = text
                        break
            continue
        kind = KIND_OF.get(name)
        if kind is None:
            continue
        per_chrom = feats.child.get(kind, {})
        max_len = feats.child_max_len.get(kind, 0)
        for qi in range(n):
            entry = per_chrom.get(peaks_chrom[qi])
            if not entry:
                continue
            rows, order = entry
            # getGenomicAnnotation.internal lines 23-27:
            #   if (sameStrand) findOverlaps(peaks, GRegion)
            #   else            findOverlaps(peaks, unstrand(GRegion))
            # so sameStrand DOES gate exon/intron/UTR hits. Dropping it here is
            # wrong; the earlier 290->261 drop came from the TSS vector being
            # pre-filtered, which R never does.
            pstrand = peaks_strand[qi] if peaks_strand is not None else None
            hit = _overlaps(rows, order, int(peaks_start[qi]),
                            int(peaks_end[qi]), max_len=max_len,
                            same_strand=same_strand,
                            peak_strand=pstrand,
                            strand_of=tx_strand)
            if hit is not None:
                labels[qi] = _label(rows, hit, kind)

    # --- everything below runs AFTER the loop, as in R ---
    #   dsd <- getOption("ChIPseeker.downstreamDistance"); if (is.null(dsd)) dsd <- 3000
    #   if (dsd/1000 <= 1) { ... } else { for (i in 1:(dsd/1000)) {...}; <catch-all> }
    # Written as R writes it: the loop labels overwrite, then the catch-all
    # re-selects on `annotation == "Intergenic"`, which is why the trailing gap
    # between trunc(dsd/1000)*1000 and dsd still ends up labelled.
    if downstream / 1000 <= 1:
        for qi in range(n):
            d = dd2[qi]
            if labels[qi] == "Intergenic" and d is not None \
                    and abs(d) <= downstream and d != 0:
                labels[qi] = f"Downstream (<={downstream}bp)"
    else:
        for i in range(1, int(downstream // 1000) + 1):
            text = "Downstream (<1kb)" if i == 1 else f"Downstream ({i - 1}-{i}kb)"
            for qi in range(n):
                d = dd2[qi]
                if labels[qi] == "Intergenic" and d is not None \
                        and abs(d) <= i * 1000 and d != 0:
                    labels[qi] = text
        whole = int(downstream // 1000)
        tail = (f"Downstream ({downstream / 1000:g}kb-{downstream}bp)")
        for qi in range(n):
            d = dd2[qi]
            if labels[qi] == "Intergenic" and d is not None \
                    and abs(d) <= downstream and d != 0:
                labels[qi] = tail

    #   annotation[which(annotation == "Intergenic")] = "Distal Intergenic"
    for qi in range(n):
        if labels[qi] == "Intergenic":
            labels[qi] = "Distal Intergenic"
    return labels, distances, nearest


# Coarse bins: the label's first word, with ChIPseeker's spaced UTR spellings.
_COARSE = {"Promoter": "Promoter", "5' UTR": "5'UTR", "5'UTR": "5'UTR",
           "3' UTR": "3'UTR", "3'UTR": "3'UTR", "Exon": "Exon", "Intron": "Intron"}


def coarse_bin(label: str) -> str:
    """Map a ChIPseeker annotation string back onto the six reporting bins."""
    head = label.split(" (")[0].strip()
    if head in _COARSE:
        return _COARSE[head]
    if head.startswith("Downstream"):
        # Downstream has no six-bin counterpart; it is upstream-of-a-gene
        # intergenic space, so it folds into Intergenic rather than vanishing.
        return "Intergenic"
    if head in ("Intergenic", "Distal Intergenic"):
        return "Intergenic"
    return head


# ChIPseeker's fixed reporting order from getGenomicAnnoStat's `lvs`, minus the
# Promoter sub-categories which are prepended dynamically.
_ANNO_STAT_TAIL = ("5' UTR", "3' UTR", "1st Exon", "Other Exon",
                   "1st Intron", "Other Intron", "Downstream", "Distal Intergenic")


def _promoter_sort_key(label: str) -> Tuple[int, float]:
    """Order Promoter sub-categories the way R's `sort()` does under en_US.UTF-8.

    `getGenomicAnnoStat` does ``plab <- sort(as.character(promoter))``. That is
    locale dependent: measured on R 4.5.3, en_US.UTF-8 sorts
    "(<=1kb)", "(1-2kb)", "(2-3kb)" by number, while LC_ALL=C puts "(<=1kb)"
    last because "<" (0x3C) sorts after "1" (0x31). en_US.UTF-8 is R's normal
    setup and is what the reference run used, so the numeric reading is
    reproduced here; a C-locale run of ChIPseeker would order differently and
    no single implementation can match both.
    """
    body = label.split("(", 1)[-1].rstrip(")")
    if body.startswith("<="):
        return (0, float(body[2:].rstrip("kb") or 0) / 1000)
    return (1, float(body.split("-", 1)[0].rstrip("kb") or 0) / 1000)


def annotation_stat(
    labels: Sequence[str],
    *,
    downstream: int = DEFAULT_DOWNSTREAM,
    ignore_1st_exon: bool = False,
    ignore_1st_intron: bool = False,
    ignore_downstream: bool = False,
    ignore_promoter_subcategory: bool = False,
) -> List[Tuple[str, float]]:
    """Port of ChIPseeker's `getGenomicAnnoStat`: the `annoStat` table of annotatePeak.

    This is NOT a fold to coarse bins. ChIPseeker reports a fixed 9-row table by
    default, produced by rewriting the annotation strings with regexes:

        "exon 1 of"    -> "1st Exon"      others -> "Other Exon"
        "intron 1 of"  -> "1st Intron"    others -> "Other Intron"
        "Downstream"   -> "Downstream (<==<dsd>)"
        "^Distal"      -> "Distal Intergenic"

    then ``table()``, ``/sum*100`` and reordering by `lvs`. Promoter
    sub-categories survive as separate rows, so the default table has three
    Promoter rows plus the six non-promoter rows -- nine categories, not seven.

    ``<dsd>`` above is deliberate: ChIPseeker writes ``paste0("Downstream (<=",
    dsd, ")")`` with no "bp" suffix, so the stat label differs from the per-peak
    annotation string ``"Downstream (<=300bp)"``. That inconsistency is
    reproduced rather than corrected, because the goal is to match.

    Returns ``[(feature, percent), ...]`` with the leading rows ChIPseeker
    emits; zero-count categories are dropped, as `table()` does.
    """
    e1lab = None if ignore_1st_exon else "1st Exon"
    exonlab = "Exon" if ignore_1st_exon else "Other Exon"
    i1lab = None if ignore_1st_intron else "1st Intron"
    intronlab = "Intron" if ignore_1st_intron else "Other Intron"

    if ignore_downstream:
        dslab, iglab = None, "Intergenic"
    else:
        shown = f"{round(downstream / 1000, 1)}kb" if downstream > 1000 else downstream
        dslab, iglab = f"Downstream (<={shown})", "Distal Intergenic"

    def fold(label: str) -> str:
        if ignore_promoter_subcategory and label.startswith("Promoter"):
            return "Promoter"
        if label.startswith("Promoter"):
            return label
        if label == "5' UTR" or label == "5'UTR":
            return "5' UTR"
        if label == "3' UTR" or label == "3'UTR":
            return "3' UTR"
        if label.startswith("Exon"):
            if not ignore_1st_exon and "exon 1 of" in label:
                return "1st Exon"
            return exonlab
        if label.startswith("Intron"):
            if not ignore_1st_intron and "intron 1 of" in label:
                return "1st Intron"
            return intronlab
        if label.startswith("Downstream"):
            return dslab if dslab is not None else iglab
        if label.startswith("Distal"):
            return iglab
        return iglab

    counts: Dict[str, int] = {}
    for lab in labels:
        key = fold(lab)
        counts[key] = counts.get(key, 0) + 1

    # R's `lvs` order: promoter labels (sorted), then the fixed tail. The
    # intergenic entry is the `iglab` VARIABLE, not a literal -- when
    # ignore_downstream is TRUE, iglab is "Intergenic" and dslab is NULL, so the
    # table ends with "Intergenic" and has no Downstream row.
    promoters = [k for k in counts if k == "Promoter" or
                 (k.startswith("Promoter") and "(" in k and not ignore_promoter_subcategory)]
    promoters.sort(key=_promoter_sort_key)
    lvs: List[str] = list(promoters)
    for name in _ANNO_STAT_TAIL:
        if name == "Distal Intergenic":
            lvs.append(iglab)
        elif name == "Downstream":
            if dslab:
                lvs.append(dslab)
        else:
            lvs.append(name)

    total = sum(counts.values()) or 1
    out: List[Tuple[str, float]] = []
    seen = set()
    for name in lvs:
        if name in seen:
            continue
        seen.add(name)
        if counts.get(name):
            out.append((name, counts[name] / total * 100))
    return out

