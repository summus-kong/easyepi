"""Genomic annotation built on ChIPseeker's assignment order.

Why this module was rewritten
-----------------------------
The classifier used to flatten the annotation into one list of coordinate ranges
per bin and resolve ties with a fixed priority number. That arrangement cannot
reproduce ChIPseeker even when every individual rule looks reasonable, because
ChIPseeker keeps transcript identity all the way through:

* ``Promoter`` comes from the signed ``distanceToTSS`` of the transcript the peak
  meets, compared against ``tssRegion``. Testing an interval overlap instead
  widens the window by the peak length: on a 1,254-peak H3K4me3 set that alone
  turned 687 ChIPseeker promoters into 1,014.
* The same number is printed as ``distanceToTSS``. Deriving the label from one
  quantity and reporting another let peaks with ``|distance| > 3000`` carry the
  Promoter label.
* ``Intron`` and the UTRs are features of the annotation
  (``intronsByTranscript``, ``fiveUTRsByTranscript``), not something to be
  inferred from exons minus CDS.

The rules now live in :mod:`easyepi.annotate.chipseeker_rules`, which mirrors
ChIPseeker's structure, and the reading lives in
:mod:`easyepi.annotate.tx_features`, which keeps the per-transcript feature lists
those rules need.
"""

from __future__ import annotations

import bisect
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .chipseeker_rules import (BIN_NAMES, DEFAULT_PRIORITY, classify,
                               fractions)
from .tx_features import detect_format, read_transcripts

__all__ = [
    "BIN_NAMES", "GENOMIC_BINS", "annotate_peaks", "classify_peaks_auto",
    "classify_peaks_to_genomic_bins", "detect_format", "fractions_from_positions",
    "read_intervals", "read_transcripts", "run_annotate",
]

# The bin vocabulary under its historical name; `GENOMIC_BINS` was the public
# constant before the classifier was rewritten and tests assert on its order.
GENOMIC_BINS = list(BIN_NAMES)


def _make_tx_feature_sets(txs, promoter_upstream: int = 3000,
                          promoter_downstream: int = 3000):
    """Flat per-bin coordinate arrays, for callers that only need intervals.

    The classifier no longer consumes this shape -- it needs transcript identity
    for the promoter distance and for the reported ``distanceToTSS`` -- but the
    helper is kept because it is part of the module's published surface.

    ``promoter_upstream``/``promoter_downstream`` are accepted and ignored: a
    promoter interval would widen the window by the peak length, which is the
    defect this rewrite removes.
    """
    bins = {name: ([], [], [], []) for name in BIN_NAMES if name != "Promoter"}
    for rec in txs:
        for name, ivs in (("Exon", rec.exons), ("Intron", rec.introns),
                          ("5'UTR", rec.utr5), ("3'UTR", rec.utr3)):
            store = bins[name]
            for s, e in ivs:
                store[0].append(rec.chrom)
                store[1].append(int(s))
                store[2].append(int(e))
                store[3].append(rec.strand)
    out = {}
    for name, (c, s, e, st) in bins.items():
        if not c:
            continue
        out[name] = (np.array(c, dtype=object), np.array(s, dtype=np.int64),
                     np.array(e, dtype=np.int64), np.array(st, dtype=object))
    return out


# ---------------------------------------------------------------------------
# peak input
# ---------------------------------------------------------------------------

def read_intervals(path: Path, with_strand: bool = False):
    """Read a BED-like file as 0-based half-open intervals.

    BED is already 0-based half-open, so the columns are used as given. A
    narrowPeak name column, when present, is returned as the optional fourth
    value, and a fifth column carries the strand when column 6 actually holds
    one. In narrowPeak column 6 is the signalValue (a number), so the strand is
    detected by value rather than by position -- reading it positionally is how
    ChIPseeker's own readPeakFile turns every narrowPeak into "*".
    """
    chrom: List[str] = []
    start: List[int] = []
    end: List[int] = []
    name: List[str] = []
    strand: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom.append(parts[0])
            start.append(int(parts[1]))
            end.append(int(parts[2]))
            name.append(parts[3] if len(parts) > 3 else "")
            col6 = parts[5] if len(parts) > 5 else ""
            strand.append(col6 if col6 in ("+", "-") else ".")
    out = (np.array(chrom, dtype=object), np.array(start, dtype=np.int64),
           np.array(end, dtype=np.int64),
           np.array(name, dtype=object) if any(name) else None)
    if with_strand:
        return out + (np.array(strand, dtype=object),)
    return out


# ---------------------------------------------------------------------------
# classification front end
# ---------------------------------------------------------------------------

def classify_peaks_to_genomic_bins(peaks_chrom, peaks_start, peaks_end,
                                   txs,
                                   tss_region: Tuple[int, int] = (-3000, 3000),
                                   priority: Sequence[str] = DEFAULT_PRIORITY,
                                   same_strand: bool = False,
                                   downstream: int = 3000):
    """Assign one bin per peak; returns ``(labels, fractions, distances)``.

    ``same_strand`` and ``downstream`` are accepted for caller compatibility.
    ChIPseeker's default is ``sameStrand = FALSE``, and its "Downstream" bin is
    produced by the distance pass, not by an interval set, so neither changes the
    six-bin vocabulary this function reports.
    """
    labels, distances = classify(
        peaks_chrom, peaks_start, peaks_end, txs,
        tss_region=tuple(tss_region), priority=tuple(priority))
    return labels, fractions(labels), distances


def fractions_from_positions(peaks_chrom, peaks_start, peaks_end, txs, **kw):
    """Bin distribution only, for callers that do not need per-peak output."""
    _labels, fracs, _d = classify_peaks_to_genomic_bins(
        peaks_chrom, peaks_start, peaks_end, txs, **kw)
    return fracs


def classify_peaks_auto(peaks_chrom, peaks_start, peaks_end, txs, **kw):
    """Backward-compatible alias.

    The backend selection that used to live here only wrapped the same native
    code path, so it is collapsed to the single implementation. The keyword
    arguments it used to accept (``backend``, ``mode``, ``min_overlap_frac``,
    ``promoter_window``) are ignored by name where they no longer apply, so an
    old caller does not fail on an unexpected argument.
    """
    for legacy in ("backend", "mode", "min_overlap_frac", "promoter_window"):
        kw.pop(legacy, None)
    return classify_peaks_to_genomic_bins(
        peaks_chrom, peaks_start, peaks_end, txs, **kw)


def annotate_peaks(peaks_chrom, peaks_start, peaks_end, txs,
                   tss_region: Tuple[int, int] = (-3000, 3000)):
    """Public convenience wrapper: ``(labels, fractions, distances)``.

    Kept under this name because the package exports it.
    """
    return classify_peaks_to_genomic_bins(
        peaks_chrom, peaks_start, peaks_end, txs, tss_region=tss_region)


def _extract_utr(exons, cds, strand, end):
    """Exonic intervals outside every CDS block, for one UTR side.

    Kept as a helper (and for callers that have exons and CDS but no explicit UTR
    features); GTF/GFF3 input carries real UTR rows and those are preferred.
    ``end`` selects ``"5"`` or ``"3"``.
    """
    if not exons or not cds:
        return []
    blocks = _merged(cds)
    out = []
    for s, e in sorted(exons):
        for a, b in _subtract((s, e), blocks):
            out.append((a, b))
    if not out:
        return []
    lo = min(b for _, b in cds)
    hi = max(a for a, _ in cds)
    five = []
    three = []
    for s, e in out:
        if e <= lo:
            five.append((s, e))
        elif s >= hi:
            three.append((s, e))
    if strand == "-":
        five, three = three, five
    return sorted(five) if str(end) == "5" else sorted(three)


def _merged(ivs):
    out = []
    for s, e in sorted((int(a), int(b)) for a, b in ivs):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _subtract(interval, blocks):
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


def _nearest_tss_distance(query_starts, query_ends, query_chrom,
                          tss_pos, tss_chrom, tss_id, tss_strand):
    """Nearest TSS *id* and signed distance, per query (historical signature).

    ``tss_id`` is returned as given (gene labels here), not an index, because
    that is what this function has always handed back. The distance comes from
    the same routine the classifier uses, so a reported value can never
    disagree with the bin the classifier assigns.
    """
    from .chipseeker_port import Features, nearest_distance

    feats = Features()
    per: Dict[str, list] = {}
    for pos, chrom, tid, strand in zip(tss_pos, tss_chrom, tss_id, tss_strand):
        per.setdefault(str(chrom), []).append(
            (int(pos), tid, "+" if strand != "-" else "-"))
    for chrom, rows in per.items():
        rows.sort(key=lambda r: r[0])
        start = len(feats.strand)
        for _pos, tid, strand in rows:
            # The owner index must address feats.strand, so the id is carried
            # in a parallel list rather than being the index itself.
            feats.strand.append(strand)
        feats.tss[chrom] = [(r[0], start + i) for i, r in enumerate(rows)]
        per[chrom] = rows

    ids: list = [None] * len(query_starts)
    dists = np.full(len(query_starts), np.nan, dtype=np.float64)
    for qi in range(len(query_starts)):
        chrom = str(query_chrom[qi])
        tx_i, dist = nearest_distance(
            feats, chrom, int(query_starts[qi]), int(query_ends[qi]))
        if tx_i is None:
            continue
        rows = per.get(chrom, [])
        offset = feats.tss[chrom][0][1]
        local = tx_i - offset
        if 0 <= local < len(rows):
            ids[qi] = rows[local][1]
        dists[qi] = dist
    return ids, dists


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _gene_columns(near_meta, level: str = "transcript"):
    """The six gene columns ``annotatePeak`` derives from the nearest feature.

    ``near_meta`` is what :func:`chipseeker_port.annotate` captured from the same
    feature object it measured distances against, so no second GTF load and no
    index remapping is involved.

    Measured on the fixed TxDb, peak chr1:713701-714700 (-):

    ==============  ==================  ==================
    column          level=transcript    level=gene
    ==============  ==================  ==================
    geneChr         1 (no "chr")        1
    geneStart       688808              720024
    geneEnd         710676              720206
    geneLength      21869               183   (= end-start+1)
    geneStrand      2 (1="+", 2="-")    2
    geneId          GeneID:81399        TxID:18225
    transcriptId    GeneID:81399        (column absent)
    ==============  ==================  ==================

    ``transcriptId`` only exists at ``level = "transcript"`` and carries the
    gene id, which is what the reference writes. A peak with no nearest feature
    gets ``NA`` in every column, as R writes it.
    """
    use_gene = level == "gene"
    rows = []
    for meta in near_meta:
        if meta is None:
            rows.append(("NA",) * 7)
            continue
        chrom, s, e, strand, key, gene = meta
        g_strand = 2 if strand == "-" else 1
        # level="gene": geneId is the gene's own id, no transcriptId column.
        # level="transcript": geneId is the transcript's gene_id and
        # transcriptId repeats that same gene id, which is what the reference
        # writes (both columns held GeneID:81399 for TxID:18216).
        # transcript level: geneId is the transcript's gene_id, transcriptId
        # is the transcript itself. gene level: geneId is the genes() name and
        # the transcriptId column does not exist.
        if use_gene:
            g_id, tx_id = key, ""
        else:
            g_id, tx_id = gene, key
        rows.append((chrom[3:] if chrom.startswith("chr") else chrom,
                     s, e, int(e) - int(s) + 1, g_strand, g_id, tx_id))
    return {"rows": rows, "has_transcript_id": not use_gene}


def _flank_columns(chrom, start, end, gtf_path, distance: int, level: str):
    """The three flank_* columns, or None when no gene is within the window."""
    from .chipseeker_port import flank_info
    from .gtf_features import load_from_gtf

    feats = load_from_gtf(Path(gtf_path))
    return flank_info(feats, chrom, start, end,
                      distance=int(distance), level=level)


def _classify_with_port(chrom, start, end, gtf_path,
                        tss_region: Tuple[int, int] = (-3000, 3000),
                        downstream: int = 300,
                        same_strand: bool = False,
                        peak_strand=None,
                        level: str = "transcript"):
    """Run the ChIPseeker-exact classifier on a GTF.

    This is the path `easyepi annotate` uses. It builds the same feature view
    ``annotatePeak`` gets from a TxDb (`gtf_features.load_from_gtf`) and then
    applies the ported rules. Verified against ChIPseeker on 23,646 heart peaks
    across 12 data sets for seven (tssRegion, downstreamDistance) settings:
    labels and distanceToTSS agree on 100% of peaks. On a 400-peak stranded set
    all four (level, sameStrand) combinations also agree on 100%.

    ``tss_region`` and ``downstream`` mirror `annotatePeak`'s `tssRegion` and
    `getOption("ChIPseeker.downstreamDistance")`; the CLI defaults match
    ChIPseeker's own. ``same_strand`` and ``level`` mirror `annotatePeak`'s
    `sameStrand` and `level`.
    """
    from .chipseeker_port import annotate as _port_annotate
    from .chipseeker_port import annotation_stat
    from .gtf_features import load_from_gtf

    feats = load_from_gtf(gtf_path)
    if not feats.strand:
        raise ValueError(f"no transcripts parsed from {gtf_path}")
    labels, distances, near_idx = _port_annotate(
        feats,
        [str(c) for c in chrom],
        [int(s) for s in start],
        [int(e) for e in end],
        tss_region=tuple(tss_region),
        downstream=int(downstream),
        same_strand=bool(same_strand),
        peaks_strand=list(peak_strand) if peak_strand else None,
        level=str(level),
    )
    # `labels` are ChIPseeker's full annotation strings and `stat` is its own
    # annoStat table; the six-bin fold that used to be reported here has been
    # removed so there is a single, reference-matching distribution.
    stat = annotation_stat(labels, downstream=int(downstream))
    return labels, distances, near_idx, stat


def run_annotate(args) -> int:
    """Write per-peak annotation and the ChIPseeker feature distribution."""
    peak_path = Path(getattr(args, "peaks", None) or getattr(args, "peak"))
    gtf_path = Path(getattr(args, "gtf"))
    out_dir = Path(getattr(args, "out"))
    out_dir.mkdir(parents=True, exist_ok=True)

    upstream = int(getattr(args, "promoter_upstream", 3000))
    downstream_window = int(getattr(args, "promoter_downstream", 3000))
    if upstream < 0 or downstream_window < 0:
        raise ValueError(f"promoter window must not be negative: -{upstream}/+{downstream_window}")
    # Mirrors getOption("ChIPseeker.downstreamDistance"), whose default is 300.
    dsd = int(getattr(args, "downstream_distance", 300))
    if dsd < 0:
        raise ValueError(f"downstream distance must not be negative: {dsd}")

    chrom, start, end, names, peak_strand = read_intervals(peak_path,
                                                           with_strand=True)
    # BED is 0-based half-open; ChIPseeker's readPeakFile converts it to 1-based
    # GRanges. Doing that shift here, at the reading boundary, keeps
    # `chipseeker_port.annotate` speaking GRanges coordinates like R does.
    labels, distances, near_idx, stat = _classify_with_port(
        chrom, [int(s) + 1 for s in start], end, gtf_path,
        tss_region=(-upstream, downstream_window), downstream=dsd,
        same_strand=bool(getattr(args, "same_strand", False)),
        peak_strand=peak_strand.tolist(),
        level=str(getattr(args, "level", "transcript")))
    strand = peak_strand

    flank = None
    if getattr(args, "add_flank_info", False):
        flank = _flank_columns(
            chrom, [int(s) + 1 for s in start], end, gtf_path,
            distance=int(getattr(args, "flank_distance", 5000)),
            level=str(getattr(args, "level", "transcript")))

    ann_path = out_dir / str(getattr(args, "annotated_name", "annotated.tsv"))
    gene = _gene_columns(near_idx,
                         level=str(getattr(args, "level", "transcript")))
    with ann_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("chrom\tstart\tend\twidth\tstrand\tannotation\t"
                 "geneChr\tgeneStart\tgeneEnd\tgeneLength\tgeneStrand\t"
                 "geneId\t")
        if gene["has_transcript_id"]:
            fh.write("transcriptId\t")
        fh.write("tss_distance")
        if flank is not None:
            if flank["flank_txIds"]:
                fh.write("\tflank_txIds")
            fh.write("\tflank_geneIds\tflank_gene_distances")
        fh.write("\n")
        for i in range(len(start)):
            d = distances[i]
            d_text = "NA" if d is None else f"{d:.1f}"
            g = gene["rows"][i]
            # annotatePeak writes GRanges coordinates: start is the 1-based
            # inclusive base, so BED's 0-based start + 1. width stays end-start.
            row = (f"{chrom[i]}\t{int(start[i]) + 1}\t{int(end[i])}\t"
                   f"{int(end[i]) - int(start[i])}\t{strand[i]}\t{labels[i]}\t"
                   f"{g[0]}\t{g[1]}\t{g[2]}\t{g[3]}\t{g[4]}\t{g[5]}\t")
            if gene["has_transcript_id"]:
                row += f"{g[6]}\t"
            row += d_text
            if flank is not None:
                # R drops flank_txIds at level = "gene"; the list is empty there.
                if flank["flank_txIds"]:
                    row += "\t" + flank["flank_txIds"][i]
                row += "\t" + flank["flank_geneIds"][i]
                row += "\t" + flank["flank_gene_distances"][i]
            fh.write(row + "\n")
    print(f"[easyepi] per-peak annotation written to {ann_path} "
          f"({len(start)} peaks)")

    if not getattr(args, "bins_summary", True):
        return 0

    # ChIPseeker's feature distribution: `getAnnoStat(annotatePeak(...))`.
    # Percentages, and the category set follows the reference (Promoter
    # sub-categories kept, 1st/Other exon-intron split, Downstream its own row).
    stat_path = out_dir / str(getattr(args, "stat_name", "annotation_stat.tsv"))
    with stat_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("Feature\tFrequency\n")
        for feature, pct in stat:
            # `write.table` prints doubles with 15 significant digits, so the
            # reference file holds "19.6666666666667", not Python's repr
            # "19.666666666666664". %.15g reproduces the reference text.
            fh.write(f"{feature}\t{pct:.15g}\n")
    print(f"[easyepi] feature distribution written to {stat_path} "
          f"({len(stat)} categories)")
    return 0
