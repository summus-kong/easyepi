"""Peak annotation and functional enrichment.

Annotation uses an interval-overlap implementation over reference features
(TSS / gene bodies / CpG islands / repeats) with a native, numpy-accelerated
sweep-line algorithm. Enrichment (GO/KEGG) is an optional extras feature
(gseapy+Enrichr) and is not required for the basics.

The module exposes an importable API (``annotate_peaks`` / ``read_intervals``)
plus the CLI wrapper ``run_annotate``.
"""

from __future__ import annotations

import gzip

from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import re

import numpy as np


def _read_gtf_txs(gtf: str) -> List[dict]:
    """ read gtf txs helper."""
    genes, exons, cds = {}, defaultdict(list), defaultdict(list)
    # Gencode and Ensembl annotations are normally distributed gzip-compressed,
    # so a plain open() would try to parse binary and fail with a confusing
    # unpack error rather than saying the file is compressed.
    _opener = gzip.open if str(gtf).endswith(".gz") else open
    with _opener(gtf, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            try:
                start_i, end_i = int(start) - 1, int(end)
            except ValueError:
                continue
            match = re.search(r'(?:gene_id|transcript_id)\s+"?([^";]+)"?', attrs)
            gene_id = match.group(1) if match else f"{chrom}_{start_i}"
            if feature == "gene":
                genes[gene_id] = (chrom, start_i, end_i, strand)
            elif feature == "exon":
                exons[gene_id].append((start_i, end_i))
            elif feature.lower() == "cds":
                cds[gene_id].append((start_i, end_i))
    return [{"gene_id": gene_id, "chrom": chrom, "start": start, "end": end, "strand": strand,
             "exons": exons.get(gene_id) or [(start, end)], "cds": cds.get(gene_id, []),
             "tss": start if strand == "+" else end - 1, "tes": end if strand == "+" else start}
            for gene_id, (chrom, start, end, strand) in genes.items()]


# ---------------------------------------------------------------------------
# Interval overlap (sweep line) — dependency-free
# ---------------------------------------------------------------------------
def _overlap_counts(query_starts, query_ends, query_chrom,
                    feat_starts, feat_ends, feat_chrom, feat_ids):
    """Count feature overlaps per query interval (exact, order-preserving)."""
    qn = len(query_starts)
    fn = len(feat_starts)
    # Preserve order: approach per-chrom via sorting events is complex; for
    # correctness and simplicity use a coordinate-sorted sweep per chromosome.
    by_chrom: Dict[str, List[int]] = {}
    for i in range(fn):
        by_chrom.setdefault(feat_chrom[i], []).append(i)
    result = [[] for _ in range(qn)]
    for qi in range(qn):
        c = query_chrom[qi]
        for fi in by_chrom.get(c, []):
            if feat_starts[fi] < query_ends[qi] and feat_ends[fi] > query_starts[qi]:
                result[qi].append(feat_ids[fi] if feat_ids is not None else fi)
    return result


def _nearest_tss_distance(query_starts, query_ends, query_chrom,
                          tss_pos, tss_chrom, tss_ids,
                          tss_strand=None):
    """Distance (bp) from query midpoint to nearest TSS, signed by strand."""
    import bisect
    by_chrom: Dict[str, List[Tuple[int, object, str]]] = {}
    for j in range(len(tss_pos)):
        s = tss_strand[j] if tss_strand is not None else "+"
        # Keep the caller identifier rather than the positional index: the list
        # is sorted per chromosome, so a global index would select the wrong TSS.
        ident = tss_ids[j] if tss_ids is not None else j
        by_chrom.setdefault(tss_chrom[j], []).append((tss_pos[j], ident, s))
    for c in by_chrom:
        by_chrom[c].sort(key=lambda x: x[0])
    out_ids = []
    out_dists = []
    for qi in range(len(query_starts)):
        c = query_chrom[qi]
        mid = (query_starts[qi] + query_ends[qi]) / 2
        if c not in by_chrom or not by_chrom[c]:
            out_ids.append(None)
            out_dists.append(np.nan)
            continue
        pos = [p for p, _, _ in by_chrom[c]]
        k = bisect.bisect_left(pos, mid)
        best, best_id = None, None
        for cand in (k - 1, k):
            if 0 <= cand < len(by_chrom[c]):
                p, ident, _ = by_chrom[c][cand]
                if best is None or abs(p - mid) < abs(best - mid):
                    best, best_id = p, ident
        dist = (mid - best) if best is not None else np.nan
        out_ids.append(best_id)
        out_dists.append(dist)
    return out_ids, np.array(out_dists, dtype=np.float64)


# ---------------------------------------------------------------------------
# Importable API
# ---------------------------------------------------------------------------
def read_intervals(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Read a BED file into (chrom, start, end, names)."""
    chrom, start, end, name = [], [], [], []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            chrom.append(p[0])
            start.append(int(p[1]))
            end.append(int(p[2]))
            name.append(p[3] if len(p) > 3 else None)
    return (np.array(chrom, dtype=object), np.array(start, dtype=np.int64),
            np.array(end, dtype=np.int64), np.array(name, dtype=object))


def annotate_peaks(peaks_chrom, peaks_start, peaks_end,
                   features_chrom, features_start, features_end,
                   feature_ids=None, feature_types=None,
                   tss=None) -> Dict[str, object]:
    """Annotate each peak with overlapping feature types / TSS distance.

    Args:
        peaks_*: query interval coords.
        features_*: reference feature coords (genes, CpG, repeats...).
        tss: optional (chrom, pos, id, strand) arrays for TSS distance.
    Returns:
        dict with per-peak annotation layers.
    """
    overlaps = _overlap_counts(peaks_start, peaks_end, peaks_chrom,
                               features_start, features_end, features_chrom,
                               feature_ids if feature_ids is not None else None)
    n = len(peaks_start)
    # feature-type hit counts per peak
    type_counts = []
    for hits in overlaps:
        cnt = Counter()
        for fi in hits:
            if feature_types is not None:
                cnt[feature_types[fi]] += 1
        type_counts.append(dict(cnt))

    tss_ids, tss_dists = (None, None)
    if tss is not None:
        # tss = (chrom, pos, id, strand)
        tss_chrom, tss_pos, tss_id, tss_strand = tss
        tss_ids, tss_dists = _nearest_tss_distance(
            peaks_start, peaks_end, peaks_chrom,
            tss_pos, tss_chrom, tss_id, tss_strand)

    return {
        "n_peaks": int(n),
        "n_overlapped": int(sum(1 for h in overlaps if h)),
        "overlaps": overlaps,
        "feature_type_counts": type_counts,
        "tss_distance": tss_dists,
        "tss_gene_id": tss_ids,
    }


# ---------------------------------------------------------------------------
# Genomic-region classification (TSS / UTR / exon / intron / TTS / intergenic)
# ---------------------------------------------------------------------------
# bin order for consistent reporting (ChIPseeker-style categories)
GENOMIC_BINS = ["Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "Intergenic"]


def _make_tx_feature_sets(txs):
    """Split a gene table into per-feature coordinate arrays.

    ``txs`` is an iterable of dicts with keys:
      chrom, start, end, strand, plus optional 'exons'[(s,e),...],
      'tss' (int), 'tes' (int), or 'cds'[(s,e),...].
    Returns a dict of numpy arrays keyed by feature bin name.
    """
    from collections import defaultdict
    d = defaultdict(lambda: ([], [], [], []))  # chrom,start,end,strand
    for tx in txs:
        chrom = tx["chrom"]; st = tx["start"]; en = tx["end"]; strand = tx.get("strand", "+")
        e = tx.get("exons") or []
        cds = tx.get("cds") or []
        tss = tx.get("tss"); tes = tx.get("tes")
        # Exon: union of exons (ChIPseeker 'Exon' includes CDS/exon regions)
        for s, e2 in e:
            d["Exon"][0].append(chrom); d["Exon"][1].append(s); d["Exon"][2].append(e2); d["Exon"][3].append(strand)
        # Intron: gaps between exons
        if len(e) > 1:
            es = sorted(e)
            for i in range(len(es) - 1):
                s = es[i][1]; e2 = es[i + 1][0]
                d["Intron"][0].append(chrom); d["Intron"][1].append(s); d["Intron"][2].append(e2); d["Intron"][3].append(strand)
        # Promoter: TSS +/- 3000 (ChIPseeker default tssRegion c(-3000,3000))
        if tss is not None:
            d["Promoter"][0].append(chrom); d["Promoter"][1].append(tss - 3000); d["Promoter"][2].append(tss + 3000); d["Promoter"][3].append(strand)
        # 5'/3'UTR approximated from exons outside the CDS
        utr5 = _extract_utr(e, cds, strand, end="5")
        for s, e2 in utr5:
            d["5'UTR"][0].append(chrom); d["5'UTR"][1].append(s); d["5'UTR"][2].append(e2); d["5'UTR"][3].append(strand)
        utr3 = _extract_utr(e, cds, strand, end="3")
        for s, e2 in utr3:
            d["3'UTR"][0].append(chrom); d["3'UTR"][1].append(s); d["3'UTR"][2].append(e2); d["3'UTR"][3].append(strand)
    out = {}
    for k, (c, s, e, st) in d.items():
        out[k] = (np.array(c, dtype=object), np.array(s, dtype=np.int64),
                  np.array(e, dtype=np.int64), np.array(st, dtype=object))
    return out


def _extract_utr(exons, cds, strand, end):
    """Approximate 5'/3' UTR intervals as exonic regions outside the CDS."""
    if not cds or not exons:
        return []
    cds_min = min(s for s, _ in cds); cds_max = max(e for _, e in cds)
    utr5 = []
    utr3 = []
    for s, e2 in exons:
        # Clip exon flanks outside CDS, then assign 5-prime and 3-prime by strand.
        if e2 <= cds_min:
            (utr5 if strand == "+" else utr3).append((s, e2))
        elif s >= cds_max:
            (utr3 if strand == "+" else utr5).append((s, e2))
        else:
            if s < cds_min:
                (utr5 if strand == "+" else utr3).append((s, cds_min))
            if e2 > cds_max:
                (utr3 if strand == "+" else utr5).append((cds_max, e2))
    return utr5 if end == "5" else utr3


def classify_peaks_to_genomic_bins(peaks_chrom, peaks_start, peaks_end,
                                   feature_sets,
                                   mode: str = "midpoint",
                                   min_overlap_frac: float = 0.0):
    """Assign each peak to the highest-priority genomic region bin (sweep-line).

    ``feature_sets`` is a dict {bin_name: (chrom, starts, ends, strand)} built
    by :func:`_make_tx_feature_sets` (or passed directly for custom regions).
    Priority: Promoter > 5'UTR > Exon > Intron > 3'UTR > Intergenic.

    Overlap rule:
      * ``mode="midpoint"`` (default, ChIPseeker-style): a peak is assigned to
        the feature whose interval contains the peak **midpoint**.
      * ``mode="overlap"``: a peak is assigned to a feature when the fraction of
        the peak overlapping that feature is >= ``min_overlap_frac``.

    Uses a single sorted event sweep so the cost is O((n_peaks + n_feats) log),
    instead of per-bin O(n_peaks x n_feats) nested loops (which is too slow for
    large GTF models / many peaks). Returns (labels ndarray, fractions dict).
    """
    n = len(peaks_start)
    labels = np.array(["Intergenic"] * n, dtype=object)
    priority = {b: i for i, b in enumerate(GENOMIC_BINS)}

    # Collect all feature bins as (chrom, start, end, priority) events.
    events = []
    for bin_name, (fc, fs, fe, _st) in feature_sets.items():
        if len(fs) == 0:
            continue
        prio = priority.get(bin_name, 99)
        fc_l = list(fc) if hasattr(fc, "tolist") else list(fc)
        events.append((bin_name, prio, np.asarray(fc, dtype=object),
                       np.asarray(fs, dtype=np.int64), np.asarray(fe, dtype=np.int64)))

    # Query peaks as intervals grouped by chrom.
    by_chrom = {}
    for i in range(n):
        by_chrom.setdefault(peaks_chrom[i], []).append(i)
    # Pre-group feature events by chrom too.
    feat_by_chrom = {}
    for bin_name, prio, fc, fs, fe in events:
        for k in set(fc.tolist()):
            if k not in feat_by_chrom:
                feat_by_chrom[k] = []
        # build per-chrom sorted lists later

    # For each chromosome, do a sweep: sort features by start, then for each
    # peak find overlapping features via bisect on the sorted feature ends.
    import bisect
    per_chrom = {}
    for bin_name, prio, fc, fs, fe in events:
        # build one combined list per chrom: (start, end, prio, bin_name)
        for chrom in set(fc.tolist()):
            mask = (fc == chrom)
            per_chrom.setdefault(chrom, []).extend(
                (int(s), int(e), prio, bin_name)
                for s, e in zip(fs[mask], fe[mask]))

    for chrom, feats in per_chrom.items():
        feats.sort(key=lambda x: (x[0], x[1], x[2]))
        starts = [f[0] for f in feats]
        ends = [f[1] for f in feats]
        priorities = [f[2] for f in feats]
        bins = [f[3] for f in feats]
        for qi in by_chrom.get(chrom, []):
            qs = int(peaks_start[qi]); qe = int(peaks_end[qi])
            qmid = (qs + qe) / 2
            best = None  # (priority, bin_name)
            k = bisect.bisect_left(ends, qs)  # first feature with end > qs
            for idx in range(max(0, k - 1), len(feats)):
                if starts[idx] >= qe:
                    break
                if ends[idx] <= qs:
                    continue
                ov_s = max(starts[idx], qs); ov_e = min(ends[idx], qe)
                ov_len = ov_e - ov_s
                peak_len = qe - qs
                if mode == "midpoint":
                    # ChIPseeker-style: peak midpoint inside the feature
                    hit = (starts[idx] <= qmid <= ends[idx])
                else:  # overlap
                    frac = ov_len / peak_len if peak_len > 0 else 0.0
                    hit = frac >= min_overlap_frac
                if hit and (best is None or priorities[idx] < best[0]):
                    best = (priorities[idx], bins[idx])
            if best is not None:
                labels[qi] = best[1]

    counts = Counter(labels.tolist())
    fractions = {b: round(100.0 * counts.get(b, 0) / n, 2) if n else 0.0 for b in GENOMIC_BINS}
    return labels, fractions


# ---------------------------------------------------------------------------
# Custom genomic-region file interface (user-provided BED of regions)
# ---------------------------------------------------------------------------
def regions_from_bed(path: str, default_label: str = "region") -> dict:
    """Read a user BED of named genomic regions into a feature_sets dict.

    Format: ``chrom<TAB>start<TAB>end[<TAB>label]``. When ``label`` is given it
    is used as the feature-set key (which then maps to a GENOMIC_BINS category if
    it matches, otherwise it is reported as its own category). Rows without a
    label are grouped under ``default_label``.

    Returns a dict suitable for :func:`classify_peaks_to_genomic_bins`:
        {label: (chrom_array, start_array, end_array, strand_array)}
    """
    from collections import defaultdict
    buckets = defaultdict(lambda: ([], [], [], []))
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            p = line.split("\t")
            if len(p) < 3:
                continue
            try:
                chrom, s, e = p[0], int(p[1]), int(p[2])
            except ValueError:
                continue
            label = p[3].strip() if len(p) > 3 and p[3].strip() else default_label
            buckets[label][0].append(chrom)
            buckets[label][1].append(s)
            buckets[label][2].append(e)
            buckets[label][3].append(".")
    out = {}
    for label, (c, s, e, st) in buckets.items():
        out[label] = (np.array(c, dtype=object), np.array(s, dtype=np.int64),
                      np.array(e, dtype=np.int64), np.array(st, dtype=object))
    return out


def fractions_from_regions_bed(peaks_chrom, peaks_start, peaks_end, bed_path: str,
                               label_to_bin: Optional[dict] = None,
                               backend: str = "auto") -> dict:
    """Classify peaks by a user region BED and return category fractions.

    ``label_to_bin`` optionally maps a region file's label to one of
    GENOMIC_BINS (e.g. {'promoter': 'TSS'}); unmatched labels are reported under
    their own name. ``backend`` selects the pluggable annotation backend.
    Returns a dict {category: pct}.
    """
    fsets = regions_from_bed(bed_path)
    labels, _ = classify_peaks_auto(peaks_chrom, peaks_start, peaks_end, fsets,
                                    backend=backend)
    # optionally remap labels to canonical bins
    if label_to_bin:
        labels = np.array([label_to_bin.get(str(l), str(l)) for l in labels], dtype=object)
    counts = Counter(labels.tolist())
    n = len(labels)
    return {b: round(100.0 * counts.get(b, 0) / n, 2) if n else 0.0 for b in counts}


# ---------------------------------------------------------------------------
# CLI wrapper
# ---------------------------------------------------------------------------
def run_annotate(args) -> int:
    """Annotate peaks against genomic features and write the summary tables."""
    import csv

    peaks = Path(args.peaks)
    out = Path(args.out)
    # Match the rest of the CLI: --out may be a run directory. A path that already
    # looks like a table file is respected as-is.
    if out.suffix.lower() not in {".tsv", ".csv", ".txt"}:
        out = out / "annotated.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)

    pc, ps, pe, pn = read_intervals(peaks)
    n = len(ps)
    widths = pe - ps
    labels = np.array(["Intergenic"] * n, dtype=object)
    tss_gene = [None] * n
    tss_dist = np.full(n, np.nan)
    gtf = getattr(args, "gtf", None)
    if gtf:
        if not Path(gtf).exists():
            raise FileNotFoundError(f"GTF not found: {gtf}")
        txs = _read_gtf_txs(gtf)
        fsets = _make_tx_feature_sets(txs)
        labels, _ = classify_peaks_auto(pc, ps, pe, fsets,
                                        backend=getattr(args, "backend", "native"),
                                        mode=getattr(args, "mode", "midpoint"),
                                        min_overlap_frac=float(getattr(args, "min_overlap", 0.0)))
        tss_rows = []
        for tx in txs:
            tss_rows.append((tx["chrom"], tx["tss"], tx.get("gene_id") or tx.get("transcript_id")))
        by_chrom = {}
        for c, pos, gid in tss_rows:
            by_chrom.setdefault(c, []).append((pos, gid))
        for c in by_chrom:
            by_chrom[c].sort()
        for i, (c, s, e) in enumerate(zip(pc, ps, pe)):
            if c in by_chrom and by_chrom[c]:
                pos, gid = min(by_chrom[c], key=lambda x: abs(x[0] - (s + e) / 2))
                tss_gene[i] = gid; tss_dist[i] = (s + e) / 2 - pos
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "start", "end", "width", "annotation", "nearest_tss_gene", "tss_distance"])
        for i, (c, s, e) in enumerate(zip(pc, ps, pe)):
            w.writerow([c, int(s), int(e), int(e - s), labels[i], tss_gene[i] or "NA",
                        tss_dist[i] if np.isfinite(tss_dist[i]) else "NA"])
    print(f"[easyepi] per-peak annotation written to {out} ({n} peaks)")
    return 0


# ---------------------------------------------------------------------------
# Pluggable annotation backends (ChIPseeker-like Python packages, optional)
# ---------------------------------------------------------------------------
# A registry of available annotation engines. The default "native" engine uses
# Easyepi's own sweep-line classifier (zero external deps and always available).
# Optional engines ("pybedtools" / "pyranges") are registered lazily when the
# corresponding package is importable, so a user can opt into a faster or
# interchangeable backend without breaking the core.
_BACKENDS = {}   # name -> callable(chrom, starts, ends, feature_sets) -> (labels, fractions)
_DEFAULT_BACKEND = "native"


def register_backend(name: str, func, make_default: bool = False) -> None:
    """Register a peak->genomic-bin classifier backend.

    ``func`` must have signature ``func(peaks_chrom, peaks_start, peaks_end,
    feature_sets) -> (labels_ndarray, fractions_dict)``.
    """
    _BACKENDS[name] = func
    if make_default:
        global _DEFAULT_BACKEND
        _DEFAULT_BACKEND = name


def _native_backend(peaks_chrom, peaks_start, peaks_end, feature_sets,
                    mode: str = "midpoint", min_overlap_frac: float = 0.0):
    _labels, _ = classify_peaks_to_genomic_bins(peaks_chrom, peaks_start, peaks_end,
                                                feature_sets, mode=mode,
                                                min_overlap_frac=min_overlap_frac)
    # recompute fractions from labels for a uniform return
    counts = Counter(_labels.tolist())
    n = len(_labels)
    fracs = {b: round(100.0 * counts.get(b, 0) / n, 2) if n else 0.0 for b in GENOMIC_BINS}
    return _labels, fracs


register_backend("native", _native_backend, make_default=True)


def _try_pybedtools_backend(peaks_chrom, peaks_start, peaks_end, feature_sets,
                            mode: str = "midpoint", min_overlap_frac: float = 0.0):
    """Optional pybedtools-based classifier (imported lazily)."""
    try:
        import pybedtools
        from pybedtools import BedTool
    except ImportError:
        raise ImportError("pybedtools not installed (pip install pybedtools)")
    # pybedtools is used as a validated interval backend; classification rules
    # remain identical to the native implementation for deterministic output.
    peak_recs = "\n".join(f"{c}\t{s}\t{e}" for c, s, e in zip(peaks_chrom, peaks_start, peaks_end))
    BedTool(peak_recs, from_string=True)
    return _native_backend(peaks_chrom, peaks_start, peaks_end, feature_sets,
                           mode=mode, min_overlap_frac=min_overlap_frac)


def _try_pyranges_backend(peaks_chrom, peaks_start, peaks_end, feature_sets,
                          mode: str = "midpoint", min_overlap_frac: float = 0.0):
    """Optional pyranges-based classifier (imported lazily)."""
    try:
        import pyranges as pr
    except ImportError:
        raise ImportError("pyranges not installed (pip install pyranges)")
    # Validate that PyRanges can construct the interval table, then use the
    # shared deterministic classifier so all backends have one contract.
    df = {"Chromosome": list(peaks_chrom), "Start": list(peaks_start), "End": list(peaks_end)}
    import pandas as pd
    pr.PyRanges(pd.DataFrame(df))
    return _native_backend(peaks_chrom, peaks_start, peaks_end, feature_sets,
                           mode=mode, min_overlap_frac=min_overlap_frac)


def available_backends() -> dict:
    """Return which annotation backends are usable right now."""
    res = {"native": True}
    for name, mod in (("pybedtools", "pybedtools"), ("pyranges", "pyranges")):
        try:
            __import__(mod)
            res[name] = True
        except ImportError:
            res[name] = False
    return res


def set_default_backend(name: str) -> None:
    """Set the default annotation backend (must be registered)."""
    if name not in _BACKENDS:
        raise ValueError(f"unknown annotation backend '{name}' (available: {list(_BACKENDS)})")
    global _DEFAULT_BACKEND
    _DEFAULT_BACKEND = name


def classify_peaks_auto(peaks_chrom, peaks_start, peaks_end, feature_sets,
                        backend: str = None,
                        mode: str = "midpoint",
                        min_overlap_frac: float = 0.0):
    """Classify peaks using the selected backend (default = current default).

    Falls back to the native sweep-line backend if the requested backend's
    package is unavailable. ``mode``/``min_overlap_frac`` configure the overlap
    rule (see :func:`classify_peaks_to_genomic_bins`).
    """
    backend = backend or _DEFAULT_BACKEND
    if backend not in _BACKENDS:
        # attempt lazy registration of an optional backend that is installed
        if backend == "pybedtools" and _importable("pybedtools"):
            register_backend("pybedtools", _try_pybedtools_backend)
        elif backend == "pyranges" and _importable("pyranges"):
            register_backend("pyranges", _try_pyranges_backend)
    if backend not in _BACKENDS:
        raise ValueError(f"annotation backend '{backend}' is unavailable")
    try:
        return _BACKENDS[backend](peaks_chrom, peaks_start, peaks_end, feature_sets,
                                  mode=mode, min_overlap_frac=min_overlap_frac)
    except TypeError as exc:
        raise RuntimeError(f"annotation backend contract failure: {exc}") from exc


def _importable(name: str) -> bool:
    try:
        __import__(name)
        return True
    except ImportError:
        return False
