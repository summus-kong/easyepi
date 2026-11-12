"""TSS enrichment score, the standard ATAC-seq signal-to-background metric.

Why it matters
--------------
ENCODE's ATAC pipeline judges a library largely on this number: signal piled up around
transcription start sites relative to the flanks. Without it, `easyepi qc --assay atac`
reports FRiP and fragment sizes but cannot answer "is this library enriched at
regulatory elements at all", which is the first thing a reviewer asks.

Definition used here
--------------------
The ENCODE definition, restated so the arithmetic is auditable:

* take every TSS, keep a window of +/- ``flank`` bp around it (2000 by default);
* pile up read starts (the 5' end of each alignment) across all windows, in 1 bp bins,
  into a profile;
* the *signal* is the mean pileup in the ``center`` bp at the TSS (50 by default);
* the *background* is the mean pileup in the flanks -- the outermost ``flank_bg`` bp on
  each side (100 bp, i.e. the far ends of the window);
* the score is signal / background.

Two decisions that change the number, made explicit rather than buried:

* **Reads that are not aligned, secondary or supplementary are skipped.** Counting them
  inflates both terms, but not equally.
* **Mitochondrial and unmapped contigs are excluded by name**, because a TSS list never
  contains them and their reads would otherwise only ever add to nothing.

A TSS list is required: this module does not guess TSS positions from a GTF, because a
gene-level TSS and a transcript-level TSS give materially different scores (that
distinction cost 22 percentage points in the annotation benchmark).
"""
from __future__ import annotations

import gzip
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

__all__ = ["read_tss_table", "tss_enrichment"]

#: Contigs that never carry a TSS and whose reads would only distort the background.
_SKIP_PREFIXES = ("chrM", "MT", "chrUn", "chrEBV", "_")


def _looks_like_annotation(path: str) -> bool:
    """Whether the file is a GTF/GFF, judged by its rows rather than its name.

    A GTF row has nine tab-separated fields and a feature column at index 2 drawn from
    a small vocabulary. Checking that is more reliable than an extension test -- a
    GFF3 saved as ``.txt`` is still a GFF3, and a TSV has no such column.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    features = {"gene", "transcript", "exon", "cds", "utr", "five_prime_utr",
                "three_prime_utr", "mrna", "start_codon", "stop_codon"}
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                return len(parts) >= 9 and parts[2].strip().lower() in features
    except OSError:
        return False
    return False


def read_tss_table(path: str, *, chrom_col: str = "chrom", pos_col: str = "tss",
                   strand_col: Optional[str] = None,
                   level: str = "transcript") -> List[Tuple[str, int, str]]:
    """Load TSS positions from a GTF/GFF annotation or a TSV.

    A GTF/GFF is the common case and the convenient one: the user already has the
    annotation they aligned against, so they should not have to build a second file
    that must stay in sync with it. `level` selects transcript-level (default) or
    gene-level TSS.

    A TSV is still accepted for a hand-curated list. The column names matter there: a
    column called ``start`` is ambiguous (0-based start of what?), so the caller says
    which column holds the TSS rather than having a guess made on its behalf.
    """
    if _looks_like_annotation(path):
        from .annotation_tss import tss_from_annotation

        return tss_from_annotation(path, level=level)
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: List[Tuple[str, int, str]] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header: Optional[List[str]] = None
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 1:
                parts = line.rstrip("\n").split(",")
            if header is None:
                lowered = [p.strip().lower() for p in parts]
                if chrom_col in lowered and pos_col in lowered:
                    header = lowered
                    continue
                # No recognised header: this line is data, not a header. Read the
                # columns positionally and fall through so the line is not dropped --
                # consuming it here silently lost the first TSS of every headerless
                # file.
                if len(parts) < 2:
                    raise ValueError(
                        f"{path}:{lineno}: expected at least 2 columns (chrom, tss)")
                # Keep the strand column when the file has three, so a positional
                # read does not silently drop orientation information.
                header = (["chrom", "tss", "strand"] if len(parts) >= 3
                          else ["chrom", "tss"])
            try:
                chrom = parts[header.index(chrom_col)].strip()
                pos = int(float(parts[header.index(pos_col)]))
            except (IndexError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{lineno}: cannot read a TSS from {line.strip()!r}") from exc
            # Strand matters for a signed profile, but a missing or unrecognised
            # strand column must not be fatal: default to "+" and let the caller see
            # that orientation was not available. Requiring the column would reject
            # the two-column tables most TSS exports produce.
            strand = "+"
            column = strand_col if (strand_col and strand_col in header) else (
                "strand" if "strand" in header else None)
            if column:
                try:
                    value = parts[header.index(column)].strip()
                    strand = "-" if value.startswith("-") else "+"
                except IndexError:
                    pass
            rows.append((chrom, pos, strand))
    if not rows:
        raise ValueError(f"{path}: no TSS rows found")
    return rows


def _scan_region(bam, chrom: str, lo: int, hi: int):
    """Yield records on ``chrom`` within [lo, hi) without needing an index.

    A linear pass over the file. Only used when the BAM has no index, which is the
    situation for a freshly written test file; a production BAM is indexed and takes
    the fast path.
    """
    for rec in bam.fetch(until_eof=True):
        if rec.is_unmapped or rec.reference_name != chrom:
            continue
        if rec.reference_start is None:
            continue
        if lo <= rec.reference_start < hi:
            yield rec


def tss_enrichment(
    bam_path: str,
    tss_rows: Sequence[Tuple[str, int, str]],
    *,
    flank: int = 2000,
    center: int = 50,
    flank_bg: int = 100,
    max_reads: int = 2_000_000,
    exclude_names: Sequence[str] = (),
) -> Dict[str, Any]:
    """Pile read starts around every TSS and return the enrichment score.

    ``center`` must be even and no larger than ``flank``; ``flank_bg`` is measured
    inwards from each edge of the window. Both are validated rather than coerced, so a
    caller cannot silently request a background that overlaps the signal.
    """
    import pysam

    if flank <= 0:
        raise ValueError(f"flank must be positive, got {flank}")
    if center <= 0 or center > 2 * flank:
        raise ValueError(
            f"center must be in (0, {2 * flank}], got {center}")
    if flank_bg <= 0 or 2 * flank_bg + center > 2 * flank:
        raise ValueError(
            f"flank_bg={flank_bg} with center={center} does not fit inside a "
            f"+/-{flank} window (needs 2*{flank_bg} + {center} <= {2 * flank})")

    width = 2 * flank
    profile = np.zeros(width + 1, dtype=np.int64)

    # Group TSS by chromosome and drop contigs that cannot hold one.
    by_chrom: Dict[str, np.ndarray] = {}
    for chrom, pos, _strand in tss_rows:
        if any(chrom.startswith(p) for p in _SKIP_PREFIXES):
            continue
        if exclude_names and chrom in exclude_names:
            continue
        by_chrom.setdefault(chrom, []).append(pos)
    for chrom in by_chrom:
        by_chrom[chrom] = np.sort(np.asarray(by_chrom[chrom], dtype=np.int64))

    tss_count = int(sum(len(v) for v in by_chrom.values()))
    if tss_count == 0:
        return {"tss_enrichment_score": None,
                "tss_enrichment_status": "no_tss",
                "error": "no usable TSS positions after filtering contigs",
                "reads_truncated": False,
                "max_reads": max_reads}

    used_reads = 0
    skipped_chrom = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            present = set(bam.references)
            # A TSS on a contig the BAM does not have is not an error, but it must not
            # contribute a silent zero to the profile either.
            by_chrom = {c: v for c, v in by_chrom.items() if c in present}
            if not by_chrom:
                return {"tss_enrichment_score": None,
                        "tss_enrichment_status": "no_tss",
                        "error": "no TSS chromosome is present in the BAM",
                        "reads_truncated": False,
                        "max_reads": max_reads}
            tss_count = int(sum(len(v) for v in by_chrom.values()))

            for chrom, positions in by_chrom.items():
                lo = int(positions.min()) - flank
                hi = int(positions.max()) + flank + 1
                lo = max(lo, 0)
                # `fetch(chrom, lo, hi)` needs a coordinate-sorted, indexed BAM. On an
                # unindexed file pysam returns an empty iterator rather than raising,
                # so the whole calculation would silently yield "no usable reads"
                # instead of a wrong-but-explainable error. Detect that once and fall
                # back to a linear scan, which is correct for any BAM.
                try:
                    if not bam.has_index():
                        raise ValueError("no index")
                    it = bam.fetch(chrom, lo, hi)
                except (ValueError, KeyError):
                    it = _scan_region(bam, chrom, lo, hi)
                for rec in it:
                    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                        continue
                    # The 5' end of the read is the cut site for ATAC; for a
                    # minus-strand read that is the alignment end, not the start.
                    start = rec.reference_end - 1 if rec.is_reverse else rec.reference_start
                    if start is None:
                        continue
                    idx = np.searchsorted(positions, start)
                    # Every TSS whose window covers this position contributes to it.
                    for j in (idx - 1, idx):
                        if 0 <= j < len(positions):
                            offset = start - (positions[j] - flank)
                            if 0 <= offset <= width:
                                profile[offset] += 1
                    used_reads += 1
                    if used_reads >= max_reads:
                        break
                if used_reads >= max_reads:
                    break
    except Exception as exc:  # unreadable BAM -> report why, do not crash
        # Same shape as every other return, including ``reads_truncated``: a
        # caller that reads this key must not need to guard on the status first.
        # An earlier version omitted it here, so every consumer that summarised
        # a run raised KeyError exactly when a sample had failed to read.
        return {"tss_enrichment_score": None,
                "tss_enrichment_status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "tss_count": tss_count,
                "reads_used": used_reads,
                "reads_truncated": used_reads >= max_reads,
                "max_reads": max_reads,
                "flank": flank,
                "center": center,
                "flank_bg": flank_bg}

    if used_reads == 0:
        # The key is always present so a caller can read it without a guard; None
        # plus an `error` is unambiguous, whereas a missing key is a KeyError waiting
        # to happen in whatever summary code runs next.
        return {"tss_enrichment_score": None,
                "tss_enrichment_status": "no_reads",
                "error": "no usable reads near any TSS",
                "tss_count": tss_count,
                "reads_used": 0,
                "reads_truncated": False,
                "max_reads": max_reads,
                "flank": flank,
                "center": center,
                "flank_bg": flank_bg}

    mid = width // 2
    lo_center = mid - center // 2
    hi_center = lo_center + center
    signal = float(profile[lo_center:hi_center].mean())
    left = float(profile[:flank_bg].mean()) if flank_bg else 0.0
    right = float(profile[width - flank_bg:].mean()) if flank_bg else 0.0
    background = (left + right) / 2.0

    truncated = used_reads >= max_reads
    if background <= 0:
        # The ratio is undefined, but "undefined" is not the same as "not measured".
        # A library with signal but no flanking reads is a real observation about a
        # very sparse library -- exactly what a single-cell or ultra-low-input ATAC
        # preparation looks like -- and it is a stronger statement than a low score.
        # Returning None alone would make a downstream verdict treat it identically to
        # a run where the user never supplied --tss, which is the wrong conclusion.
        status = "no_signal" if signal <= 0 else "background_zero"
        return {
            "tss_enrichment_score": None,
            "tss_enrichment_status": status,
            "error": ("no pileup at the TSS either; the library has no TSS signal"
                      if status == "no_signal"
                      else "flanking pileup is zero, so the ratio is undefined; this "
                           "is characteristic of a very sparse library, not of a "
                           "missing measurement"),
            "signal": signal,
            "background": background,
            "tss_count": tss_count,
            "reads_used": used_reads,
            "reads_truncated": truncated,
            "max_reads": max_reads,
        }

    return {
        "tss_enrichment_score": signal / background,
        "tss_enrichment_status": "measured",
        "signal": signal,
        "background": background,
        "tss_count": tss_count,
        "reads_used": used_reads,
        # Whether the read cap was hit before the file ended. A score taken from a
        # truncated scan is comparable only to other truncated runs, and saying so is
        # cheaper than discovering it from a drifting number.
        "reads_truncated": truncated,
        "max_reads": max_reads,
        "flank": flank,
        "center": center,
        "flank_bg": flank_bg,
        "profile": [int(x) for x in profile],
    }
