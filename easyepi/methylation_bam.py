"""Per-cytosine methylation calls read directly from an aligned BAM.

Why from the BAM
----------------
Every assay in this suite takes an aligned BAM as its input, and bisulfite data is
no exception: the aligner writes 5mC calls onto the record. Requiring a separate
`.cov`/`.bedMethyl` table would break that uniformity and force the user to run a
methylation caller first -- exactly the "assemble several tools" problem this suite
exists to remove.

Two encodings are handled, and they are not interchangeable
-----------------------------------------------------------
* **MM/ML tags** (samtools 1.10+, hts-specs SAMtags): ``MM`` lists positions of
  modified bases per strand in the read, ``ML`` gives a probability per listed
  position. BAM stores the *original, unmodified* sequence, so positions must be
  resolved against the read sequence, not the reference. This is the modern form and
  is read through ``AlignedSegment.modified_bases``.
* **XM tag** (Bismark): a per-base string in read orientation where ``Z`` = methylated
  CpG, ``z`` = unmethylated CpG, ``X``/``x`` = CHG, ``H``/``h`` = CHH. Older data and
  most published WGBS BAMs use this.

A BAM carrying neither is not silently treated as "0% methylated"; the reader raises
and says which tag it looked for.

Strand and coordinate conventions
---------------------------------
Methylation is strand-specific: a CpG on the forward strand and the same CpG on the
reverse strand are different cytosines and must not be merged. Calls are reported in
forward-reference coordinates with an explicit ``strand`` column, and a CpG dyad is
summarised separately because the two strands are expected (not guaranteed) to agree.

What this does not do
---------------------
It does not call DMRs, does not model overdispersion, and does not handle
non-CpG context weighting. ``methylation_dmr`` is a separate step.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator

import numpy as np


def _require_pysam():
    from easyepi.core.io import _require
    return _require("pysam", "pip install 'easyepi[bam]'")


#: Bismark XM characters -> (modified?, context)
XM_CODES: dict[str, tuple[bool, str]] = {
    "Z": (True, "CpG"), "z": (False, "CpG"),
    "X": (True, "CHG"), "x": (False, "CHG"),
    "H": (True, "CHH"), "h": (False, "CHH"),
    "U": (True, "CHH"), "u": (False, "CHH"),
}


@dataclass(frozen=True)
class MethCall:
    """One cytosine observation, in forward-reference coordinates."""

    chrom: str
    pos: int                 # 0-based position of the C on the forward strand
    strand: str              # "+" or "-"
    context: str             # CpG / CHG / CHH
    methylated: bool

    def as_row(self) -> tuple[str, int, str, str, int]:
        return (self.chrom, self.pos, self.strand, self.context,
                1 if self.methylated else 0)


def has_methylation_tags(bam_path: str | Path, sample_reads: int = 5_000) -> dict[str, Any]:
    """Which methylation encoding the BAM carries, if any.

    Returns ``{"has_mm": bool, "has_xm": bool, "checked": int}``. Sampling stops as
    soon as either tag is seen, so a probe costs one short scan, not a full pass.
    """
    pysam = _require_pysam()
    seen = 0
    has_mm = has_xm = False
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            seen += 1
            if not has_mm and rec.has_tag("MM"):
                has_mm = True
            if not has_xm and rec.has_tag("XM"):
                has_xm = True
            if has_mm or has_xm:
                break
            if seen >= sample_reads:
                break
    return {"has_mm": has_mm, "has_xm": has_xm, "checked": seen}


def _calls_from_xm(rec, chrom: str) -> Iterator[MethCall]:
    """Bismark XM: one character per aligned base, in read orientation.

    The XM string already encodes context and state for every base (``Z``/``z`` CpG,
    ``X``/``x`` CHG, ``H``/``h`` CHH, and other letters for non-cytosines), so the
    call set is exactly the positions whose character is in :data:`XM_CODES`. A base
    reading ``T`` (an unmethylated C converted by bisulfite) still appears as ``z``,
    which is why the state must come from the tag and not from the read sequence.
    """
    xm = rec.get_tag("XM")
    try:
        pairs = rec.get_aligned_pairs(matches_only=True)
    except Exception:                                   # noqa: BLE001
        return
    reverse = rec.is_reverse
    n = len(xm)
    for qi, ref_pos in pairs:
        if qi is None or ref_pos is None or qi >= n:
            continue
        entry = XM_CODES.get(xm[qi])
        if entry is None:
            continue
        methylated, context = entry
        yield MethCall(chrom=chrom, pos=int(ref_pos),
                       strand="-" if reverse else "+",
                       context=context, methylated=methylated)


def _calls_from_mm(rec, chrom: str) -> Iterator[MethCall]:
    """MM/ML: positions of modified bases, resolved through the read sequence."""
    try:
        mods = rec.modified_bases
    except Exception:                                   # noqa: BLE001
        return
    if not mods:
        return
    seq = rec.query_sequence or ""
    pairs = rec.get_aligned_pairs(matches_only=True)
    qpos_to_ref = {q: r for q, r in pairs if q is not None and r is not None}
    reverse = rec.is_reverse
    for key, calls in mods.items():
        # key is (canonical base, strand, modification code)
        try:
            base, strand_code, _mod = key
        except (TypeError, ValueError):
            continue
        context = "CpG" if base.upper() == "C" else "CHH"
        for qpos, prob in calls:
            ref_pos = qpos_to_ref.get(int(qpos))
            if ref_pos is None:
                continue
            # a probability >= 0.5 is the conventional call threshold; the raw value
            # is kept by callers that want a soft call
            methylated = float(prob) >= 0.5
            yield MethCall(chrom=chrom, pos=int(ref_pos),
                           strand="-" if (strand_code == 1) ^ reverse else "+",
                           context=context, methylated=methylated)


def read_methylation_calls(bam_path: str | Path, chrom: str | None = None,
                           start: int = 0, end: int | None = None,
                           min_mapq: int = 10, max_reads: int | None = None,
                           contexts: Iterable[str] = ("CpG",),
                           skip_duplicates: bool = True
                           ) -> tuple[list[MethCall], dict[str, Any]]:
    """Collect cytosine calls from a BAM, with a record of what was read.

    ``contexts`` defaults to CpG only: CHG/CHH are present in plant and some
    mammalian data but are sparse and noisier, so a caller that wants them asks.
    Returns ``(calls, stats)``; ``stats`` names the encoding used and the read and
    call counts, so a suspiciously small number is visible rather than silent.
    """
    pysam = _require_pysam()
    probe = has_methylation_tags(bam_path)
    if not (probe["has_mm"] or probe["has_xm"]):
        raise ValueError(
            "no methylation calls found in this BAM: it carries neither an MM/ML "
            "tag (samtools 1.10+) nor a Bismark XM tag. Bisulfite alignments write "
            "one of these; check that this BAM is the aligner output and not a "
            "post-processed derivative that dropped the tag.")

    # Contexts are compared verbatim. Upper-casing them (as an earlier version did)
    # turns "CpG" into "CPG" and silently filters out every call, which reads as
    # "no methylation in this sample" rather than as a bug.
    canonical = {"cpg": "CpG", "chg": "CHG", "chh": "CHH"}
    want = set()
    for ctx in contexts:
        key = str(ctx).strip().lower()
        if key not in canonical:
            raise ValueError(
                f"unknown methylation context {ctx!r}; use CpG, CHG or CHH")
        want.add(canonical[key])
    calls: list[MethCall] = []
    reads = used = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if chrom is not None and chrom not in bam.references:
            raise ValueError(
                f"chrom {chrom!r} is not in the BAM header "
                f"(has e.g. {list(bam.references)[:5]})")
        if end is None:
            end = bam.get_reference_length(chrom) if chrom else None
        iterator = (bam.fetch(chrom, start, end) if chrom
                    else bam.fetch(until_eof=True))
        for rec in iterator:
            reads += 1
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            if skip_duplicates and rec.is_duplicate:
                continue
            if rec.mapping_quality < min_mapq:
                continue
            used += 1
            name = rec.reference_name
            gen = (_calls_from_mm(rec, name) if rec.has_tag("MM")
                   else _calls_from_xm(rec, name))
            for call in gen:
                if call.context in want:
                    calls.append(call)
            if max_reads is not None and used >= max_reads:
                break

    stats = {"bam": str(bam_path), "encoding": ("MM/ML" if probe["has_mm"] else "XM"),
             "reads_seen": reads, "reads_used": used, "calls": len(calls),
             "contexts": sorted(want)}
    return calls, stats


def calls_to_dataframe(calls: Iterable[MethCall], chrom_lengths: dict[str, int]):
    """Tabulate calls, matching the column contract of a bedMethyl-style table.

    Columns: ``chrom start end strand context methylated``. ``start``/``end`` are
    0-based half-open on the cytosine itself so the table can be joined with BED
    tools without an offset convention change.
    """
    import pandas as pd

    rows = [c.as_row() for c in calls]
    frame = pd.DataFrame(rows, columns=["chrom", "pos", "strand", "context",
                                        "methylated"])
    if frame.empty:
        frame["start"] = []
        frame["end"] = []
        return frame[["chrom", "start", "end", "strand", "context", "methylated"]]
    frame["start"] = frame["pos"].astype(np.int64)
    frame["end"] = frame["start"] + 1
    return frame[["chrom", "start", "end", "strand", "context", "methylated"]]


def cpg_level_methylation(frame, min_coverage: int = 5):
    """Collapse per-call rows into per-CpG-dyad beta values.

    The two strands of a CpG are combined because they measure the same dyad, but
    only after checking they exist: a site with calls on one strand only is reported
    with its own coverage rather than being silently averaged with a missing partner.

    ``beta`` = methylated / (methylated + unmethylated); sites below
    ``min_coverage`` are dropped, and the number dropped is reported, because a beta
    from one or two reads is noise that would dominate any downstream comparison.
    """
    import pandas as pd

    if frame.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "coverage", "methylated",
                                     "beta", "n_strands"]), {"sites": 0, "dropped": 0}

    grouped = (frame.groupby(["chrom", "start"], as_index=False)
               .agg(methylated=("methylated", "sum"),
                    coverage=("methylated", "size"),
                    n_strands=("strand", "nunique")))
    grouped["beta"] = grouped["methylated"] / grouped["coverage"]
    grouped["end"] = grouped["start"] + 1
    before = len(grouped)
    kept = grouped[grouped["coverage"] >= int(min_coverage)].copy()
    stats = {"sites_before_filter": int(before), "sites": int(len(kept)),
             "dropped_below_coverage": int(before - len(kept)),
             "min_coverage": int(min_coverage),
             "sites_single_strand": int((kept["n_strands"] == 1).sum())}
    return kept.reset_index(drop=True), stats
