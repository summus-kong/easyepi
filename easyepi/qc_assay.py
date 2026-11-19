"""One QC entry point that asks the user what they measured.

Why one command instead of five
-------------------------------
``atac-qc``, ``chip-qc``, ``histone-qc`` and ``mnase-qc`` were four names for one
implementation: three of them were a single ``dispatch(...)`` call, and the only
real difference is which metric set applies (``ASSAY_METRIC_BUILDERS``).

This module answers the two questions that decide a QC: *which assay is this*, and
*can the inputs it was given support that assay's metrics*. The vocabulary itself
lives in :mod:`easyepi.assays`, so a spelling accepted here is accepted by
``bam2bw`` and ``peaks`` as well.

``--assay`` is required when the inputs cannot settle it, and is recorded in the
report so the thresholds applied are auditable. Nothing is inferred from a filename
without saying so.

The metric set per assay is the shipped table, not a new one:

===========  ==============================================
assay        metrics
===========  ==============================================
atac         fragments, FRiP, signal enrichment
chip         FRiP, signal enrichment
histone      FRiP, signal enrichment
mnase        periodicity, fragments
methylation  beta-value distribution
===========  ==============================================
"""

from __future__ import annotations

from typing import Any

from easyepi.assays import (ASSAY_ALIASES, CANONICAL_ASSAYS, AssayError,
                            canonical_assay, describe_spellings, is_known)


#: What each assay's QC can consume. Used to reject impossible requests early
#: instead of producing a report full of "not evaluated".
ASSAY_INPUTS: dict[str, tuple[str, ...]] = {
    "atac": ("bam", "bigwig", "peak"),
    "chip": ("bam", "bigwig", "peak"),
    "histone": ("bam", "bigwig", "peak"),
    "mnase": ("bam",),
    # A bisulfite BAM is the documented and preferred input: the aligner writes 5mC
    # onto the record (MM/ML or Bismark XM), so `--bam` alone is enough. A
    # `.cov`/`.bedMethyl` table stays accepted for data that arrives already called,
    # which is why both names are listed.
    "methylation": ("bam", "methylation"),
}


def resolve_assay(
    assay: str | None,
    available: dict[str, dict[str, str]],
) -> str:
    """Return the canonical assay label, or explain why it cannot be chosen.

    ``assay`` is the user's own statement of what they measured and always wins.
    Only when it is absent does the command try to infer, and the inference is
    returned as a plain label so the caller can print it; an unresolved case is an
    error that lists the choices rather than a silent default.
    """
    if assay:
        return canonical_assay(assay)

    kinds = {kind for tracks in available.values() for kind in tracks}
    if not kinds:
        raise AssayError(
            "no input tracks were resolved, and no --assay was given; "
            "pass --assay to state what you measured")
    if kinds == {"methylation"}:
        return "methylation"
    # BAM/peak/BigWig alone cannot distinguish ATAC from ChIP from histone: the files
    # look identical. Say so instead of guessing a threshold table.
    raise AssayError(
        "cannot infer the assay from "
        + ", ".join(sorted(kinds))
        + f"; pass --assay {describe_spellings()}. "
        "The assay selects which QC metrics and thresholds apply, so it is not "
        "guessed.")


def check_inputs(assay: str, available: dict[str, dict[str, str]]) -> list[str]:
    """Return human-readable notes about inputs this assay cannot use.

    Not an error: a sample may legitimately lack a BigWig, in which case the
    BigWig-only metrics are simply not computed. The note exists so the user is not
    surprised by missing columns.
    """
    key = canonical_assay(assay)
    usable = set(ASSAY_INPUTS[key])
    present = {kind for tracks in available.values() for kind in tracks}
    unusable = present - usable
    notes = []
    if unusable:
        notes.append(
            f"{key} QC does not use: {', '.join(sorted(unusable))}; "
            "these files are ignored")
    for sample, tracks in available.items():
        if "bam" in usable and "bam" not in tracks and key != "methylation":
            notes.append(f"{sample}: no BAM, read-based metrics are skipped")
        if key == "methylation" and not ({"bam", "methylation"} & set(tracks)):
            # Either input carries beta: a bisulfite BAM through its MM/ML or XM
            # tag, or a table that was already called. Only their absence stops
            # the metric, so the note names both rather than demanding a table.
            notes.append(
                f"{sample}: no BAM with MM/ML or XM tags and no methylation "
                "table, beta metrics are skipped")
    return notes


def describe() -> str:
    """Table of assay -> metrics, for --help and for an error message."""
    from easyepi.assay_metrics import ASSAY_METRIC_BUILDERS

    lines = ["assay        metrics"]
    for assay in CANONICAL_ASSAYS:
        metrics = ASSAY_METRIC_BUILDERS.get(f"{assay}-qc")
        lines.append(f"{assay:12s} {', '.join(metrics) if metrics else 'n/a'}")
    return "\n".join(lines)
