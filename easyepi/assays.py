"""One definition of what an assay is, and of every spelling a user may type.

Why this module exists
----------------------
The assay name settled whether a track was shifted, which MACS3 format was used,
which QC metric set applied, and which inputs were acceptable. Each of those lived in
its own table with its own spellings, and they had already drifted:

* ``bam2bw`` keyed on ``cut&tag``/``cut&run`` while ``peaks`` keyed on
  ``cuttag``/``cutrun``, so ``--assay cuttag`` worked in one command and was rejected
  by the other.
* ``assay_metrics._alias`` repeated the same vocabulary a third time, in a fourth
  normalisation style (strip separators, map ``&`` to ``and``).
* ``qc --assay`` was passed straight to the payload without being checked at all, so
  a typo produced a report with an empty ``assay_qc`` instead of an error.

Everything that needs to know "which assay is this" goes through
:func:`canonical_assay`, and everything that needs to offer the choice to argparse
uses :func:`accepted_spellings`. A spelling added here is accepted everywhere at once.
"""

from __future__ import annotations

#: The five assays the suite models. Values are the labels used in reports and in
#: the per-assay tables below.
CANONICAL_ASSAYS: tuple[str, ...] = ("atac", "chip", "histone", "mnase", "methylation")

#: Every spelling a user may reasonably type -> canonical label.
#:
#: ``cut&tag`` and ``cut&run`` are included because ``&`` is a shell operator, which
#: forces the unquoted abbreviations ``cuttag``/``cutrun``; both forms have to work so
#: that quoting a value does not change whether it is accepted. ``tf`` means a
#: transcription-factor ChIP, which uses the ChIP metric set. ``wgbs``/``rrbs``/
#: ``bisulfite``/``methylseq`` all describe bisulfite sequencing.
ASSAY_ALIASES: dict[str, str] = {
    "atac": "atac", "atac-seq": "atac", "atacseq": "atac",
    "chip": "chip", "chip-seq": "chip", "chipseq": "chip",
    "chip-gc": "chip", "tf": "chip", "tfchip": "chip",
    "cut&tag": "chip", "cuttag": "chip", "cut-tag": "chip", "cut_tag": "chip",
    "cut&run": "chip", "cutrun": "chip", "cut-run": "chip", "cut_run": "chip",
    "histone": "histone", "histone-mark": "histone", "histonemark": "histone",
    "marks": "histone",
    "mnase": "mnase", "mnase-seq": "mnase", "mnaseseq": "mnase",
    "methylation": "methylation", "methyl": "methylation",
    "methylseq": "methylation", "wgbs": "methylation", "rrbs": "methylation",
    "bisulfite": "methylation", "bisulfite-seq": "methylation",
}


class AssayError(ValueError):
    """Raised when an assay is missing, unknown, or incompatible with its inputs."""


def _normalise(assay: object) -> str:
    """Lower-case, trim, and reduce separators so equivalent spellings coincide."""
    text = str(assay).strip().lower()
    for sep in ("_", " ", "."):
        text = text.replace(sep, "-")
    return text


def canonical_assay(assay: object) -> str:
    """Return the canonical label for any accepted spelling.

    Raises :class:`AssayError` naming the alternatives rather than a bare KeyError,
    because the caller is usually a command line and the user needs the choices.
    """
    key = ASSAY_ALIASES.get(_normalise(assay))
    if key is None:
        raise AssayError(
            f"unknown assay {assay!r}; choose one of " + ", ".join(CANONICAL_ASSAYS))
    return key


def is_known(assay: object) -> bool:
    """Whether ``assay`` is a spelling this module understands (no exception)."""
    return _normalise(assay) in ASSAY_ALIASES


def accepted_spellings(only: tuple[str, ...] | None = None) -> tuple[str, ...]:
    """Spellings to offer in argparse ``choices``, canonical labels first.

    ``only`` restricts the result to assays drawn from those canonical labels, which
    is how ``peaks`` excludes ``methylation`` (MACS3 has no bisulfite mode) without
    keeping a second copy of the vocabulary.
    """
    allowed = set(only) if only is not None else set(CANONICAL_ASSAYS)
    canonical = [name for name in CANONICAL_ASSAYS if name in allowed]
    # every alias that resolves into the allowed set, minus the canonical ones which
    # are already listed first
    aliases = sorted(
        spelling for spelling, target in ASSAY_ALIASES.items()
        if target in allowed and spelling not in canonical)
    return tuple(canonical + aliases)


def describe_spellings() -> str:
    """One-line list of canonical labels, for help text and error messages."""
    return ", ".join(CANONICAL_ASSAYS)
