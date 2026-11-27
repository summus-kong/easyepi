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

from dataclasses import dataclass

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


#: What each canonical assay *measures*, as distinct from which group it belongs
#: to. Two assays can share a grouping and differ here -- which is exactly the
#: CUT&Tag/CUT&RUN case that a single `canonical_assay` answer could not express.
#:
#: ``tn5_insertion``  strand-resolved transposase insertion; needs the +4/-5 offset
#: ``mnase_footprint`` MNase-protected fragment; no Tn5 offset exists to apply
#: ``fragment_depth``  read-body depth over a ChIP fragment
#: ``density``         accessibility or occupancy density, no strand offset
#: ``methylation``     per-cytosine calls with coverage
SIGNAL_CLASSES: tuple[str, ...] = (
    "tn5_insertion", "mnase_footprint", "fragment_depth", "density", "methylation",
)

#: Canonical assay -> the signal it measures.
#:
#: `cuttag` and `cutrun` share the `chip` *group* (MACS3 calls peaks on them alike)
#: but must not share a signal model, so they carry distinct entries here. This is
#: the table that the old spelling-matching in `bam2bw_cli` was standing in for.
ASSAY_SIGNAL: dict[str, str] = {
    "atac": "density",
    "chip": "fragment_depth",
    "histone": "fragment_depth",
    "mnase": "mnase_footprint",
    "methylation": "methylation",
    # Filled per spelling below, because `canonical_assay` maps both to `chip`.
}

#: Spellings that must NOT inherit their group's signal model.
#:
#: This is deliberately keyed by normalised spelling rather than by canonical
#: label: two different assays collide on `chip`, and the collision is the whole
#: reason this table exists.
SIGNAL_OVERRIDES: dict[str, str] = {
    "cut&tag": "tn5_insertion",
    "cuttag": "tn5_insertion",
    "cut-tag": "tn5_insertion",
    "cut_tag": "tn5_insertion",
    "cut&run": "mnase_footprint",
    "cutrun": "mnase_footprint",
    "cut-run": "mnase_footprint",
    "cut_run": "mnase_footprint",
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


def signal_class(assay: object) -> str:
    """What ``assay`` measures, independent of which group it peaks with.

    Raises for an unknown spelling, like :func:`canonical_assay`, so a caller cannot
    get a silently wrong signal model from a typo.
    """
    key = _normalise(assay)
    if key in SIGNAL_OVERRIDES:
        return SIGNAL_OVERRIDES[key]
    canonical = ASSAY_ALIASES.get(key)
    if canonical is None:
        raise AssayError(
            f"unknown assay {assay!r}; choose one of " + ", ".join(CANONICAL_ASSAYS))
    return ASSAY_SIGNAL[canonical]


def uses_tn5_shift(assay: object) -> bool:
    """Whether the +4/-5 Tn5 offset applies.

    Answered from the signal class, so it survives normalisation:
    ``uses_tn5_shift("cuttag")`` and ``uses_tn5_shift(canonical_assay("cuttag"))``
    differ only if the signal model is read from the group, which it no longer is.
    """
    return signal_class(assay) == "tn5_insertion"


def is_cutsite_assay(assay: object) -> bool:
    """Whether the assay resolves a cut site, which wants a fine bin.

    True for Tn5 insertion and for MNase footprints: both are positional in a way
    read-body depth is not.
    """
    return signal_class(assay) in {"tn5_insertion", "mnase_footprint"}


@dataclass(frozen=True)
class Assay:
    """Everything a processing decision needs to know about an assay.

    Three fields, three questions, answered once:

    ``spelling``  what the user typed, normalised -- keeps the audit trail honest
    ``group``     which analysis family it belongs to, for shared rules (MACS3
                  format, QC metric set)
    ``signal``    what is actually measured, which the group does not determine
    """

    spelling: str
    group: str
    signal: str

    @property
    def uses_tn5_shift(self) -> bool:
        return self.signal == "tn5_insertion"

    @property
    def is_cutsite(self) -> bool:
        return self.signal in {"tn5_insertion", "mnase_footprint"}


def describe_assay(assay: object) -> Assay:
    """Resolve any accepted spelling into spelling, group and signal class.

    This is the entry point for anything that decides how to *process* data. Use
    :func:`canonical_assay` only when nothing but the group matters.

    Named `describe_` rather than `resolve_` because `easyepi.qc_assay` already
    exports a `resolve_assay(assay, available)` that answers a different
    question (which metric set applies to a sample set); two functions with one
    name and different arity is how a silent mis-import happens.

    Raises :class:`AssayError` for an unknown spelling, so a typo cannot silently
    yield a default signal model.
    """
    key = _normalise(assay)
    group = ASSAY_ALIASES.get(key)
    if group is None:
        raise AssayError(
            f"unknown assay {assay!r}; choose one of " + ", ".join(CANONICAL_ASSAYS))
    signal = SIGNAL_OVERRIDES.get(key, ASSAY_SIGNAL[group])
    return Assay(spelling=key, group=group, signal=signal)


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
