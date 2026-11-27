"""The assay vocabulary has one definition, and every command reads that one.

These tests exist because the same vocabulary used to live in four places, with two
different spellings of CUT&Tag, and had already produced a real bug: `--assay cuttag`
worked in `bam2bw` and was rejected by `peaks`. A single table is only useful if it
stays single, so this file asserts the consumers agree and that the old duplicates
are gone.
"""
from __future__ import annotations

import pytest

from easyepi.assays import (ASSAY_ALIASES, CANONICAL_ASSAYS, AssayError,
                            accepted_spellings, canonical_assay, is_known,
                            describe_assay, signal_class, uses_tn5_shift)
from easyepi.bam2bw_cli import ASSAY_DEFAULTS, resolve_defaults
from easyepi.cli import ASSAY_CHOICES_FOR_CLI
from easyepi.peaks_cli import ASSAY_CHOICES as PEAK_CHOICES
from easyepi.qc_assay import ASSAY_INPUTS, check_inputs, resolve_assay


# --------------------------------------------------------------------------- table
def test_the_five_canonical_assays_are_the_keys_everywhere():
    assert set(CANONICAL_ASSAYS) == set(ASSAY_DEFAULTS) == set(ASSAY_INPUTS)
    # peaks deliberately omits methylation: MACS3 has no bisulfite mode
    assert set(PEAK_CHOICES) == set(CANONICAL_ASSAYS) - {"methylation"}


def test_every_alias_resolves_to_a_canonical_label():
    for spelling, target in ASSAY_ALIASES.items():
        assert target in CANONICAL_ASSAYS, f"{spelling} -> {target}"
        assert canonical_assay(spelling) == target
        # case, padding and separator variants of the same spelling agree
        assert canonical_assay(f"  {spelling.upper()} ") == target


@pytest.mark.parametrize("spelling", ["cut&tag", "cuttag", "cut-tag", "cut_tag"])
def test_cuttag_spellings_agree_across_commands(spelling):
    """The bug that motivated this module: one spelling worked, the other did not."""
    assert canonical_assay(spelling) == "chip"
    assert resolve_defaults(spelling)["shift"] == "tn5"
    assert is_known(spelling)
    assert spelling in ASSAY_CHOICES_FOR_CLI


def test_unknown_assay_raises_with_the_choices():
    with pytest.raises(AssayError, match="choose one of"):
        canonical_assay("rna-seq")
    with pytest.raises(ValueError):
        resolve_defaults("rna-seq")


def test_accepted_spellings_can_be_restricted():
    everything = accepted_spellings()
    assert set(CANONICAL_ASSAYS) <= set(everything)
    no_methylation = accepted_spellings(("atac", "chip", "histone", "mnase"))
    assert "methylation" not in no_methylation
    assert "wgbs" not in no_methylation
    assert canonical_assay("wgbs") == "methylation"  # still valid, just not offered


# ------------------------------------------------------------------------ consumers
def test_qc_accepts_every_spelling_bam2bw_accepts():
    available = {"S1": {"bam": "a.bam"}}
    for spelling in ASSAY_CHOICES_FOR_CLI:
        assert canonical_assay(spelling) in ASSAY_INPUTS
        # the table-indexing check must not raise for anything argparse accepts
        assert isinstance(check_inputs(spelling, available), list)


def test_resolve_assay_returns_what_check_inputs_accepts():
    for spelling in ("ATAC", "atac-seq", "cutrun", "wgbs", "bisulfite"):
        resolved = resolve_assay(spelling, {"S1": {"bam": "a.bam"}})
        assert check_inputs(resolved, {"S1": {"bam": "a.bam"}}) == check_inputs(
            spelling, {"S1": {"bam": "a.bam"}})


# ------------------------------------------------------------------ no duplicates left
def test_no_module_keeps_its_own_alias_table():
    """A second copy is what drifted last time; assert the known ones are gone."""
    import easyepi.assay_metrics as metrics
    import easyepi.bam2bw_cli as b2b
    import easyepi.peaks_cli as peaks

    for module in (metrics, b2b, peaks):
        assert not hasattr(module, "_ASSAY_SPELLINGS"), module.__name__


# --------------------------------------------------------------------- peaks specifics
def test_peaks_choices_come_from_the_shared_table():
    """`peaks` must offer the shared spellings minus the assays MACS3 cannot do."""
    from easyepi.cli import build_parser

    parser = build_parser()
    subs = {}
    for action in parser._actions:
        if hasattr(action, "choices") and hasattr(action.choices, "keys"):
            subs = action.choices
            break
    peaks = subs["peaks"]
    choices = []
    for action in peaks._actions:
        if action.dest == "assay":
            choices = list(action.choices or [])

    assert choices, "peaks --assay has no choices"
    assert "methylation" not in choices
    assert "wgbs" not in choices
    # every offered spelling must resolve into the peaks table
    from easyepi.peaks_cli import ASSAY_CHOICES as PEAK_TABLE
    for spelling in choices:
        assert canonical_assay(spelling) in PEAK_TABLE, spelling

# --------------------------------------------------------------------------
# Group and signal class are independent axes
# --------------------------------------------------------------------------
def test_cuttag_and_cutrun_share_a_group_but_not_a_signal():
    """Both peak with MACS3 as ChIP; only CUT&Tag measures a Tn5 insertion.

    This is the boundary a single `canonical_assay` answer could not express, and
    the reason the Tn5 decision used to depend on whether a caller had normalised.
    """
    tag = describe_assay("CUT&Tag")
    run = describe_assay("CUT&RUN")

    assert tag.group == run.group == "chip"
    assert tag.signal == "tn5_insertion"
    assert run.signal == "mnase_footprint"
    assert tag.uses_tn5_shift is True
    assert run.uses_tn5_shift is False


def test_cutrun_never_gets_the_tn5_shift_in_any_spelling():
    """MNase has no Tn5 footprint; shifting it moves every cut site off the truth."""
    for spelling in ("cutrun", "CUT&RUN", "cut-run", "cut_run", "CUT&Run"):
        resolved = describe_assay(spelling)
        assert resolved.signal == "mnase_footprint", spelling
        assert resolved.uses_tn5_shift is False, spelling


def test_cuttag_gets_the_shift_in_every_spelling():
    for spelling in ("cuttag", "CUT&Tag", "cut-tag", "cut_tag", "Cut&Tag"):
        assert describe_assay(spelling).uses_tn5_shift is True, spelling


def test_bam2bw_defaults_agree_with_describe_assay():
    """The defaults table and the vocabulary must answer the same thing."""
    from easyepi.bam2bw_cli import resolve_defaults

    for spelling, expected in (("cuttag", "tn5"), ("CUT&RUN", "none"),
                               ("atac", "tn5"), ("chip", "none")):
        assert resolve_defaults(spelling)["shift"] == expected, spelling


def test_describe_assay_rejects_an_unknown_spelling():
    """A typo must not yield a default signal model."""
    with pytest.raises(AssayError):
        describe_assay("not-an-assay")


def test_describe_assay_reports_the_spelling_it_was_given():
    """The audit trail keeps what the user typed, in normalised form."""
    assert describe_assay("CUT&Tag").spelling == "cut&tag"


def test_atac_is_a_density_assay_and_still_shifted():
    """ATAC shifts because of Tn5 as well, but it is not a cut-site-per-fragment assay."""
    atac = describe_assay("atac")
    assert atac.group == "atac"
    assert atac.signal == "density"
    assert atac.uses_tn5_shift is False
    assert atac.is_cutsite is False


def test_methylation_is_its_own_signal_class():
    m = describe_assay("wgbs")
    assert m.group == "methylation"
    assert m.signal == "methylation"
    assert m.is_cutsite is False
