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
                            accepted_spellings, canonical_assay, is_known)
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
