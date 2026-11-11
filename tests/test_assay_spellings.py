"""The `--assay` spellings a user can actually type.

Two spellings exist for the same assay because `&` is a shell operator: `cut&tag` has
to be quoted, so every page tells the user to write `cuttag` instead. Each command has
to honour both, and `peaks` has to keep rejecting `methylation`, which has no MACS3
mode. These are one-line assertions against the tables the commands read.
"""
from __future__ import annotations

import pytest

from easyepi.bam2bw_cli import ASSAY_CHOICES as BAM2BW_CHOICES
from easyepi.bam2bw_cli import resolve_defaults
from easyepi.cli import ASSAY_CHOICES_FOR_CLI, build_parser
from easyepi.peaks_cli import ASSAY_CHOICES as PEAKS_CHOICES
from easyepi.peaks_cli import ASSAY_MACS_OPTIONS, _ASSAY_SPELLINGS


@pytest.mark.parametrize("spelling", ["cut&tag", "cuttag", "cut-tag"])
def test_bam2bw_resolves_cuttag_spellings_to_one_table(spelling):
    got = resolve_defaults(spelling)
    assert got == resolve_defaults("cut&tag")
    assert got["shift"] == "tn5" and got["bin_size"] == 10


def test_bam2bw_still_rejects_an_unknown_assay():
    with pytest.raises(ValueError, match="unknown assay"):
        resolve_defaults("not-an-assay")


def test_cli_choices_include_both_spellings():
    for spelling in ("cut&tag", "cuttag", "cut&run", "cutrun"):
        assert spelling in ASSAY_CHOICES_FOR_CLI
    # the canonical table keeps one key per assay; the aliases live in the lookup
    assert "cuttag" not in BAM2BW_CHOICES


def test_peaks_normalises_the_shell_safe_spellings():
    for spelling, key in _ASSAY_SPELLINGS.items():
        assert key in ASSAY_MACS_OPTIONS, f"{spelling} -> {key} is not a table key"


def test_peaks_does_not_advertise_methylation():
    parser = build_parser()
    peaks = None
    for action in parser._actions:
        if hasattr(action, "choices") and hasattr(action.choices, "keys"):
            peaks = action.choices.get("peaks")
            break
    assert peaks is not None
    choices = []
    for action in peaks._actions:
        if action.dest == "assay":
            choices = list(action.choices or [])
    assert "methylation" not in choices
    assert "cut&tag" in choices and "cuttag" in choices
    assert set(choices) <= set(PEAKS_CHOICES) | set(_ASSAY_SPELLINGS)
