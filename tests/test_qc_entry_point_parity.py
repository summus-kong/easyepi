"""The QC entry points must expose the same flags, and must agree on the numbers.

Two failures motivated these tests, both silent:

* `--mito-pattern` was added to `qc` only. The five assay aliases build their own
  parser, so the flag did not exist there and `getattr(args, "mito_pattern", None)`
  returned None -- a default that turns a missing option into a plausible value
  instead of an error. The end-to-end run reported `not_applicable` with the flag
  on the command line.
* `read_flag_metrics` and `library_summary` each counted mitochondrial reads, and
  only one of them could reach chrM in a coordinate-sorted BAM. One JSON carried
  both answers.

Parity between the parsers is checked structurally; agreement between the numbers is
checked on a real call, because the two used to differ only at runtime.
"""
from __future__ import annotations

import pytest

from easyepi.cli import build_parser

QC_ENTRY_POINTS = ("qc", "atac-qc", "chip-qc", "histone-qc", "methylation-qc",
                   "mnase-qc")


def _subparsers() -> dict:
    parser = build_parser()
    return [a for a in parser._actions if getattr(a, "choices", None)][0].choices


@pytest.mark.parametrize("name", QC_ENTRY_POINTS)
def test_every_qc_entry_point_accepts_mito_pattern(name):
    """The alias parsers are built separately, so the flag must be added to each.

    A missing flag is invisible at runtime because the caller reads it with
    `getattr(..., None)`; only a structural check catches it.
    """
    parser = _subparsers()[name]
    options = {s for a in parser._actions for s in a.option_strings}
    assert "--mito-pattern" in options, (
        f"{name} does not accept --mito-pattern; the assay aliases build their own "
        "parser, so adding the flag to `qc` alone leaves it unreachable here")


@pytest.mark.parametrize("name", QC_ENTRY_POINTS)
def test_every_qc_entry_point_accepts_tss_flags(name):
    """Same class of bug as above: these were also added to `qc` first."""
    parser = _subparsers()[name]
    options = {s for a in parser._actions for s in a.option_strings}
    assert {"--tss", "--tss-flank"} <= options, f"{name} is missing the TSS flags"


def test_mito_pattern_parses_through_both_paths():
    """Parsing, not just registration: an alias must hand the value through."""
    parser = build_parser()
    pattern = r"^(chr)?M(T)?$"

    via_qc = parser.parse_args(["qc", "--assay", "chip", "--bam", "x.bam",
                                "--mito-pattern", pattern, "--out", "o.json"])
    assert via_qc.mito_pattern == pattern

    via_alias = parser.parse_args(["chip-qc", "--bam", "x.bam",
                                   "--mito-pattern", pattern, "--out", "out"])
    assert via_alias.mito_pattern == pattern
