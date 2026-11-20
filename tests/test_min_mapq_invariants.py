"""`dmr --min-mapq` defaults to 10 on purpose; that must not drift to 30.

Every other command filters at MAPQ 30. `dmr` reads bisulfite alignments, where
aligners assign lower MAPQ (a bisulfite-converted read matches fewer reference
k-mers, so BWA-meth/bwa-meth-style mappers and Bismark both emit many 20-30
alignments that are perfectly usable). Raising the floor to 30 would silently
drop real coverage at CpG sites and change every DMR call.

A help string alone cannot hold that: whoever changes `default=10` to `default=30`
will also happily change the text. These tests fail in that case, and they also
require the reason to stay written down.
"""
from __future__ import annotations

import re

from easyepi.cli import build_parser


def _subparsers():
    p = build_parser()
    for action in p._actions:
        if getattr(action, "choices", None):
            return action.choices
    raise AssertionError("no subcommand group found")


SUBPARSERS = _subparsers()


def _min_mapq(name):
    for a in SUBPARSERS[name]._actions:
        if "--min-mapq" in a.option_strings:
            return a
    raise AssertionError(f"{name} has no --min-mapq option")


# the commands that read ordinary DNA alignments, where MAPQ 30 is the convention
DNA_COMMANDS = ("bam2bw", "differential-peaks", "nucleosome-positioning")


def test_dmr_defaults_lower_than_the_dna_commands():
    dmr = _min_mapq("dmr").default
    higher = {n: _min_mapq(n).default for n in DNA_COMMANDS
              if n in SUBPARSERS and _min_mapq(n).default > dmr}
    assert higher, (
        f"dmr --min-mapq is {dmr}; it is meant to stay below the DNA commands "
        f"({ {n: _min_mapq(n).default for n in DNA_COMMANDS if n in SUBPARSERS} }) "
        "because bisulfite aligners report lower MAPQ. If raising it is intended, "
        "update this test and the reason in the help text together."
    )


def test_dmr_keeps_the_bisulfite_reason_in_help():
    text = str(_min_mapq("dmr").help or "").lower()
    assert "bisulfite" in text, (
        "the dmr --min-mapq help must keep the bisulfite explanation; it is the "
        "only place the difference from the other commands is recorded"
    )
    # the stated number must agree with the actual default
    stated = re.search(r"default (\d+)", text)
    assert stated, "dmr --min-mapq help must state its default"
    assert int(stated.group(1)) == _min_mapq("dmr").default, (
        f"help says default {stated.group(1)} but the parser uses "
        f"{_min_mapq('dmr').default}"
    )


def test_dna_commands_agree_on_30():
    """The DNA-alignment commands share one convention; keep them aligned."""
    values = {}
    for name in DNA_COMMANDS:
        if name in SUBPARSERS:
            values[name] = _min_mapq(name).default
    assert set(values.values()) == {30}, (
        f"the DNA-alignment commands should all default to 30, got {values}"
    )


def test_every_stated_default_matches_the_parser():
    """A help string that names a default must name the real one.

    This is the general form of the dmr guard: it catches any option whose text
    and behaviour drift apart, not just --min-mapq.
    """
    wrong = []
    for name, sp in SUBPARSERS.items():
        for a in sp._actions:
            if not a.option_strings:
                continue
            text = str(a.help or "")
            # 1e-4 and 0.0001 must both parse; an exponent-blind pattern reads
            # "default 1e-4" as 1.0 and reports a false mismatch.
            m = re.search(r"default (\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", text)
            if not m:
                continue
            stated, actual = float(m.group(1)), a.default
            if isinstance(actual, (int, float)) and float(actual) != stated:
                wrong.append(f"{name} {a.option_strings[0]}: "
                             f"help says {stated}, parser uses {actual}")
    assert not wrong, "help text disagrees with the parser:\n  " + "\n  ".join(wrong)
