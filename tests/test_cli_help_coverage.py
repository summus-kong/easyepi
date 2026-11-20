"""No option may ship without help text, on any subcommand.

`nucleosome-positioning --help` once printed eleven bare option names and `--bam`
never said single-end input is refused, so the constraint only surfaced after a
run failed. Covering two commands fixed those two; this covers all of them,
because the same gap was sitting in histone-states, annotate, normalize and nine
others (28 bare options across 12 subcommands).

Two things here are semantic rather than structural, and both guard a real defect
that a coverage count would miss:

  * --bam must say single-end input is refused, not merely "paired-end";
  * the coordinate options must state the 0-based half-open convention, which is
    what callers get wrong.
"""
from __future__ import annotations

import pytest

from easyepi.cli import build_parser


def _subparsers():
    p = build_parser()
    for action in p._actions:
        if getattr(action, "choices", None):
            return action.choices
    raise AssertionError("no subcommand group found")


SUBPARSERS = _subparsers()
NAMES = sorted(SUBPARSERS)


def _options(parser):
    """(option_strings, help) for every user-facing option."""
    out = []
    for a in parser._actions:
        if not a.option_strings:
            continue
        if set(a.option_strings) <= {"-h", "--help"}:
            continue
        out.append((a.option_strings, a.help))
    return out


def test_there_is_a_real_set_of_subcommands():
    # a guard against a refactor that silently empties the group
    assert len(NAMES) >= 25, f"only {len(NAMES)} subcommands found"


@pytest.mark.parametrize("name", NAMES)
def test_every_option_has_help(name):
    missing = [opts[0] for opts, h in _options(SUBPARSERS[name])
               if not h or not str(h).strip()]
    assert not missing, f"{name}: options without help text: {missing}"


@pytest.mark.parametrize("name", NAMES)
def test_help_is_not_a_placeholder(name):
    """A one-word or restated-flag help is worse than none: it looks covered."""
    weak = []
    for opts, h in _options(SUBPARSERS[name]):
        text = str(h or "").strip()
        if not text:
            continue
        flag = opts[0].lstrip("-").replace("-", " ")
        if len(text) < 12 or text.lower().rstrip(".") == flag.lower():
            weak.append(f"{opts[0]}: {text!r}")
    assert not weak, f"{name}: help text too thin: {weak}"


def _bam_help():
    helps = {opts[0]: str(h or "") for opts, h in
             _options(SUBPARSERS["nucleosome-positioning"])}
    return helps["--bam"]


def test_bam_help_states_that_single_end_is_refused():
    text = _bam_help()
    assert "paired-end" in text
    assert "refused" in text, "--bam must say single-end input is refused"
    # the refusal message points at the working alternative; keep them aligned
    assert "easyepi nucleosome" in text


@pytest.mark.parametrize("name", ["nucleosome-positioning", "nucleosome",
                                  "histone-states"])
def test_coordinates_are_documented_as_zero_based(name):
    helps = {opts[0]: str(h or "") for opts, h in _options(SUBPARSERS[name])}
    if name == "histone-states":
        # start has a default, so it is documented but not required
        assert "0-based" in helps["--start"]
        assert "0-based" in helps["--end"]
        assert "exclusive" in helps["--end"]
        return
    assert "0-based" in helps["--start"]
    assert "0-based" in helps["--end"]
    # the half-open convention is the part people get wrong
    assert "exclusive" in helps["--end"]
