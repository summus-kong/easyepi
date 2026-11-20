"""No option may ship without help text.

`nucleosome-positioning --help` used to print nine bare option names, and
`--bam` did not say single-end input is refused, so the constraint only
appeared after a run failed. Both are easy to regress while adding options.
"""
from __future__ import annotations

from easyepi.cli import build_parser


def _subparser(name):
    p = build_parser()
    for action in p._actions:
        if getattr(action, "choices", None) and name in action.choices:
            return action.choices[name]
    raise AssertionError(f"subcommand {name} not registered")


def _options(parser):
    """(option_strings, help) for every non-default action."""
    out = []
    for a in parser._actions:
        if not a.option_strings:
            continue
        if set(a.option_strings) <= {"-h", "--help"}:
            continue
        out.append((a.option_strings, a.help))
    return out


def test_every_option_has_help():
    for name in ("nucleosome-positioning", "nucleosome"):
        missing = [opts[0] for opts, h in _options(_subparser(name))
                   if not h or not str(h).strip()]
        assert not missing, f"{name}: options without help text: {missing}"


def test_bam_help_states_that_single_end_is_refused():
    helps = {opts[0]: str(h or "") for opts, h in
             _options(_subparser("nucleosome-positioning"))}
    text = helps["--bam"]
    assert "paired-end" in text
    assert "refused" in text, "--bam must say single-end input is refused"
    # the refusal message points at the working alternative; keep them aligned
    assert "easyepi nucleosome" in text


def test_coordinates_are_documented_as_zero_based():
    helps = {opts[0]: str(h or "") for opts, h in
             _options(_subparser("nucleosome-positioning"))}
    assert "0-based" in helps["--start"]
    assert "0-based" in helps["--end"]
    # the half-open convention is the part people get wrong
    assert "exclusive" in helps["--end"]
