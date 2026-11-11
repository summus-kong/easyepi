"""External-tool wrapping: defaults first, verbatim passthrough last."""
from __future__ import annotations

import argparse

import pytest

from easyepi.tool_wrapper import (ArgvPassthrough, describe,
                                  validate_against_parser)


def _parser_with(*args, **kwargs) -> argparse.ArgumentParser:
    """Build a parser; returned as a factory because that is the API contract."""
    return argparse.ArgumentParser(prog="tool")


def _factory(configure) -> "callable":
    """Wrap a parser-builder so ``validate_against_parser`` can call it.

    The helper takes a *factory*, not a parser instance: a fresh parser has to be
    built per call, because ``parse_args`` on a reused one would accumulate state.
    """
    def build() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(prog="tool")
        configure(parser)
        return parser
    return build


# --------------------------------------------------------------------------
# Argv assembly
# --------------------------------------------------------------------------
def test_modelled_options_come_before_passthrough():
    """argparse takes the last occurrence, so passthrough must come last to win."""
    pw = ArgvPassthrough("macs3")
    pw.add("callpeak").add_flag("-q", 0.01).add_flag("--broad")
    pw.add_extra(["--qvalue 0.001"])
    argv = pw.argv
    assert argv[0] == "macs3"
    assert argv.index("-q") < argv.index("--qvalue")


def test_add_flag_skips_none_values():
    pw = ArgvPassthrough("tool").add_flag("--extsize", None).add_flag("--shift", 0)
    assert pw.modelled == ["--extsize", "--shift", "0"]


def test_passthrough_accepts_one_quoted_string():
    """A user writes --macs-arg="--nomodel --extsize 200" first."""
    pw = ArgvPassthrough("macs3").add_extra(["--nomodel --extsize 200"])
    assert pw.passthrough == ["--nomodel", "--extsize", "200"]


def test_passthrough_accepts_a_list_and_repeats():
    pw = ArgvPassthrough("macs3").add_extra(["--nomodel", "--extsize 200", "--bw=300"])
    assert pw.passthrough == ["--nomodel", "--extsize", "200", "--bw=300"]


def test_passthrough_splits_quoted_values_with_spaces():
    pw = ArgvPassthrough("tool").add_extra(['--label "two words"'])
    assert pw.passthrough == ["--label", "two words"]


def test_empty_passthrough_leaves_argv_unchanged():
    pw = ArgvPassthrough("tool").add("run")
    assert pw.argv == ["tool", "run"]
    assert pw.as_dict()["passthrough"] == []


def test_as_dict_records_what_was_run():
    pw = ArgvPassthrough("macs3").add("callpeak").add_extra(["--broad"])
    payload = pw.as_dict()
    assert payload["tool"] == "macs3"
    assert payload["argv"] == ["macs3", "callpeak", "--broad"]
    assert payload["passthrough"] == ["--broad"]


# --------------------------------------------------------------------------
# Validation against the tool's own parser
# --------------------------------------------------------------------------
def test_validation_accepts_known_options():
    def configure(p):
        p.add_argument("--broad", action="store_true")
        p.add_argument("--extsize", type=int)

    ns, err = validate_against_parser(["tool", "--broad", "--extsize", "200"],
                                      _factory(configure))
    assert err is None
    assert ns.broad is True and ns.extsize == 200


def test_unknown_option_raises_with_the_offending_tokens():
    def configure(p):
        p.add_argument("--broad", action="store_true")

    with pytest.raises(ValueError, match="rejected the options"):
        validate_against_parser(["tool", "--nonsense"], _factory(configure))


def test_passthrough_overrides_a_modelled_default():
    """The whole point: the user's later flag wins."""
    def configure(p):
        p.add_argument("-q", type=float, default=0.01)

    pw = ArgvPassthrough("tool").add_flag("-q", 0.01).add_extra(["-q 0.001"])
    ns, err = validate_against_parser(pw.argv, _factory(configure))
    assert err is None
    assert ns.q == 0.001


def test_leaf_parser_skips_the_subcommand_token():
    """A leaf parser does not know the sub-command name, so it is dropped."""
    def configure(p):
        p.add_argument("--broad", action="store_true")

    ns, err = validate_against_parser(["tool", "callpeak", "--broad"],
                                      _factory(configure))
    assert err is None and ns.broad is True


def test_missing_factory_reports_instead_of_failing():
    ns, err = validate_against_parser(["tool", "--x"], None)
    assert ns is None
    assert "no parser factory" in err


def test_broken_factory_is_reported_not_raised():
    def boom():
        raise RuntimeError("no console script")

    ns, err = validate_against_parser(["tool"], boom)
    assert ns is None
    assert "no console script" in err


# --------------------------------------------------------------------------
# Summary output
# --------------------------------------------------------------------------
def test_describe_lists_the_passthrough():
    pw = ArgvPassthrough("macs3").add("callpeak").add_extra(["--nomodel"])
    text = describe(pw)
    assert "macs3" in text and "--nomodel" in text


def test_describe_reports_validation_state():
    pw = ArgvPassthrough("macs3")
    pw.validated_by = "macs3-argparse"
    assert "macs3-argparse" in describe(pw)


def test_describe_reports_a_validation_note():
    pw = ArgvPassthrough("macs3")
    pw.validation_error = "console script missing"
    assert "console script missing" in describe(pw)
