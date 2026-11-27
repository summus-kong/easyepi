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

def test_homer_helper_not_found_is_not_reported_as_success(tmp_path, monkeypatch):
    """HOMER exits 0 when its helper programs are missing; that is a failure.

    Reproduced on the analysis host with HOMER installed as a conda package and
    the environment not activated:

        sh: fasta2tab.pl: command not found
        sh: homer2: command not found
        sh: findKnownMotifs.pl: command not found
        sh: compareMotifs.pl: command not found
        Job finished                              <- exit status 0
        homerMotifs.all.motifs  (0 bytes)

    `easyepi motifs --backend homer` accepted that and wrote `completed` beside an
    empty motif file. Read as a result, it says "this peak set has no motifs",
    which is the opposite of what happened. The wrapper must raise instead.
    """
    from easyepi.motifs import external_tools

    fake_bin = tmp_path / "env" / "bin"
    fake_bin.mkdir(parents=True)
    exe = fake_bin / "findMotifs.pl"
    exe.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    exe.chmod(0o755)

    target = tmp_path / "t.fa"
    bg = tmp_path / "b.fa"
    target.write_text(">a\nACGTACGT\n", encoding="utf-8")
    bg.write_text(">b\nACGTACGT\n", encoding="utf-8")
    out = tmp_path / "out"

    class _Proc:
        returncode = 0
        stdout = ("Selected Options:\n\tInput file = t.fa\n"
                  "sh: homer2: command not found\n"
                  "sh: findKnownMotifs.pl: command not found\n"
                  "\tJob finished\n")
        stderr = ""

    monkeypatch.setattr(external_tools, "which_tool", lambda *a, **k: str(exe))
    monkeypatch.setattr(external_tools, "_run", lambda *a, **k: _Proc())

    (out / "homerResults").mkdir(parents=True)
    (out / "homerMotifs.all.motifs").write_text("", encoding="utf-8")

    with pytest.raises(RuntimeError, match="helper programs were not found"):
        external_tools.run_homer_motifs(str(target), str(bg), str(out))


def test_homer_run_reports_the_motifs_it_found(tmp_path, monkeypatch):
    """A real run must report how many motifs it produced, not just `available`."""
    from easyepi.motifs import external_tools

    fake_bin = tmp_path / "env" / "bin"
    fake_bin.mkdir(parents=True)
    exe = fake_bin / "findMotifs.pl"
    exe.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    exe.chmod(0o755)
    target = tmp_path / "t.fa"
    bg = tmp_path / "b.fa"
    target.write_text(">a\nACGTACGT\n", encoding="utf-8")
    bg.write_text(">b\nACGTACGT\n", encoding="utf-8")
    out = tmp_path / "out"

    class _Proc:
        returncode = 0
        stdout = "\tJob finished\n"
        stderr = ""

    monkeypatch.setattr(external_tools, "which_tool", lambda *a, **k: str(exe))
    monkeypatch.setattr(external_tools, "_run", lambda *a, **k: _Proc())

    (out / "homerResults").mkdir(parents=True)
    (out / "homerMotifs.all.motifs").write_text(
        ">MOTIF1\nline\n>MOTIF2\nline\n>MOTIF3\nline\n", encoding="utf-8")

    res = external_tools.run_homer_motifs(str(target), str(bg), str(out))
    assert res["n_denovo_motifs"] == 3, res


def _fasta(tmp_path, name, seqs):
    path = tmp_path / name
    path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)),
                    encoding="utf-8")
    return path


def test_package_and_scan_work_with_no_homer_at_all(tmp_path, monkeypatch):
    """HOMER is not a dependency of the package or of PWM scanning.

    This is the property the user asked for: not installing HOMER must not break
    the install, the import, or a scan. Only de novo discovery needs it.
    """
    import easyepi
    from easyepi.motifs import external_tools

    monkeypatch.setattr(external_tools, "which_tool", lambda *a, **k: None)
    monkeypatch.setattr(external_tools, "tool_available", lambda *a, **k: False)

    # import and the public surface still work
    assert callable(easyepi.scan_sequences)
    assert not hasattr(easyepi, "homer_like"), "homer_like must be gone"

    # a scan runs with no external program
    from easyepi.motifs import read_pwm, scan_sequences

    pwm_file = tmp_path / "m.pfm"
    pwm_file.write_text("10 0 0 0\n0 10 0 0\n0 0 0 0\n0 0 10 10\n", encoding="utf-8")
    res = scan_sequences(read_pwm(pwm_file), ["ACGT", "GGGG"])
    assert len(res["scores"]) == 2


def test_known_enrichment_needs_no_homer(tmp_path, monkeypatch):
    """Known-motif enrichment is numpy/scipy only; it must not probe HOMER."""
    from easyepi.motifs import external_tools, homer_workflow

    def _boom(*a, **k):
        raise AssertionError("known enrichment must not call HOMER")

    monkeypatch.setattr(external_tools, "run_homer_motifs", _boom)
    target = ["ACGTACGTAC"] * 4
    bg = ["GGGGCCCCGG"] * 4

    # `.motifs` is the concatenated-collection reader: a tab-separated header,
    # then one row of four probabilities per position. (`.motif` dispatches to a
    # different parser, which expects base letters on the rows.)
    motif_file = tmp_path / "k.motifs"
    motif_file.write_text(
        ">ACGT\tMOTIF1/Test\t6.0\t0\t0\tT:1.0,B:0.1\n"
        "0.9 0.1 0.0 0.0\n0.0 0.9 0.1 0.0\n0.0 0.1 0.9 0.0\n0.1 0.0 0.1 0.8\n",
        encoding="utf-8")
    motifs = homer_workflow.load_known_motifs(extra_pfm_files=[str(motif_file)])
    assert motifs, "fixture did not parse"

    rows = homer_workflow.known_motif_enrichment_homer(target, bg, motifs)
    assert isinstance(rows, list)


def test_de_novo_without_homer_reports_what_is_missing(tmp_path, monkeypatch):
    """De novo must fail clearly and point at the alternatives, not substitute."""
    from easyepi.motifs import external_tools, homer_workflow

    def _absent(*a, **k):
        raise external_tools.ToolNotAvailableError(
            "External tool 'findMotifs.pl' is not installed.")

    monkeypatch.setattr(external_tools, "run_homer_motifs", _absent)
    target = _fasta(tmp_path, "t.fa", ["ACGTACGTAC"] * 3)
    bg = _fasta(tmp_path, "b.fa", ["GGGGCCCCGG"] * 3)

    try:
        homer_workflow.run_homer_workflow(
            str(target), str(bg), str(tmp_path / "o"), k=4, backend="homer")
    except external_tools.ToolNotAvailableError as exc:
        assert "not installed" in str(exc)
    else:
        raise AssertionError("de novo discovery must require HOMER")


def test_native_backend_is_refused_with_an_explanation(tmp_path):
    """--backend native must not silently do nothing after the removal."""
    from easyepi.motifs import homer_workflow

    target = _fasta(tmp_path, "t.fa", ["ACGTACGTAC"] * 3)
    bg = _fasta(tmp_path, "b.fa", ["GGGGCCCCGG"] * 3)

    try:
        homer_workflow.run_homer_workflow(
            str(target), str(bg), str(tmp_path / "o"), k=4, backend="native")
    except ValueError as exc:
        msg = str(exc)
        assert "no built-in engine" in msg
        # The message has to be actionable, not just a refusal.
        assert "--pwm" in msg or "--db" in msg
    else:
        raise AssertionError("--backend native should be refused now")


def test_homer_failure_still_raises(tmp_path, monkeypatch):
    """An installed-but-failing HOMER must not be reported as a result."""
    from easyepi.motifs import external_tools, homer_workflow

    monkeypatch.setattr(external_tools, "run_homer_motifs",
                        lambda *a, **k: {"returncode": 3, "stderr": "boom"})
    target = _fasta(tmp_path, "t.fa", ["ACGTACGTAC"] * 3)
    bg = _fasta(tmp_path, "b.fa", ["GGGGCCCCGG"] * 3)

    try:
        homer_workflow.run_homer_workflow(
            str(target), str(bg), str(tmp_path / "o"), k=4, backend="auto")
    except RuntimeError as exc:
        assert "HOMER findMotifs.pl failed" in str(exc)
    else:
        raise AssertionError("a failing HOMER run must raise")
