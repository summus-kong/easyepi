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

import re

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


# A help string may not be a bare format name. "GTF/GFF annotation" and
# "output directory" tell a caller only what they already read off the flag.
_FORMAT_NOUNS = {
    "bam", "bigwig", "bw", "bed", "tsv", "csv", "json", "txt", "text",
    "fasta", "fa", "fastq", "fq", "gtf", "gff", "gff3", "narrowpeak",
    "broadpeak", "chromsizes", "bigbed", "bb", "vcf", "pdf", "svg", "png",
    "jpeg", "jpg", "tiff", "html", "toml", "yaml", "yml",
}
_FILLER = {
    "a", "an", "the", "of", "for", "to", "in", "from", "into", "with",
    "input", "output", "path", "file", "files", "directory", "dir", "folder",
    "name", "optional", "default", "and", "or", "annotation",
}


def _words(text):
    return [w.strip("().,;:'\"") for w in re.split(r"[^A-Za-z0-9_]+", text.lower())
            if w.strip("().,;:'\"")]


@pytest.mark.parametrize("name", NAMES)
def test_help_is_not_only_a_format_name(name):
    """`--annotation 'GTF/GFF annotation'` is not help.

    This is the check that a length floor cannot make: the string is long enough
    to pass test_help_is_not_a_placeholder and still says nothing a caller could
    not read off the flag. It caught seven -QC commands sharing one such string.
    """
    offenders = []
    for opts, h in _options(SUBPARSERS[name]):
        text = str(h or "").strip()
        if not text:
            continue
        flag_words = set(_words(opts[0].lstrip("-")))
        rest = [w for w in _words(text)
                if w not in _FILLER and w not in flag_words]
        if rest and all(w in _FORMAT_NOUNS for w in rest):
            offenders.append(f"{opts[0]}: {text!r}")
    assert not offenders, (
        f"{name}: help names only a format, not a purpose: {offenders}"
    )


@pytest.mark.parametrize("name", NAMES)
def test_out_help_says_what_lands_there(name):
    """`--out 'output directory'` leaves the reader guessing what appears inside.

    Every --out on every command must name at least one artifact or say what the
    path is used for, because the 13 commands had drifted into 13 phrasings of
    "output directory".
    """
    for opts, h in _options(SUBPARSERS[name]):
        if not set(opts) & {"--out", "-o"}:
            continue
        text = str(h or "").strip()
        assert text, f"{name}: --out has no help text"
        words = _words(text)
        content = [w for w in words if w not in _FILLER and w not in
                   ("out", "o", "file", "path")]
        assert content, (
            f"{name}: --out help says nothing about its contents: {text!r}"
        )


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
