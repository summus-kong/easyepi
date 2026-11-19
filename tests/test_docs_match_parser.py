"""The command pages must match the parser, and the nav must match the files.

`scripts/gen_docs.py` builds `docs/commands/*.md` from the live parser, which is
why the `--no-smoothing` flag had to be regenerated before it appeared. Nothing
was checking that anyone ran it, so a flag could ship documented nowhere, or a
page could describe a flag that no longer exists.

Two cheap invariants, no mkdocs needed:

1. every option the parser accepts for a subcommand appears in that subcommand's
   page, and every option the page lists exists in the parser;
2. every ``.md`` target in ``mkdocs.yml``'s nav exists on disk.

Both read the shipped files, so they fail on the same commit that introduces the
drift.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DOCS = REPO / "docs"
MKDOCS = REPO / "mkdocs.yml"

#: Subcommand -> the page gen_docs.py writes for it.
COMMAND_PAGES = {
    "dmr": "commands/dmr.md",
    "annotate": "commands/annotate.md",
    "qc": "commands/qc.md",
    "peaks": "commands/peaks.md",
    "bam2bw": "commands/bam2bw.md",
    "diff": "commands/diff.md",
}


def _parser_actions(subcommand: str):
    """The option strings the CLI actually accepts for one subcommand."""
    from easyepi.cli import build_parser

    parser = build_parser()
    for action in parser._actions:
        if getattr(action, "choices", None) and subcommand in action.choices:
            sub = action.choices[subcommand]
            return [opt for a in sub._actions for opt in a.option_strings]
    pytest.skip(f"subcommand {subcommand!r} not present in the parser")


def _documented_flags(page: Path) -> set[str]:
    """Every flag in the "Options" table, short and long.

    gen_docs.py writes one row per action and puts *all* of its option strings in
    the first cell, as ``| `-c`, `--control` | ...``, so a match anchored on a
    leading ``--`` would silently skip every short option.
    """
    text = page.read_text(encoding="utf-8")
    flags: set[str] = set()
    for line in text.splitlines():
        if not line.startswith("| `"):
            continue
        cell = line.split("|")[1]
        flags.update(re.findall(r"`(-{1,2}[A-Za-z0-9-]+)`", cell))
    return flags


@pytest.mark.parametrize("subcommand,rel", sorted(COMMAND_PAGES.items()))
def test_every_documented_flag_exists_in_the_parser(subcommand, rel):
    page = DOCS / rel
    if not page.exists():
        pytest.skip(f"{rel} not generated in this checkout")
    accepted = set(_parser_actions(subcommand))
    documented = _documented_flags(page)
    assert documented, f"{rel} documents no flags at all"
    stale = documented - accepted
    assert not stale, (
        f"{rel} documents flag(s) the parser does not accept: {sorted(stale)}. "
        "Re-run scripts/gen_docs.py.")


@pytest.mark.parametrize("subcommand,rel", sorted(COMMAND_PAGES.items()))
def test_every_parser_flag_is_documented(subcommand, rel):
    page = DOCS / rel
    if not page.exists():
        pytest.skip(f"{rel} not generated in this checkout")
    accepted = set(_parser_actions(subcommand))
    documented = _documented_flags(page)
    # argparse adds -h/--help itself; the generated pages do not list them.
    missing = accepted - documented - {"-h", "--help"}
    assert not missing, (
        f"{rel} does not document {sorted(missing)}. "
        "Re-run scripts/gen_docs.py.")


def test_the_no_smoothing_flag_is_documented():
    """The flag this round added, pinned by name.

    The two tests above are generic and would pass if *both* the flag and its
    doc row were removed. This one ties the page to the flag existing at all.
    """
    from easyepi.cli import build_parser

    parser = build_parser()
    sub = parser._actions[[i for i, a in enumerate(parser._actions)
                          if getattr(a, "choices", None)].pop()].choices["dmr"]
    assert "--no-smoothing" in {opt for a in sub._actions for opt in a.option_strings}

    page = DOCS / "commands" / "dmr.md"
    if not page.exists():
        pytest.skip("commands/dmr.md not generated in this checkout")
    assert "--no-smoothing" in _documented_flags(page)


def test_nav_targets_exist():
    """Every ``*.md`` in the nav must be a real file under docs/."""
    if not MKDOCS.exists():
        pytest.skip("no mkdocs.yml")
    text = MKDOCS.read_text(encoding="utf-8")
    targets = re.findall(r"^\s*-\s+[^:\n]+:\s*(\S+\.md)\s*$", text, re.MULTILINE)
    assert targets, "mkdocs.yml lists no .md pages"
    missing = [t for t in targets if not (DOCS / t).exists()]
    assert not missing, f"mkdocs nav points at missing pages: {missing}"


def test_the_dev_page_is_in_the_nav():
    """The replacement notes must be reachable, not just present on disk."""
    if not MKDOCS.exists():
        pytest.skip("no mkdocs.yml")
    text = MKDOCS.read_text(encoding="utf-8")
    assert "dev/DSS_REPLACEMENT_PLAN.md" in text
    assert (DOCS / "dev" / "DSS_REPLACEMENT_PLAN.md").exists()

def test_no_relative_link_escapes_the_docs_tree():
    """A ``../`` link from a docs page is broken on the published site.

    MkDocs only serves ``docs/``, so a page linking to ``../<something>`` outside
    that tree renders fine in the editor and 404s on the site. The build warns
    rather than fails, which is exactly the kind of thing that ships unnoticed.
    """
    offenders = []
    for page in DOCS.rglob("*.md"):
        text = page.read_text(encoding="utf-8")
        for target in re.findall(r"\]\((\.\./[^)\s#]+)", text):
            resolved = (page.parent / target).resolve()
            if not resolved.exists():
                offenders.append(f"{page.relative_to(REPO)} -> {target} (missing)")
            elif DOCS.resolve() not in resolved.parents and resolved != DOCS.resolve():
                offenders.append(f"{page.relative_to(REPO)} -> {target} (outside docs/)")
    assert not offenders, "links that break on the published site:\n" + "\n".join(offenders)


def test_docs_do_not_reference_untracked_files():
    """A link to a gitignored file 404s on GitHub as well as on the site."""
    out = subprocess.run(["git", "ls-files"], cwd=REPO, capture_output=True, text=True)
    tracked = set(out.stdout.split("\n"))

    offenders = []
    for page in [REPO / "README.md", *DOCS.rglob("*.md")]:
        text = page.read_text(encoding="utf-8")
        for target in re.findall(r"\]\(([^)\s#]+\.md)\)", text):
            if "://" in target or target.startswith("#"):
                continue
            resolved = (page.parent / target).resolve()
            if resolved.exists():
                rel = resolved.relative_to(REPO.resolve()).as_posix()
                if rel not in tracked:
                    offenders.append(f"{page.relative_to(REPO)} -> {target} (ignored by git)")
    assert not offenders, "\n".join(offenders)


def test_the_published_docs_state_which_machine_does_what():
    """The two roles must be in the shipped docs, not only in a local file.

    The rules file names an internal host, an account and this project's network
    situation, so it stays local; the published pages carry the substance.
    """
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    development = (DOCS / "DEVELOPMENT.md").read_text(encoding="utf-8")

    for text, where in ((readme, "README.md"), (development, "docs/DEVELOPMENT.md")):
        assert "relay" in text, f"{where} must say the local box is only a relay"
        assert "analysis host" in text, f"{where} must name the analysis host"
