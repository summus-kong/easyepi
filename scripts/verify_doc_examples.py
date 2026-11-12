"""Parse every documented `easyepi ...` example with the live argparse tree.

A doc example that argparse rejects is a defect a reader hits immediately, and it is
invisible to `mkdocs build --strict`, which only checks links. Shell continuation lines
are joined first, so a multi-line invocation is checked as the one command a reader
would paste. Examples containing shell expansions (`${s}`, pipes, redirects) are skipped
because they are not a single argv.

Run from anywhere: ``python scripts/verify_doc_examples.py``
"""
from __future__ import annotations

import contextlib
import io
import re
import shlex
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from easyepi.cli import build_parser  # noqa: E402

SHELLY = ("${", "$(", "|", ">", "<", "&>", "&&")


def subparsers(parser) -> dict:
    for action in parser._actions:
        if hasattr(action, "choices") and hasattr(action.choices, "keys"):
            return action.choices
    return {}


def logical_lines(text: str) -> list[tuple[int, str]]:
    """Join backslash-continued shell lines, keeping the first line number."""
    out: list[tuple[int, str]] = []
    pending = ""
    start = 0
    for lineno, line in enumerate(text.splitlines(), 1):
        stripped = line.rstrip()
        if not pending:
            start = lineno
        if stripped.endswith("\\"):
            pending += stripped[:-1] + " "
            continue
        out.append((start, (pending + stripped) if pending else stripped))
        pending = ""
    if pending:
        out.append((start, pending))
    return out


#: Flags whose value is a *produced* artefact this suite writes. A page may not read
#: one of these unless an earlier step on that page created it. Raw user inputs (BAM,
#: BigWig, FASTA, GTF, motif files) are deliberately excluded: a quick start is not
#: expected to create the sequencing data it analyses.
PRODUCED_FLAGS = {"--peak", "--peaks"}

#: Flags whose value the example creates.
OUTPUT_FLAGS = {"--out", "-o"}


def _values(tokens: list[str], flags: set[str]) -> set[str]:
    """Collect the values of ``flags`` from a token list (comma lists split)."""
    found, want = set(), False
    for token in tokens:
        if want:
            found.update(piece for piece in token.split(",") if piece)
            want = False
            continue
        if token in flags:
            want = True
        elif "=" in token and token.split("=", 1)[0] in flags:
            found.update(p for p in token.split("=", 1)[1].split(",") if p)
    return found


def check_documented_paths() -> list[str]:
    """Flag a page that reads a peak set no earlier step on that page produces.

    The parser check cannot see this. `--peak peaks/sample.narrowPeak` parses whether
    or not that file is ever written, so the quick start used to tell a reader to QC a
    file that step 2 writes under a different name, and all 81 examples "passed".

    Only peak paths are checked, because they are the one intermediate this suite both
    produces and consumes within a page. A peak path counts as produced when an
    ``--out`` on the page is that path, is its parent directory, or is an ancestor:
    `peaks` writes ``<out>/<stem>/<stem>_peaks.narrowPeak``, so the leaf name is the
    writer's choice. Pages that create nothing are fragments and are skipped.
    """
    problems: list[str] = []
    for page in sorted((REPO / "docs").rglob("*.md")):
        # `commands/<name>.md` documents one invocation in isolation; the inputs it
        # names are whatever the user already has. Only pages that present a *sequence*
        # can be inconsistent about what they produce.
        if page.parent.name == "commands":
            continue
        text = page.read_text(encoding="utf-8", errors="replace")
        outs: set[str] = set()
        peak_reads: list[tuple[int, str]] = []

        for lineno, line in logical_lines(text):
            match = re.match(r"^\s*easyepi\s+([a-z][a-z0-9-]*)\s*(.*)$", line)
            if not match:
                continue
            rest = match.group(2).split(" #")[0].strip()
            try:
                tokens = shlex.split(rest)
            except ValueError:
                continue
            outs.update(v.rstrip("/") for v in _values(tokens, OUTPUT_FLAGS))
            for v in _values(tokens, PRODUCED_FLAGS):
                v = v.rstrip("/")
                # Only paths are checked. A bare token (`--peaks ctrl1,ctrl2`) is a
                # sample label that `differential-peaks` resolves against --peak-dir,
                # not a file this page must create.
                if "/" not in v:
                    continue
                if v.startswith(("/", "~", "$")) or "${" in v:
                    continue
                peak_reads.append((lineno, v))

        if not outs:
            continue

        def produced(path: str) -> bool:
            """Whether an --out on this page plausibly writes ``path``.

            Being merely *under* an --out is not enough: `peaks -o peaks/` writes
            ``peaks/<stem>/<stem>_peaks.narrowPeak``, so ``peaks/sample.narrowPeak``
            sits inside that directory yet is never created. The leaf therefore has
            to look like something the writer names: the stem appears in it, or it
            carries a peak-style suffix the writer appends.
            """
            leaf = path.rsplit("/", 1)[-1]
            stem = leaf.split(".")[0]
            if stem.endswith("_peaks"):
                stem = stem[: -len("_peaks")]
            for out in outs:
                if path == out:
                    return True
                if not path.startswith(out + "/"):
                    continue
                parent = path[len(out) + 1:].rsplit("/", 1)[0] if "/" in path[len(out) + 1:] else ""
                # `peaks` names the per-sample directory and prefixes the file with
                # the same stem, so both carry it.
                if parent and parent in leaf:
                    return True
                if stem and stem in parent:
                    return True
            return False

        for lineno, path in peak_reads:
            if not produced(path):
                problems.append(
                    f"{page.relative_to(REPO)}:{lineno}: reads peak set {path!r}, "
                    f"which no step on this page writes (writes: "
                    f"{', '.join(sorted(outs))})")
    return problems


SITE_URL = "https://summus-kong.github.io/easyepi/"
DOC_LINK = re.compile(re.escape(SITE_URL) + r"[A-Za-z0-9_./-]*")


def nav_pages() -> set[str]:
    """Absolute URLs of every page listed in mkdocs.yml, derived the same way."""
    config = (REPO / "mkdocs.yml").read_text(encoding="utf-8")
    pages = set()
    for md in re.findall(r":\s+([A-Za-z0-9_./-]+\.md)\s*$", config, flags=re.M):
        stem = md[:-3]
        if stem == "index":
            pages.add(SITE_URL)
        elif stem.endswith("/index"):
            pages.add(SITE_URL + stem[:-6] + "/")
        else:
            pages.add(SITE_URL + stem + "/")
    return pages


def check_readme_links() -> list[str]:
    """README.md links to published pages; every one must exist in the nav.

    Checked against ``mkdocs.yml`` rather than by fetching the site: a network check
    would fail for reasons unrelated to the repository, and the nav is the same source
    the build uses to decide what exists. A link that is not in the nav is either a
    typo or a page that was removed, and both should fail before publishing.
    """
    readme = REPO / "README.md"
    if not readme.exists():
        return []
    known = nav_pages()
    problems = []
    text = readme.read_text(encoding="utf-8", errors="replace")
    for lineno, line in enumerate(text.splitlines(), 1):
        for url in DOC_LINK.findall(line):
            if url not in known:
                problems.append(f"README.md:{lineno}: {url} is not a page in mkdocs.yml")
    return problems


def main() -> int:
    parser = build_parser()
    known = subparsers(parser)
    sources = sorted((REPO / "docs").rglob("*.md")) + [REPO / "README.md"]
    pattern = re.compile(r"^\s*easyepi\s+([a-z][a-z0-9-]*)\s*(.*)$")

    checked = 0
    rejected: list[tuple[str, int, str]] = []
    skipped_shell: list[str] = []
    unknown: set[str] = set()

    for path in sources:
        text = path.read_text(encoding="utf-8", errors="replace")
        for lineno, line in logical_lines(text):
            match = pattern.match(line)
            if not match:
                continue
            command, rest = match.group(1), match.group(2)
            if command not in known:
                unknown.add(command)
                continue
            rest = rest.split(" #")[0].strip()
            if any(token in rest for token in SHELLY):
                skipped_shell.append(f"{path.relative_to(REPO)}:{lineno}")
                continue
            try:
                tokens = shlex.split(rest)
            except ValueError:
                continue
            checked += 1
            try:
                with contextlib.redirect_stderr(io.StringIO()):
                    parser.parse_args([command, *tokens])
            except SystemExit:
                rejected.append((str(path.relative_to(REPO)), lineno, line.strip()))

    print(f"examples parsed : {checked}")
    if unknown:
        print(f"unknown commands: {sorted(unknown)}")
    if skipped_shell:
        print(f"skipped (shell expansions), {len(skipped_shell)}: "
              + ", ".join(skipped_shell))
    if rejected:
        print(f"\nREJECTED by the live parser: {len(rejected)}")
        for rel, lineno, line in rejected:
            print(f"  {rel}:{lineno}: {line}")
        return 1
    print("every example is accepted by the live parser")

    path_problems = check_documented_paths()
    if path_problems:
        print(f"\nPATHS a page reads but never produces: {len(path_problems)}")
        for problem in path_problems:
            print(f"  {problem}")
        return 1
    print("every example reads paths the same page produces")

    link_problems = check_readme_links()
    if link_problems:
        print(f"\nREADME doc links that match no nav page: {len(link_problems)}")
        for problem in link_problems:
            print(f"  {problem}")
        return 1
    print("every documentation link in README.md resolves to a nav page")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
