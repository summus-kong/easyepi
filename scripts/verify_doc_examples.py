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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
