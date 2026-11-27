#!/usr/bin/env python3
"""Check an *installed* easyepi from outside its source tree.

`docs/PUBLISHING.md` step 4 asks for two things that a source-tree test run
cannot answer:

  * every subcommand must print `--help` from an installation
  * importing `easyepi` must not be shadowed by the source directory

A checkout always satisfies both by accident (the tree is on `sys.path`), so
they only mean something when run against a wheel installed elsewhere.  This
script therefore refuses to run from inside the source tree: it exits non-zero
if `easyepi` resolves next to a `pyproject.toml`.

Usage::

    python scripts/verify_installed.py            # from any directory
    python scripts/verify_installed.py --quiet    # only failures

Exit status is 0 only when every subcommand helps and every module imports.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import pkgutil
import subprocess
import sys
from pathlib import Path

#: Paths that exist in old layouts.  Their presence in a wheel means the build
#: swept in stale directories rather than only the current package.
LEGACY_MODULES = ("easyepi.python", "easyepi.backends", "easyepi.scripts")


def _installed_root() -> Path:
    import easyepi

    return Path(easyepi.__file__).resolve().parent


def _source_tree_root() -> Path | None:
    """Walk up from the import location looking for a checkout marker."""
    for parent in _installed_root().parents:
        if (parent / "pyproject.toml").is_file():
            return parent
    return None


def check_not_shadowed(quiet: bool) -> list[str]:
    root = _installed_root()
    problems = []
    if _source_tree_root() is not None:
        problems.append(
            f"easyepi imported from {root}, which sits inside a source tree; "
            "run this from an installation that is not the checkout"
        )
    if not quiet:
        print(f"imported from : {root}")
    return problems


def check_legacy_absent(quiet: bool) -> list[str]:
    problems = []
    for name in LEGACY_MODULES:
        try:
            spec = importlib.util.find_spec(name)
        except ModuleNotFoundError:
            spec = None
        if spec is not None:
            problems.append(f"legacy module shipped: {name}")
        elif not quiet:
            print(f"absent        : {name}")
    return problems


def subcommands() -> list[str]:
    from easyepi.cli import build_parser

    parser = build_parser()
    for action in parser._actions:
        names = getattr(action, "_name_parser_map", None)
        if names:
            return sorted(names)
    return []


def check_subcommands(quiet: bool) -> list[str]:
    names = subcommands()
    problems = []
    if not names:
        return ["no subcommands found on the parser"]
    for name in names:
        proc = subprocess.run(
            [sys.executable, "-m", "easyepi.cli", name, "--help"],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            first = (proc.stderr or proc.stdout).strip().splitlines()
            problems.append(f"{name} --help failed: {first[0] if first else '?'}")
        elif not quiet:
            print(f"ok            : {name} --help")
    if not quiet:
        print(f"\nsubcommands checked: {len(names)}")
    return problems


def check_module_imports(quiet: bool) -> list[str]:
    import easyepi

    problems = []
    walked = 0
    for mod in pkgutil.walk_packages(easyepi.__path__, prefix="easyepi."):
        walked += 1
        try:
            importlib.import_module(mod.name)
        except Exception as exc:  # noqa: BLE001 - report, never abort the sweep
            problems.append(f"{mod.name}: {type(exc).__name__}: {exc}")
    if not quiet:
        print(f"modules imported: {walked}")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true", help="print only failures")
    args = ap.parse_args()

    problems: list[str] = []
    for check in (
        check_not_shadowed,
        check_legacy_absent,
        check_subcommands,
        check_module_imports,
    ):
        problems.extend(check(args.quiet))
        if not args.quiet:
            print()

    if problems:
        print("FAILURES")
        for p in problems:
            print("  " + p)
        return 1
    print("installed layout verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
