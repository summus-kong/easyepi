"""Actionable error messages.

Why a dedicated module
----------------------
Three classes of failure account for almost every support question, and each has a
different fix:

1. **A path does not exist.** The user needs to know which path was looked for and
   what the nearest existing alternative is.
2. **An optional dependency is missing.** The user needs the exact install command
   for *their* interpreter, not a generic pip line.
3. **A required argument is missing or contradictory.** The user needs a working
   example, not the argparse usage wall.

None of these are silent, but the default Python messages are hard to act on:
``ModuleNotFoundError: No module named 'pysam'`` does not say how to install it into
the environment that is actually running, and ``FileNotFoundError: [Errno 2] ...``
does not suggest the typo.

The helpers here raise the same exception types as before (so callers and tests that
catch ``FileNotFoundError`` keep working) but with a message that ends in an
action. They never swallow an error and never substitute a default for a missing
input.
"""

from __future__ import annotations

import difflib
import importlib.util
import sys
from pathlib import Path
from typing import Iterable, Sequence


#: Install hint per optional dependency. The extras name matches pyproject.toml.
DEPENDENCY_HINTS: dict[str, tuple[str, str]] = {
    # module: (extras group, what it is for)
    "pysam": ("bam", "reading BAM alignments"),
    "pyBigWig": ("bam", "reading and writing BigWig signal"),
    "numpy": ("core", "numeric arrays"),
    "scipy": ("core", "statistics and clustering"),
    "pandas": ("core", "tables"),
    "matplotlib": ("viz", "figures"),
    "hmmlearn": ("states", "chromatin state learning"),
    "sklearn": ("viz", "clustering in figures"),
}


def missing_dependency(module: str, purpose: str | None = None) -> ImportError:
    """Build the ImportError for a missing optional dependency.

    The message names the module, what it is needed for, and the exact install
    command for the running interpreter -- including the interpreter path, because
    a user with several environments otherwise installs into the wrong one.
    """
    extra, default_purpose = DEPENDENCY_HINTS.get(module, (None, None))
    why = purpose or default_purpose
    lines = [f"{module} is required" + (f" for {why}." if why else ".")]
    if extra:
        lines.append(
            f'  install: "{sys.executable}" -m pip install "easyepi[{extra}]"')
    else:
        lines.append(f'  install: "{sys.executable}" -m pip install {module}')
    lines.append(f"  python in use: {sys.executable}")
    return ImportError("\n".join(lines))


def require_dependency(module: str, purpose: str | None = None):
    """Import ``module`` or raise with an actionable message."""
    try:
        return importlib.import_module(module)
    except ImportError:
        pass
    # a package can exist but fail to import (broken build); re-raise that as-is so a
    # real installation problem is not reported as "not installed"
    if importlib.util.find_spec(module) is not None:
        raise
    raise missing_dependency(module, purpose) from None


def _suggest(name: str, candidates: Iterable[str], limit: int = 3) -> str:
    matches = difflib.get_close_matches(str(name), list(candidates), n=limit,
                                       cutoff=0.6)
    if not matches:
        return ""
    return "  did you mean: " + ", ".join(matches)


def missing_file(path: str | Path, purpose: str | None = None,
                 search_dir: str | Path | None = None,
                 suffixes: Sequence[str] | None = None) -> FileNotFoundError:
    """Build the FileNotFoundError for an absent input.

    Includes the resolved absolute path (so a relative-path mistake is obvious), the
    nearest existing filename in the same folder when one is close, and the files
    that *do* exist when a suffix filter was expected -- which is what turns "file
    not found" into "you passed a .bam where a .bw was expected".
    """
    target = Path(path).expanduser()
    lines = [f"input not found: {path}"
             + (f" (needed for {purpose})" if purpose else "")]
    if target.is_absolute():
        lines.append(f"  looked for: {target}")
    else:
        lines.append(f"  looked for: {target.resolve()}  (relative to {Path.cwd()})")

    parent = target.parent if str(target.parent) not in ("", ".") else Path.cwd()
    if search_dir is not None:
        parent = Path(search_dir).expanduser()
    if parent.is_dir():
        if suffixes:
            present = sorted(p.name for p in parent.iterdir()
                             if p.is_file() and any(p.name.endswith(s)
                                                    for s in suffixes))
            if present:
                lines.append(f"  found these instead in {parent}: "
                             + ", ".join(present[:6])
                             + (" ..." if len(present) > 6 else ""))
            else:
                lines.append(f"  no {', '.join(suffixes)} files in {parent}")
        else:
            names = [p.name for p in parent.iterdir() if p.is_file()]
            hint = _suggest(target.name, names)
            if hint:
                lines.append(hint)
    else:
        lines.append(f"  the directory does not exist either: {parent}")
    return FileNotFoundError("\n".join(lines))


def require_file(path: str | Path, purpose: str | None = None,
                 search_dir: str | Path | None = None,
                 suffixes: Sequence[str] | None = None) -> Path:
    """Return the path if it exists, else raise an actionable error."""
    target = Path(path).expanduser()
    if not target.exists():
        raise missing_file(path, purpose, search_dir, suffixes)
    if not target.is_file():
        raise IsADirectoryError(
            f"{path} is a directory but a file was expected"
            + (f" (needed for {purpose})" if purpose else "")
            + "\n  pass the file inside it, not the folder")
    return target


def missing_executable(name: str, purpose: str, alternatives: Sequence[str] = ()) -> FileNotFoundError:
    """Build the error for an external program that is not on PATH."""
    lines = [f"{name} not found on PATH" + (f" (needed for {purpose})" if purpose else "")]
    if alternatives:
        lines.append("  alternatives: " + ", ".join(alternatives))
    lines.append(f"  PATH checked: {__import__('os').environ.get('PATH', '')[:200]}")
    return FileNotFoundError("\n".join(lines))


def require_executable(name: str, purpose: str,
                       alternatives: Sequence[str] = ()) -> str:
    """Locate an external program or raise with the options."""
    import shutil

    found = shutil.which(name)
    if found:
        return found
    for alt in alternatives:
        found = shutil.which(alt)
        if found:
            return found
    raise missing_executable(name, purpose, alternatives)


class UsageError(ValueError):
    """A contradictory or insufficient set of arguments, with a working example."""

    def __init__(self, message: str, example: str | None = None):
        text = message
        if example:
            text += "\n  working example:\n    " + example.replace("\n", "\n    ")
        super().__init__(text)
        self.example = example
