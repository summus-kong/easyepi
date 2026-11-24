"""Record what went into a run so a later step can check it still does.

BACKLOG item 10 asked for four things that were missing: input hashes,
output hashes, direct dependency versions and external tool versions.
The design constraints here come from how the runs are actually used:

  * A BAM can be tens of gigabytes. Hashing it is minutes of IO, so the
    caller decides whether to pay for it (`digest_files`), unlike the
    dependency versions which are free.
  * Everything must survive a missing file, a missing tool or a missing
    package. A run that already produced correct output should not fail
    because `macs3 --version` is slow or the file was moved afterwards,
    but the record has to say which part could not be established rather
    than silently returning an empty dict.
  * Paths are recorded relative to a root when possible, so a manifest
    stays readable when the run directory moves.

Nothing here writes files; `RunRecorder` owns that. This module only
builds the dictionaries.
"""

from __future__ import annotations

import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
from typing import Any, Iterable, Sequence

# Packages whose version can change a number that ends up in a report.
# `easyepi` itself is included because a bug fix changes results too.
DEFAULT_PACKAGES = (
    "easyepi",
    "numpy",
    "scipy",
    "pandas",
    "pysam",
    "pyBigWig",
    "pydeseq2",
)

# External binaries the package shells out to. Each entry is
# `(executable, arguments that print a version)`.
DEFAULT_TOOLS = {
    "samtools": ("--version",),
    "bedtools": ("--version",),
    "bowtie2": ("--version",),
    "macs3": ("--version",),
    "homer": ("findPeaks", "-h"),
    "bedGraphToBigWig": (),
}

_CHUNK = 1024 * 1024


def _sha256(path: Path) -> str:
    """Stream a file through SHA-256 so memory stays flat on big BAMs."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _relative(path: Path, root: Path | None) -> str:
    """Prefer a path relative to `root`; fall back to the absolute one."""
    if root is not None:
        try:
            return str(path.resolve().relative_to(root.resolve()))
        except (ValueError, OSError):
            pass
    return str(path)


def digest_files(
    paths: Iterable[str | os.PathLike[str]],
    *,
    root: str | os.PathLike[str] | None = None,
    limit_bytes: int | None = None,
) -> dict[str, dict[str, Any]]:
    """SHA-256 and size for each existing file, keyed by its relative path.

    Returns a per-file record rather than a flat mapping so a caller can
    see *why* an entry is missing:

        {"sample.bam": {"sha256": "...", "size": 123}}
        {"gone.bam":   {"error": "not_found"}}

    A directory or a vanished path becomes an `error` entry instead of
    raising, because a manifest written after the fact should still be
    produced. `limit_bytes` skips files larger than the given size and
    records `skipped_too_large`, which keeps an accidental call on a
    multi-terabyte tree from stalling a run.
    """
    root_path = Path(root) if root is not None else None
    out: dict[str, dict[str, Any]] = {}
    for raw in paths:
        path = Path(raw)
        key = _relative(path, root_path)
        if key in out:
            continue
        if not path.exists():
            out[key] = {"error": "not_found"}
            continue
        if not path.is_file():
            out[key] = {"error": "not_a_file"}
            continue
        try:
            size = path.stat().st_size
            if limit_bytes is not None and size > limit_bytes:
                out[key] = {"error": "skipped_too_large", "size": size}
                continue
            out[key] = {"sha256": _sha256(path), "size": size}
        except OSError as exc:
            out[key] = {"error": f"unreadable: {exc.strerror or exc}"}
    return out


def package_versions(names: Sequence[str] = DEFAULT_PACKAGES) -> dict[str, str]:
    """Installed version of each package, or `not installed`.

    The literal string "not installed" is used instead of `None` because
    it reads correctly in `status.json` and because existing consumers of
    `RunRecorder.status["versions"]` already carry that convention.
    """
    out: dict[str, str] = {"python": platform.python_version()}
    for name in names:
        try:
            out[name] = version(name)
        except PackageNotFoundError:
            out[name] = "not installed"
        except Exception:  # a broken dist-info must not abort the manifest
            out[name] = "unknown"
    return out


_VERSION_RE = re.compile(r"\d+\.\d+(?:\.\d+)?")


def _first_version_token(text: str) -> str | None:
    """First dotted version number in `text`.

    Tool banners vary wildly ("samtools 1.19", "bedtools v2.31.0",
    "bowtie2-align-s version 2.5.1"). Extracting the first number is
    crude but stable, and the raw banner is kept alongside it so a
    surprising value can be traced back.
    """
    match = _VERSION_RE.search(text)
    return match.group(0) if match else None


def resolve(executable: str, search_paths: Sequence[str] = ()) -> str | None:
    """Absolute path to `executable`, or None.

    Tried in this order: the value itself if it is an existing file,
    then `PATH` via `shutil.which`, then each directory in
    `search_paths`. The last step is what makes a conda env usable
    without activating it, which is how the analysis host runs.
    """
    candidate = Path(executable)
    if candidate.is_file():
        return str(candidate)
    found = shutil.which(executable)
    if found is not None:
        return found
    for directory in search_paths:
        trial = Path(directory) / executable
        if trial.is_file():
            return str(trial)
    return None


def _banner(executable: str, args: Sequence[str],
            search_paths: Sequence[str] = ()) -> tuple[str | None, str, str | None]:
    """Run a version command and return `(banner, status, resolved_path)`.

    `status` is "ok", "not_found" or "failed"; the banner is the first
    line of stdout, else the first line of stderr, because some tools
    print their banner to stderr. The resolved path is reported so a
    reader can tell which build produced the number.
    """
    exe = resolve(executable, search_paths)
    if exe is None:
        return None, "not_found", None
    try:
        proc = subprocess.run(
            [exe, *args],
            capture_output=True,
            text=True,
            timeout=20,
        )
    except (OSError, subprocess.SubprocessError):
        return None, "failed", exe
    for stream in (proc.stdout, proc.stderr):
        if stream and stream.strip():
            return stream.strip().splitlines()[0], "ok", exe
    return None, "failed", exe


def external_tools(
    tools: dict[str, Sequence[str]] | None = None,
    search_paths: Sequence[str] = (),
) -> dict[str, dict[str, Any]]:
    """Version banner for each external binary this package may call.

    An unavailable tool is recorded as `not_found` rather than omitted:
    "we looked and it was absent" is the useful statement, and it makes a
    manifest diffable across machines.

    `search_paths` lets a caller point at directories outside `PATH`,
    which is the normal situation when the tools live in a conda
    environment that the run does not activate.
    """
    wanted = tools if tools is not None else DEFAULT_TOOLS
    out: dict[str, dict[str, Any]] = {}
    for name, args in wanted.items():
        banner, status, path = _banner(name, args, search_paths)
        entry: dict[str, Any] = {"status": status}
        if path is not None:
            entry["path"] = path
        if banner is not None:
            entry["banner"] = banner
            found = _first_version_token(banner)
            if found:
                entry["version"] = found
        out[name] = entry
    return out


def environment() -> dict[str, Any]:
    """Platform facts that change nothing but explain a lot.

    Deliberately excludes hostname, user and absolute prefixes: a
    manifest ends up in supplementary material, and those three are
    machine identity rather than reproducibility.
    """
    return {
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor_count": os.cpu_count(),
        "argv0": Path(sys.argv[0]).name,
    }


def build(
    *,
    inputs: Iterable[str | os.PathLike[str]] = (),
    outputs: Iterable[str | os.PathLike[str]] = (),
    root: str | os.PathLike[str] | None = None,
    hash_inputs: bool = True,
    hash_outputs: bool = False,
    packages: Sequence[str] = DEFAULT_PACKAGES,
    tools: dict[str, Sequence[str]] | None = None,
    tool_search_paths: Sequence[str] = (),
    limit_bytes: int | None = None,
) -> dict[str, Any]:
    """Assemble the whole record.

    Output hashing defaults to off because it is normally done *after* a
    run's files exist, i.e. by a second call; `finish` in the recorder
    passes `hash_outputs=True` at the end. Input hashing defaults to on
    because that is the check a reader will want.
    """
    record: dict[str, Any] = {
        "environment": environment(),
        "packages": package_versions(packages),
        "external_tools": external_tools(tools, tool_search_paths),
    }
    if hash_inputs:
        record["inputs"] = digest_files(inputs, root=root, limit_bytes=limit_bytes)
    if hash_outputs:
        record["outputs"] = digest_files(outputs, root=root, limit_bytes=limit_bytes)
    return record


def describe_gaps(record: dict[str, Any]) -> list[str]:
    """Human-readable list of what the record could not establish.

    A manifest that quietly says "not_found" for everything looks the
    same as one with nothing to check. This turns the silent cases into
    sentences a user can act on, and the recorder logs them.
    """
    gaps: list[str] = []
    for section in ("inputs", "outputs"):
        for path, entry in (record.get(section) or {}).items():
            error = entry.get("error")
            if error:
                gaps.append(f"{section}: {path} -> {error}")
    for name, entry in (record.get("external_tools") or {}).items():
        if entry.get("status") != "ok":
            gaps.append(f"external tool: {name} -> {entry.get('status')}")
    for name, value in (record.get("packages") or {}).items():
        if value in ("not installed", "unknown"):
            gaps.append(f"package: {name} -> {value}")
    return gaps
