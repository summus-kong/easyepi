"""Shared runner for the independent ``easyepi.universal.*`` script entry points.

Each of those modules used to repeat the same five steps by hand:

    output = Path(args.out)
    directory = output.parent if output.suffix == ".tsv" else output
    logger = setup_logging(NAME, directory, verbose=..., quiet=..., log_file=...)
    with RunRecorder(NAME, directory, logger, figure_formats=selected_formats(args)) as run:
        do_work(args)
        run.set_metrics(...)

That boilerplate is what made the entry points look like duplicated logic. The
analysis itself already lives in shared modules (``easyepi.normalize``,
``easyepi.integration``, ``easyepi.motifs.external_tools``, …); only the wrapper
was copied. ``run_entry`` is that wrapper, in one place.

Kept deliberately dumb: it does not parse arguments and does not know what the
work does, so each entry point keeps its own distinct CLI surface.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

from ._plots import selected_formats
from .logging_utils import RunRecorder, setup_logging

# Suffixes that mean "args.out is a file, not a directory".
_FILE_SUFFIXES = {".tsv", ".csv", ".json", ".txt"}


def output_directory(out: str | Path) -> Path:
    """Return the directory that should hold logs/figures for this ``--out``."""
    path = Path(out)
    return path.parent if path.suffix.lower() in _FILE_SUFFIXES else path


@contextmanager
def run_entry(module: str, args: Any, **metrics: Any) -> Iterator[RunRecorder]:
    """Open a logged, figure-aware run and record summary metrics on success.

    Usage::

        with run_entry("NORMALIZE", args, method=args.method) as run:
            run_normalize(args)

    ``render_figures`` follows the entry point's ``--figure-format`` flags, the
    same as the main CLI.
    """
    directory = output_directory(args.out)
    directory.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(
        module,
        directory,
        verbose=getattr(args, "verbose", False),
        quiet=getattr(args, "quiet", False),
        log_file=getattr(args, "log_file", None),
    )
    formats = selected_formats(args)
    with RunRecorder(module, directory, logger, figure_formats=formats) as run:
        yield run
        if metrics:
            run.set_metrics(**metrics)


def redirect_out_to_file(args: Any, filename: str) -> None:
    """Turn a directory ``--out`` into a specific file path, in place.

    Some entry points accept either a file or a directory for ``--out``; this
    keeps that behaviour in one place instead of in every module.
    """
    path = Path(args.out)
    if path.suffix.lower() not in _FILE_SUFFIXES:
        args.out = str(path / filename)
