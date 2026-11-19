"""Dispatch independent script entry points to shared CLI functions."""

from __future__ import annotations

import sys
from typing import Sequence


def dispatch(command: str, argv: Sequence[str] | None = None) -> int:
    """Forward an independent script entry point to the shared CLI."""
    from easyepi.cli import main
    return main([command, *(list(argv) if argv is not None else sys.argv[1:])])

