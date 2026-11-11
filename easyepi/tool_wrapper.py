"""A shared way to wrap an external tool: sensible defaults plus verbatim passthrough.

Why a shared helper
-------------------
Several commands drive an external program that has more options than the suite will
ever model: MACS3 for peak calling, and potentially samtools or bedtools elsewhere.
Two rules must hold for every one of them:

1. **The suite supplies assay-appropriate defaults.** A user types the inputs and an
   assay; the wrapper chooses formats, shifts and cut-offs.
2. **Anything the wrapper does not model can be passed through verbatim**, using the
   upstream tool's own spelling, so a user is never blocked waiting for a wrapper
   release. The passthrough is recorded in the run manifest, so a result stays
   reproducible even for options the suite knows nothing about.

The mechanism
-------------
``ArgvPassthrough`` collects the modelled options as an argv list, appends the user's
extra tokens, and hands the whole thing to a ``build`` callable. Because argparse
takes the last occurrence of a repeated flag, appending the passthrough means the user
overrides the suite's defaults -- the behaviour a caller expects from "extra options".

Validation happens against the upstream parser, not against a copy of its rules, so a
typo raises immediately with the tool's own message instead of silently doing nothing.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence


@dataclass
class ArgvPassthrough:
    """Accumulate argv for an external tool, with a recorded passthrough."""

    tool: str
    #: Tokens the wrapper modelled itself.
    modelled: list[str] = field(default_factory=list)
    #: Tokens the user supplied via ``--<tool>-arg``.
    passthrough: list[str] = field(default_factory=list)
    #: How the final argv was validated, for the run manifest.
    validated_by: str | None = None
    #: A non-fatal explanation when validation could not run.
    validation_error: str | None = None

    def add(self, *tokens: str) -> "ArgvPassthrough":
        self.modelled.extend(str(t) for t in tokens)
        return self

    def add_flag(self, flag: str, value: Any | None = None) -> "ArgvPassthrough":
        """Add ``flag`` and, when not ``None``, its value."""
        self.modelled.append(flag)
        if value is not None:
            self.modelled.append(str(value))
        return self

    def add_extra(self, extra: Sequence[str] | None) -> "ArgvPassthrough":
        """Append user passthrough, splitting shell-quoted strings.

        ``--macs-arg="--nomodel --extsize 200"`` is what a user writes first, so the
        single-string form has to work; ``shlex.split`` also handles a quoted value
        that contains spaces.
        """
        for raw in extra or []:
            self.passthrough.extend(shlex.split(str(raw)))
        return self

    @property
    def argv(self) -> list[str]:
        """The full argv: modelled options first, passthrough last."""
        return [self.tool, *self.modelled, *self.passthrough]

    def as_dict(self) -> dict[str, Any]:
        return {"tool": self.tool,
                "argv": self.argv,
                "modelled": list(self.modelled),
                "passthrough": list(self.passthrough),
                "validated_by": self.validated_by,
                "validation_error": self.validation_error}


def validate_against_parser(
    argv: Sequence[str],
    parser_factory: "Callable[[], Any] | None",
    *,
    has_program_name: bool = True,
    parser_is_top_level: bool = False,
) -> tuple[Any | None, str | None]:
    """Parse ``argv`` with the tool's own parser.

    The two flags name two independent facts, because conflating them mis-parses:

    ``has_program_name``
        ``argv`` starts with the executable (as ``ArgvPassthrough.argv`` does) and that
        token must be dropped before parsing.
    ``parser_is_top_level``
        the parser owns sub-commands, so the token after the program name is a
        sub-command and must be **kept**. For a leaf parser (one command's options only)
        the sub-command token must be dropped.

    Returns ``(namespace, error)``; ``namespace`` is ``None`` when no factory was given
    or parsing failed.

    A ``SystemExit`` from argparse means the user's passthrough was rejected. That is a
    usage error, raised as ``ValueError`` naming the tokens, because silently falling
    back to the suite's defaults would hide the typo.
    """
    if parser_factory is None:
        return None, "no parser factory supplied"
    try:
        parser = parser_factory()
    except Exception as exc:                        # noqa: BLE001
        return None, f"could not build the tool's parser: {type(exc).__name__}: {exc}"

    tokens = list(argv)
    if has_program_name and tokens:
        tokens = tokens[1:]
    if tokens and not parser_is_top_level and not tokens[0].startswith("-"):
        # a leaf parser does not know the sub-command name
        tokens = tokens[1:]

    try:
        return parser.parse_args(tokens), None
    except SystemExit:
        raise ValueError(
            f"{argv[0] if argv else 'tool'} rejected the options: {' '.join(tokens)}\n"
            f"  check the values passed through against `"
            f"{argv[0] if argv else 'tool'} --help`")


def describe(passthrough: ArgvPassthrough, indent: str = "  ") -> str:
    """One-line-per-fact summary for the CLI to print before running."""
    lines = [f"{indent}caller    : {passthrough.tool}"]
    if passthrough.passthrough:
        lines.append(f"{indent}passthrough: {' '.join(passthrough.passthrough)}")
    if passthrough.validated_by:
        lines.append(f"{indent}validated : {passthrough.validated_by}")
    if passthrough.validation_error:
        lines.append(f"{indent}note      : {passthrough.validation_error}")
    return "\n".join(lines)
