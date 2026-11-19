"""High-level defaults with a passthrough for the underlying tool's own options.

The rule this module implements
------------------------------
A user should be able to type nothing but the essential inputs and get sensible,
assay-appropriate behaviour; and a user who needs an option the suite does not model
should be able to pass it straight through, using the upstream tool's own spelling,
without waiting for a wrapper release.

    easyepi peaks sample.bam -o peaks/ --assay atac
    easyepi peaks sample.bam -o peaks/ --assay atac --macs-arg="--nomodel --extsize 200"

The passthrough is recorded verbatim in the run manifest, so a result can be
reproduced even when the option is one this suite knows nothing about.

How the Namespace is built
--------------------------
MACS3's console script has ``prepare_argparser()``, which builds the full argparse tree
including every default. Loading the script as a module normally runs ``main()`` and
calls ``sys.exit``, so the factory is reached by executing the script's source with the
``main`` call guarded: the module-level code defines the functions, and
``prepare_argparser`` is then called directly with ``sys.argv`` set to the intended
command line.

That produces exactly the Namespace the command line would have produced -- every
default included -- which is why this module no longer hand-builds one. Hand-building
missed attributes that MACS3 compares against numbers, surfacing as
``TypeError: '<' not supported between 'NoneType' and 'int'`` from inside MACS3.

If the console script cannot be read, the fallback is explicit and recorded.
"""

from __future__ import annotations

import shlex
import sys
from argparse import Namespace
from pathlib import Path
from typing import Any, Sequence


def macs3_script() -> Path | None:
    """Path of the ``macs3`` console script for the running interpreter."""
    import shutil

    found = shutil.which("macs3")
    if found:
        return Path(found)
    candidate = Path(sys.executable).parent / "macs3"
    return candidate if candidate.exists() else None


def _call_argparser_factory(argv: Sequence[str]):
    """Run MACS3's own ``prepare_argparser`` with ``argv`` in place.

    The script is executed with ``__name__`` set to something other than
    ``"__main__"``, so its module-level ``main()`` guard does not fire, and its source
    is read rather than imported so no ``.py`` file is required on disk.

    Returns ``(parser_or_namespace, error_message)``. On failure the message explains
    which step broke, because a silent fallback to the hand-built Namespace is how the
    original TypeError reached the user.
    """
    script = macs3_script()
    if script is None:
        return None, "no macs3 console script next to the interpreter"

    source = script.read_text(encoding="utf-8", errors="replace")
    # the script ends with a bare `if __name__ == "__main__": main()`; making __name__
    # different prevents main() from parsing the real process argv.
    namespace = {"__name__": "_easyepi_macs3_loader", "__file__": str(script)}

    saved_argv = sys.argv
    try:
        sys.argv = list(argv)
        exec(compile(source, str(script), "exec"), namespace)  # noqa: S102
        factory = namespace.get("prepare_argparser")
        if factory is None:
            return None, "prepare_argparser not found in the macs3 script"
        parser = factory()
        return parser, None
    except SystemExit as exc:
        return None, f"macs3 script exited while building the parser (code {exc.code})"
    except Exception as exc:                      # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"
    finally:
        sys.argv = saved_argv


def split_passthrough(extra: Sequence[str] | None) -> list[str]:
    """Split ``--macs-arg`` values into individual tokens.

    Accepts both ``--macs-arg="--nomodel --extsize 200"`` and a list, because quoting
    one string is what a user reaches for first and ``shlex`` handles it correctly
    including quoted values with spaces.
    """
    tokens: list[str] = []
    for raw in extra or []:
        tokens.extend(shlex.split(str(raw)))
    return tokens


def build_namespace(*, treatment: list[str], out_dir: Path, name: str, genome: str,
                    fmt: str, broad: bool, control: list[str] | None, qvalue: float,
                    extsize: int | None, shift: int, keep_dup: str,
                    extra_args: Sequence[str] | None = None
                    ) -> tuple[Namespace, dict[str, Any]]:
    """Return ``(namespace, provenance)`` ready for ``callpeak_cmd.run``.

    Suite-modelled options come first, then the passthrough, so a user can override
    anything this module chose -- argparse takes the last occurrence of a repeated
    flag, which is the behaviour a caller expects from "extra options".
    """
    argv = ["macs3", "callpeak", "-t", *treatment]
    if control:
        argv += ["-c", *control]
    argv += ["-f", fmt, "-g", genome, "-n", name, "--outdir", str(out_dir),
             "-q", str(qvalue), "--keep-dup", keep_dup]
    if broad:
        argv.append("--broad")
    if extsize is not None:
        argv += ["--extsize", str(extsize)]
    if shift:
        argv += ["--shift", str(shift)]

    passthrough = split_passthrough(extra_args)
    argv += passthrough

    parser, error = _call_argparser_factory(argv[1:])
    provenance: dict[str, Any] = {
        "argv": argv,
        "passthrough": passthrough,
        "namespace_from": None,
        "namespace_error": error,
    }

    if parser is not None:
        try:
            options = parser.parse_args(argv[1:])
        except SystemExit:
            # a passthrough token was rejected by MACS3's own parser; re-raise as a
            # usage error so the message names the offending argument
            raise ValueError(
                "MACS3 rejected the options: " + " ".join(argv[1:])
                + "\n  check the --macs-arg value against `macs3 callpeak --help`")
        provenance["namespace_from"] = "macs3-argparse"
        return options, provenance

    provenance["namespace_from"] = "fallback"
    return _fallback_namespace(treatment=treatment, out_dir=out_dir, name=name,
                               genome=genome, fmt=fmt, broad=broad,
                               control=control, qvalue=qvalue, extsize=extsize,
                               shift=shift, keep_dup=keep_dup), provenance


#: Values MACS3's own defaults produce, used only when its parser cannot be loaded.
#: Attribute names come from what ``opt_validate_callpeak`` reads.
_FALLBACK: dict[str, Any] = {
    "bdg_control": None, "bdg_treat": None, "broadcutoff": None, "bw": 200,
    "call_summits": False, "cutoff_analysis": False, "cutoff_analysis_file": None,
    "d_min": 20, "do_SPMR": False, "downsample": False, "fecutoff": 0.1,
    "gzip_flag": False, "largelocal": 10000, "lmfold": 1.0, "maxcount": None,
    "maxgap": 30, "mfold": 5.0, "minlen": 200, "nomodel": False, "pvalue": None,
    "ratio": 1.0, "refine_peaks": False, "scaleto": None, "seed": -1,
    "smalllocal": 1000, "store_bdg": False, "tsize": 147, "umfold": 100.0,
    "verbose": 2,
}


def _fallback_namespace(**kwargs) -> Namespace:
    """Minimal Namespace for an install whose console script cannot be read."""
    fmt = kwargs["fmt"]
    values = dict(_FALLBACK)
    values.update({
        "tfile": list(kwargs["treatment"]),
        "cfile": list(kwargs["control"] or []),
        "format": fmt,
        "gsize": kwargs["genome"],
        "name": kwargs["name"],
        "outdir": str(kwargs["out_dir"]),
        "broad": bool(kwargs["broad"]),
        "qvalue": float(kwargs["qvalue"]),
        "extsize": kwargs["extsize"],
        "shift": int(kwargs["shift"]),
        "keepduplicates": kwargs["keep_dup"],
        "PE_MODE": fmt in ("BAMPE", "BEDPE", "FRAG"),
    })
    return Namespace(**values)
