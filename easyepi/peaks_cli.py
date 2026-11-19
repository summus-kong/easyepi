"""Peak calling through MACS3, invoked as a library rather than a subprocess.

Why import instead of shelling out
----------------------------------
MACS3 ships ``MACS3.Commands.callpeak_cmd.run(args)``, which takes an argparse
``Namespace`` and runs the whole pipeline in-process. Importing it:

* removes the requirement that a ``macs3`` executable be on ``PATH`` (a library
  install is enough, and conda/pip layouts differ in where the binary lands);
* lets the caller see the return value and exceptions instead of parsing stderr;
* keeps one process, so a cluster job does not spend its accounting on a wrapper.

The previous version of ``easyepi peaks`` had two backends: this and a hand-written
Poisson caller. The Poisson caller is being retired for the reason the concordance
run in ``scripts/verify_peak_concordance.py`` demonstrated: with a global null it
found only 37% of MACS3's peaks and its own peaks were 1000 bp wide against MACS3's
191 bp. Keeping a second, worse caller behind a flag invites a user to pick it.

Backend policy
--------------
MACS3 is the only backend. If it is not importable the command fails with the exact
install line, rather than silently falling back to something less accurate.
"""

from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path
from typing import Any, Sequence

from easyepi.assays import AssayError, canonical_assay

#: Per-assay MACS3 options, keyed by the canonical label from
#: :mod:`easyepi.assays`. Stated explicitly so an ATAC run cannot inherit ChIP
#: assumptions; ``BAMPE`` is required for paired-end ATAC and CUT&Tag. CUT&Tag
#: and CUT&RUN resolve to ``chip`` because MACS3 needs the same format for them;
#: the shift is a `bam2bw` concern, not a peak-calling one.
ASSAY_MACS_OPTIONS: dict[str, dict[str, Any]] = {
    "atac": {"format": "BAMPE", "broad": False, "shift": 0, "extsize": None},
    "chip": {"format": "BAM", "broad": False, "shift": 0, "extsize": None},
    "histone": {"format": "BAM", "broad": True, "shift": 0, "extsize": None},
    "mnase": {"format": "BAM", "broad": False, "shift": 0, "extsize": None},
}

#: Canonical assays peak calling supports: MACS3 has no bisulfite mode.
ASSAY_CHOICES = tuple(ASSAY_MACS_OPTIONS)

_BAM_EXTS = (".bam",)


def expand_paths(values: Sequence[str] | None) -> list[str]:
    """Flatten repeated and comma-separated path arguments."""
    out: list[str] = []
    for raw in values or []:
        for piece in str(raw).split(","):
            path = piece.strip()
            if path:
                out.append(path)
    return out


def macs3_available() -> bool:
    """Whether MACS3 can be imported in this interpreter."""
    import importlib.util

    return importlib.util.find_spec("MACS3") is not None


def require_macs3():
    """Import the MACS3 entry point or raise with the interpreter-specific hint."""
    from easyepi.errors import missing_dependency, require_dependency

    try:
        from MACS3.Commands import callpeak_cmd
        return callpeak_cmd
    except ImportError as exc:
        raise missing_dependency(
            "MACS3", "peak calling in `easyepi peaks`") from exc


def run_peaks_command(args) -> int:
    """CLI driver for ``easyepi peaks``: MACS3 only, one output per input BAM."""
    from easyepi.errors import UsageError, require_file

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    bams = [require_file(p, purpose="peak calling input", suffixes=_BAM_EXTS)
            for p in expand_paths(getattr(args, "bam", None))]
    if not bams:
        raise UsageError(
            "no input: peak calling needs a BAM",
            example="easyepi peaks sample.bam -o peaks/ --assay atac")

    # `cuttag`, `cut&tag` and the rest resolve through the shared vocabulary, so a
    # label accepted by `bam2bw` one step earlier is accepted here too. MACS3 has
    # no bisulfite mode, so `methylation` is rejected with the supported list.
    raw_assay = getattr(args, "assay", None)
    if raw_assay:
        try:
            assay = canonical_assay(raw_assay)
        except AssayError:
            assay = None
        if assay not in ASSAY_MACS_OPTIONS:
            raise UsageError(
                f"unknown assay {raw_assay!r}; choose one of "
                f"{', '.join(ASSAY_CHOICES)}",
                example="easyepi peaks sample.bam -o peaks/ --assay atac")
    else:
        assay = ""
    opts = ASSAY_MACS_OPTIONS.get(assay, {"format": "BAM", "broad": False,
                                          "shift": 0, "extsize": None})

    broad_arg = getattr(args, "broad", None)
    broad = bool(broad_arg) if broad_arg is not None else bool(opts["broad"])
    fmt = (getattr(args, "format", None) or opts["format"]).upper()
    genome = getattr(args, "genome", None) or "hs"
    qvalue = float(getattr(args, "qvalue", None) or 0.01)
    extsize = getattr(args, "extsize", None) or opts["extsize"]
    shift = getattr(args, "shift_size", None)
    shift = int(shift) if shift is not None else int(opts["shift"])
    keep_dup = "auto" if not getattr(args, "keep_dup_all", False) else "all"
    extra = getattr(args, "macs_arg", None)
    if extra:
        print(f"  passthrough: {extra}")

    controls = [str(require_file(p, purpose="control BAM", suffixes=_BAM_EXTS))
                for p in expand_paths(getattr(args, "control", None))]

    print(f"[easyepi] peaks (MACS3): {len(bams)} treatment BAM(s)")
    print(f"  assay  : {assay or 'not stated'} -> format={fmt}, "
          f"{'broad' if broad else 'narrow'}")
    print(f"  genome : {genome}   q<={qvalue}")
    if controls:
        print(f"  control: {len(controls)} file(s)")

    from easyepi.macs_options import build_namespace as build_macs_namespace

    callpeak_cmd = require_macs3()

    written: list[dict[str, Any]] = []
    failures: list[str] = []
    for bam in bams:
        stem = bam.name[: -len(".bam")] if bam.name.endswith(".bam") else bam.stem
        try:
            # MACS3 writes <name>_peaks.<type>Peak into outdir; a per-sample
            # subdirectory keeps two runs with the same stem from colliding.
            sample_dir = out_dir / stem
            sample_dir.mkdir(parents=True, exist_ok=True)
            # The Namespace comes from MACS3's own argparse definition, read out of
            # its console script. A hand-written Namespace missed attributes that
            # opt_validate_callpeak compares to numbers, which surfaces as a
            # TypeError deep inside MACS3 instead of an argument error.
            options, provenance = build_macs_namespace(
                treatment=[str(bam)], out_dir=sample_dir, name=stem,
                genome=genome, fmt=fmt, broad=broad, control=controls,
                qvalue=qvalue, extsize=extsize, shift=shift, keep_dup=keep_dup,
                extra_args=getattr(args, "macs_arg", None))
            callpeak_cmd.run(options)
        except SystemExit as exc:                     # MACS3 calls sys.exit on error
            failures.append(f"{bam.name}: MACS3 aborted (exit {exc.code})")
            print(f"  [FAIL] {bam.name}: MACS3 aborted with exit code {exc.code}; "
                  "see the MACS3 log above")
            continue
        except Exception as exc:                      # noqa: BLE001
            failures.append(f"{bam.name}: {type(exc).__name__}: {exc}")
            print(f"  [FAIL] {bam.name}: {exc}")
            continue

        produced = sorted(sample_dir.glob(f"{stem}_peaks.*Peak"))
        if not produced:
            failures.append(f"{bam.name}: MACS3 produced no peak file")
            print(f"  [FAIL] {bam.name}: no peak file in {sample_dir}")
            continue
        target = produced[0]
        n = sum(1 for _ in target.open(encoding="utf-8", errors="replace"))
        written.append({"input": str(bam), "peaks": str(target), "n_peaks": n,
                        "format": fmt, "broad": broad,
                        "options_from": provenance["namespace_from"],
                        "macs_argv": provenance["argv"],
                        "macs_passthrough": provenance["passthrough"]})
        print(f"  [OK]   {bam.name} -> {target.relative_to(out_dir)} ({n} peaks)")

    manifest = out_dir / "peaks_manifest.json"
    manifest.write_text(json.dumps({
        "caller": "MACS3 (imported as a library)",
        "assay": assay or None, "format": fmt, "genome": genome,
        "broad": broad, "qvalue": qvalue, "extsize": extsize, "shift": shift,
        "outputs": written, "failures": failures,
    }, indent=2), encoding="utf-8")

    if failures:
        print(f"[easyepi] {len(failures)} BAM(s) failed")
        return 1
    if written:
        print(f"[easyepi] next: easyepi annotate {written[0]['peaks']} "
              f"--gtf genes.gtf -o annot/")
    return 0
