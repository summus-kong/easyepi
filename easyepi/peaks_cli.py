"""`easyepi peaks` -- peak calling from a BAM or a BigWig.

The subcommand used to raise unconditionally:

    raise ValueError("Use python -m easyepi.specific.atac.peaks or
    easyepi.specific.chip.narrow_peaks with MACS3; the legacy stub has been disabled")

while being the third entry in ``--help``, so the most visible command in the suite
failed every time and pointed at a different set of module paths. The README and the
per-command pages also advertised ``--bigwig``, ``--assay`` and ``--broad``, none of
which the parser accepted.

This module supplies what those docs promised, on top of the functions that exist:

    call_peaks_from_bam(bam_path, chrom, bin_size=200, scale=1.0,
                        threshold=5.0, max_qvalue=0.01, t5_shift=True) -> PeakSet
    poisson_peaks(signal, chrom, start, bin_size, ...) -> PeakSet
    peaks_to_bed(peaks, out_path, format="narrowPeak")
    bigwig_to_signal(path, chrom, start, end, bin_size=50) -> np.ndarray

The caller is per-chromosome and has no broad/narrow switch, so ``--broad`` selects
the output BED format and the manifest field, not the statistical test. That is
stated in the printed summary rather than implied.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Sequence

#: Per-assay calling defaults, kept explicit so a ChIP run cannot inherit ATAC's Tn5
#: assumption.
ASSAY_PEAK_DEFAULTS: dict[str, dict[str, Any]] = {
    "atac": {"shift": "tn5", "genome": "hs", "broad": False},
    "cuttag": {"shift": "tn5", "genome": "hs", "broad": False},
    "cutrun": {"shift": "tn5", "genome": "hs", "broad": False},
    "chip": {"shift": "none", "genome": "hs", "broad": False},
    "histone": {"shift": "none", "genome": "hs", "broad": True},
    "mnase": {"shift": "none", "genome": "hs", "broad": False},
}

ASSAY_CHOICES = tuple(ASSAY_PEAK_DEFAULTS)

_BAM_EXTS = (".bam",)
_BW_EXTS = (".bw", ".bigWig", ".bigwig")


def expand_paths(values: Sequence[str] | None) -> list[str]:
    """Flatten repeated and comma-separated path arguments."""
    out: list[str] = []
    for raw in values or []:
        for piece in str(raw).split(","):
            path = piece.strip()
            if path:
                out.append(path)
    return out


def resolve_inputs(args) -> tuple[list[tuple[str, Path]], list[str]]:
    """Return ``[(kind, path)]`` plus notes, rejecting missing paths early."""
    from easyepi.errors import UsageError, require_file

    notes: list[str] = []
    items: list[tuple[str, Path]] = []

    for path in expand_paths(getattr(args, "bam", None)):
        items.append(("bam", require_file(path, purpose="peak calling input",
                                          suffixes=_BAM_EXTS)))
    for path in expand_paths(getattr(args, "bigwig", None)):
        items.append(("bigwig", require_file(path, purpose="peak calling input",
                                             suffixes=_BW_EXTS)))

    if not items:
        raise UsageError(
            "no input: pass --bam or --bigwig",
            example="easyepi peaks track.bw -o peaks/ --assay atac")

    seen: set[str] = set()
    unique: list[tuple[str, Path]] = []
    for kind, path in items:
        key = str(path)
        if key in seen:
            notes.append(f"duplicate input ignored: {path.name}")
            continue
        seen.add(key)
        unique.append((kind, path))
    return unique, notes


def caller_available() -> str | None:
    for name in ("macs3", "macs2"):
        found = shutil.which(name)
        if found:
            return found
    return None


def chromosome_names(kind: str, path: Path) -> list[str]:
    if kind == "bam":
        import pysam

        with pysam.AlignmentFile(str(path), "rb") as bam:
            return list(bam.references)
    import pyBigWig

    bw = pyBigWig.open(str(path))
    try:
        return list(bw.chroms())
    finally:
        bw.close()


def chromosome_length(kind: str, path: Path, chrom: str) -> int:
    if kind == "bam":
        import pysam

        with pysam.AlignmentFile(str(path), "rb") as bam:
            return int(bam.get_reference_length(chrom))
    import pyBigWig

    bw = pyBigWig.open(str(path))
    try:
        return int(bw.chroms()[chrom])
    finally:
        bw.close()


def call_native(kind: str, path: Path, target: Path, *, assay_defaults: dict,
                chroms: list[str] | None, bin_size: int, threshold: float,
                max_qvalue: float, broad: bool) -> tuple[int, list[str]]:
    """Run the built-in caller over the requested chromosomes, one BED out.

    Per-chromosome results are written with the shipped ``peaks_to_bed``; the first
    chromosome creates the file and the rest are appended by writing to a temporary
    path and concatenating, because ``peaks_to_bed`` takes a path rather than a
    handle.
    """
    from easyepi.peaks import (bigwig_to_signal, call_peaks_from_bam,
                              peaks_to_bed, poisson_peaks)

    notes: list[str] = []
    available = chromosome_names(kind, path)
    requested = list(chroms or available)
    missing = [c for c in requested if c not in available]
    if missing:
        notes.append(f"absent from {path.name}, skipped: {', '.join(missing)}")
    targets = [c for c in requested if c in available]
    if not targets:
        raise ValueError(f"{path.name} contains none of the requested chromosomes")

    use_t5 = assay_defaults.get("shift") == "tn5"
    if kind != "bam" and use_t5:
        notes.append("BigWig input carries no read ends, so the Tn5 shift is not applied")

    fmt = "broadPeak" if broad else "narrowPeak"
    total = 0
    parts: list[Path] = []
    for chrom in targets:
        if kind == "bam":
            peak_set = call_peaks_from_bam(str(path), chrom, bin_size=bin_size,
                                           threshold=threshold,
                                           max_qvalue=max_qvalue, t5_shift=use_t5)
        else:
            length = chromosome_length(kind, path, chrom)
            signal = bigwig_to_signal(path, chrom, 0, length, bin_size=bin_size)
            peak_set = poisson_peaks(signal, chrom, 0, bin_size,
                                     threshold=threshold, max_qvalue=max_qvalue)
        total += int(getattr(peak_set, "n_peaks", 0))
        part = target.with_name(f".{target.name}.{chrom}.part")
        peaks_to_bed(peak_set, part, format=fmt)
        parts.append(part)

    with target.open("w", encoding="utf-8") as out:
        for part in parts:
            out.write(part.read_text(encoding="utf-8"))
            part.unlink()
    return total, notes


def run_peaks_command(args) -> int:
    """CLI driver for ``easyepi peaks``."""
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    items, notes = resolve_inputs(args)
    assay = (getattr(args, "assay", None) or "").strip().lower()
    if assay and assay not in ASSAY_PEAK_DEFAULTS:
        from easyepi.errors import UsageError

        raise UsageError(
            f"unknown assay {assay!r}; choose one of {', '.join(ASSAY_CHOICES)}",
            example="easyepi peaks track.bw -o peaks/ --assay atac")
    defaults = ASSAY_PEAK_DEFAULTS.get(assay, {"shift": "none", "genome": "hs",
                                               "broad": False})

    broad_arg = getattr(args, "broad", None)
    broad = bool(broad_arg) if broad_arg is not None else bool(defaults["broad"])

    backend = (getattr(args, "backend", None) or "auto").lower()
    macs = caller_available()
    if backend == "auto":
        backend = "macs3" if macs else "native"
    if backend == "macs3" and not macs:
        from easyepi.errors import UsageError

        raise UsageError(
            "--backend macs3 was requested but neither macs3 nor macs2 is on PATH",
            example="easyepi peaks track.bw -o peaks/ --backend native")

    chroms = None
    if getattr(args, "chrom", None):
        chroms = [c.strip() for c in str(args.chrom).split(",") if c.strip()]

    bin_size = int(getattr(args, "bin_size", None) or 200)
    threshold = float(getattr(args, "threshold", None) or 5.0)
    max_qvalue = float(getattr(args, "max_qvalue", None) or 0.01)

    print(f"[easyepi] peaks: {len(items)} input(s)")
    print(f"  assay    : {assay or 'not stated'}")
    print(f"  output    : {'broadPeak' if broad else 'narrowPeak'}"
          + (" (from the assay default)" if broad_arg is None and broad else ""))
    print(f"  backend   : {backend}"
          + (f" ({macs})" if macs and backend == "macs3" else ""))
    if backend == "native":
        print(f"  parameters: bin={bin_size} threshold={threshold} q<={max_qvalue}")
    for note in notes:
        print(f"  note: {note}")

    written: list[dict[str, Any]] = []
    failures: list[str] = []
    for kind, path in items:
        stem = path.name
        for ext in _BAM_EXTS + _BW_EXTS:
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        target = out_dir / f"{stem}.{'broadPeak' if broad else 'narrowPeak'}"
        try:
            if backend == "macs3" and kind == "bam" and macs:
                import logging
                from types import SimpleNamespace

                from easyepi._peak_calling import macs3_command
                from easyepi.logging_utils import run_external

                sample = SimpleNamespace(sample_id=stem, bam=str(path), bigwig=None)
                command = macs3_command(
                    macs, sample, out_dir, genome_size=defaults["genome"],
                    assay="ATAC" if assay in ("atac", "cuttag", "cutrun") else "CHIP",
                    broad=broad)
                run_external(command, logging.getLogger("easyepi.peaks"))
                produced = sorted(out_dir.glob(f"{stem}_peaks.*Peak"))
                if produced:
                    produced[0].rename(target)
                n = sum(1 for _ in target.open()) if target.exists() else 0
            else:
                if backend == "macs3" and kind != "bam":
                    print("  note: MACS3 needs a BAM; using the native caller")
                n, extra = call_native(kind, path, target, assay_defaults=defaults,
                                       chroms=chroms, bin_size=bin_size,
                                       threshold=threshold, max_qvalue=max_qvalue,
                                       broad=broad)
                for note in extra:
                    print(f"  note: {note}")
        except Exception as exc:                       # noqa: BLE001
            failures.append(f"{path.name}: {type(exc).__name__}: {exc}")
            print(f"  [FAIL] {path.name}: {exc}")
            if getattr(args, "verbose", False):
                import traceback

                traceback.print_exc()
            continue
        written.append({"input": str(path), "kind": kind, "peaks": str(target),
                        "n_peaks": n})
        print(f"  [OK]   {path.name} -> {target.name} ({n} peaks)")

    manifest = out_dir / "peaks_manifest.json"
    manifest.write_text(json.dumps({
        "assay": assay or None, "format": "broadPeak" if broad else "narrowPeak",
        "backend": backend, "shift": defaults["shift"], "chroms": chroms,
        "bin_size": bin_size, "threshold": threshold, "max_qvalue": max_qvalue,
        "outputs": written, "failures": failures,
        "note": "the built-in caller scores significance per bin; --broad selects the "
                "BED format, not a different test",
    }, indent=2), encoding="utf-8")

    if failures:
        print(f"[easyepi] {len(failures)} input(s) failed")
        return 1
    if written:
        print(f"[easyepi] next: easyepi annotate {written[0]['peaks']} "
              f"--gtf genes.gtf -o annot/")
    return 0
