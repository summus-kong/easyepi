"""`easyepi bam2bw` : turn BAM files into BigWig tracks.

This is the first step of every analysis in the suite: peaks, QC, tracks and
states all read the BigWig (or the BAM) produced here. The CLI is explicit-argument
first by design:

    easyepi bam2bw sample.bam -o tracks/sample.bw --assay atac

Multiple BAMs in one call are supported, because a study normally has several:

    easyepi bam2bw rep1.bam rep2.bam -o tracks/ --assay atac

Sample names default to the file stem and can be overridden with
``--name BAM=NAME`` when the stem is not what the report should say.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from easyepi.assays import CANONICAL_ASSAYS, canonical_assay

#: Sensible defaults per assay, keyed by the canonical label from
#: :mod:`easyepi.assays`. Printed by the command and recorded in the metadata,
#: so a user can always see which shift and bin size were applied.
ASSAY_DEFAULTS: dict[str, dict[str, Any]] = {
    "atac": {"shift": "tn5", "bin_size": 10, "normalize": "cpm"},
    "chip": {"shift": "none", "bin_size": 25, "normalize": "cpm"},
    "histone": {"shift": "none", "bin_size": 25, "normalize": "cpm"},
    "mnase": {"shift": "none", "bin_size": 10, "normalize": "cpm"},
    "methylation": {"shift": "none", "bin_size": 100, "normalize": "raw"},
}

#: Overrides applied on top of ASSAY_DEFAULTS for Tn5 cut-site assays. CUT&Tag
#: and CUT&RUN map to `chip` in the shared vocabulary because MACS3 treats them
#: identically, but they measure a cut site, not read-body depth: handling them
#: as plain `chip` silently dropped the +4/-5 shift that `--assay cut&tag` used
#: to apply. Matching is by canonical label plus a separator-insensitive test, so
#: `cut_tag` and `Cut Tag` are covered without listing every variant.
_CUTSITE_DEFAULTS: dict[str, Any] = {"shift": "tn5", "bin_size": 10}


def _is_cutsite_assay(assay: str) -> bool:
    """Whether this assay measures a Tn5/MNase cut site rather than read depth."""
    flat = str(assay).strip().lower()
    for ch in ("_", " ", "-", "&"):
        flat = flat.replace(ch, "")
    return flat in {"cuttag", "cutrun"}

ASSAY_CHOICES = CANONICAL_ASSAYS


def resolve_defaults(assay: str | None) -> dict[str, Any]:
    """Per-assay defaults, with a neutral fallback when no assay is stated.

    The fallback is deliberately conservative (no shift, raw units): it changes
    nothing about the signal, so a wrong assumption here cannot silently alter the
    values. The command prints which defaults it used.
    """
    if assay is None:
        return {"shift": "none", "bin_size": 25, "normalize": "raw",
                "source": "fallback"}
    # The vocabulary lives in easyepi.assays, so `cuttag`, `cut-tag`, `cut&tag`
    # and every other accepted spelling resolve to one of the five keys here.
    # AssayError is a ValueError, which is what callers of this function catch.
    key = canonical_assay(assay)
    defaults = dict(ASSAY_DEFAULTS[key])
    if _is_cutsite_assay(assay):
        # CUT&Tag/CUT&RUN share the `chip` metric set but not its signal handling.
        defaults.update(_CUTSITE_DEFAULTS)
    return {**defaults, "source": "assay"}


def _sample_names(bams: Sequence[str], overrides: Sequence[str] | None) -> dict[str, str]:
    """Map BAM path -> sample name, honouring ``--name BAM=NAME`` overrides."""
    names = {bam: Path(bam).stem for bam in bams}
    for spec in overrides or []:
        if "=" not in spec:
            raise ValueError(f"--name expects BAM=NAME, got {spec!r}")
        bam, name = spec.split("=", 1)
        bam, name = bam.strip(), name.strip()
        if bam not in names:
            raise ValueError(
                f"--name refers to {bam!r}, which is not among the BAM arguments")
        if not name:
            raise ValueError(f"--name gave an empty name for {bam!r}")
        names[bam] = name
    dupes = sorted({n for n in names.values() if list(names.values()).count(n) > 1})
    if dupes:
        raise ValueError(
            "duplicate sample name(s): " + ", ".join(dupes)
            + "; disambiguate with --name BAM=NAME")
    return names


def run_bam2bw(args) -> int:
    """CLI driver: one or more BAMs -> one BigWig each, plus metadata."""
    import json

    from easyepi.bam2bigwig import bam_to_bigwig

    bams = [str(b) for b in (args.bam or [])]
    if not bams:
        raise ValueError("at least one BAM is required")

    from easyepi.errors import require_file

    # Checked here so a typo fails before any work starts, with the resolved path
    # and the nearest existing filename in the message.
    for bam in bams:
        require_file(bam, purpose="BAM input", suffixes=(".bam",))

    out_dir = Path(args.out)
    if args.output_file and len(bams) > 1:
        raise ValueError(
            "--output-file names a single BigWig but several BAMs were given; "
            "use --out DIR to write one track per sample")

    defaults = resolve_defaults(getattr(args, "assay", None))
    bin_size = int(args.bin_size if args.bin_size is not None else defaults["bin_size"])
    shift = args.shift if args.shift is not None else defaults["shift"]
    normalize = args.normalize if args.normalize is not None else defaults["normalize"]

    names = _sample_names(bams, getattr(args, "name", None))

    from easyepi.coverage_parallel import available_workers
    jobs = available_workers(getattr(args, "jobs", None))

    print(f"[easyepi] bam2bw: {len(bams)} BAM(s)")
    print(f"  assay      : {getattr(args, 'assay', None) or 'not stated'}"
          f" (defaults from {'assay table' if defaults['source'] == 'assay' else 'fallback'})")
    print(f"  shift      : {shift}"
          + ("  (+4/-5 Tn5 cut sites)" if shift == "tn5" else "  (read body)"))
    print(f"  bin size   : {bin_size} bp")
    print(f"  jobs       : {jobs} process(es)")
    print(f"  normalize  : {normalize}"
          + (f"  (scale applied to the values)" if normalize != "raw" else ""))

    written = []
    failures = []
    for bam in bams:
        name = names[bam]
        if args.output_file:
            target = Path(args.output_file)
        else:
            out_dir.mkdir(parents=True, exist_ok=True)
            target = out_dir / f"{name}.bw"
        try:
            meta = bam_to_bigwig(
                bam, target,
                chrom=getattr(args, "chrom", None),
                bin_size=bin_size, normalize=normalize, shift=shift,
                min_mapq=int(args.min_mapq), max_depth=args.max_depth,
                jobs=jobs)
        except Exception as exc:                        # noqa: BLE001
            failures.append(f"{name}: {type(exc).__name__}: {exc}")
            print(f"  [FAIL] {name}: {exc}")
            continue
        written.append({"sample": name, "bigwig": str(target), "metadata": meta})
        mb = target.stat().st_size / 1e6
        print(f"  [OK]   {name} -> {target.name} ({mb:.1f} MB, "
              f"{meta['mapped_reads']['mapped']:,} mapped reads)")

    if args.manifest:
        Path(args.manifest).write_text(
            json.dumps({"tracks": written, "failures": failures},
                       indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"  manifest   : {args.manifest}")
    if failures:
        print(f"[easyepi] {len(failures)} BAM(s) failed; see the messages above")
        return 1
    if written:
        print(f"[easyepi] next: easyepi peaks {' '.join(w['bigwig'] for w in written[:2])}"
              f" -o peaks/ --assay {getattr(args, 'assay', None) or 'chip'}")
    return 0
