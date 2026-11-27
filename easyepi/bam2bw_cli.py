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

from easyepi.assays import CANONICAL_ASSAYS, canonical_assay, describe_assay

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

#: Overrides for assays that measure a cut site rather than read-body depth: a
#: fine bin is wanted in both cases. The *shift* is not shared, because the two
#: enzymes cut differently:
#:
#: * **CUT&Tag** uses hyperactive Tn5, whose dimer inserts at the two ends of its
#:   footprint -- the standard +4/-5 strand offset is what recovers the true
#:   insertion site.
#: * **CUT&RUN** uses a Protein A/G-MNase fusion. MNase cuts the linker DNA around
#:   the targeted factor and leaves a protected footprint; there is no Tn5
#:   footprint and therefore no +4/-5 offset to correct. Applying one moved every
#:   cut site by 4-5 bp and was justified only by the fact that MACS3 treats the
#:   two assays the same way, which is a peak-calling detail, not a signal-model
#:   one.
#:
#: Both map to `chip` in the shared vocabulary (MACS3 argues they call peaks
#: alike), so the distinction has to be made here rather than inherited.
_CUTSITE_BIN_SIZE = 10

#: Kept as a name other code may import, but no longer the decision: the signal
#: model lives in `easyepi.assays` now, so a caller cannot lose the shift by
#: normalising the assay first. See `assays.SIGNAL_OVERRIDES`.
_TN5_SHIFT_ASSAYS = frozenset({"cuttag"})


def _assay_flat(assay: str) -> str:
    """Lower-case an assay label with every separator removed.

    So `cut_tag`, `Cut Tag` and `CUT&TAG` compare equal without listing each
    variant.
    """
    flat = str(assay).strip().lower()
    for ch in ("_", " ", "-", "&", "."):
        flat = flat.replace(ch, "")
    return flat


def _is_cutsite_assay(assay: str) -> bool:
    """Whether this assay measures a cut site rather than read-body depth.

    Both CUT&Tag and CUT&RUN qualify for a fine bin. This only decides the bin
    size; the strand shift is decided by :func:`_uses_tn5_shift`.

    Resolves the assay rather than matching its spelling, so the answer cannot
    drift from `assays.py`. Note that this still depends on `assay` carrying a
    non-group spelling: `_is_cutsite_assay("chip")` is False, because a bare group
    label genuinely does not say whether the signal is a cut site.
    """
    from easyepi.assays import describe_assay
    return describe_assay(assay).is_cutsite


def _uses_tn5_shift(assay: str) -> bool:
    """Whether the +4/-5 Tn5 offset applies to this assay.

    True for CUT&Tag, false for CUT&RUN: MNase has no Tn5 footprint to correct,
    so shifting a CUT&RUN library moves its cut sites off the truth.

    Delegates to the shared vocabulary. This used to match the raw spelling here,
    which meant `_uses_tn5_shift(canonical_assay("cuttag"))` returned False -- the
    shift silently disappeared for any caller that normalised first. Callers that
    need both the group and the signal must use `resolve_assay` once and read both
    fields from it, rather than round-tripping through `canonical_assay`.
    """
    from easyepi.assays import describe_assay
    return describe_assay(assay).uses_tn5_shift

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
    # Resolve once: the group picks the table, the signal class picks the shift.
    # Resolving twice from a round-tripped group label is what lost the Tn5 offset.
    resolved = describe_assay(assay)
    defaults = dict(ASSAY_DEFAULTS[resolved.group])
    if resolved.is_cutsite:
        # A cut-site assay wants fine bins. The shift is decided separately: only
        # a Tn5-based library has the +4/-5 footprint to correct.
        defaults["bin_size"] = _CUTSITE_BIN_SIZE
        defaults["shift"] = "tn5" if _uses_tn5_shift(assay) else "none"
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
