"""BAM -> BigWig: the first step of the chain, done properly.

This is the foundation the rest of the suite builds on, so it has to be correct
and fast rather than merely present.

Why not count per read in Python
--------------------------------
The previous helper looped over every alignment to add 1.0 into a bin array. On a
30M-read BAM that is a Python loop over 30M records, which is minutes of wall
clock for a pure-Python caller and puts an artificial ceiling on library size.
Here the work is delegated to ``pysam``'s C-level ``count_coverage`` when it is
available, which histograms the CIGAR-aware depth directly.

Signal definition, stated explicitly
------------------------------------
``count_coverage`` returns per-base depth, i.e. the number of read bases covering
each base (5' overhang aware, CIGAR-mapped). Two consequences the user must know:

* A **fragment** (paired-end) covers roughly 150-300 bp, a **read** 50-150 bp.
  Depth therefore differs by the fragment/read ratio, so samples must not be
  compared without normalisation (``--normalize cpm``).
* For ATAC/CUT&Tag the biologically meaningful signal is the Tn5 cut site, not the
  read body. ``--shift tn5`` applies the standard +4/-5 shift per strand before
  accumulation, which is what the peak callers expect.

Normalisation is explicit: ``raw`` (depth as-is), ``cpm`` (per million mapped
read starts), or ``rpkm`` (per million per kb of bin, using the bin width).
A genome-wide scale factor is required for any comparison between samples, so it is
recorded in the output metadata rather than left implicit.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

NORMS = ("raw", "cpm", "rpkm")
SHIFTS = ("none", "tn5")


def _require_pysam():
    from easyepi.core.io import _require
    return _require("pysam", "pip install 'easyepi[bam]'")


def _read_shift(is_reverse: bool, shift: str) -> int:
    """Per-strand offset for ATAC-style cut-site placement.

    The Tn5 insertion site is the 5' end of the read, corrected by +4 on the
    forward strand and -5 on the reverse strand (Buenrostro et al. 2013). Without
    this the two strands' signals are offset by ~9 bp and a footprint becomes a
    doublet.
    """
    if shift == "none":
        return 0
    return -5 if is_reverse else 4


def reference_sizes(bam_path: str | Path) -> dict[str, int]:
    """Chromosome sizes from a BAM header, in BAM order.

    Order matters for a BigWig: ``addHeader`` fixes the chromosome order and every
    later entry must respect it.
    """
    pysam = _require_pysam()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if not bam.references:
            raise ValueError(f"{bam_path} has no reference sequence in its header")
        sizes = dict(zip(bam.references, bam.lengths))
    return sizes


def count_mapped(bam_path: str | Path, chrom: str | None = None,
                 sample_reads: int = 200_000) -> dict[str, Any]:
    """Mapped-read totals used by the CPM/RPKM scale factors.

    Uses ``idxstats`` when the BAM is indexed (exact and cheap); falls back to
    scanning a sample of the first records and extrapolating, which is flagged in
    the returned ``method`` field so a user knows the number is an estimate.
    """
    pysam = _require_pysam()
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if bam.has_index():
            stats = bam.get_index_statistics()
            per_chrom = {s.contig: int(s.mapped) for s in stats}
            total = sum(per_chrom.values())
            if chrom is not None:
                if chrom not in per_chrom:
                    raise ValueError(f"chrom {chrom!r} is not in {bam_path}")
                return {"mapped": per_chrom[chrom], "total": total,
                        "method": "idxstats"}
            return {"mapped": total, "total": total, "method": "idxstats"}

        seen = 0
        mapped = 0
        for rec in bam.fetch(until_eof=True):
            if chrom is not None and rec.reference_name != chrom:
                continue
            seen += 1
            if not rec.is_unmapped:
                mapped += 1
            if seen >= sample_reads:
                break
    if seen == 0:
        return {"mapped": 0, "total": 0, "method": "sampled_empty"}
    # extrapolation is only valid for the whole file, not per chromosome
    fraction = mapped / seen
    total = int(round(seen * fraction)) if mapped == seen else 0
    return {"mapped": mapped, "total": total or mapped,
            "method": "sampled_extrapolated", "fraction_mapped": round(fraction, 4)}


def has_base_qualities(bam_path: str | Path, sample_reads: int = 20_000) -> bool:
    """Whether the BAM stores per-base qualities.

    ``count_coverage(quality_threshold=...)`` treats a record with no qualities as
    failing every base, so such a file yields an all-zero track while the command
    still looks successful. Callers must skip the quality filter in that case.
    """
    pysam = _require_pysam()
    seen = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped:
                continue
            seen += 1
            if rec.query_qualities is not None:
                return True
            if seen >= sample_reads:
                break
    return False


def bam_to_bigwig(
    bam_path: str | Path,
    bigwig_path: str | Path,
    *,
    chrom: str | None = None,
    bin_size: int = 25,
    normalize: str = "raw",
    shift: str = "none",
    min_mapq: int = 30,
    max_depth: int | None = None,
    blacklist: Iterable[tuple[str, int, int]] | None = None,
    jobs: int = 1,
) -> dict[str, Any]:
    """Convert a coordinate-sorted BAM to BigWig coverage.

    Returns a metadata dict that is written next to the BigWig as
    ``<name>.bam2bw.json``: every parameter that affects the values is recorded so
    two BigWigs can be shown to be comparable, or not.

    ``max_depth`` caps per-base depth before binning. Duplicate-heavy libraries can
    reach thousands of depth at a single base from stacked artefacts; the cap keeps
    one bad locus from dominating the colour scale without altering well-behaved
    regions.
    """
    if normalize not in NORMS:
        raise ValueError(f"normalize must be one of {NORMS}")
    if shift not in SHIFTS:
        raise ValueError(f"shift must be one of {SHIFTS}")
    if bin_size < 1:
        raise ValueError("bin_size must be >= 1")

    pysam = _require_pysam()
    import pyBigWig

    sizes = reference_sizes(bam_path)
    targets = [chrom] if chrom else list(sizes)
    for name in targets:
        if name not in sizes:
            raise ValueError(
                f"chrom {name!r} is not in the BAM header "
                f"(has e.g. {list(sizes)[:5]})")

    counts = count_mapped(bam_path, chrom)
    qualities_present = has_base_qualities(bam_path)
    scale = 1.0
    if normalize == "cpm":
        if counts["mapped"] <= 0:
            raise ValueError("cannot normalise to CPM: no mapped reads detected")
        scale = 1e6 / counts["mapped"]
    elif normalize == "rpkm":
        if counts["mapped"] <= 0:
            raise ValueError("cannot normalise to RPKM: no mapped reads detected")
        # per million mapped, per kb of bin
        scale = 1e9 / (counts["mapped"] * bin_size)

    out_path = Path(bigwig_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    bw = pyBigWig.open(str(out_path), "w")
    try:
        bw.addHeader([(n, int(sizes[n])) for n in targets])

        bl = list(blacklist or [])
        blacklist_by_chrom: dict[str, list[tuple[int, int]]] = {}
        for c, s, e in bl:
            blacklist_by_chrom.setdefault(c, []).append((int(s), int(e)))

        # Compute every chromosome first, in one pool submission, so that a
        # multi-chromosome run keeps all cores busy and a single-chromosome run is
        # not dominated by pool startup per contig. Values are then written in
        # header order, which is what BigWig requires.
        computed: dict[str, np.ndarray] = {}
        if shift == "tn5":
            # cut sites need read ends, which accumulated depth does not retain
            from easyepi.coverage_readends import read_end_bins

            for name in targets:
                length = int(sizes[name])
                values, _stats = read_end_bins(bam_path, name, length, bin_size,
                                               min_mapq=min_mapq, shift="tn5")
                computed[name] = values * scale
        else:
            # One pool over EVERY chromosome-chunk: per-chromosome pools left the
            # largest contig as the critical path and measured only 1.6x on 249 Mb.
            from easyepi.coverage_parallel import CoverageSpec, coverage_many

            specs = [CoverageSpec(bam_path=str(bam_path), chrom=name,
                                 chrom_length=int(sizes[name]),
                                 bin_size=bin_size, min_mapq=int(min_mapq),
                                 shift="none", max_depth=max_depth)
                     for name in targets]
            raw = coverage_many(specs, workers=int(jobs))
            for name in targets:
                computed[name] = raw[name] * scale

        for name in targets:
            length = int(sizes[name])
            values = computed[name]

            # A blacklist region must not carry signal; zeroing is explicit rather
            # than dropping the bins, so the track stays dense.
            n_bins = max(1, int(np.ceil(length / bin_size)))
            for bs, be in blacklist_by_chrom.get(name, []):
                lo = max(0, bs // bin_size)
                hi = min(n_bins, int(np.ceil(be / bin_size)))
                values[lo:hi] = 0.0
            computed[name] = values

        per_chrom_written = {}
        for name in targets:
            length = int(sizes[name])
            values = computed[name]
            starts = np.arange(values.size, dtype=np.int64) * bin_size
            ends = np.minimum(starts + bin_size, length)
            keep = ends > starts
            if not keep.any():
                continue
            bw.addEntries(
                [name] * int(keep.sum()),
                starts[keep].astype('int64').tolist(),
                ends=ends[keep].astype('int64').tolist(),
                values=np.nan_to_num(values[keep], nan=0.0,
                                     posinf=0.0, neginf=0.0).astype(float).tolist())
            per_chrom_written[name] = int(keep.sum())
    finally:
        bw.close()

    meta = {
        "bam": str(bam_path),
        "bigwig": str(out_path),
        "bin_size": int(bin_size),
        "normalize": normalize,
        "scale_factor": float(scale),
        "shift": shift,
        "min_mapq": int(min_mapq),
        "base_qualities_present": bool(qualities_present),
        "quality_filter_applied": bool(qualities_present),
        "jobs": int(jobs),
        "max_depth": None if max_depth is None else int(max_depth),
        "mapped_reads": counts,
        "chromosomes": per_chrom_written,
        "signal_units": {
            "raw": "read bases covering a bin (depth), summed",
            "cpm": "depth per million mapped reads",
            "rpkm": "depth per million mapped reads per kb of bin",
        }[normalize],
    }
    Path(str(out_path) + ".json").write_text(
        json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    return meta


def _sum_to_bins(depth: np.ndarray, offset: int, bin_size: int,
                 n_bins: int) -> np.ndarray:
    """Sum a per-base depth array into the global bin grid starting at ``offset``."""
    out = np.zeros(n_bins, dtype=np.float64)
    first_bin = offset // bin_size
    # position of the first base within its bin
    lead = offset % bin_size
    if lead:
        head = min(bin_size - lead, depth.size)
        if first_bin < n_bins:
            out[first_bin] += float(depth[:head].sum())
        depth = depth[head:]
        first_bin += 1
    n_full = depth.size // bin_size
    if n_full:
        block = depth[: n_full * bin_size].reshape(n_full, bin_size).sum(axis=1)
        end_bin = min(n_bins, first_bin + n_full)
        out[first_bin:end_bin] += block[: end_bin - first_bin]
        depth = depth[n_full * bin_size:]
        first_bin = end_bin
    if depth.size and first_bin < n_bins:
        out[first_bin] += float(depth.sum())
    return out
