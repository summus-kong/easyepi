"""Nucleosome positioning from MNase-seq fragment ends.

Why fragment ends, not read starts
----------------------------------
MNase-digested chromatin yields fragments whose *ends* are the micrococcal nuclease
cut sites. Nucleosome dyads sit between paired cut sites and are protected from
cutting, so the informative signal is (a) the cut-site density and (b) the
inferred fragment midpoints. Using pileup read starts would put the signal at the
fragment edge and shift every call by half a fragment.

Method (published, no fitted parameters)
----------------------------------------
1. **Cut sites** — each fragment contributes one cut at each end. Forward reads give
   ``start`` and reverse reads give ``end``; a paired record gives both plus the
   dyad at the fragment midpoint.
2. **Dyad density** — fragment midpoints binned at ``bin_size``; this is the
   standard nucleosome-centre track used by NucleoATAC/  
   MINCE-style analyses (Schep et al. 2015, Genome Res 25:1757; Chereji et al.
   2018, Bioinformatics 34:i296 for the dyad-density framing).
3. **Calls** — local maxima of the smoothed dyad density that exceed a robust
   background (median + k * MAD) and are separated by at least ``min_spacing`` bp.
   MAD is used instead of SD because dyad density is dominated by a few deep
   nucleosomes; an SD threshold would either admit noise or reject real arrays.
4. **Fuzziness** — per call, the bp window that contains a given fraction of the
   local dyad density; bulk MNase-seq gives well-positioned calls (narrow) while
   soft-digested samples give fuzzier ones. Reported per call, not interpreted.

What this does not do
---------------------
It does not resolve overlapping nucleosomes, does not model the digestion
kinetics, and does not call positions from ATAC-seq (which reports accessibility,
not protection). Those need NucleoATAC or a template-based algorithm with a
fitted linker-length model; see the module docstring in
``easyepi/specific/mnase/positioning.py``.

References
----------
* Schep AN et al. (2015) NucleoATAC. Genome Res 25:1757-1768.
* Chereji RV, Bryson TD, Henikoff S (2018) Quantitative MNase-seq. Bioinformatics 34:i296.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Curve primitives
# ---------------------------------------------------------------------------
def _smooth(values: np.ndarray, window: int) -> np.ndarray:
    """Moving average with edge-safe padding.

    Smoothing is applied before peak calling because dyad density is a sparse
    count signal: a raw argmax would lock onto a single fragment, not a positioned
    nucleosome. The window is in bins, so the caller controls the bp scale.
    """
    sig = np.asarray(values, dtype=np.float64)
    window = max(1, int(window))
    if window == 1 or sig.size == 0:
        return sig.copy()
    if window % 2 == 0:
        window += 1
    pad = window // 2
    padded = np.pad(sig, pad, mode="edge")
    kernel = np.ones(window, dtype=np.float64) / window
    return np.convolve(padded, kernel, mode="valid")[: sig.size]


def robust_background(values: np.ndarray) -> tuple[float, float]:
    """Median and scaled MAD of a non-negative signal.

    Returns ``(median, sigma)`` where ``sigma`` is the MAD scaled by 1.4826 so it
    is comparable to a standard deviation. MAD is used because a handful of deeply
    occupied nucleosomes otherwise inflate the SD and hide real calls.
    """
    sig = np.asarray(values, dtype=np.float64)
    if sig.size == 0:
        return 0.0, 0.0
    median = float(np.median(sig))
    mad = float(np.median(np.abs(sig - median)))
    return median, 1.4826 * mad


def call_positions(dyad_density: np.ndarray, bin_size: int, smooth_bins: int = 3,
                   mad_k: float = 2.0, min_spacing: int = 100,
                   min_height: float = 0.0,
                   min_relative_height: float = 0.05) -> list[dict[str, Any]]:
    """Call nucleosome positions as local maxima above a robust background.

    ``dyad_density`` is the per-bin fragment-midpoint count. A bin is called when
    it is a strict local maximum of the smoothed curve, sits above
    ``median + mad_k * MAD`` (and above ``min_height`` when that is larger), and is
    at least ``min_spacing`` bp from the previously accepted call. Calls are
    returned left to right sorted by coordinate.
    """
    sig = np.asarray(dyad_density, dtype=np.float64)
    if sig.size == 0:
        raise ValueError("dyad_density is empty")
    if np.any(sig < 0):
        raise ValueError("dyad_density must be non-negative")
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    if min_spacing < 0:
        raise ValueError("min_spacing must be non-negative")
    if not (0.0 <= min_relative_height < 1.0):
        raise ValueError("min_relative_height must be in [0, 1)")

    smoothed = _smooth(sig, smooth_bins)
    median, sigma = robust_background(smoothed)
    peak = float(smoothed.max()) if smoothed.size else 0.0
    relative_floor = min_relative_height * peak if min_relative_height > 0 else 0.0
    threshold = max(median + mad_k * sigma, float(min_height), relative_floor)

    candidates: list[int] = []
    for i in range(smoothed.size):
        val = smoothed[i]
        if val <= threshold:
            continue
        left = smoothed[i - 1] if i > 0 else -np.inf
        right = smoothed[i + 1] if i + 1 < smoothed.size else -np.inf
        # strict on the right so a flat plateau yields one call, not many
        if val > left and val >= right:
            candidates.append(i)

    # Strongest first, then enforce spacing so a broad plateau cannot emit
    # adjacent calls that are really one nucleosome.
    candidates.sort(key=lambda i: (-smoothed[i], i))
    accepted: list[int] = []
    for idx in candidates:
        if all(abs(idx - other) * bin_size >= min_spacing for other in accepted):
            accepted.append(idx)
    accepted.sort()

    total = float(sig.sum())
    calls = []
    for idx in accepted:
        calls.append({
            "bin": int(idx),
            "position": int(idx * bin_size),
            "dyad_density": round(float(sig[idx]), 6),
            "smoothed_density": round(float(smoothed[idx]), 6),
            "score": round(float(smoothed[idx] - median), 6),
            "frac_fragments": round(float(sig[idx] / total), 8) if total > 0 else 0.0,
        })
    return calls


def fuzziness(dyad_density: np.ndarray, bin_size: int, position: int,
              fraction: float = 0.5, max_radius: int = 500) -> dict[str, Any]:
    """Breadth of the dyad cluster around one call, in bp.

    Grows outward from ``position`` until ``fraction`` of the local signal within
    ``max_radius`` has been enclosed. A narrow value means the fragments agree on
    one dyad (well positioned); a wide value means a fuzzy or delocalised
    nucleosome. The number is descriptive: no position-quality cutoff is asserted.
    """
    sig = np.asarray(dyad_density, dtype=np.float64)
    if sig.size == 0:
        raise ValueError("dyad_density is empty")
    if not (0.0 < fraction <= 1.0):
        raise ValueError("fraction must be in (0, 1]")
    centre = int(position) // int(bin_size)
    if not (0 <= centre < sig.size):
        raise ValueError("position lies outside the supplied track")
    radius = max(1, int(max_radius) // int(bin_size))
    lo = max(0, centre - radius)
    hi = min(sig.size, centre + radius + 1)
    window = sig[lo:hi]
    target = fraction * float(window.sum())
    if target <= 0:
        return {"position": int(position), "width_bp": 0, "frac_covered": 0.0}
    order = np.argsort(-window)
    cum = np.cumsum(window[order])
    keep = order[: int(np.searchsorted(cum, target)) + 1]
    span = (int(keep.max()) - int(keep.min()) + 1) * int(bin_size)
    return {"position": int(position), "width_bp": int(span),
            "frac_covered": round(float(cum[len(keep) - 1] / window.sum()), 4)}


def spacing_summary(calls: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Neighbour spacing statistics across a sorted call set."""
    if len(calls) < 2:
        return {"n_pairs": 0, "median_spacing": None, "mean_spacing": None}
    positions = sorted(int(c["position"]) for c in calls)
    gaps = np.diff(np.asarray(positions, dtype=np.int64))
    gaps = gaps[gaps > 0]
    if gaps.size == 0:
        return {"n_pairs": 0, "median_spacing": None, "mean_spacing": None}
    return {"n_pairs": int(gaps.size),
            "median_spacing": float(np.median(gaps)),
            "mean_spacing": round(float(gaps.mean()), 2),
            "min_spacing": int(gaps.min()),
            "max_spacing": int(gaps.max())}


# ---------------------------------------------------------------------------
# Fragment extraction from alignments
# ---------------------------------------------------------------------------
def fragments_from_bam(bam_path: str | Path, chrom: str, start: int, end: int,
                       max_fragments: int = 2_000_000,
                       min_mapq: int = 30) -> tuple[np.ndarray, bool]:
    """Collect fragment (dyad) intervals overlapping a region.

    Each record contributes its midpoint and length. Forward reads use
    ``[start, end)``; reads on the minus strand are placed by ``end``, because a
    reverse alignment's ``start`` is the right edge in reference coordinates.
    Paired records are used once (first-in-pair) to avoid double counting.
    Unmapped, secondary, supplementary and duplicate alignments are skipped, as
    are alignments below ``min_mapq``.
    """
    from easyepi.core.io import _require
    pysam = _require("pysam", "pip install 'easyepi[bam]'")

    buf_m: list[int] = []
    buf_l: list[int] = []
    paired_records = 0
    total_records = 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        if chrom not in bam.references:
            available = list(bam.references[:5])
            raise ValueError(
                f"chromosome '{chrom}' absent from {bam_path} (has e.g. {available})")
        for rec in bam.fetch(chrom, int(start), int(end)):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            if rec.is_duplicate:
                continue
            if rec.mapping_quality < min_mapq:
                continue
            total_records += 1
            if rec.is_paired:
                paired_records += 1
                if not rec.is_read1:
                    continue
            f_start = int(rec.reference_start)
            f_end = int(rec.reference_end) if rec.reference_end is not None else f_start
            if f_end <= f_start:
                continue
            buf_m.append((f_start + f_end) // 2)
            buf_l.append(f_end - f_start)
            if len(buf_m) >= max_fragments:
                break
    # A single-end library has no fragment span, so the dyad is unknowable.
    paired = total_records > 0 and paired_records > 0
    if not buf_m:
        return np.empty((0, 2), dtype=np.int64), paired
    mids = np.asarray(buf_m, dtype=np.int64)
    lengths = np.asarray(buf_l, dtype=np.int64)
    return np.column_stack([mids, lengths]), paired


def bam_to_dyad_density(bam_path: str | Path, chrom: str, start: int, end: int,
                        bin_size: int = 10, min_length: int = 120,
                        max_length: int = 200,
                        center: bool = True,
                        max_fragments: int = 2_000_000,
                        min_mapq: int = 30) -> dict[str, Any]:
    """Bin fragment midpoints into a dyad-density track.

    ``min_length``/``max_length`` default to the mononucleosome window (about
    147 bp protected plus linker), so sub-nucleosomal and dinucleosomal fragments
    do not smear the dyad track. ``center=True`` places fragment midpoints at the
    bin centre; with ``center=False`` the midpoint is truncated into its bin, which
    is what a bedGraph-style pileup does.
    """
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    if end <= start:
        raise ValueError("end must exceed start")
    if max_length <= min_length:
        raise ValueError("max_length must exceed min_length")

    frags, paired = fragments_from_bam(bam_path, chrom, start, end,
                                       max_fragments=max_fragments, min_mapq=min_mapq)
    if not paired:
        raise ValueError(
            "nucleosome positioning needs paired-end fragments: a single-end read "
            "reports one cut site but no fragment span, so the dyad cannot be "
            "located. Supply paired-end MNase-seq/ATAC-seq BAM, or use "
            "`easyepi nucleosome` with a measured occupancy signal instead.")
    n_bins = max(1, int(np.ceil((end - start) / bin_size)))
    density = np.zeros(n_bins, dtype=np.float64)
    if frags.size == 0:
        return {"density": density, "n_fragments": 0, "n_used": 0,
                "n_bins": n_bins, "bin_size": int(bin_size),
                "fragment_length_median": None}

    mids, lengths = frags[:, 0], frags[:, 1]
    keep = (lengths >= min_length) & (lengths <= max_length)
    used = mids[keep]
    idx = (used - start) // bin_size
    idx = idx[(idx >= 0) & (idx < n_bins)]
    np.add.at(density, idx, 1.0)
    return {"density": density, "n_fragments": int(frags.shape[0]),
            "n_used": int(used.size), "n_bins": n_bins,
            "bin_size": int(bin_size),
            "fragment_length_median": int(np.median(lengths)) if lengths.size else None}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run_positioning(args) -> int:
    """CLI driver: BAM -> dyad density -> calls -> report + figure."""
    import json

    from easyepi._plots import selected_formats
    from easyepi.footprint import fig_nucleosome_positions
    from easyepi.nucleosome import periodicity

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    track = bam_to_dyad_density(
        args.bam, args.chrom, args.start, args.end,
        bin_size=int(args.bin_size),
        min_length=int(args.min_length), max_length=int(args.max_length),
        max_fragments=int(args.max_fragments), min_mapq=int(args.min_mapq))
    density = track["density"]
    if track["n_used"] == 0:
        raise ValueError(
            f"no fragments between {args.min_length}-{args.max_length} bp in "
            f"{args.chrom}:{args.start}-{args.end}; check the region, the length "
            f"window and that the BAM is coordinate-sorted")

    calls = call_positions(density, int(args.bin_size),
                           smooth_bins=int(args.smooth_bins),
                           mad_k=float(args.mad_k),
                           min_spacing=int(args.min_spacing),
                           min_height=float(args.min_height),
                           min_relative_height=float(args.min_relative_height))
    for call in calls:
        call["genomic_position"] = int(args.start + call["position"])
        call["fuzziness"] = fuzziness(density, int(args.bin_size), call["position"],
                                      fraction=float(args.fuzziness_fraction),
                                      max_radius=int(args.fuzziness_radius))

    spacing = spacing_summary(calls)
    period = periodicity(density, fs_hint=1.0 / max(1, int(args.bin_size)))

    bed_path = out / "nucleosome_positions.bed"
    with bed_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for i, call in enumerate(calls, start=1):
            gpos = call["genomic_position"]
            half = int(args.bin_size)
            writer.writerow([args.chrom, max(0, gpos - half), gpos + half,
                             f"nuc_{i}", call["score"],
                             "+" if i % 2 else "-"])

    report = {
        "bam": str(args.bam),
        "region": f"{args.chrom}:{args.start}-{args.end}",
        "bin_size": int(args.bin_size),
        "length_window": [int(args.min_length), int(args.max_length)],
        "n_fragments_seen": track["n_fragments"],
        "n_fragments_used": track["n_used"],
        "fragment_length_median": track["fragment_length_median"],
        "n_calls": len(calls),
        "spacing": spacing,
        "periodicity": period,
        "parameters": {"smooth_bins": int(args.smooth_bins), "mad_k": float(args.mad_k),
                       "min_spacing": int(args.min_spacing),
                       "min_height": float(args.min_height),
                       "min_relative_height": float(args.min_relative_height)},
        "calls": calls,
        "method": "dyad density (fragment midpoints) + robust local maxima",
    }
    (out / "nucleosome_positions.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    formats = selected_formats(args)
    if formats:
        fig_nucleosome_positions(
            density, calls, bin_size=int(args.bin_size), region=report["region"],
            title=f"Nucleosome dyad density, {report['region']}",
            save=str(out / "nucleosome_positions"), formats=list(formats))

    print(f"[easyepi] nucleosome-positioning: {track['n_used']} fragments "
          f"(median {track['fragment_length_median']} bp) -> {len(calls)} calls "
          f"in {report['region']}")
    if spacing["n_pairs"]:
        print(f"[easyepi] median spacing {spacing['median_spacing']:.0f} bp "
              f"over {spacing['n_pairs']} neighbouring pairs")
    print(f"[easyepi] wrote {bed_path.name}, nucleosome_positions.json")
    return 0
