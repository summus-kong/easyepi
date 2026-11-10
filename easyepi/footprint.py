"""Aggregate motif-site signal profiles for descriptive visualization.

Profile-bin tests in legacy helpers do not model biological replication and
must not be interpreted as evidence of TF binding or differential occupancy.
For bias-corrected footprint analysis use an established caller and appropriate
replicate-aware inference. Signal aggregation is related to deepTools profile
visualization (Ramirez et al. 2016, doi:10.1093/nar/gkw257); equivalence has not
been benchmarked.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np


def extract_profiles(signal: np.ndarray, sites: Sequence[int], center: int = 0,
                     flank: int = 100) -> np.ndarray:
    """Extract centered profiles for a list of motif center positions.

    ``signal`` is the per-base accessibility/cut array. ``sites`` are motif center
    positions **within the signal array** (absolute indices). ``center`` is the
    reference coordinate that the returned profiles are centered on (usually 0,
    meaning profiles are indexed from ``-flank`` to ``+flank`` relative to each
    site). Returns (n_sites, 2*flank+1) array or empty if no sites fit.
    """
    sig = np.asarray(signal, dtype=np.float64)
    n = len(sites)
    L = 2 * flank + 1
    out = np.full((n, L), np.nan, dtype=np.float64)
    for i, pos in enumerate(sites):
        # pos is an absolute index into sig; center the window on pos
        lo = pos - flank
        hi = pos + flank + 1
        if lo < 0 or hi > len(sig):
            continue
        out[i] = sig[lo:hi]
    return out[~np.isnan(out).all(axis=1)]


def footprint_profile(signal: np.ndarray, sites: Sequence[int], center: int = 0,
                      flank: int = 100, background_sites: Optional[Sequence[int]] = None,
                      background: Optional[np.ndarray] = None) -> dict:
    """Compute an aggregate binding-site signal profile (+ optional background).

    ``sites`` are absolute indices into ``signal`` (0-based, in the array returned
    by :func:`local_signal`). Returns a dict with:
      * ``target_profile``   — mean per-base signal across all sites (the curve),
      * ``background_profile``— background mean profile (or shoulder baseline),
      * ``score``            — (shoulder_baseline - centre) / shoulder_baseline.
        A *positive* value = centre LOWER than shoulders (a real footprint dip);
        *negative* = centre HIGHER than shoulders (a binding-site signal peak,
        the typical CUT&Tag case). This is a descriptive contrast, NOT itself a
        significance test — see :func:`footprint_significance`.
      * ``n_sites``          — number of sites used.
    """
    target = extract_profiles(signal, sites, center, flank)
    if target.size == 0:
        return {"target_profile": np.array([]), "background_profile": np.array([]),
                "score": np.nan, "n_sites": 0}
    if background is not None:
        bg_profile = np.asarray(background, dtype=np.float64)
    elif background_sites is not None:
        bg = extract_profiles(signal, background_sites, center, flank)
        bg_profile = np.nanmean(bg, axis=0) if bg.size else np.nanmean(target, axis=0)
    else:
        # No explicit background control. Build a flat "shoulder baseline" from the
        # average of the left/right flanking windows of the target profile itself,
        # so the centre-vs-shoulder contrast is well defined without a control set.
        L = len(np.nanmean(target, axis=0))
        shoulder = int(round(0.25 * L))      # ~25% into each flank
        tgt_mean0 = np.nanmean(target, axis=0)
        left = np.nanmean(tgt_mean0[:shoulder]) if shoulder > 0 else np.nan
        right = np.nanmean(tgt_mean0[-shoulder:]) if shoulder > 0 else np.nan
        base = np.nanmean([left, right])
        if not np.isfinite(base):
            base = np.nanmean(signal)
        bg_profile = np.full(L, base, dtype=np.float64)
    tgt_mean = np.nanmean(target, axis=0)
    mid = len(tgt_mean) // 2
    # Footprint score: how much the centre is depleted relative to the local
    # shoulder baseline (a positive score = a dip below the flanks).
    base = float(bg_profile[mid])
    tgt_center = float(tgt_mean[mid])
    if base and base > 0:
        score = (base - tgt_center) / base
    else:
        score = 0.0
    return {"target_profile": tgt_mean, "background_profile": bg_profile,
            "score": float(score), "n_sites": int(target.shape[0]),
            "flank": int(flank)}


def footprint_score_for_sites(signal: np.ndarray, sites: Sequence[int], center: int,
                              flank: int = 100,
                              background_sites: Optional[Sequence[int]] = None) -> float:
    """Convenience: return just the center-dip footprint score for a locus set."""
    return footprint_profile(signal, sites, center, flank,
                             background_sites=background_sites)["score"]


# ---------------------------------------------------------------------------
# Significance: target vs random-background signal contrast
# ---------------------------------------------------------------------------
def _shoulder_depth_per_site(profiles: np.ndarray, core_frac: float = 0.2) -> np.ndarray:
    """Per-locus centre-vs-shoulder contrast statistic.

    ``profiles`` is (n_sites, 2*flank+1). For each site we take the mean signal
    over a small central window (``core``) and the mean over symmetric left/right
    shoulder windows, and return ``(shoulder - core)``. A larger (positive) value
    means the centre is DEPLETED relative to its shoulders (a footprint dip); a
    negative value means the centre is a signal peak (typical CUT&Tag binding
    site). NaN-aware per row.
    """
    prof = np.asarray(profiles, dtype=np.float64)
    if prof.ndim != 2 or prof.shape[1] < 5:
        return np.array([])
    L = prof.shape[1]
    mid = L // 2
    core = max(1, int(round(core_frac * L / 2)))
    # shoulder windows: 20%..40% into each flank
    lo_shoulder = int(round(0.15 * L))
    hi_shoulder = int(round(0.35 * L))
    out = []
    for row in prof:
        c = np.nanmean(row[mid - core:mid + core + 1]) if (mid - core >= 0 and mid + core + 1 <= L) else np.nan
        left = np.nanmean(row[lo_shoulder:hi_shoulder]) if hi_shoulder <= mid else np.nan
        right = np.nanmean(row[L - hi_shoulder:L - lo_shoulder]) if (L - hi_shoulder) >= mid else np.nan
        shoulder = np.nanmean([left, right])
        if np.isfinite(c) and np.isfinite(shoulder):
            out.append(shoulder - c)
        else:
            out.append(np.nan)
    return np.array(out, dtype=np.float64)


def footprint_significance(target_profiles: np.ndarray,
                           bg_profiles: Optional[np.ndarray] = None,
                           core_frac: float = 0.2) -> dict:
    """Test whether target loci differ from (random) background in signal shape.

    The per-locus contrast statistic is ``shoulder - centre`` (from
    :func:`_shoulder_depth_per_site`). We report:
      * the mean contrast for target and (if given) background,
      * a Mann-Whitney U test comparing target vs background (two-sided),
      * a one-sample Wilcoxon signed-rank test that the target contrast is > 0
        (i.e. centre depleted vs its own shoulders — a footprint dip).

    Returns a dict with ``target_mean_dip``, ``background_mean_dip``, ``n_target``,
    ``n_bg``, ``mannwhitney_u``, ``p_mannwhitney`` and ``p_wilcoxon``. Note that
    ``*_dip`` keys keep their historical name; a negative value means the centre
    is a signal PEAK (typical CUT&Tag binding site), not a footprint dip.
    """
    t_dip = _shoulder_depth_per_site(target_profiles, core_frac)
    t_dip = t_dip[np.isfinite(t_dip)]
    res = {"target_mean_dip": float(np.mean(t_dip)) if t_dip.size else np.nan,
           "n_target": int(t_dip.size),
           "target_dips": t_dip.tolist()}
    if bg_profiles is not None and np.asarray(bg_profiles).size:
        b_dip = _shoulder_depth_per_site(bg_profiles, core_frac)
        b_dip = b_dip[np.isfinite(b_dip)]
        res["background_mean_dip"] = float(np.mean(b_dip)) if b_dip.size else np.nan
        res["n_bg"] = int(b_dip.size)
        if t_dip.size > 2 and b_dip.size > 2:
            from scipy.stats import mannwhitneyu
            mw = mannwhitneyu(t_dip, b_dip, alternative="two-sided")
            res["mannwhitney_u"] = float(mw.statistic)
            res["p_mannwhitney"] = float(mw.pvalue)
    if t_dip.size > 2:
        from scipy.stats import wilcoxon
        try:
            w = wilcoxon(t_dip, alternative="greater")  # dip > 0
            res["p_wilcoxon"] = float(w.pvalue)
        except Exception as exc:
            # e.g. all pairwise differences are zero, so the test is undefined.
            res["p_wilcoxon"] = np.nan
            res["p_wilcoxon_error"] = f"{type(exc).__name__}: {exc}"
    return res


def local_signal(access_path: str, chrom: str, start: int, end: int,
                 bins: Optional[int] = None) -> np.ndarray:
    """Read per-base accessibility/cut signal from a bigWig into an array."""
    import pyBigWig
    bw = pyBigWig.open(access_path)
    try:
        vals = bw.values(chrom, int(start), int(end))
        if vals is None:
            return np.zeros(end - start, dtype=np.float64)
        arr = np.array([0.0 if v is None or np.isnan(v) else float(v) for v in vals])
        if bins is not None and bins != len(arr):
            # resample to a fixed bin count
            idx = np.linspace(0, len(arr) - 1, bins)
            arr = np.interp(idx, np.arange(len(arr)), arr)
        return arr
    finally:
        bw.close()


# ---------------------------------------------------------------------------
# Figure: footprint profile
# ---------------------------------------------------------------------------
def fig_footprint(profile: dict, title: str = "TF footprint", save: Optional[str] = None,
                  formats: Optional[List[str]] = None):
    """Draw the binding-signal profile figure; returns the Matplotlib Figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from easyepi.pubfig import finalize_figure
    tgt = np.asarray(profile.get("target_profile", []))
    bg = np.asarray(profile.get("background_profile", []))
    fig, ax = plt.subplots(figsize=(4, 3))
    flank = int(profile.get("flank") or (len(tgt) // 2 if tgt.size else 0))
    # Profiles are indexed from -flank to +flank around the reference position, so
    # the x axis must show signed bp offsets. Plotting them against their array
    # index instead produces a 0..2*flank ramp that looks like a gradient along the
    # region and mislabels the axis.
    x = np.arange(-flank, flank + 1) if tgt.size == 2 * flank + 1 else None
    if bg.size:
        ax.plot(x if x is not None and bg.size == x.size else np.arange(bg.size),
                bg, color="#767676", lw=1, label="background")
    if tgt.size:
        ax.plot(x if x is not None else np.arange(tgt.size), tgt,
                color="#0F4D92", lw=1.6, label="target")
        mid = len(tgt) // 2
        ax.axvspan(-5, 5, color="#B64342", alpha=0.2, lw=0)
        ax.axvline(0, color="#B64342", lw=0.6, ls="--")
    ax.set_title(title)
    ax.set_xlabel(f"distance from reference position (bp, {flank} bp flank)")
    ax.set_ylabel("accessibility / cut")
    ax.legend(frameon=False, fontsize=6)
    if save:
        return finalize_figure(fig, save, formats=["png"] if formats is None else formats, dpi=300)
    return fig


def fig_nucleosome_positions(density, calls, bin_size: int = 10, region: str = "",
                             max_calls: int = 40,
                             title: str = "Nucleosome dyad density",
                             save: Optional[str] = None,
                             formats: Optional[List[str]] = None):
    """Dyad-density track with the called positions marked.

    Only the first ``max_calls`` calls are annotated; labelling every call on a
    100 kb window would overplot the axis and make the figure unreadable.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from easyepi.pubfig import finalize_figure

    dens = np.asarray(density, dtype=np.float64)
    x = np.arange(dens.size) * int(bin_size)
    fig, ax = plt.subplots(figsize=(7, 2.8))
    ax.fill_between(x, dens, color="#0F4D92", alpha=0.28, lw=0)
    ax.plot(x, dens, color="#0F4D92", lw=0.8)
    if calls:
        shown = calls[:max_calls]
        xc = np.array([c["position"] for c in shown], dtype=float)
        yc = np.array([dens[min(dens.size - 1, int(c["position"]) // int(bin_size))]
                       for c in shown], dtype=float)
        ax.scatter(xc, yc, s=14, color="#B64342", zorder=3,
                   label=f"{len(calls)} calls")
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.set_xlabel("position in region (bp)")
    ax.set_ylabel("dyad density (fragments/bin)")
    ax.set_xlim(0, max(1.0, float(x[-1]) if x.size else 1.0))
    ax.set_title(f"{title}  [{region}]" if region else title, fontsize=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    if save:
        return finalize_figure(fig, save, formats=["png"] if formats is None else formats,
                               dpi=300)
    return fig


def fig_tf_annotation(summary: Sequence[dict], n_peaks: int,
                      title: str = "Known-motif annotation",
                      top: int = 12, save: Optional[str] = None,
                      formats: Optional[List[str]] = None):
    """Horizontal bar chart of the fraction of peaks carrying a motif site.

    Only the ``top`` most frequent TFs are drawn; the table holds the rest. Sorted
    the same way as the written summary so the figure and ``tf_summary.tsv`` agree.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from easyepi.pubfig import finalize_figure

    rows = [r for r in summary if r.get("n_peaks_with_site", 0) > 0][:top]
    fig, ax = plt.subplots(figsize=(4.6, max(2.0, 0.26 * len(rows) + 1.0)))
    if not rows:
        ax.text(0.5, 0.5, "no motif site above threshold", ha="center", va="center")
        ax.set_axis_off()
    else:
        names = [r["tf"] for r in rows][::-1]
        pct = [r["pct_peaks_with_site"] for r in rows][::-1]
        ax.barh(names, pct, color="#0F4D92", height=0.7)
        ax.set_xlabel(f"% of {n_peaks} peaks with a site")
        ax.set_xlim(0, max(1.0, min(100.0, max(pct) * 1.25)))
        ax.tick_params(axis="y", labelsize=7)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    ax.set_title(title, fontsize=8)
    fig.tight_layout()
    if save:
        return finalize_figure(fig, save, formats=["png"] if formats is None else formats, dpi=300)
    return fig


# ---------------------------------------------------------------------------
# High-level import entry point (mirrors the CLI subcommand).
# ---------------------------------------------------------------------------
def read_site_bed(path: str) -> List[tuple]:
    """Return [(chrom, center), ...] from a BED file (chrom start end [name])."""
    centers = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            centers.append((parts[0].strip(), (start + end) // 2))
    return centers


def run_footprint(bigwig: str, chrom: str, start: int, end: int,
                  target_bed: str, background_bed: Optional[str] = None,
                  flank: int = 100, title: str = "TF footprint",
                  out_dir: Optional[str] = None, save_fig: bool = True) -> dict:
    """Import-style entry: aggregate TF binding-site signal profile + significance.

    Reads the per-base signal from a bigWig and averages it over motif-centre loci
    (and optional random-background sites), producing the aggregate signal curve
    (same as ``deepTools computeMatrix + plotProfile``) plus a significance test
    of target vs background. See the module docstring for terminology — this is a
    binding-signal profile, NOT a claim of an ATAC "footprint" dip.

    Args:
        bigwig: per-base accessibility/cut bigWig path.
        chrom: chromosome for the probe region.
        start, end: probe region (bp), in which ``target_bed`` centers must fall.
        target_bed: BED of target motif sites (chrom start end).
        background_bed: optional BED of random background/control sites.
        flank: flank (bp) each side of the motif center.
        out_dir: if given, write ``footprint.json`` and the figure here.
        save_fig: whether to render the figure (requires matplotlib).

    Returns:
        dict with ``score`` (centre-vs-shoulder contrast; negative = signal peak),
        ``n_sites``, ``target_profile``, ``background_profile``, ``significance``.
    """
    signal = local_signal(bigwig, chrom, start, end)
    offset = -start
    N = len(signal)
    t_idx = [pos + offset for chrom_, pos in read_site_bed(target_bed) if chrom_ == chrom]
    t_idx = [i for i in t_idx if 0 <= i < N]
    b_idx = None
    if background_bed:
        b_idx = [pos + offset for chrom_, pos in read_site_bed(background_bed) if chrom_ == chrom]
        b_idx = [i for i in b_idx if 0 <= i < N]
    prof = footprint_profile(signal, t_idx, flank=flank, background_sites=b_idx)
    prof["chrom"], prof["start"], prof["end"] = chrom, start, end

    # --- significance: target vs random-background signal contrast ---
    target_mat = extract_profiles(signal, t_idx, 0, flank)
    bg_mat = extract_profiles(signal, b_idx, 0, flank) if b_idx else None
    sig = footprint_significance(target_mat, bg_mat)
    prof["significance"] = sig

    if out_dir:
        import json
        from pathlib import Path
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        # annotate the title with the key significance result
        pval = sig.get("p_mannwhitney")
        p_str = f"(Mann-Whitney p={pval:.2e})" if pval is not None and np.isfinite(pval) else ""
        full_title = f"{title}\ncentre-vs-shoulder={prof['score']:.3f} {p_str}".strip()
        with open(Path(out_dir) / "footprint.json", "w", encoding="utf-8") as fh:
            json.dump({"score": prof["score"], "n_sites": prof["n_sites"],
                       "target_profile": np.asarray(prof["target_profile"]).tolist(),
                       "background_profile": np.asarray(prof["background_profile"]).tolist(),
                       "significance": sig}, fh, indent=2)
        if save_fig:
            fig_footprint(prof, title=full_title, save=str(Path(out_dir) / "footprint"))
    return prof
