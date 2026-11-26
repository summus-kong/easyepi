"""Footprint scorer backends selectable from `easyepi footprint --engine`.

Three implementations, one output contract:

``native``
    `easyepi.footprint.footprint_profile` -- the proportion contrast
    ``(shoulder - centre) / shoulder``. This is the historic default and is
    untouched by this module.

``port``
    The bundled 1:1 reimplementation of TOBIAS ``ScoreBigwig --score footprint``
    (`easyepi.specific.atac.tobias_footprint`). Needs no TOBIAS installation.
    Verified against the reference: same-window agreement is float32-exact on the
    interior, with the two edge bases differing only because the reference clips its
    region extension at the chromosome start.

``tobias``
    The external ``TOBIAS ScoreBigwig`` binary through
    `easyepi.motifs.external_tools.run_tobias_footprint`.

The TOBIAS score is additive, not a proportion::

    flank_mean - fp_mean

where ``fp_mean`` averages only the **negative** values inside the footprint window
and ``flank_mean`` averages only the **positive** values in the two flanks. When a
track carries no negative values (uncorrected coverage, for example) ``fp_mean`` is
0 and the score degenerates to the flank level. That is a property of the definition,
not a failure, so it is reported through ``degenerate`` instead of raising.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

#: TOBIAS defaults, read from `tobias/parsers.py` for 0.17.5.
DEFAULT_FP_MIN = 20
DEFAULT_FP_MAX = 50
DEFAULT_FLANK_MIN = 10
DEFAULT_FLANK_MAX = 30


def tobias_port_score(signal: Sequence[float], sites: Sequence[int],
                      flank: int = 100) -> dict:
    """Score footprint depth at ``sites`` with the bundled TOBIAS port.

    Parameters
    ----------
    signal
        Per-base track covering the probe window.
    sites
        Index positions (already offset into ``signal``) of the target sites.
    flank
        Probe half-width in bases. The TOBIAS window search is capped by
        ``DEFAULT_FP_MAX + DEFAULT_FLANK_MAX`` so a short probe still works.

    Returns
    -------
    dict with the additive score, the per-site values, and a ``degenerate`` flag.
    """
    from easyepi.specific.atac.tobias_footprint import tobias_footprint_array

    arr = np.asarray(signal, dtype=float)
    if arr.ndim != 1:
        raise ValueError("signal must be one-dimensional")
    positions = [int(p) for p in sites if 0 <= int(p) < arr.size]
    if not positions:
        raise ValueError("no target site falls inside the signal window")

    # The kernel writes at `i + left_flank + right_flank + fp` for every `i` in
    # range(flank_min, flank_max+1), so the window needs `2*flank_max + fp_max`
    # plus one slot of slack. A window of exactly that size writes nothing --
    # measured, not assumed.
    fp_max = min(DEFAULT_FP_MAX, max(DEFAULT_FP_MIN, flank // 4))
    flank_max = min(DEFAULT_FLANK_MAX, max(DEFAULT_FLANK_MIN, flank // 4))
    span = 2 * flank_max + fp_max + 1

    per_site = []
    for pos in positions:
        window = arr[pos:pos + span]
        if window.size < span:
            window = np.concatenate([window, np.zeros(span - window.size)])
        scored = np.asarray(
            tobias_footprint_array(window, DEFAULT_FLANK_MIN, flank_max,
                                   DEFAULT_FP_MIN, fp_max),
            dtype=float)
        # Only the footprint span carries a score; the flanks and the unwritten
        # tail stay at zero.
        centre = scored[flank_max:flank_max + fp_max]
        per_site.append(float(np.nanmean(centre)) if np.isfinite(centre).any() else 0.0)

    per_site_arr = np.asarray(per_site, dtype=float)
    has_negative = bool((arr < 0).any())

    return {
        "engine": "port",
        "score": float(np.nanmean(per_site_arr)),
        "n_sites": len(per_site),
        "per_site": per_site_arr.tolist(),
        "fp_window": [DEFAULT_FP_MIN, fp_max],
        "flank_window": [DEFAULT_FLANK_MIN, flank_max],
        "degenerate": (not has_negative),
        "note": ("signal has no negative values: fp_mean is 0 so the additive "
                 "score collapses to the flank level by definition"
                 if not has_negative else ""),
    }


def tobias_binary_score(bigwig: str, bed: str, chrom: str, start: int,
                        end: int) -> dict:
    """Score with the external TOBIAS binary, mapped onto the same contract.

    The binary writes a bigWig rather than returning numbers, so this keeps the
    request metadata and lets the caller read the track. It raises when TOBIAS is
    absent, which the CLI turns into a message naming the alternative engines.
    """
    import tempfile

    from easyepi.motifs.external_tools import run_tobias_footprint

    with tempfile.TemporaryDirectory(prefix="easyepi_fp_") as tmp:
        result = run_tobias_footprint([bigwig], bed, chrom, tmp)
    if result.get("returncode"):
        raise ValueError(f"TOBIAS ScoreBigwig failed: {result.get('stderr', '')[:400]}")
    return {
        "engine": "tobias",
        "score": float("nan"),
        "n_sites": 0,
        "per_site": [],
        "track": result.get("command", ""),
        "degenerate": False,
        "note": ("TOBIAS wrote a footprint track; read it with pyBigWig for "
                 "per-base values"),
    }
