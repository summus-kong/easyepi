"""Report how well a background actually matches its target set.

Why this exists as its own module
---------------------------------
`gc_background` picks matched controls, but matching is a claim that has to be
checkable. "GC-matched background" in an output file is unfalsifiable: the reader
cannot tell whether the tolerance was met for 99% of targets or 40%, and a
background that silently loses the hard cases biases the comparison in the
direction of whatever survived. So every run reports the two distributions it
used, side by side, with the test that says whether they differ.

What is reported
----------------
* GC fraction and length: n, mean, sd, min, quartiles, max, for both sets.
* A Kolmogorov-Smirnov statistic and p-value for GC and for length, because
  "matched within tolerance 0.05" and "the distributions agree" are different
  statements and only the second one is what an enrichment test assumes.
* The number of targets that received no control, with reasons.

What is deliberately not reported
---------------------------------
No claim about mappability, repeat content, chromatin state or accessibility.
Those are not measured here, so a background matched on GC and length is
composition-matched only, and the caller is told that in the payload.
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np


def _summarize(values: Sequence[float]) -> dict[str, Any]:
    """n/mean/sd/min/quartiles/max, with NaN rather than a fabricated zero."""
    arr = np.asarray([v for v in values if v is not None and np.isfinite(v)],
                     dtype=np.float64)
    if arr.size == 0:
        return {"n": 0, "mean": None, "sd": None, "min": None,
                "q25": None, "median": None, "q75": None, "max": None}
    return {
        "n": int(arr.size),
        "mean": round(float(arr.mean()), 6),
        "sd": round(float(arr.std(ddof=1)), 6) if arr.size > 1 else None,
        "min": round(float(arr.min()), 6),
        "q25": round(float(np.percentile(arr, 25)), 6),
        "median": round(float(np.percentile(arr, 50)), 6),
        "q75": round(float(np.percentile(arr, 75)), 6),
        "max": round(float(arr.max()), 6),
    }


def _ks(a: Sequence[float], b: Sequence[float]) -> dict[str, Any]:
    """Two-sample KS test, or a stated reason it could not run.

    Returning `None` without a reason is how a report ends up silent about the
    one check that would have caught an unmatched background.
    """
    xa = np.asarray([v for v in a if v is not None and np.isfinite(v)], dtype=np.float64)
    xb = np.asarray([v for v in b if v is not None and np.isfinite(v)], dtype=np.float64)
    if xa.size == 0 or xb.size == 0:
        return {"statistic": None, "pvalue": None,
                "note": "empty sample; KS not run"}
    try:
        from scipy.stats import ks_2samp
    except Exception:  # scipy is a core dependency, but never crash a report for it
        return {"statistic": None, "pvalue": None, "note": "scipy unavailable"}
    res = ks_2samp(xa, xb)
    return {"statistic": round(float(res.statistic), 6),
            "pvalue": float(res.pvalue)}


def compare_distributions(target_gc: Sequence[float],
                          bg_gc: Sequence[float],
                          target_len: Sequence[float],
                          bg_len: Sequence[float]) -> dict[str, Any]:
    """Side-by-side GC and length distributions with a KS test on each."""
    return {
        "gc": {
            "target": _summarize(target_gc),
            "background": _summarize(bg_gc),
            "ks": _ks(target_gc, bg_gc),
            "units": "fraction of A+C+G+T",
        },
        "length": {
            "target": _summarize(target_len),
            "background": _summarize(bg_len),
            "ks": _ks(target_len, bg_len),
            "units": "bp",
        },
    }


def background_report(targets: Sequence[tuple[str, int, int]],
                      controls: Sequence[tuple[str, int, int]],
                      target_seqs: Sequence[str] | None = None,
                      control_seqs: Sequence[str] | None = None,
                      match_summary: dict[str, Any] | None = None,
                      matched_on: Sequence[str] = ("gc", "length")) -> dict[str, Any]:
    """Assemble the report a caller must publish alongside an enrichment result.

    ``target_seqs``/``control_seqs`` are optional: when supplied, GC is measured
    on the actual sequences (which is what a motif scan sees); otherwise it falls
    back to the interval lengths only and says so, rather than reporting a GC
    comparison it never made.
    """
    from easyepi.gc_background import gc_fraction

    out: dict[str, Any] = {
        "matched_on": list(matched_on),
        "n_target": len(targets),
        "n_control": len(controls),
        "composition_matched_only": True,
        "not_matched": ["mappability", "repeat content", "chromatin state",
                        "accessibility", "expression"],
    }

    if target_seqs is not None and control_seqs is not None:
        t_gc = [gc_fraction(s) for s in target_seqs]
        c_gc = [gc_fraction(s) for s in control_seqs]
        t_len = [len(s) for s in target_seqs]
        c_len = [len(s) for s in control_seqs]
        out["distributions"] = compare_distributions(t_gc, c_gc, t_len, c_len)
        out["gc_source"] = "measured on sequences"
    else:
        t_len = [int(e) - int(s) for _c, s, e in targets]
        c_len = [int(e) - int(s) for _c, s, e in controls]
        out["distributions"] = {
            "gc": {"status": "not measured",
                   "note": "no sequences supplied, so GC was not compared"},
            "length": {"target": _summarize(t_len), "background": _summarize(c_len),
                       "ks": _ks(t_len, c_len), "units": "bp"},
        }
        out["gc_source"] = "not measured"

    if match_summary is not None:
        out["matching"] = {
            "n_matched": match_summary.get("n_matched"),
            "n_unmatched": match_summary.get("n_unmatched"),
            "frac_matched": match_summary.get("frac_matched"),
            "tolerance": match_summary.get("tolerance"),
            "min_distance": match_summary.get("min_distance"),
        }
        # A low match rate is the failure mode that biases a result, so it is
        # surfaced rather than left for the reader to notice in the numbers.
        frac = match_summary.get("frac_matched")
        if isinstance(frac, (int, float)) and frac < 0.8:
            out["warning"] = (
                f"only {frac:.0%} of targets received a matched control; the "
                "unmatched targets are not a random subset, so the comparison "
                "covers a biased subset of the peaks")

    return out
