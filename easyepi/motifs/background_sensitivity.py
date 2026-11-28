"""Measure how much the choice of background changes a motif enrichment ranking.

This closes the evidence gap the roadmap calls M3. The background construction
exists (composition-matched, and accessibility-matched since `--accessible`), but
nothing yet shows the construction *matters*: a ranking that is identical under
three different nulls could mean the biology is strong, or that the test is
insensitive and the ranking is meaningless.

What it does, and why each part is needed:

  configs   One target sequence set, several background sets. Each is a different
            claim about "what would this look like with no enrichment".
  enrich    The same PWM scan and the same test for every configuration, so the
            only thing that varies is the null.
  compare   Spearman rho between configurations, the overlap of each top-N, and
            the largest rank move per motif. A single "stable / unstable" verdict
            would hide which motifs moved and by how much.

The comparison deliberately does not decide what a large move means. A motif whose
rank collapses under a matched background was being carried by composition; one
that survives is not thereby proven to be bound. That distinction is the report's
job to expose, not to adjudicate.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence

import numpy as np


def _rank_of(rows: Sequence[Dict], name: str) -> int | None:
    """1-based position of ``name`` in ``rows``, or None when absent."""
    for i, row in enumerate(rows, 1):
        if row.get("name") == name:
            return i
    return None


def _spearman(a: Sequence[float], b: Sequence[float]) -> float | None:
    """Spearman rho for two aligned score vectors, or None when undefined."""
    x = np.asarray(a, dtype=np.float64)
    y = np.asarray(b, dtype=np.float64)
    if x.size < 3 or y.size != x.size:
        return None
    if np.all(x == x[0]) or np.all(y == y[0]):
        # A constant vector has no rank correlation; returning 0.0 would read as
        # "uncorrelated" when the truth is "not measurable".
        return None
    rx = np.argsort(np.argsort(x)).astype(np.float64)
    ry = np.argsort(np.argsort(y)).astype(np.float64)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.sqrt((rx ** 2).sum() * (ry ** 2).sum()))
    return None if denom == 0 else float((rx * ry).sum() / denom)


def compare_rankings(rankings: Dict[str, List[Dict]],
                     top_n: int = 10,
                     score_key: str = "log_p") -> Dict[str, Any]:
    """Compare enrichment rankings produced under different backgrounds.

    ``rankings`` maps a configuration label to the rows returned by
    :func:`easyepi.motifs.homer_workflow.known_motif_enrichment_homer`.

    Every configuration is compared against the **first** one, which the caller
    should make explicit in its own labels; this function does not guess which
    background is "correct".
    """
    labels = list(rankings)
    if len(labels) < 2:
        raise ValueError("need at least two configurations to compare")
    base = labels[0]
    base_rows = rankings[base]

    # Motifs common to every configuration: a motif missing from one ranking has no
    # rank to compare, and treating it as last would invent an ordering.
    common: set[str] | None = None
    for label in labels:
        names = {r.get("name") for r in rankings[label]}
        common = names if common is None else (common & names)
    common = sorted(common or set())

    base_by_name = {r.get("name"): r for r in base_rows}
    pairs: Dict[str, Any] = {}
    for label in labels[1:]:
        other_by_name = {r.get("name"): r for r in rankings[label]}
        shared = [n for n in common if n in base_by_name and n in other_by_name]
        rho = _spearman(
            [float(base_by_name[n].get(score_key, 0.0)) for n in shared],
            [float(other_by_name[n].get(score_key, 0.0)) for n in shared],
        )
        base_top = [r.get("name") for r in base_rows[:top_n]]
        other_top = [r.get("name") for r in rankings[label][:top_n]]
        overlap = sorted(set(base_top) & set(other_top))

        moves = []
        for name in shared:
            r_base = _rank_of(base_rows, name)
            r_other = _rank_of(rankings[label], name)
            if r_base is not None and r_other is not None:
                moves.append({"name": name, "rank_in_base": r_base,
                              "rank_here": r_other,
                              "delta": r_other - r_base})
        moves.sort(key=lambda m: -abs(m["delta"]))

        pairs[label] = {
            "spearman_rho": None if rho is None else round(rho, 4),
            "spearman_undefined_because": (
                None if rho is not None else
                "fewer than three shared motifs, or a constant score vector"),
            "n_shared_motifs": len(shared),
            "top_overlap": {"n": len(overlap), f"top_{top_n}": overlap},
            "largest_rank_moves": moves[:10],
        }

    return {
        "base_configuration": base,
        "configurations": labels,
        "top_n": top_n,
        "score_key": score_key,
        "n_common_motifs": len(common),
        "comparisons": pairs,
        "how_to_read": (
            "Each comparison is against the base configuration, which is the first "
            "one supplied. A low rho or a large rank move means the ranking depends "
            "on the null, not that either ranking is wrong; a high rho means the "
            "ranking is stable, which is necessary but not sufficient for it to be "
            "biologically meaningful."),
    }
