"""Background sensitivity: does the ranking depend on the null?

The comparison function is the part that must be right before any real-data number
means anything, so it is tested against rankings whose answer is known by
construction.
"""
from __future__ import annotations

import pytest

from easyepi.motifs.background_sensitivity import compare_rankings


def _rows(pairs):
    """Build enrichment rows from (name, log_p) pairs, already ordered."""
    return [{"name": n, "log_p": lp} for n, lp in pairs]


def test_identical_rankings_are_perfectly_correlated():
    base = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0), ("D", 6.0)])
    same = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0), ("D", 6.0)])
    out = compare_rankings({"genomic": base, "matched": same}, top_n=4)
    cmp = out["comparisons"]["matched"]
    assert cmp["spearman_rho"] == pytest.approx(1.0)
    assert cmp["top_overlap"]["n"] == 4
    assert cmp["largest_rank_moves"][0]["delta"] == 0


def test_reversed_ranking_is_perfectly_anti_correlated():
    base = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0), ("D", 6.0)])
    flipped = _rows([("D", 9.0), ("C", 8.0), ("B", 7.0), ("A", 6.0)])
    out = compare_rankings({"genomic": base, "matched": flipped}, top_n=2)
    cmp = out["comparisons"]["matched"]
    assert cmp["spearman_rho"] == pytest.approx(-1.0)
    assert cmp["top_overlap"]["n"] == 0


def test_rank_moves_are_reported_with_direction():
    base = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0), ("D", 6.0)])
    moved = _rows([("C", 9.0), ("A", 8.0), ("B", 7.0), ("D", 6.0)])
    out = compare_rankings({"genomic": base, "matched": moved}, top_n=4)
    moves = {m["name"]: m for m in out["comparisons"]["matched"]["largest_rank_moves"]}
    assert moves["C"]["delta"] == -2      # up two places
    assert moves["A"]["delta"] == 1       # down one


def test_a_constant_score_vector_reports_undefined_not_zero():
    """"Not measurable" and "uncorrelated" are different findings."""
    base = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0)])
    flat = _rows([("A", 5.0), ("B", 5.0), ("C", 5.0)])
    out = compare_rankings({"genomic": base, "matched": flat}, top_n=3)
    cmp = out["comparisons"]["matched"]
    assert cmp["spearman_rho"] is None
    assert "constant" in cmp["spearman_undefined_because"]


def test_motifs_missing_from_one_ranking_are_excluded_not_ranked_last():
    """A motif absent from a ranking has no rank; inventing one is worse than omitting."""
    base = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0), ("D", 6.0)])
    partial = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0)])
    out = compare_rankings({"genomic": base, "matched": partial}, top_n=4)
    cmp = out["comparisons"]["matched"]
    assert cmp["n_shared_motifs"] == 3
    assert out["n_common_motifs"] == 3


def test_too_few_shared_motifs_reports_undefined_correlation():
    base = _rows([("A", 9.0), ("B", 8.0)])
    other = _rows([("A", 9.0), ("B", 8.0)])
    out = compare_rankings({"genomic": base, "matched": other}, top_n=2)
    assert out["comparisons"]["matched"]["spearman_rho"] is None


def test_one_configuration_is_rejected():
    with pytest.raises(ValueError, match="at least two"):
        compare_rankings({"genomic": _rows([("A", 1.0)])})


def test_the_base_configuration_is_the_first_supplied():
    """The function must not guess which null is "correct"."""
    a = _rows([("A", 9.0), ("B", 8.0), ("C", 7.0)])
    b = _rows([("C", 9.0), ("B", 8.0), ("A", 7.0)])
    out = compare_rankings({"first": a, "second": b})
    assert out["base_configuration"] == "first"
    assert list(out["comparisons"]) == ["second"]
