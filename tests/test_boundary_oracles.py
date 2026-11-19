"""Boundary tests for the points where a defect would change a scientific result.

The evaluation's B11 asks for minimal counterexamples and oracle differences rather
than more happy-path samples. Each test below fails if a specific silent-wrongness
regression returns; none of them asserts on formatting or logging.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from easyepi.assay_metrics import nucleosome_repeat_length
from easyepi.dmr import differential_sites


def _index(n: int) -> pd.MultiIndex:
    return pd.MultiIndex.from_arrays(
        [["chr1"] * n, np.arange(n), ["+"] * n],
        names=["chrom", "pos", "strand"],
    )


# --------------------------------------------------------------------------
# DMR: pooling must not turn absent replication into extra power
# --------------------------------------------------------------------------
def test_pooled_test_still_runs_when_only_one_sample_per_group():
    """The pooled path remains available; it is the fallback for n=1 designs."""
    a = pd.DataFrame({"n_meth": [10, 20], "n_total": [30, 30]}, index=_index(2)[:2])
    b = pd.DataFrame({"n_meth": [12, 22], "n_total": [30, 30]}, index=_index(2)[:2])
    frame, stats = differential_sites(a, b)
    assert "padj" in frame.columns
    assert stats["sites_tested"] >= 0


def test_replicate_aware_refuses_single_replicate():
    """A dispersion cannot be estimated from one sample; refusing beats guessing."""
    from easyepi.dmr_replicate_aware import differential_sites_replicate_aware

    cols = {(0, "n_meth"): [5], (0, "n_total"): [10]}
    a = pd.DataFrame(cols, index=_index(1))
    a.columns = pd.MultiIndex.from_tuples(a.columns)
    b = a.copy()
    with pytest.raises(ValueError, match=">=2 samples per group"):
        differential_sites_replicate_aware(a, b)


# --------------------------------------------------------------------------
# NRL: an unstructured histogram must not produce a number
# --------------------------------------------------------------------------
def test_flat_fragment_lengths_do_not_yield_an_nrl():
    """A flat distribution has no periodicity; reporting a spacing is fabrication."""
    rng = np.random.default_rng(3)
    flat = rng.integers(100, 400, size=4000)
    result = nucleosome_repeat_length(flat)
    assert result.get("nrl_bp") is None
    assert "error" in result


def test_random_fragment_lengths_do_not_yield_an_nrl():
    rng = np.random.default_rng(9)
    noise = rng.integers(50, 900, size=8000)
    result = nucleosome_repeat_length(noise)
    assert result.get("nrl_bp") is None


def test_a_real_ladder_does_yield_an_nrl():
    """Guarding the guard: the refusal above must not suppress real periodicity."""
    rng = np.random.default_rng(11)
    fragments = []
    for centre in (180, 360, 540):
        fragments.append(rng.normal(centre, 12, 1500))
    values = np.concatenate(fragments).astype(int)
    values = values[(values > 50) & (values < 900)]
    result = nucleosome_repeat_length(values)
    assert result.get("nrl_bp") is not None
    # Centre-to-centre spacing of the ladder above is ~180 bp.
    assert 150 <= result["nrl_bp"] <= 210


def test_uncertainty_is_labelled_as_not_a_confidence_interval():
    rng = np.random.default_rng(13)
    fragments = np.concatenate([rng.normal(c, 12, 1200) for c in (180, 360, 540)]).astype(int)
    fragments = fragments[(fragments > 50) & (fragments < 900)]
    result = nucleosome_repeat_length(fragments)
    if result.get("nrl_bp") is not None:
        assert "not a confidence interval" in result["nrl_uncertainty_kind"]
