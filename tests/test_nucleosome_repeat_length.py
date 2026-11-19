"""Nucleosome repeat length from a fragment-length distribution.

The NRL is the spacing of the mono/di/tri-nucleosome ladder, so the test builds a
library whose spacing is known exactly and checks that it is recovered. A test
that only asserted a plausible number would pass on a method that returns a
constant.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.assay_metrics import nucleosome_repeat_length

pytest.importorskip("numpy")


def _ladder(nrl, *, n=6000, mono=147.0, width=12.0, seed=0,
            weights=(0.72, 0.22, 0.06)):
    """Fragments at mono, mono+nrl, mono+2*nrl with Gaussian spread."""
    rng = np.random.default_rng(seed)
    values = []
    for k, w in enumerate(weights):
        centre = mono + k * nrl
        count = int(n * w)
        values.append(rng.normal(loc=centre, scale=width, size=count))
    return np.concatenate(values)


@pytest.mark.parametrize("nrl", [165, 175, 185, 195])
def test_known_repeat_length_is_recovered(nrl):
    res = nucleosome_repeat_length(_ladder(nrl))
    assert "nrl_bp" in res, res
    # Within 5 bp: the parabola refinement on a 1 bp grid is finer than that,
    # but a real library has a shoulder the synthetic one does not.
    assert abs(res["nrl_bp"] - nrl) <= 5.0, res


def test_nrl_is_the_spacing_not_the_mononucleosome_size():
    """A method that returns the mononucleosome peak would fail this.

    Mononucleosome fragments are ~147 bp; the repeat adds the linker. Reporting
    the former as the NRL is the classic 30-50 bp error, so the two numbers must
    come out different here.
    """
    res = nucleosome_repeat_length(_ladder(190, mono=145.0))
    assert res["nrl_bp"] > res["mononucleosome_length_bp"] + 20
    assert abs(res["mononucleosome_length_bp"] - 145) <= 6


def test_a_flat_distribution_refuses_to_invent_a_spacing():
    """No ladder, no NRL: a uniform length histogram has no spacing to find."""
    rng = np.random.default_rng(3)
    flat = rng.uniform(100, 1200, size=4000)
    res = nucleosome_repeat_length(flat)
    # Either it declines, or it reports a weak autocorrelation -- but it must not
    # claim a confident spacing.
    if "nrl_bp" in res:
        assert res["nrl_autocorrelation"] < 0.5 or res["nrl_uncertainty_bp"] > 20
    else:
        assert "error" in res


def test_too_few_fragments_is_refused():
    res = nucleosome_repeat_length([150.0, 180.0, 210.0])
    assert "error" in res and "50" in res["error"]


def test_out_of_range_fragments_are_excluded():
    """A ladder outside the search window must not be fitted from noise."""
    res = nucleosome_repeat_length(_ladder(180), min_length=400, max_length=1200)
    assert "error" in res


def test_uncertainty_is_reported_with_the_estimate():
    res = nucleosome_repeat_length(_ladder(180))
    assert "nrl_uncertainty_bp" in res
    assert res["nrl_uncertainty_bp"] >= 0
    assert res["nrl_search_range"] == [130, 300]


def test_the_note_names_the_common_confusion():
    res = nucleosome_repeat_length(_ladder(180))
    assert "not the mononucleosome fragment length" in res["note"]
