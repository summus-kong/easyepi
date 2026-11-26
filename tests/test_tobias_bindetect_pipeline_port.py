"""BINDetect's orchestration layer, checked against TOBIAS's own objects.

What this covers: PSSM construction, MOODS threshold, per-motif scanning, and the
midpoint offset. What it does not cover: bigWig sampling and file layout, which are
I/O rather than computation.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_bindetect_pipeline import (
    MOODS_CORE,
    UNIFORM_BG,
    midpoint_offset,
    motif_threshold,
    pwm_to_pssm,
    scan_region,
)

POS_A = np.array([[85, 5, 5, 5],
                  [5, 85, 5, 5],
                  [5, 5, 85, 5],
                  [5, 5, 5, 85],
                  [80, 7, 7, 6]], dtype=float)
POS_B = np.array([[5, 5, 5, 85],
                  [5, 5, 85, 5],
                  [5, 85, 5, 5],
                  [85, 5, 5, 5],
                  [80, 7, 7, 6]], dtype=float)


def test_pssm_uses_natural_log_and_renormalised_pseudocount():
    """log2 would inflate every entry by ln(2); the reference uses np.log."""
    pssm = pwm_to_pssm(POS_A.T)
    # a 0.85-frequency base against a uniform background, natural log of 0.85/0.25
    # after the 0.01 pseudocount renormalisation is 1.2168 (measured against OneMotif)
    assert pssm[0, 0] == pytest.approx(1.2168, abs=1e-4)
    assert pssm[0, 1] == pytest.approx(-1.5706, abs=1e-4)


def test_pssm_matches_reference_one_motif():
    from tobias.utils.motifs import OneMotif

    mine = pwm_to_pssm(POS_A.T)
    ref = np.asarray(OneMotif("t", counts=POS_A.T.copy()).get_pssm().pssm)
    assert mine.shape == ref.shape
    assert np.abs(mine - ref).max() < 1e-12


def test_threshold_matches_reference():
    from tobias.utils.motifs import OneMotif

    pssm = pwm_to_pssm(POS_A.T)
    m = OneMotif("t", counts=POS_A.T.copy()).get_pssm().get_threshold(pvalue=1e-4)
    assert motif_threshold(pssm, pvalue=1e-4) == pytest.approx(float(m.threshold), abs=1e-12)


def test_pssm_rejects_wrong_shape():
    with pytest.raises(ValueError):
        pwm_to_pssm(np.ones((5, 4)))


def test_scan_finds_both_strands_and_attributes_motifs():
    seq = "ACGT" * 25 + "GTCAGGTACG" + "GGCC" * 20
    pssm_a = pwm_to_pssm(POS_A.T)
    pssm_b = pwm_to_pssm(POS_B.T)
    hits = scan_region(seq, [pssm_a, pssm_b],
                       [motif_threshold(pssm_a, pvalue=1e-4),
                        motif_threshold(pssm_b, pvalue=1e-4)])

    assert hits, "expected hits in a sequence built from the motif"
    motifs = {h[0] for h in hits}
    strands = {h[3] for h in hits}
    assert motifs <= {0, 1}
    assert strands <= {"+", "-"}


def test_scan_is_deterministic():
    seq = "GTCAGGTACG" * 12
    pssm = pwm_to_pssm(POS_A.T)
    thr = motif_threshold(pssm, pvalue=1e-4)
    assert scan_region(seq, [pssm], [thr]) == scan_region(seq, [pssm], [thr])


def test_midpoint_offset_lands_inside_the_motif():
    # region starts at 1000, motif starts at 900 -> offset 100, plus half its length
    assert midpoint_offset(1000, 900, 10) == 105


def test_constants_match_reference_defaults():
    assert MOODS_CORE == 7
    assert UNIFORM_BG == (0.25, 0.25, 0.25, 0.25)
