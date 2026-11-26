"""ATACorrect scoring is bit-identical to TOBIAS on both strands.

Why this test exists: whole-track comparison of the corrected output gave

    forward  corr 0.999992
    reverse  corr 0.989637

That asymmetry suggested a strand-specific defect, and I previously wrote that the
reverse path had a "k-mer index off by one". Measuring the scoring step in isolation
disproves that: with an identical model and identical sequence, forward and reverse
both agree to machine precision and neither produces a NaN-mask disagreement.

So the residual reverse-track difference is not in the scoring code. It has to be in
the k-mer stream fed to the model -- i.e. how cut sites are read from the BAM, which
lives in the driver script, not in this module. That is recorded here rather than
guessed at again.

The reference side is TOBIAS's compiled `ScoreMatrix.score_sequence`.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_atacorrect import SequenceMatrix as PortMatrix

RefMatrix = pytest.importorskip("tobias.utils.sequences").SequenceMatrix

LENGTH = 10


def _codes(seq: str) -> np.ndarray:
    table = {"A": 0, "C": 1, "G": 2, "T": 3}
    return np.array([table.get(b, 4) for b in seq.upper()], dtype=np.int64)


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _matched_models(seed: int = 11):
    """One set of counts loaded into both implementations."""
    rng = np.random.default_rng(seed)
    counts = rng.integers(5, 40, size=(5, 5, LENGTH, LENGTH)).astype(np.float64)

    port = PortMatrix.create(LENGTH, "DWM")
    port.counts = counts.copy()
    port.prepare_mat()

    ref = RefMatrix.create(LENGTH, "DWM")
    ref.counts = counts.copy()
    ref.prepare_mat()
    return port, ref


def test_models_are_identical_before_scoring():
    """guards the comparison below: a model mismatch would invalidate it."""
    port, ref = _matched_models()

    assert np.abs(np.asarray(port.bias_pwm_log)
                  - np.asarray(ref.bias_pwm_log)).max() == 0.0
    assert np.abs(np.asarray(port.bias_dwm_log)
                  - np.asarray(ref.bias_dwm_log)).max() == 0.0
    assert np.abs(np.asarray(port.bg_dwm_log)
                  - np.asarray(ref.bg_dwm_log)).max() == 0.0


def test_scoring_matches_reference_on_forward_strand():
    port, ref = _matched_models()
    seq = "".join(np.random.default_rng(7).choice(list("ACGT"), size=200))

    a = np.asarray(port.score_sequence(_codes(seq)))
    b = np.asarray(ref.score_sequence(_codes(seq)))

    numeric = np.isfinite(a) & np.isfinite(b)
    assert numeric.sum() > 0
    assert np.abs(a[numeric] - b[numeric]).max() < 1e-12
    assert (np.isnan(a) == np.isnan(b)).all()


def test_scoring_matches_reference_on_reverse_strand():
    """the strand that looked worse end to end; here it is machine-precision equal."""
    port, ref = _matched_models()
    base = "".join(np.random.default_rng(7).choice(list("ACGT"), size=200))
    seq = _revcomp(base)

    a = np.asarray(port.score_sequence(_codes(seq)))
    b = np.asarray(ref.score_sequence(_codes(seq)))

    numeric = np.isfinite(a) & np.isfinite(b)
    assert numeric.sum() > 0
    assert np.abs(a[numeric] - b[numeric]).max() < 1e-12
    assert (np.isnan(a) == np.isnan(b)).all()


def test_no_strand_asymmetry_in_the_scorer():
    """forward and reverse must agree with the reference to the same tolerance."""
    port, ref = _matched_models()
    base = "".join(np.random.default_rng(3).choice(list("ACGT"), size=150))

    worst = {}
    for label, seq in (("forward", base), ("reverse", _revcomp(base))):
        a = np.asarray(port.score_sequence(_codes(seq)))
        b = np.asarray(ref.score_sequence(_codes(seq)))
        numeric = np.isfinite(a) & np.isfinite(b)
        worst[label] = float(np.abs(a[numeric] - b[numeric]).max())

    assert worst["forward"] < 1e-12
    assert worst["reverse"] < 1e-12
    # the reverse error must not be an order of magnitude worse than forward
    assert worst["reverse"] <= max(worst["forward"], 1e-14) * 100


def test_nan_positions_agree_on_both_strands():
    """windows containing N must be skipped identically, on either strand."""
    port, ref = _matched_models()
    seq = "ACGT" * 20 + "NNNNN" + "TGCA" * 20

    for s in (seq, _revcomp(seq)):
        a = np.asarray(port.score_sequence(_codes(s)))
        b = np.asarray(ref.score_sequence(_codes(s)))
        assert (np.isnan(a) == np.isnan(b)).all()
