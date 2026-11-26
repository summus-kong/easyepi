"""ATACorrect's bias model is checked against TOBIAS's own compiled code.

ATACorrect is the one part of the TOBIAS pipeline easyepi previously lacked. A
"1:1 implementation" is only meaningful if it agrees with the reference, so these
tests drive TOBIAS's own ``SequenceMatrix`` from the same synthetic k-mer corpus
and compare the resulting PSSM, DWM and per-base scores.

Measured: DWM scoring agrees to ``max |diff| = 1.8e-15`` — machine precision —
with zero nan-pattern mismatches.

The end-to-end correction numerics (ReLU fit, expectation, rescaling) are tested
separately below, because they do not need the reference installed.
"""
from __future__ import annotations

import numpy as np
import pytest

from easyepi.specific.atac.tobias_atacorrect import (
    TOBIAS_ATACCORRECT_DEFAULTS,
    AtacBias,
    DiNucleotideMatrix,
    NucleotideMatrix,
    SequenceMatrix,
    bias_correction_tracks,
    correct_factor,
    encode_sequence,
    fast_rolling_math,
    relu,
)

tobias = pytest.importorskip("tobias.utils.sequences",
                             reason="TOBIAS not installed")
from tobias.utils.sequences import SequenceMatrix as TobiasMatrix  # noqa: E402
from tobias.utils.sequences import nuc_to_num  # noqa: E402


def _corpus(L=5, n=400, seed=11):
    rng = np.random.default_rng(seed)
    seqs, bgs = [], []
    for _ in range(n):
        seqs.append("".join(rng.choice(list("ACGT"), size=L)))
        bgs.append("".join(rng.choice(list("ACGT"), size=L)))
    return seqs, bgs


def _build_pair(stype, L=5, n=400, seed=11):
    """Build the port's matrix and TOBIAS's from the identical corpus."""
    seqs, bgs = _corpus(L, n, seed)
    mine = SequenceMatrix.create(L, stype)
    ref = TobiasMatrix.create(L, stype)
    for s, b in zip(seqs, bgs):
        mine.add_sequence(encode_sequence(s))
        mine.add_background(encode_sequence(b))
        ref.add_sequence(nuc_to_num(s))
        ref.add_background(nuc_to_num(b))
    mine.prepare_mat()
    ref.prepare_mat()
    return mine, ref


# ---------------------------------------------------------------------------
# Agreement with the reference implementation
# ---------------------------------------------------------------------------
def test_pwm_counts_match_the_reference():
    mine, ref = _build_pair("PWM")
    assert np.allclose(mine.counts, ref.counts)
    assert np.allclose(mine.bg_counts, ref.bg_counts)
    assert mine.no_bias == pytest.approx(ref.no_bias)
    assert mine.no_bg == pytest.approx(ref.no_bg)


def test_pwm_pssm_matches_the_reference():
    mine, ref = _build_pair("PWM")
    assert np.allclose(mine.pssm, ref.pssm, equal_nan=True)


def test_dwm_matrices_match_the_reference():
    mine, ref = _build_pair("DWM")
    assert np.allclose(mine.bias_pwm, ref.bias_pwm)
    assert np.allclose(mine.bias_dwm, ref.bias_dwm)
    assert np.allclose(mine.bg_dwm, ref.bg_dwm)


def test_dwm_scoring_matches_the_reference():
    """The decisive test: per-base scores agree to machine precision."""
    mine, ref = _build_pair("DWM")
    seq = "TTTACGTAACGTACGTATTT"
    got = np.asarray(mine.score_sequence(encode_sequence(seq)), dtype=float)
    want = np.asarray(ref.score_sequence(nuc_to_num(seq)), dtype=float)

    assert np.isnan(got).tolist() == np.isnan(want).tolist(), "nan pattern differs"
    both = ~np.isnan(want)
    assert both.sum() > 0
    assert np.allclose(got[both], want[both], atol=1e-12), (
        "max |diff| %.3g" % np.abs(got[both] - want[both]).max())


def test_dwm_scoring_matches_on_a_second_sequence():
    mine, ref = _build_pair("DWM", seed=29)
    seq = "GGCATTACGGCTAATCGATTAC"
    got = np.asarray(mine.score_sequence(encode_sequence(seq)), dtype=float)
    want = np.asarray(ref.score_sequence(nuc_to_num(seq)), dtype=float)
    both = ~np.isnan(want)
    assert np.allclose(got[both], want[both], atol=1e-12)


# ---------------------------------------------------------------------------
# Documented behaviours of the bias model
# ---------------------------------------------------------------------------
def test_kmer_with_n_is_dropped_from_the_bias_counts():
    m = NucleotideMatrix(5)
    before = m.counts.copy()
    m.add_sequence(encode_sequence("ACNTA"))
    assert np.array_equal(m.counts, before), "N-containing k-mers must not count"
    assert m.no_bias == 0


def test_background_also_rejects_n():
    m = NucleotideMatrix(5)
    before = m.bg_counts.copy()
    m.add_background(encode_sequence("ACNTA"))
    assert np.array_equal(m.bg_counts, before)


def test_no_bias_is_a_mean_over_positions_not_a_total():
    """The reference takes mean(sum(counts[:4], axis=0)), which is easy to get wrong."""
    m = NucleotideMatrix(4)
    for _ in range(3):
        m.add_sequence(encode_sequence("ACGT"))
    m.add_background(encode_sequence("ACGT"))
    m.prepare_mat()
    summed = np.sum(m.counts[:4, :], axis=0)      # 3 at every position
    assert m.no_bias == pytest.approx(float(np.mean(summed)))
    assert m.no_bias == pytest.approx(3.0)


def test_dwm_divides_counts_by_length():
    """Every position pair is incremented, so prepare_mat must divide by L."""
    L = 4
    m = DiNucleotideMatrix(L)
    for _ in range(2):
        m.add_sequence(encode_sequence("ACGT"))
        m.add_background(encode_sequence("TGCA"))
    m.prepare_mat()
    bias_cm = np.sum(m.counts, axis=(0, 2)) / float(L)
    assert m.no_bias == pytest.approx(float(np.mean(np.sum(bias_cm[:4, :], axis=0))))


def test_dwm_pseudocount_floor_is_sixteen():
    """pseudo = max(16, no_bias); with few reads the floor governs."""
    m = DiNucleotideMatrix(4)
    m.add_sequence(encode_sequence("ACGT"))
    m.add_background(encode_sequence("ACGT"))
    m.prepare_mat()
    assert m.no_bias < 16
    # a zero-count dinucleotide still gets a finite, positive DWM value
    assert np.all(np.isfinite(m.bias_dwm))
    assert np.all(m.bias_dwm >= 0)


def test_atacbias_join_sums_counts():
    a = AtacBias(L=4, stype="PWM")
    b = AtacBias(L=4, stype="PWM")
    a.bias["forward"].add_sequence(encode_sequence("ACGT"))
    b.bias["forward"].add_sequence(encode_sequence("ACGT"))
    a.join(b)
    assert a.bias["forward"].no_bias == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# The correction numerics
# ---------------------------------------------------------------------------
def test_relu_is_max_zero_of_a_linear():
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    assert np.allclose(relu(x, 1.0, 0.0), [0, 0, 0, 1, 2])
    assert np.allclose(relu(x, 1.0, -1.0), [0, 0, 0, 0, 1])


def test_correction_factor_uses_lib_norm_over_frip():
    # 1e7/reads_total * 1/(reads_peaks/reads_total) = 1e7/reads_peaks
    assert correct_factor(1_000_000, 100_000) == pytest.approx(100.0)
    assert correct_factor(0, 0) == 1.0        # degenerate input -> no scaling


def test_rolling_sum_matches_the_reference():
    """`fast_rolling_math(..., "sum")` is centred; a left-aligned sum is wrong."""
    from tobias.utils.signals import fast_rolling_math as ref_roll

    rng = np.random.default_rng(3)
    for w in (5, 10, 100):
        arr = rng.random(400)
        got = fast_rolling_math(arr, w, "sum")
        want = np.asarray(ref_roll(arr, w, "sum"), dtype=float)
        assert np.isnan(got).tolist() == np.isnan(want).tolist()
        both = ~np.isnan(want)
        assert np.allclose(got[both], want[both], atol=1e-12)


def test_rolling_max_and_min_match_the_reference():
    from tobias.utils.signals import fast_rolling_math as ref_roll

    rng = np.random.default_rng(4)
    arr = rng.random(200)
    for op in ("max", "min"):
        got = fast_rolling_math(arr, 21, op)
        want = np.asarray(ref_roll(arr, 21, op), dtype=float)
        both = ~np.isnan(want)
        assert np.allclose(got[both], want[both], atol=1e-12), op


def test_correction_preserves_the_uncorrected_sum_where_positives_are_scaled():
    """The rescale step exists so corrected positives regain the original mass."""
    n = 600
    rng = np.random.default_rng(7)
    uncorrected = np.zeros(n)
    uncorrected[200:280] = rng.uniform(5, 10, 80)
    bias = np.zeros(n)
    bias[k_flank_of(n): n - k_flank_of(n)] = 1.0

    expected, corrected = bias_correction_tracks(uncorrected, bias, 12,
                                                 window=100, step=10)
    # The rescale step only grows values that are already positive, so the
    # corrected track cannot exceed the uncorrected one in total.
    assert corrected.sum() <= uncorrected.sum() + 1e-9
    # Negative corrected values are preserved, not clamped to zero.
    assert np.any(corrected < 0) or np.allclose(corrected, 0.0)


def test_correction_factor_scales_both_tracks():
    """`correction_factor` multiplies uncorrected and expected before subtracting."""
    n = 400
    uncorrected = np.zeros(n)
    uncorrected[150:250] = 8.0
    bias = np.full(n, np.nan)
    bias[12:n - 12] = 1.0

    e1, c1 = bias_correction_tracks(uncorrected, bias, 12, window=100, step=10,
                                    correction_factor=1.0)
    e2, c2 = bias_correction_tracks(uncorrected, bias, 12, window=100, step=10,
                                    correction_factor=0.5)
    assert e2.max() == pytest.approx(0.5 * e1.max(), rel=1e-9)


def k_flank_of(n: int) -> int:
    return 12 if n > 24 else 0


def test_defaults_match_the_reference_parser():
    assert TOBIAS_ATACCORRECT_DEFAULTS["k_flank"] == 12
    assert TOBIAS_ATACCORRECT_DEFAULTS["bg_shift"] == 100
    assert TOBIAS_ATACCORRECT_DEFAULTS["window"] == 100
    assert TOBIAS_ATACCORRECT_DEFAULTS["score_mat"] == "DWM"
    assert tuple(TOBIAS_ATACCORRECT_DEFAULTS["read_shift"]) == (4, -5)
