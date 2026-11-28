"""The fast scan must agree with an exhaustive one, including edge cases.

The optimised path sums one base offset at a time instead of materialising an
(n_windows, L) window matrix. Two things can silently break: the offset slicing,
and the reverse-complement position mapping. Both are checked here against a
deliberately naive reference.
"""
import numpy as np
import pytest

from easyepi.motifs import scan_sequences

_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _brute_best(pwm, seq, reverse_complement=True):
    """Best window score by explicit loops; the reference implementation."""
    L = pwm.shape[1]
    strands = [seq]
    if reverse_complement:
        strands.append("".join(_COMP.get(c, "N") for c in reversed(seq)))
    best = -np.inf
    for s in strands:
        for i in range(len(s) - L + 1):
            w = s[i:i + L]
            if any(c not in _IDX for c in w):
                continue
            score = sum(pwm[_IDX[w[j]], j] for j in range(L))
            if score > best:
                best = score
    return best


def _pwm(seed, length):
    rng = np.random.default_rng(seed)
    p = rng.random((4, length))
    return np.log2(p / p.sum(axis=0, keepdims=True) + 1e-9)


@pytest.mark.parametrize("length", [4, 8, 10, 16])
@pytest.mark.parametrize("seqlen", [20, 60, 200])
def test_best_score_matches_exhaustive(length, seqlen):
    pwm = _pwm(length, length)
    rng = np.random.default_rng(seqlen + length)
    seqs = ["".join(rng.choice(list("ACGT"), seqlen)) for _ in range(5)]
    got = scan_sequences(pwm, seqs)["scores"]
    for s, g in zip(seqs, got):
        assert abs(float(g) - _brute_best(pwm, s)) < 1e-9, (length, seqlen, s)


def test_reverse_complement_can_win():
    """A sequence whose only good match is on the minus strand must be found."""
    pwm = _pwm(0, 8)
    core = "".join("ACGT"[int(np.argmax(pwm[:, j]))] for j in range(8))
    rc = "".join(_COMP[c] for c in reversed(core))
    flank = "TTTTTTTTTT"
    seq = flank + rc + flank
    plus = _brute_best(pwm, seq, reverse_complement=False)
    both = float(scan_sequences(pwm, [seq])["scores"][0])
    assert both > plus, "minus-strand hit not found"
    assert abs(both - _brute_best(pwm, seq)) < 1e-9


def test_strand_and_position_are_consistent():
    pwm = _pwm(1, 8)
    core = "".join("ACGT"[int(np.argmax(pwm[:, j]))] for j in range(8))
    seq = "TTTTT" + core + "TTTTT"
    res = scan_sequences(pwm, [seq])
    assert res["strands"][0] == "+"
    assert res["positions"][0] == 5


def test_ambiguous_base_is_never_a_hit():
    """A window containing N must not be scored, even if the rest matches well."""
    pwm = _pwm(2, 6)
    core = "".join("ACGT"[int(np.argmax(pwm[:, j]))] for j in range(6))
    broken = core[:3] + "N" + core[4:]
    assert not np.isfinite(float(scan_sequences(pwm, [broken])["scores"][0]))


def test_lowercase_is_treated_as_the_base():
    pwm = _pwm(3, 6)
    core = "".join("ACGT"[int(np.argmax(pwm[:, j]))] for j in range(6))
    upper = float(scan_sequences(pwm, [core])["scores"][0])
    lower = float(scan_sequences(pwm, [core.lower()])["scores"][0])
    assert abs(upper - lower) < 1e-12


def test_shorter_than_motif_is_minus_infinity():
    pwm = _pwm(4, 10)
    assert not np.isfinite(float(scan_sequences(pwm, ["ACGT"])["scores"][0]))


def test_empty_sequence_list():
    pwm = _pwm(5, 8)
    assert scan_sequences(pwm, [])["scores"].size == 0
