"""1:1 reimplementation of TOBIAS ``ATACorrect``.

Tn5 inserts preferentially into certain sequences, so raw cut-site tracks carry
sequence bias on top of the biology. ATACorrect learns a k-mer preference model
from the data itself and subtracts the part of each track that the preference
predicts; the residual is the corrected signal.

Reference, read line by line:

* ``tobias/utils/sequences.pyx`` — ``NucleotideMatrix`` / ``DiNucleotideMatrix``:
  ``add_sequence``, ``add_background``, ``prepare_mat``, ``score_sequence``
* ``tobias/tools/atacorrect_functions.py`` — ``AtacBias``, ``bias_estimation``,
  ``bias_correction`` (incl. ``relu``)
* ``tobias/tools/atacorrect.py`` — driver, ``--norm-off`` / correction factor
* ``tobias/parsers.py`` — defaults: ``k_flank=12``, ``bg_shift=100``,
  ``read_shift=[4,-5]``, ``window=100``, ``score_mat=DWM``

Behaviours that make or break agreement:

1. **Two count matrices per strand.** ``counts`` is incremented at a read's cut
   site; ``bg_counts`` at the same k-mer shifted ``bg_shift`` bp upstream so the
   background is never inside the fragment. A k-mer containing ``N`` is dropped
   from **both** (the reference early-returns in each).
2. **The default model is a dinucleotide weight matrix (DWM)**, not a PWM.
   Counts are 4-D ``(bias_base, bg_base, bias_pos, bg_pos)`` and every position
   pair is incremented, so after dividing by ``L`` each position holds the
   per-position nucleotide profile.
3. **PSSM uses a mean, not a total:**
   ``no_bias = mean(sum(counts[:4], axis=0))``, then
   ``pwm = (counts+1)/(no_bias+4)`` and ``pssm = log2(pwm_bias / pwm_bg)``.
4. **The DWM is pseudocounted toward the PWM with a floor of 16:**
   ``pseudo = max(16, no_bias)`` and
   ``dwm[Sb,Sg,m,n] = (counts[Sb,Sg,m,n] + pseudo*pwm[Sa,n]*pwm[Sb,m]) / (no_bias+pseudo)``.
5. **DWM scoring is a conditional-probability ratio**, not a sum of log-odds:
   for each position ``n`` the code normalises over the four possible bases
   (``Cn``) and accumulates ``proba = cond_log + pwm[Sn,n] - log2(Cn)``, then the
   score is ``bias_log - bg_log``. Scores land at the k-mer centre
   (``scores[i + L//2]``), flanks stay ``nan``, any k-mer with ``N`` is skipped.
6. **Per-window ReLU fit.** The uncorrected signal is regressed on the bias score
   through ``relu(x,a,b) = max(0, a*x + b)``; if the fit does not converge the
   fallback shifts the bias so its minimum at cut positions is 0.
7. **Expectation is a share of the local total:** ``expected = signal_sum *
   (bias / bias_sum)`` over ``window`` bp.
8. **Correction and rescaling:** ``corrected = uncorrected - expected``, then
   positive values are scaled up so the corrected mass matches the original
   rolling sum.

This module is the numerical core. Reading cut sites from a BAM, writing bigWig
tracks and the multiprocessing driver stay with the existing easyepi code, so the
algorithm can be verified against the reference without a genome.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

__all__ = [
    "NucleotideMatrix",
    "DiNucleotideMatrix",
    "SequenceMatrix",
    "AtacBias",
    "relu",
    "bias_correction_tracks",
    "correct_factor",
    "encode_sequence",
    "TOBIAS_ATACCORRECT_DEFAULTS",
]

#: Reference defaults from ``tobias/parsers.py``.
TOBIAS_ATACCORRECT_DEFAULTS = {
    "k_flank": 12,
    "bg_shift": 100,
    "read_shift": (4, -5),
    "window": 100,
    "score_mat": "DWM",
    "extend": 100,
}

#: base encoding used throughout the reference: A=0, T=1, C=2, G=3, N=4
_BASE_CODE = {"A": 0, "T": 1, "C": 2, "G": 3, "N": 4}


def encode_sequence(seq: str) -> np.ndarray:
    """Encode DNA to the reference alphabet (A=0, T=1, C=2, G=3, N=4)."""
    return np.array([_BASE_CODE.get(ch, 4) for ch in str(seq).upper()], dtype=np.int64)


class SequenceMatrix:
    """Factory mirroring the reference's ``SequenceMatrix.create``."""

    @staticmethod
    def create(length: int, stype: str):
        if stype == "PWM":
            return NucleotideMatrix(length)
        if stype == "DWM":
            return DiNucleotideMatrix(length)
        raise ValueError("unknown score type %r (expected PWM or DWM)" % (stype,))


class NucleotideMatrix:
    """Reference ``NucleotideMatrix`` (``score_mat=PWM``)."""

    def __init__(self, length: int):
        self.length = int(length)
        self.counts = np.zeros((5, self.length))
        self.neg_counts = np.zeros((5, self.length))
        self.bg_counts = np.zeros((5, self.length))
        self.no_bias = 0.0
        self.no_bg = 0.0
        self.pssm = np.zeros((5, self.length))

    def add_sequence(self, sequence, amount: float = 1.0) -> None:
        """Add a k-mer to the bias counts; drop it entirely if it contains N."""
        seq = np.asarray(sequence, dtype=np.int64)
        if seq.shape[0] != self.length:
            raise ValueError("sequence length %d != matrix length %d"
                             % (seq.shape[0], self.length))
        if np.any(seq == 4):
            return
        target = self.counts if amount > 0 else self.neg_counts
        target[seq, np.arange(self.length)] += amount
        self.no_bias += amount

    def add_background(self, sequence, amount: float = 1.0) -> None:
        seq = np.asarray(sequence, dtype=np.int64)
        if np.any(seq == 4):
            return
        self.bg_counts[seq, np.arange(self.length)] += amount
        self.no_bg += amount

    def add_counts(self, other) -> None:
        self.counts += other.counts
        self.bg_counts += other.bg_counts
        self.neg_counts += other.neg_counts
        self.no_bias += other.no_bias
        self.no_bg += other.no_bg

    def prepare_mat(self) -> "NucleotideMatrix":
        # mean over positions of the summed A/T/C/G counts -- a mean, not a total
        self.no_bias = float(np.mean(np.sum(self.counts[:4, :], axis=0)))
        self.no_bg = float(np.mean(np.sum(self.bg_counts[:4, :], axis=0)))
        self.bias_pwm = np.true_divide(self.counts + 1, self.no_bias + 4)
        self.bg_pwm = np.true_divide(self.bg_counts + 1, self.no_bg + 4)
        with np.errstate(divide="ignore", invalid="ignore"):
            self.pssm = np.log2(np.true_divide(self.bias_pwm, self.bg_pwm))
        return self

    def score_sequence(self, sequence) -> np.ndarray:
        """Sum ``pssm[base, m]`` per window; place it at the window centre."""
        seq = np.asarray(sequence, dtype=np.int64)
        reglen = seq.shape[0]
        L = self.length
        F = int(L / 2.0)
        scores = np.full(reglen, np.nan, dtype=np.float64)
        for i in range(reglen - L + 1):
            window = seq[i:i + L]
            if np.any(window > 3):
                continue
            scores[i + F] = float(self.pssm[window, np.arange(L)].sum())
        return scores


class DiNucleotideMatrix(NucleotideMatrix):
    """Reference ``DiNucleotideMatrix`` (``score_mat=DWM``, the default).

    Counts every ordered position pair, so each position ends up holding the
    per-position nucleotide profile once the 4-D array is summed and divided by
    ``L``. Scoring uses conditional probabilities against a DWM pseudocounted
    toward the PWM.
    """

    def __init__(self, length: int):
        super().__init__(length)
        self.counts = np.zeros((5, 5, self.length, self.length))
        # the reference keeps neg_counts 2-D here purely so add_counts can join
        self.neg_counts = np.zeros((5, self.length))
        self.bg_counts = np.zeros((5, 5, self.length, self.length))
        self.bias_pwm = np.zeros((5, self.length))
        self.bg_pwm = np.zeros((5, self.length))
        self.bias_dwm = np.zeros((5, 5, self.length, self.length))
        self.bg_dwm = np.zeros((5, 5, self.length, self.length))

    def add_sequence(self, sequence, amount: float = 1.0) -> None:
        seq = np.asarray(sequence, dtype=np.int64)
        if np.any(seq == 4):
            return
        L = self.length
        for m in range(L):
            for n in range(L):
                self.counts[seq[m], seq[n], m, n] += amount
        self.no_bias += amount

    def add_background(self, sequence, amount: float = 1.0) -> None:
        seq = np.asarray(sequence, dtype=np.int64)
        if np.any(seq == 4):
            return
        L = self.length
        for m in range(L):
            for n in range(L):
                self.bg_counts[seq[m], seq[n], m, n] += amount
        self.no_bg += amount

    def prepare_mat(self) -> "DiNucleotideMatrix":
        """Build PWM and DWM exactly as the reference does."""
        L = self.length

        # divide by L because every position was added L times
        bias_cm = np.sum(self.counts, axis=(0, 2)) / float(L)
        self.no_bias = float(np.mean(np.sum(bias_cm[:4, :], axis=0)))
        self.bias_pwm = np.true_divide(bias_cm + 1, self.no_bias + 4)

        bg_cm = np.sum(self.bg_counts, axis=(0, 2)) / float(L)
        self.no_bg = float(np.mean(np.sum(bg_cm[:4, :], axis=0)))
        self.bg_pwm = np.true_divide(bg_cm + 1, self.no_bg + 4)

        with np.errstate(divide="ignore", invalid="ignore"):
            self.pssm = np.log2(np.true_divide(self.bias_pwm, self.bg_pwm))

        pseudo_bias = max(16.0, self.no_bias)
        pseudo_bg = max(16.0, self.no_bg)

        # dwm[Sb,Sg,m,n] = (counts + pseudo*pwm[Sg,n]*pwm[Sb,m]) / (no+pseudo)
        Sb = np.arange(5)[:, None, None, None]
        Sg = np.arange(5)[None, :, None, None]
        M = np.arange(L)[None, None, :, None]
        N = np.arange(L)[None, None, None, :]
        self.bias_dwm = ((self.counts + pseudo_bias
                          * self.bias_pwm[Sg, N] * self.bias_pwm[Sb, M])
                         / float(self.no_bias + pseudo_bias))
        self.bg_dwm = ((self.bg_counts + pseudo_bg
                        * self.bg_pwm[Sg, N] * self.bg_pwm[Sb, M])
                       / float(self.no_bg + pseudo_bg))

        with np.errstate(divide="ignore"):
            self.bias_pwm_log = np.log2(self.bias_pwm)
            self.bg_pwm_log = np.log2(self.bg_pwm)
            self.bias_dwm_log = np.log2(self.bias_dwm)
            self.bg_dwm_log = np.log2(self.bg_dwm)
        return self

    def _score_one(self, window: np.ndarray, dwm_log, pwm_log) -> float:
        """Conditional-probability score for a single k-mer (reference inner loop)."""
        L = self.length
        total = 0.0
        for n in range(L):
            Sn = window[n]
            if Sn > 3:
                return np.nan

            # normalise over the four possible bases at position n
            Cn = 0.0
            for a in range(4):
                prod_log = 0.0
                for m in range(L):
                    if m != n:
                        prod_log += dwm_log[window[m], a, m, n] - pwm_log[a, n]
                Cn += 2.0 ** (prod_log + pwm_log[a, n])

            cond_log = 0.0
            for m in range(L):
                if m != n:
                    cond_log += dwm_log[window[m], Sn, m, n] - pwm_log[Sn, n]

            total += cond_log + pwm_log[Sn, n] - np.log2(Cn)
        return total

    def score_sequence(self, sequence) -> np.ndarray:
        """Score each k-mer as ``bias_log - bg_log`` at the window centre."""
        seq = np.asarray(sequence, dtype=np.int64)
        reglen = seq.shape[0]
        L = self.length
        F = int(L / 2.0)
        scores = np.full(reglen, np.nan, dtype=np.float64)

        for i in range(reglen - L + 1):
            window = seq[i:i + L]
            if np.any(window > 3):
                continue
            bias_log = self._score_one(window, self.bias_dwm_log, self.bias_pwm_log)
            bg_log = self._score_one(window, self.bg_dwm_log, self.bg_pwm_log)
            if np.isnan(bias_log) or np.isnan(bg_log):
                continue
            scores[i + F] = bias_log - bg_log
        return scores


class AtacBias:
    """Container for the estimated bias of both strands (reference ``AtacBias``)."""

    def __init__(self, L: int = 10, stype: str = "DWM"):
        self.stype = stype
        self.bias = {s: SequenceMatrix.create(L, stype)
                     for s in ("forward", "reverse", "both")}
        self.no_reads = 0.0
        self.correction_factor = 1.0

    def join(self, other) -> "AtacBias":
        for s in ("forward", "reverse", "both"):
            self.bias[s].add_counts(other.bias[s])
        self.no_reads += other.no_reads
        return self


def relu(x, a, b):
    """Reference ``relu``: ``max(0, a*x + b)``."""
    return np.maximum(0.0, a * np.asarray(x, dtype=float) + b)


def correct_factor(reads_total: float, reads_peaks: float) -> float:
    """Reference normalisation ``lib_norm * (1/FRiP)``; ``--norm-off`` uses 1.0."""
    if reads_total <= 0:
        return 1.0
    lib_norm = 1e7 / float(reads_total)
    frip = float(reads_peaks) / float(reads_total)
    if frip <= 0:
        return 1.0
    return lib_norm * (1.0 / frip)


def fast_rolling_math(arr, w: int, operation: str = "sum") -> np.ndarray:
    """``fast_rolling_math`` from the reference: centred, NaN padded.

    The reference writes ``roll_arr[i + lf]`` where ``lf = floor(w/2)``, so the
    window sum is attributed to the CENTRE of the window, and the first ``lf`` and
    last ``ceil(w/2)`` entries stay NaN.
    """
    v = np.asarray(arr, dtype=np.float64)
    n = v.shape[0]
    w = int(w)
    out = np.full(n, np.nan, dtype=np.float64)
    lf = int(np.floor(w / 2.0))
    if w <= 0 or n < w:
        return out
    if operation == "sum":
        for i in range(n - w + 1):
            out[i + lf] = np.sum(v[i:i + w])
    elif operation == "mean":
        for i in range(n - w + 1):
            out[i + lf] = np.sum(v[i:i + w]) / w
    elif operation == "max":
        for i in range(n):
            lo = max(0, i - lf)
            hi = min(n, i - lf + w)
            out[i] = np.max(v[lo:hi])
    elif operation == "min":
        for i in range(n):
            lo = max(0, i - lf)
            hi = min(n, i - lf + w)
            out[i] = np.min(v[lo:hi])
    else:
        raise ValueError("unknown operation: %s" % operation)
    return out


_rolling_sum = fast_rolling_math


def _rolling_sum_window(values, window: int) -> np.ndarray:
    return fast_rolling_math(values, window, "sum")


def _fit_window_bias(signal_w: np.ndarray, bias_w: np.ndarray):
    """Fit ``signal ~ relu(bias)`` for one window, with the reference's fallback."""
    if not np.max(signal_w) > 0:
        return None
    try:
        popt, _ = curve_fit(relu, bias_w, signal_w)
        bias_predict = relu(bias_w, *popt)
    except (OptimizeWarning, RuntimeError, ValueError):
        cut_positions = np.logical_not(np.isclose(signal_w, 0))
        if not np.any(cut_positions):
            return None
        bias_predict = bias_w - np.min(bias_w[cut_positions])
        bias_predict[bias_predict < 0] = 0

    peak = np.max(bias_predict)
    return bias_predict / peak if peak > 0 else np.zeros(len(bias_w))


def bias_correction_tracks(uncorrected: np.ndarray, bias: np.ndarray,
                           k_flank: int, window: int = 100, step: int = 10,
                           correction_factor: float = 1.0):
    """Compute ``(expected, corrected)`` for one strand of one region.

    Transcribed from ``atacorrect_functions.bias_correction``. ``bias_prediction``
    is recomputed with ``nanmean`` inside the window loop, so positions not yet
    covered by any window are NaN and get no expected share -- which is exactly why
    ``expected`` is zero at the leading bases.
    """
    u = np.asarray(uncorrected, dtype=np.float64).copy()
    b = np.asarray(bias, dtype=np.float64).copy()
    reg_len = u.shape[0]
    k_flank = int(k_flank)
    window = int(window)
    step = int(step)

    reg_end = reg_len - k_flank
    overlaps = max(1, int(window / step))
    window_starts = list(range(k_flank, reg_end - window, step))
    window_ends = list(range(k_flank + window, reg_end, step))
    windows = list(zip(window_starts, window_ends))

    bias_predictions = np.zeros((overlaps, reg_len))
    # The reference writes `bias_predictions[k_flank:reg_end] = np.nan`, which slices
    # the ROW axis. With overlaps (window/step = 10) <= k_flank (12) that selects no
    # rows, so no cell is marked NaN and unwritten cells keep their 0.0. `nanmean`
    # then averages those zeros in, scaling the prediction by roughly 1/overlaps.
    # Reproduced here verbatim rather than "fixed", because the goal is 1:1.
    bias_predictions[k_flank:reg_end] = np.nan
    row = 0
    bias_prediction = np.full(reg_len, np.nan)
    for w0, w1 in windows:
        pred = _fit_window_bias(u[w0:w1], b[w0:w1])
        if pred is not None:
            bias_predictions[row, w0:w1] = pred
            row = row + 1 if row < overlaps - 1 else 0
            with np.errstate(invalid="ignore"):
                bias_prediction = np.nanmean(bias_predictions, axis=0)

    bias = bias_prediction

    # ---- expected signal ----
    signal_sum = fast_rolling_math(u, window, "sum")
    signal_sum[np.isnan(signal_sum)] = 0

    bias_sum = fast_rolling_math(bias, window, "sum")
    nulls = np.logical_or(np.isclose(bias_sum, 0), np.isnan(bias_sum))
    bias_sum[nulls] = 1
    bias_probas = bias / bias_sum
    bias_probas[nulls] = 0

    with np.errstate(invalid="ignore"):
        expected = signal_sum * bias_probas
    expected = np.where(np.isnan(expected), 0.0, expected)

    # ---- corrected signal ----
    u = u * correction_factor
    expected = expected * correction_factor
    corrected = u - expected

    # ---- rescale positive cuts to restore the uncorrected sum ----
    uncorrected_sum = fast_rolling_math(u, window, "sum")
    uncorrected_sum[np.isnan(uncorrected_sum)] = 0
    corrected_sum = fast_rolling_math(np.abs(corrected), window, "sum")
    corrected_sum[np.isnan(corrected_sum)] = 0

    corrected_pos = np.copy(corrected)
    corrected_pos[corrected_pos < 0] = 0
    corrected_pos_sum = fast_rolling_math(corrected_pos, window, "sum")
    corrected_pos_sum[np.isnan(corrected_pos_sum)] = 0
    corrected_neg_sum = corrected_sum - corrected_pos_sum

    zero_sum = corrected_pos_sum == 0
    corrected_pos_sum[zero_sum] = np.nan
    scale_factor = (uncorrected_sum - corrected_neg_sum) / corrected_pos_sum
    scale_factor[zero_sum] = 1
    scale_factor[scale_factor < 1] = 1
    pos_bool = corrected > 0
    corrected[pos_bool] *= scale_factor[pos_bool]

    return expected, corrected
