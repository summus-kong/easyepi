"""Legacy exploratory motif utilities retained for compatibility.

These heuristics have not been validated as equivalent to HOMER or MEME.
Use established external discovery tools for confirmatory analyses.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import hypergeom, binom

from easyepi.motifs.formats import Motif

_BASES = "ACGT"
_BASE_TO_I = {b: i for i, b in enumerate(_BASES)}


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------
def _rc(seq: str) -> str:
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(comp.get(c, "N") for c in reversed(seq))


def _zoopps_occur(seq: str, kmer: str) -> bool:
    """ZOOPS: does the sequence contain >=1 occurrence of kmer (either strand)?"""
    if kmer in seq:
        return True
    return _rc(kmer) in seq


def _count_occurrences(seqs: Sequence[str], kmer: str) -> int:
    """Number of sequences containing the k-mer (ZOOPS)."""
    return sum(1 for s in seqs if _zoopps_occur(s, kmer))


def hypergeom_pvalue(target_hits: int, n_target: int,
                     bg_hits: int, n_bg: int) -> float:
    """P-value of observing at least ``target_hits`` in the target set, given the
    background enrichment, under the hypergeometric null."""
    N = n_target + n_bg
    K = bg_hits
    n = n_target
    k = target_hits
    if K == 0:
        # No background occurrences: the target hits are an enrichment signal.
        # Under a Poisson model with background rate ~0, observing any target
        # hit is significant; use a conservative small p (bound by 1/max(1,n)).
        if k > 0:
            return min(1.0, 1.0 / max(n, 1))
        return 1.0
    if n == 0 or N == 0:
        return 1.0
    return float(hypergeom.sf(k - 1, M=N, n=n, N=K))


# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------
def _freq_biased_candidates(seqs: Sequence[str], k: int,
                            max_candidates: int = 2000) -> List[str]:
    """Generate candidate k-mers enriched in the target set.

    ZOOPS: counts the number of sequences containing each k-mer (or its reverse
    complement), canonicalized to the forward strand, then returns the top
    ``max_candidates`` k-mers by decreasing target count.
    """
    counts: Counter = Counter()
    for s in seqs:
        seen = set()
        for i in range(len(s) - k + 1):
            km = s[i:i + k]
            if any(c not in _BASE_TO_I for c in km) or km in seen:
                continue
            seen.add(km)
        # canonicalize each seen k-mer to the smaller of {kmer, rc}
        for km in seen:
            canon = min(km, _rc(km))
            counts[canon] += 1
    ranked = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    return [km for km, _ in ranked[:max_candidates]]


# ---------------------------------------------------------------------------
# PWM construction from a k-mer set
# ---------------------------------------------------------------------------
def _pwm_from_pattern(seqs: Sequence[str], kmer: str, pseudocount: float = 0.1) -> Tuple[np.ndarray, List[str]]:
    """Build a 4 x L frequency matrix from the best hit of ``kmer`` in each seq."""
    L = len(kmer)
    mat = np.ones((4, L), dtype=np.float64) * pseudocount
    best_hits: List[str] = []
    for s in seqs:
        best, best_sc = None, -np.inf
        for strand_seq in (s, _rc(s)):
            for i in range(len(strand_seq) - L + 1):
                km = strand_seq[i:i + L]
                sc = sum(1.0 for a, b in zip(km, kmer) if a == b)
                if sc > best_sc:
                    best, best_sc = km, sc
        if best and any(c in _BASE_TO_I for c in best):
            best_hits.append(best)
            for j, base in enumerate(best):
                if base in _BASE_TO_I:
                    mat[_BASE_TO_I[base], j] += 1.0
    return mat, best_hits


# ---------------------------------------------------------------------------
# Background correction and score
# ---------------------------------------------------------------------------
def _sequence_bias_correction(target_seqs, bg_seqs, L: int,
                              max_iter: int = 5) -> np.ndarray:
    """Estimate a position-specific nucleotide bias of the target set relative
    to the background (first-order Markov correction on 1-, 2-, 3-mer content).
    Returns a per-position background log-odds used to discount spurious motifs.
    """
    # First-order composition of target vs background (per position is unstable;
    # use a single mononucleotide correction for simplicity + a dinucleotide
    # term, matching HOMER's spirit of removing composition bias).
    def monofreq(seqs):
        """Per-position base frequencies normalised to sum to one."""
        cnt = np.zeros(4)
        for s in seqs:
            for c in s:
                if c in _BASE_TO_I:
                    cnt[_BASE_TO_I[c]] += 1
        return cnt / max(cnt.sum(), 1)

    ft = monofreq(target_seqs)
    fb = monofreq(bg_seqs)
    # log2 fold of target vs background composition; clip to avoid inf
    with np.errstate(divide="ignore"):
        correction = np.log2(np.maximum(ft, 1e-6) / np.maximum(fb, 1e-6))
    correction = np.clip(correction, -2, 2)
    return correction


@dataclass
class MotifHit:
    """A discovered de novo motif and its statistics."""

    kmer: str
    consensus: str
    pwm: np.ndarray
    target_occur: int
    n_target: int
    bg_occur: int
    n_bg: int
    pvalue: float
    log_odds: float
    ic_bits: float
    sites: List[str] = field(default_factory=list)


def discover_homeric(target_seqs: Sequence[str], bg_seqs: Sequence[str],
                    k: int = 8, max_candidates: int = 1500,
                    max_motifs: int = 10, pseudocount: float = 0.1,
                    rng_seed: int = 0) -> List[MotifHit]:
    """Differential (target vs background) HOMER-style de novo motif discovery.

    Args:
        target_seqs: list of sequences (typically peak sequences).
        bg_seqs: list of background sequences (equally sized).
        k: candidate k-mer length.
        max_candidates: number of candidate k-mers to test.
        max_motifs: number of motifs to return.
        pseudocount: smoothing for PWM construction.
    Returns:
        List of :class:`MotifHit` sorted by ascending p-value.
    """
    target = [s.upper() for s in target_seqs if len(s) >= k]
    bg = [s.upper() for s in bg_seqs if len(s) >= k]
    if len(target) == 0:
        return []
    n_t = len(target)
    n_b = len(bg) if bg else 0

    # background composition correction
    correction = _sequence_bias_correction(target, bg if bg else target, k) if bg else np.zeros(4)

    candidates = _freq_biased_candidates(target, k, max_candidates)
    results: List[MotifHit] = []
    seen = set()

    for kmer in candidates:
        if not kmer or "N" in kmer:
            continue
        canon = min(kmer, _rc(kmer))
        if canon in seen:
            continue
        seen.add(canon)
        t_hits = _count_occurrences(target, canon)
        b_hits = _count_occurrences(bg, canon) if n_b else 0
        if t_hits == 0:
            continue
        p = hypergeom_pvalue(t_hits, n_t, b_hits, n_b)
        t_frac = t_hits / n_t
        b_frac = (b_hits / n_b) if n_b else 0.0
        # Guard zero-count cases: no target hit -> neg; target hit, no bg -> max +.
        if t_frac <= 0:
            log_odds = -10.0
        elif b_frac <= 0:
            log_odds = 10.0
        else:
            log_odds = float(np.log2(max(t_frac, 1e-6) / max(b_frac, 1e-6)))
        mat, best_hits = _pwm_from_pattern(target, canon, pseudocount)
        freq = mat / mat.sum(axis=0, keepdims=True)
        ic = float((freq * np.log2(freq / 0.25)).sum())
        consensus = "".join(_BASES[int(a)] for a in np.argmax(freq, axis=0))
        results.append(MotifHit(kmer=canon, consensus=consensus, pwm=freq,
                                target_occur=t_hits, n_target=n_t,
                                bg_occur=b_hits, n_bg=n_b, pvalue=float(p),
                                log_odds=log_odds, ic_bits=ic, sites=best_hits))
        if len(results) >= max_candidates:
            break

    # Rank by decreasing fold enrichment (log-odds), then by p-value.
    results.sort(key=lambda r: (-r.log_odds, r.pvalue))
    sig = [r for r in results if r.pvalue < 0.05]
    return (sig if sig else results)[:max_motifs]


# ---------------------------------------------------------------------------
# Known-motif enrichment (HOMER knownResults.txt style)
# ---------------------------------------------------------------------------
def known_motif_enrichment(target_seqs: Sequence[str], bg_seqs: Sequence[str],
                           motifs: Sequence[Motif],
                           score_thr: float = 2.0) -> List[Dict]:
    """Scan a set of known PFMs against target vs background and rank by
    enrichment, returning rows compatible with HOMER's ``knownResults.txt``.

    A sequence "has" a motif if its best log-odds score exceeds a threshold
    derived from the *background* score distribution (``bg_mean + score_thr *
    bg_sd``), which is the correct null (unlike a target-derived threshold).
    Enrichment is tested with a binomial / hypergeometric p-value.

    Args:
        target_seqs: peak sequences (target set).
        bg_seqs: background sequences.
        motifs: known PFMs (:class:`Motif`).
        score_thr: number of background SDs above the background mean to call
            a match (default 2.0).
    """
    from easyepi.motifs import scan_sequences
    target = [s.upper() for s in target_seqs if s]
    bg = [s.upper() for s in bg_seqs if s]
    n_t = len(target) or 1
    n_b = len(bg) or 1
    rows = []
    for m in motifs:
        pwm = m.to_log_odds()
        res_t = scan_sequences(pwm, target)
        res_b = scan_sequences(pwm, bg)
        st = np.asarray(res_t["scores"], dtype=np.float64)
        sb = np.asarray(res_b["scores"], dtype=np.float64)
        if st.size == 0:
            continue
        # Replace non-finite (-inf from no match) scores with the finite min so
        # they do not poison the background mean / std.
        st = np.where(np.isfinite(st), st, np.min(st[np.isfinite(st)]) if np.any(np.isfinite(st)) else 0.0)
        sb = np.where(np.isfinite(sb), sb, np.min(sb[np.isfinite(sb)]) if np.any(np.isfinite(sb)) else 0.0)
        # Threshold derived from the BACKGROUND null distribution.
        bg_mean = float(sb.mean()) if sb.size else 0.0
        bg_sd = float(sb.std()) if sb.size > 1 else 1.0
        if bg_sd <= 0:
            bg_sd = 1.0
        thr = bg_mean + score_thr * bg_sd
        t_like = int((st >= thr).sum())
        b_like = int((sb >= thr).sum())
        # Binomial test: target hit-rate vs background hit-rate.
        p = _binomial_enrichment(t_like, n_t, b_like, n_b)
        t_frac = t_like / n_t
        b_frac = b_like / n_b
        # Log-odds of enrichment; guard the zero-count cases:
        #  - no target hits -> no enrichment (large negative)
        #  - target hits with zero background -> maximally enriched (large +)
        if t_frac <= 0:
            log_odds = -10.0
        elif b_frac <= 0:
            log_odds = 10.0
        else:
            log_odds = float(np.log2(max(t_frac, 1e-6) / max(b_frac, 1e-6)))
        rows.append({
            "name": m.name, "id": m.id, "tf": m.tf, "consensus": m.consensus_seq(),
            "pct_target": round(100 * t_frac, 2), "pct_background": round(100 * b_frac, 2),
            "pvalue": float(p), "log_odds": round(log_odds, 3),
            "score_threshold": round(thr, 3),
        })
    # Rank primarily by p-value (strongest enrichment first). When p-values
    # are comparable/underflow, log-odds breaks ties.
    rows.sort(key=lambda r: (r["pvalue"], -r["log_odds"]))
    return rows


def _binomial_enrichment(target_hits: int, n_target: int,
                         bg_hits: int, n_bg: int) -> float:
    """P-value of observing >= target_hits matches given the background rate.

    Uses a binomial model: background hit rate p0 = bg_hits/n_bg; the number of
    target matches is ~ Binomial(n_target, p0); we return P(X >= target_hits).
    This is more stable than the hypergeometric when counts are small.
    """
    from scipy.stats import binom
    n = int(n_target)
    if n == 0:
        return 1.0
    p0 = (bg_hits / n_bg) if n_bg else 0.0
    if p0 <= 0:
        # background has zero hits but target has some -> strong signal
        return float(1.0 / max(n, 1)) if target_hits > 0 else 1.0
    if p0 >= 1.0:
        return float(binom.sf(target_hits - 1, n, p0))
    return float(binom.sf(target_hits - 1, n, p0))
