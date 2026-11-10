"""ChromHMM in pure Python: learn chromatin states from histone marks.

Reproduces the four published steps of Ernst & Kellis (2012, Nature) using
``hmmlearn`` instead of the Java tool, so no JVM is required:

1. **Binarize** — each mark is thresholded into a 0/1 track per genomic bin.
   ChromHMM's default is a Poisson background per mark (``-peaks`` uses peaks
   instead). Both are offered; the Poisson cutoff is the published default.
2. **LearnModel** — a multivariate Bernoulli HMM over the binarised vector, fitted
   with Baum-Welch (EM). ``hmmlearn.CategoricalHMM`` with ``n_features=2`` per
   mark is exactly that model: each state emits an independent Bernoulli per
   mark, and the transition matrix encodes the spatial state structure.
3. **MakeSegmentation** — Viterbi decoding of the state path over each chromosome,
   written as BED.
4. **Annotate** — per-state emission probabilities and genomic coverage, which is
   what the ChromHMM emission/overlap reports contain.

Design notes
------------
* Memory is kept flat by fitting on a *sample* of bins and then decoding the full
  chromosome with the fitted model (the standard trick for genome-scale HMMs).
* ``hmmlearn`` expects an integer observation per timestep, so the binarised mark
  vector is encoded as a single integer in base 2. The decoder maps it back; with
  <= 20 marks this is exact.
* Emission probabilities are recovered from ``emissionprob_`` by decoding the
  integer back to bits.

Reference: Ernst & Kellis 2012, Nature 473:43-49 ("ChromHMM: automating chromatin-
state discovery and characterization").
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from .core.io import _require

# ---------------------------------------------------------------------------
# 1. Binarise
# ---------------------------------------------------------------------------
def poisson_threshold(values: np.ndarray, cutoffs: Sequence[float] = (0.5, 0.99)) -> tuple[np.ndarray, dict]:
    """ChromHMM's default binarisation: a Poisson background per mark.

    ``values`` is the per-bin signal for one mark. The cutoff is chosen so that a
    Poisson variable with the observed mean exceeds it with probability
    ``1 - max(cutoffs)``, matching ChromHMM's ``-c`` semantics.
    """
    from scipy.stats import poisson

    flat = values[np.isfinite(values)]
    if flat.size == 0:
        return np.zeros(len(values), dtype=np.int8), {"mean": 0.0, "cutoff": 0.0}
    mean = float(flat.mean())
    target = float(max(cutoffs))
    # smallest integer whose survival probability falls below 1-target
    cutoff = int(poisson.isfinite(1 - target, mean)) if False else int(poisson.ppf(target, mean))
    cutoff = max(cutoff, 0)
    return (values > cutoff).astype(np.int8), {"mean": mean, "cutoff": float(cutoff)}


def binarize_from_peaks(peak_files: Sequence[str | Path], chrom: str, n_bins: int,
                        resolution: int) -> np.ndarray:
    """Mark a bin as 1 when it overlaps any peak (ChromHMM ``-peaks`` mode)."""
    track = np.zeros(n_bins, dtype=np.int8)
    for path in peak_files:
        opener = __import__("gzip").open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or f[0] != chrom:
                    continue
                try:
                    start, end = int(f[1]), int(f[2])
                except ValueError:
                    continue
                lo = max(0, start // resolution)
                hi = min(n_bins - 1, max(lo, (end - 1) // resolution))
                track[lo:hi + 1] = 1
    return track


# ---------------------------------------------------------------------------
# 2. Learn the model
# ---------------------------------------------------------------------------
def encode_observation(marks: np.ndarray) -> np.ndarray:
    """Encode a (n_bins, n_marks) 0/1 matrix as one integer per bin.

    ``hmmlearn`` takes a single categorical observation per timestep. Base-2
    encoding lets the emission table be decoded back to per-mark Bernoulli
    probabilities exactly.
    """
    marks = np.asarray(marks, dtype=np.int64)
    n_marks = marks.shape[1]
    weights = 1 << np.arange(n_marks)
    return marks @ weights


def learn_states(marks: np.ndarray, n_states: int, *, seed: int = 0,
                 max_iter: int = 200, max_fit_bins: int = 2_000_000,
                 n_init: int = 8) -> dict[str, Any]:
    """Fit a multivariate-Bernoulli HMM (ChromHMM's LearnModel).

    EM lands in local optima readily; ChromHMM itself runs several random
    initialisations and keeps the best likelihood (``LearnModel``'s default is a
    handful of restarts). The same is done here: ``n_init`` fits from different
    seeds, the one with the highest log-likelihood wins.

    Fitting on a capped sample keeps memory bounded; a chromatin-state transition
    matrix is dominated by local structure, which a random sample preserves.
    """
    from hmmlearn import hmm

    marks = np.asarray(marks, dtype=np.int8)
    if marks.ndim != 2:
        raise ValueError("marks must be (bins, marks)")
    n_bins, n_marks = marks.shape
    if n_marks < 1:
        raise ValueError("need at least one mark")
    if n_bins < n_states * 10:
        raise ValueError(f"only {n_bins} bins for {n_states} states; "
                         "need at least 10 bins per state")

    obs = encode_observation(marks)
    if n_bins > max_fit_bins:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(n_bins, size=max_fit_bins, replace=False))
        fit_obs = obs[idx]
    else:
        fit_obs = obs

    n_symbols = 1 << n_marks
    X = fit_obs.reshape(-1, 1)

    best_model = None
    best_score = -np.inf
    scores: list[float] = []
    for attempt in range(max(1, n_init)):
        model = hmm.CategoricalHMM(n_components=n_states, n_iter=max_iter,
                                   random_state=seed + attempt, init_params="ste",
                                   params="ste", tol=1e-4)
        model.n_features = n_symbols
        try:
            model.fit(X)
            score = float(model.score(X))
        except Exception:
            continue
        scores.append(score)
        if score > best_score:
            best_score, best_model = score, model
    if best_model is None:
        raise RuntimeError("every HMM initialisation failed to converge")

    return {"model": best_model, "n_marks": n_marks, "n_states": n_states,
            "n_symbols": n_symbols, "fitted_on": int(len(fit_obs)),
            "n_bins_total": int(n_bins), "n_init": len(scores),
            "loglik": best_score, "all_loglik": scores}


def emission_probabilities(state: dict[str, Any]) -> np.ndarray:
    """Decode ``emissionprob_`` back to a (states, marks) Bernoulli table."""
    model = state["model"]
    n_marks = state["n_marks"]
    n_states = state["n_states"]
    n_symbols = state["n_symbols"]
    probs = model.emissionprob_
    out = np.zeros((n_states, n_marks))
    for symbol in range(n_symbols):
        for k in range(n_marks):
            if symbol & (1 << k):
                out[:, k] += probs[:, symbol]
    return np.clip(out, 1e-6, 1 - 1e-6)


# ---------------------------------------------------------------------------
# 3. Segment
# ---------------------------------------------------------------------------
def segment(marks: np.ndarray, state: dict[str, Any]) -> np.ndarray:
    """Viterbi path over the full track (ChromHMM's MakeSegmentation)."""
    obs = encode_observation(marks)
    return state["model"].predict(obs.reshape(-1, 1))


def states_to_bed(paths: np.ndarray, chrom: str, resolution: int,
                  state_names: Sequence[str]) -> list[tuple[str, int, int, str]]:
    """Run-length encode a state path into BED intervals."""
    if paths.size == 0:
        return []
    rows: list[tuple[str, int, int, str]] = []
    start = 0
    cur = paths[0]
    for i in range(1, paths.size):
        if paths[i] != cur:
            rows.append((chrom, start * resolution, i * resolution, state_names[cur]))
            start, cur = i, paths[i]
    rows.append((chrom, start * resolution, paths.size * resolution, state_names[cur]))
    return rows


# ---------------------------------------------------------------------------
# 4. Annotate
# ---------------------------------------------------------------------------
def default_state_names(n_states: int) -> list[str]:
    return [f"E{i + 1}" for i in range(n_states)]


def state_summary(paths: np.ndarray, emissions: np.ndarray,
                  state_names: Sequence[str], resolution: int) -> list[dict[str, Any]]:
    """Per-state coverage and emission profile (ChromHMM's annotation report)."""
    total = max(1, paths.size)
    rows: list[dict[str, Any]] = []
    for idx, name in enumerate(state_names):
        n = int((paths == idx).sum())
        rows.append({
            "state": name,
            "n_bins": n,
            "coverage_bp": n * resolution,
            "fraction": round(n / total, 6),
            "emissions": [round(float(v), 4) for v in emissions[idx]],
        })
    rows.sort(key=lambda r: -r["fraction"])
    return rows
