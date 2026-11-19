"""PWM parsing and scanning, optional external motif tools and legacy utilities.

Native discovery helpers are exploratory and have no validated equivalence
to HOMER or MEME. They remain available as low-level compatibility functions.
"""

from __future__ import annotations

from collections import Counter
from itertools import product
from pathlib import Path

import numpy as np

from easyepi.motifs import formats  # noqa: F401
from easyepi.motifs import homer_like  # noqa: F401
from easyepi.motifs import logo  # noqa: F401

# Relative information-content (bits) for base frequencies
_BASES = "ACGT"
_BASE_TO_I = {b: i for i, b in enumerate(_BASES)}


def read_pwm(path: Path, pseudocount: float = 0.01) -> np.ndarray:
    """Parse a JASPAR-style position frequency matrix (4 x L) into log-odds.

    The file should have 4 rows (A, C, G, T) of counts. Rows are mapped ACGT.
    Returns a (4, L) log-odds matrix.
    """
    rows = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith((">", "#")):
                continue
            fields = line.replace("[", " ").replace("]", " ").replace(",", " ").split()
            label = fields.pop(0).upper() if fields[0].upper() in _BASES else _BASES[len(rows)] if len(rows) < 4 else "extra"
            if label in rows or label == "extra":
                raise ValueError("Expected exactly one row for each of A, C, G, T")
            rows[label] = [float(value) for value in fields]
    if set(rows) != set(_BASES):
        raise ValueError("Expected four PWM rows (A,C,G,T)")
    counts = np.array([rows[base] for base in _BASES], dtype=float)
    if counts.ndim != 2 or counts.shape[1] == 0 or not np.isfinite(counts).all() or (counts < 0).any() or pseudocount <= 0:
        raise ValueError("PWM counts must be finite, nonnegative and nonempty; pseudocount must be positive")
    freq = counts + pseudocount
    freq = freq / freq.sum(axis=0, keepdims=True)
    bg = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
    with np.errstate(divide="ignore"):
        log_odds = np.log2(freq / bg[:, None])
    return log_odds  # shape (4, L)


def scan_sequences(pwm: np.ndarray, seqs, reverse_complement: bool = True,
                   max_hits: int = 5) -> dict:
    """Scan sequences (list of uppercase ACGT strings) against a PWM.

    Vectorized with numpy sliding-window so it runs tens of times faster than a
    per-position Python loop (the bottleneck in known-motif enrichment).

    Returns a dict with best score per sequence and position, plus a running
    list of matched motifs (score, pos, strand).
    """
    pwm = np.asarray(pwm, dtype=np.float64)
    if pwm.ndim != 2 or pwm.shape[0] != 4 or pwm.shape[1] == 0 or not np.isfinite(pwm).all():
        raise ValueError("PWM must be a finite (4, length) matrix")
    L = pwm.shape[1]
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rc = lambda s: "".join(comp.get(c, "N") for c in reversed(s))  # noqa: E731

    scores = []
    pos_list = []
    strand_list = []
    for seq in seqs:
        seq = str(seq).upper()
        best = -np.inf
        best_p, best_strand = -1, "+"
        for strand, s in ((("+", seq), ("-", rc(seq))) if reverse_complement else (("+", seq),)):
            # encode to int indices; mask invalid (N) positions
            idx = np.array([_BASE_TO_I.get(c, -1) for c in s], dtype=np.int64)
            if len(idx) < L:
                continue
            valid = idx != -1
            # sliding windows of length L → (n_windows, L)
            windows = _sliding_windows(idx, L)
            # gather PWM scores for each window: (n_windows, L) → sum
            # map idx -> pwm value; invalid (-1) contributes 0 (masked below)
            wm = windows.copy()
            wm[wm == -1] = 0
            window_scores = pwm[wm, np.arange(L)].sum(axis=1)
            # mask windows containing an invalid base
            window_valid = _sliding_windows(valid.astype(np.int64), L).min(axis=1).astype(bool)
            window_scores[~window_valid] = -np.inf
            if window_scores.size:
                pos = int(np.argmax(window_scores))
                sc = float(window_scores[pos])
                if sc > best:
                    best, best_p, best_strand = sc, (pos if strand == "+" else len(seq) - pos - L), strand
        scores.append(best)
        pos_list.append(best_p)
        strand_list.append(best_strand)
    return {"scores": np.array(scores), "positions": pos_list,
            "strands": strand_list, "motif_len": int(L)}


def _sliding_windows(arr: np.ndarray, L: int) -> np.ndarray:
    """Return a (n_windows, L) view of all length-L windows of a 1-D array."""
    if arr.size < L:
        return np.empty((0, L), dtype=arr.dtype)
    view = np.lib.stride_tricks.sliding_window_view(arr, L)
    return view.reshape(-1, L)


def discover_motifs(seqs, k: int = 8, max_motifs: int = 5, iters: int = 40,
                    pseudocount: float = 0.1, reverse_complement: bool = True) -> dict:
    """De novo motif discovery with a MEME-style EM algorithm.

    Runs a k-mer-seeded Expectation-Maximization: pick the most enriched k-mer
    as a seed, then iterate (assign best site per sequence -> update the PWM ->
    add pseudocount / renormalize) until the information content stabilizes.

    Returns:
        dict with per-motif PWM (4 x k log-odds), information content (bits),
        number of sites, seed k-mer, and a consensus string.
    """
    seqs = [str(s).upper() for s in seqs if str(s).strip()]
    from itertools import product
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    rcmap = lambda s: "".join(comp.get(c, "N") for c in reversed(s))  # noqa: E731

    # count k-mers (both strands) as seeds
    counts: Counter = Counter()
    for s in seqs:
        for strand_seq in (s, rcmap(s) if reverse_complement else s):
            for i in range(len(strand_seq) - k + 1):
                km = strand_seq[i:i + k]
                if any(c not in _BASE_TO_I for c in km):
                    continue
                counts[km] += 1
    seeds = [m for m, _ in counts.most_common(max_motifs)]

    results = []
    for seed in seeds:
        # init PWM (as counts)
        pw = np.ones((4, k), dtype=np.float64) * pseudocount
        for j, base in enumerate(seed):
            pw[_BASE_TO_I[base], j] += 1.0
        seqs_clean = [s for s in seqs if len(s) >= k]
        if len(seqs_clean) < 3:
            continue
        last_ic = 0.0
        for _ in range(iters):
            freq = pw + pseudocount
            freq = freq / freq.sum(axis=0, keepdims=True)
            # E-step: for each seq find the best-position k-mer for this PWM
            assigned = []
            for s in seqs_clean:
                best_s, best_sc, best_strand = None, -np.inf, "+"
                for strand_flag, ss in ((1, s), (-1, rcmap(s) if reverse_complement else s)):
                    for i in range(len(ss) - k + 1):
                        km = ss[i:i + k]
                        if any(c not in _BASE_TO_I for c in km):
                            continue
                        mat = np.array([_BASE_TO_I[c] for c in km])
                        sc = float(np.log(freq[mat, np.arange(k)]).sum())
                        if sc > best_sc:
                            best_s, best_sc, best_strand = km, sc, strand_flag
                if best_s:
                    assigned.append(best_s)
            # M-step: rebuild PWM from assigned sites
            new_pw = np.ones((4, k), dtype=np.float64) * pseudocount
            for km in assigned:
                for j, base in enumerate(km):
                    new_pw[_BASE_TO_I[base], j] += 1.0
            pw = new_pw
            freq = pw / pw.sum(axis=0, keepdims=True)
            ic = float((freq * np.log2(freq / 0.25)).sum())
            if abs(ic - last_ic) < 1e-4:
                break
            last_ic = ic
        log_odds = np.log2((freq + pseudocount) / 0.25)
        consensus = "".join(_BASES[int(a)] for a in freq.argmax(axis=0))
        results.append({
            "seed": seed, "consensus": consensus, "information_content_bits": round(ic, 3),
            "n_sites": len(assigned), "pwm": log_odds.tolist(),
        })
    return {"k": int(k), "n_seqs": len(seqs), "motifs": results}


def run_motifs(args) -> int:
    """Run motif discovery or known-motif enrichment on the supplied sequences."""
    from easyepi.motifs.homer_workflow import run_homer_workflow
    import json
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # collect sequences from a FASTA (or a plain single-column file)
    fasta = getattr(args, "fasta", None)
    background = getattr(args, "background", None)
    pwm_path = getattr(args, "pwm", None)
    motif_db = getattr(args, "db", None)  # a HOMER/JASPAR/MEME motif file
    k = getattr(args, "k", 8)
    backend = getattr(args, "backend", "auto")
    homer_dir = getattr(args, "homer_dir", None)

    # An explicitly supplied PWM always selects scanning.
    if fasta and not pwm_path and (background or backend in ("homer", "auto")):
        # If bg missing, we still allow HOMER's own background (or native with
        # target-as-null); pass bg_fasta=None and let the workflow handle it.
        bg_fasta = background if background else None
        res = run_homer_workflow(
            peak_fasta=fasta, bg_fasta=bg_fasta, out_dir=str(out),
            k=k, backend=backend, homer_dir=homer_dir,
            homer_known=motif_db)  # motif_db used as known library if given
        print(f"[easyepi] motif workflow ({res.get('backend')}) written to {out}")
        return 0

    # --- legacy/simple paths (single PWM scan, or simple known enrichment) ---
    seqs = _read_sequences(fasta) if fasta else []
    bg_seqs = _read_sequences(background) if background else []

    res = {"status": "ok"}
    if pwm_path and seqs:
        pwm = read_pwm(Path(pwm_path))
        res.update(scan_sequences(pwm, seqs))
        res["n_seqs"] = len(seqs)
        res["mode"] = "scan"
    elif motif_db and seqs:
        motifs = formats.read_motifs(Path(motif_db))
        if not bg_seqs:
            raise ValueError("Known motif enrichment requires an explicit background FASTA")
        bg = bg_seqs
        rows = homer_like.known_motif_enrichment(seqs, bg or seqs, motifs)
        res.update({"mode": "known_enrichment", "n_seqs": len(seqs), "motifs": rows})
    elif bg_seqs and seqs:
        hits = homer_like.discover_homeric(seqs, bg_seqs, k=k)
        res.update({"mode": "de_novo", "k": int(k), "n_target": len(seqs),
                    "n_background": len(bg_seqs),
                    "motifs": [{"kmer": h.kmer, "consensus": h.consensus,
                                "pvalue": h.pvalue, "log_odds": h.log_odds,
                                "ic_bits": h.ic_bits, "pct_target": 100 * h.target_occur / h.n_target,
                                "pct_background": 100 * h.bg_occur / h.n_bg,
                                "pwm": h.pwm.tolist()} for h in hits]})
    else:
        raise ValueError("Provide target FASTA and a PWM or an explicit supported motif backend")
    out.write_text(json.dumps(res, indent=2, default=lambda x: x.tolist() if hasattr(x, "tolist") else x), encoding="utf-8")
    print(f"[easyepi] motif analysis written to {out}")
    return 0


def _read_sequences(path) -> list:
    """Read sequences from FASTA or a plain one-sequence-per-line file."""
    seqs = []
    with open(path) as fh:
        buf = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf))
                buf = []
            else:
                buf.append(line.upper())
        if buf:
            seqs.append("".join(buf))
    return seqs
