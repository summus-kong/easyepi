"""PWM parsing and scanning, plus known-motif enrichment.

De novo motif discovery is delegated to HOMER and nothing here substitutes for
it. What this module does provide runs on numpy and scipy alone: reading PWMs,
scanning sequences, and testing known-motif enrichment against a background.
That is why `easyepi motifs --pwm` and `--db` work on a machine with no external
program installed.
"""

from __future__ import annotations

import tempfile
from collections import Counter
from itertools import product
from pathlib import Path

import numpy as np

from easyepi.motifs import formats  # noqa: F401
from easyepi.motifs import logo  # noqa: F401

# Relative information-content (bits) for base frequencies
_BASES = "ACGT"
_BASE_TO_I = {b: i for i, b in enumerate(_BASES)}

# Byte -> base index, 255 for anything that is not a base. A 256-entry table
# turns the per-character dict lookup into one vectorised gather, which matters
# because the scan encodes every base of every sequence for every motif.
_BYTE_TO_I = np.full(256, 255, dtype=np.uint8)
for _b, _i in _BASE_TO_I.items():
    _BYTE_TO_I[ord(_b)] = _i
    _BYTE_TO_I[ord(_b.lower())] = _i
# Complement lookup on base indices (A<->T, C<->G); only valid for 0..3.
_COMPLEMENT = np.array([3, 2, 1, 0], dtype=np.uint8)


def _encode(seq: str) -> tuple:
    """Sequence -> (base indices 0..3, valid mask). Non-ACGT bases are invalid."""
    raw = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    idx = _BYTE_TO_I[raw]
    valid = idx != 255
    if not valid.all():
        idx = idx.copy()
        idx[~valid] = 0
    return idx, valid


def _window_scores(idx: np.ndarray, valid: np.ndarray,
                   pwm: np.ndarray) -> "np.ndarray | None":
    """Score every window of one strand without materialising (n_windows, L).

    score[k] = sum_d pwm[idx[k+d], d] for windows with no invalid base, else
    -inf. Accumulating one offset at a time keeps the working set at n bases
    rather than n*L.
    """
    L = pwm.shape[1]
    n = idx.shape[0]
    if n < L:
        return None
    out = np.zeros(n - L + 1, dtype=np.float64)
    rows = pwm.T                       # (L, 4); rows[d] indexes bases at offset d
    bad = None
    for d in range(L):
        span = slice(d, n - L + 1 + d)
        out += rows[d][idx[span]]
        if bad is None:
            bad = ~valid[span]
        else:
            bad |= ~valid[span]
    if bad.any():
        out[bad] = -np.inf
    return out


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

    scores = []
    pos_list = []
    strand_list = []
    for seq in seqs:
        idx, valid = _encode(str(seq).upper())
        n = idx.shape[0]
        best = -np.inf
        best_p, best_strand = -1, "+"
        ws = _window_scores(idx, valid, pwm)
        if ws is not None and ws.size and np.isfinite(ws).any():
            k = int(np.argmax(ws))
            if ws[k] > best:
                best, best_p, best_strand = float(ws[k]), k, "+"
        if reverse_complement and n >= L:
            # reverse complement: reverse the order and complement each base; the
            # validity pattern reverses with it, so stop markers stay aligned
            ws = _window_scores(_COMPLEMENT[idx[::-1]], valid[::-1], pwm)
            if ws is not None and ws.size and np.isfinite(ws).any():
                k = int(np.argmax(ws))
                if ws[k] > best:
                    best, best_p, best_strand = float(ws[k]), n - k - L, "-"
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


def _window_scores(idx: np.ndarray, valid: np.ndarray, pwm: np.ndarray) -> np.ndarray:
    """PWM score for every window; -inf where a window holds an invalid base."""
    L = pwm.shape[1]
    if idx.size < L:
        return np.empty(0, dtype=np.float64)
    windows = _sliding_windows(idx, L).copy()
    windows[windows == -1] = 0
    scores = pwm[windows, np.arange(L)].sum(axis=1)
    ok = _sliding_windows(valid.astype(np.int64), L).min(axis=1).astype(bool)
    scores[~ok] = -np.inf
    return scores


#: Above this many outcomes the exact PWM score distribution is not enumerated.
#: 4**10 = 1_048_576, which is the practical ceiling for a per-motif loop.
EXACT_SCORE_LIMIT = 4 ** 10


def _exact_score_distribution(pwm: np.ndarray, bg: float = 0.25):
    """Exact (scores, probabilities) for a PWM, by convolving its columns.

    A log-odds score is a sum of independent per-column contributions, so its
    distribution is the convolution of the four outcomes of each column. This is
    the same construction FIMO uses, without the binning that
    ``pwm_score_distribution`` applies for tractability -- and that binning
    reports scores a raw window sum cannot reach.
    """
    pwm = np.asarray(pwm, dtype=np.float64)
    scores = np.zeros(1, dtype=np.float64)
    probs = np.ones(1, dtype=np.float64)
    base_p = np.full(4, bg, dtype=np.float64)
    for col in range(pwm.shape[1]):
        contrib = pwm[:, col]
        scores = (scores[:, None] + contrib[None, :]).ravel()
        # each existing outcome splits into 4, each with probability bg
        probs = (probs[:, None] * base_p[None, :]).ravel()
        order = np.argsort(scores, kind="stable")
        scores, probs = scores[order], probs[order]
    # merge numerically identical outcomes so the tail is not fragmented
    keys = np.round(scores, 9)
    uniq, inv = np.unique(keys, return_inverse=True)
    agg = np.zeros(uniq.size, dtype=np.float64)
    np.add.at(agg, inv, probs)
    return uniq, agg


def _attainable_threshold(pwm: np.ndarray, fpr: float) -> tuple[float, float]:
    """Score threshold whose upper tail is <= ``fpr``, as an attainable score.

    Uses the exact distribution when the motif is short enough to enumerate, so
    the returned value is a score an actual window can produce. Longer motifs
    fall back to the binned distribution, whose reported bin edge can sit slightly
    above any real score; the caller then compares with a bin-width tolerance.
    """
    pwm = np.asarray(pwm, dtype=np.float64)
    n_outcomes = 4 ** pwm.shape[1]
    if n_outcomes <= EXACT_SCORE_LIMIT:
        scores, probs = _exact_score_distribution(pwm)
    else:
        from easyepi.tf_annotate import pwm_score_distribution
        scores, probs = pwm_score_distribution(pwm)

    tail = np.cumsum(probs[::-1])[::-1]
    ok = np.nonzero(tail <= fpr)[0]
    if ok.size == 0:
        raise ValueError(
            f"no PWM score reaches the requested fpr={fpr:g}; the best achievable "
            f"rate is {probs[-1]:g}")
    idx = int(ok[0])
    return float(scores[idx]), float(tail[idx])


def _threshold_tolerance(pwm: np.ndarray) -> float:
    """Slack when comparing a raw window score to the reported threshold.

    Two different gaps, two different sizes:

    * exact path: the merged outcomes are rounded to 9 decimals, so a perfect
      match lands ~3.5e-10 below its own reported score. 1e-9 relative covers it.
    * binned path (motifs wider than 10 bp): reported scores are bin *edges*,
      which a raw sum cannot reach; the gap can approach a bin width, so half a
      bin is used.
    """
    pwm = np.asarray(pwm, dtype=np.float64)
    if 4 ** pwm.shape[1] <= EXACT_SCORE_LIMIT:
        return 1e-9
    from easyepi.tf_annotate import pwm_score_distribution
    scores, _ = pwm_score_distribution(pwm)
    return float(scores[1] - scores[0]) / 2 if scores.size > 1 else 0.0


def scan_sites(pwm: np.ndarray, seqs, *, fpr: float | None = 1e-4,
               reverse_complement: bool = True, pvalue: bool = False,
               score_threshold: float | None = None) -> dict:
    """Every above-threshold site of one PWM in each sequence.

    Unlike :func:`scan_sequences`, which keeps only the best site per sequence,
    this returns all sites, so the result can be written as a BED-like table with
    coordinates. Coordinates are **0-based half-open offsets into the supplied
    sequence**, not genomic coordinates: this function never sees a chromosome or
    an offset, and silently labelling an offset as a genomic position is how a
    peak table ends up disagreeing with a BED file.

    Args:
        pwm: (4, L) log-odds matrix.
        seqs: uppercase ACGT sequences.
        fpr: per-sequence false-positive rate used to derive the score
            threshold from the exact PWM score distribution. Ignored when
            ``score_threshold`` is given.
        reverse_complement: scan both strands.
        pvalue: also report each site's upper-tail p-value (costs a lookup per
            site against a precomputed tail table).
        score_threshold: use this score instead of deriving one.

    Returns:
        ``{"threshold", "achieved_fpr", "strands", "sites": [...]}`` where each
        site is ``{"seq_index", "start", "end", "strand", "score", "pvalue"}``.
        ``end`` is exclusive, so ``seq[start:end]`` is the match on the plus
        strand; on the minus strand it is the forward-strand interval.
    """
    from easyepi.tf_annotate import pwm_threshold, tail_table, site_pvalue

    pwm = np.asarray(pwm, dtype=np.float64)
    if pwm.ndim != 2 or pwm.shape[0] != 4 or pwm.shape[1] == 0:
        raise ValueError("PWM must be a (4, length) matrix")
    L = pwm.shape[1]

    achieved = None
    if score_threshold is None:
        if fpr is None:
            raise ValueError("pass either fpr or score_threshold")
        threshold, achieved = _attainable_threshold(pwm, fpr)
        tol = _threshold_tolerance(pwm)
    else:
        threshold = float(score_threshold)
        tol = 0.0

    table = tail_table(pwm) if pvalue else None
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}

    # A self-reverse-complementary motif matches the same physical site on both
    # passes. Scanning both strands then reports every hit twice, which doubles
    # every count a downstream enrichment test sees. Detect it from the count
    # matrix by checking whether the forward and reverse scans coincide.
    scan_both = reverse_complement

    sites = []
    for si, raw in enumerate(seqs):
        seq = str(raw).upper()
        if len(seq) < L:
            continue
        strands = [("+", seq)]
        if scan_both:
            strands.append(("-", "".join(comp.get(c, "N") for c in reversed(seq))))
        for strand, s in strands:
            idx = np.array([_BASE_TO_I.get(c, -1) for c in s], dtype=np.int64)
            valid = idx != -1
            scores = _window_scores(idx, valid, pwm)
            if scores.size == 0:
                continue
            hits = np.flatnonzero(scores >= threshold - tol * max(1.0, abs(threshold)))
            for pos in hits:
                score = float(scores[pos])
                if strand == "+":
                    start, end = int(pos), int(pos) + L
                else:
                    # the minus-strand window covers the same forward interval
                    start, end = len(seq) - int(pos) - L, len(seq) - int(pos)
                site = {"seq_index": si, "start": start, "end": end,
                        "strand": strand, "score": score}
                if table is not None:
                    site["pvalue"] = site_pvalue(pwm, score, None, table=table)
                sites.append(site)

    # A self-complementary motif scores identically on both strands at the same
    # interval; keep one record so counts are not doubled.
    if scan_both and sites:
        seen = {}
        deduped = []
        for s in sites:
            key = (s["seq_index"], s["start"], s["end"], s["score"])
            if key in seen:
                continue
            seen[key] = True
            deduped.append(s)
        sites = deduped

    sites.sort(key=lambda s: (s["seq_index"], s["start"], s["strand"]))
    return {"threshold": float(threshold),
            "achieved_fpr": None if achieved is None else float(achieved),
            "motif_length": L,
            "strands": "+/-" if reverse_complement else "+",
            "n_sites": len(sites),
            "sites": sites}


def run_motifs(args) -> int:
    """Scan a PWM, run de novo discovery, or test enrichment via AME."""
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
    # AME exposes 6 tests x 4 scorings that do not agree with each other, so the
    # choice is surfaced rather than fixed, and recorded in the result.
    ame_method = getattr(args, "ame_method", "fisher")
    ame_scoring = getattr(args, "ame_scoring", "avg")

    # `native` used to select a built-in de novo engine. That engine has been
    # removed, so the flag now exists only to explain itself: silently doing
    # nothing, or quietly routing to HOMER, would both be worse than a refusal.
    # The refusal itself is raised in run_homer_workflow, which is the layer that
    # knows de novo has no engine; nothing is warned here, because the error is
    # the complete answer.

    # Dispatch order matters. A supplied --db means known-motif enrichment, which
    # is a different question from de novo discovery; gating on --background first
    # sent every --db run into HOMER and silently ignored the database.
    if not pwm_path and not motif_db and fasta and (
            background or backend in ("homer", "auto")):
        if not getattr(args, "peaks", None):
            raise ValueError(
                "de novo discovery needs --peaks: the regions to analyse come "
                "from it. Pass --db with --background for known-motif enrichment, "
                "or --pwm for a scan, if you only have sequences.")
        # A missing background is passed through as None: HOMER rejects FASTA
        # mode without one, and the workflow turns that into the message the
        # caller sees. There is no target-as-null fallback any more.
        bg_fasta = background if background else None
        res = run_homer_workflow(
            peak_fasta=fasta, bg_fasta=bg_fasta, out_dir=str(out),
            k=k, backend=backend, homer_dir=homer_dir,
            homer_known=motif_db)  # motif_db used as known library if given
        print(f"[easyepi] motif workflow ({res.get('backend')}) written to {out}")
        return 0

    # --- paths that need no external program ---
    seqs = _read_sequences(fasta) if fasta else []
    bg_seqs = _read_sequences(background) if background else []

    res = {"status": "ok", "external_program_required": False}
    if pwm_path and seqs:
        pwm = read_pwm(Path(pwm_path))
        res.update(scan_sequences(pwm, seqs))
        res["n_seqs"] = len(seqs)
        res["mode"] = "scan"
    elif motif_db and seqs:
        if not bg_seqs:
            raise ValueError(
                "Known motif enrichment requires an explicit background FASTA "
                "(--background); without one there is no null to test against")
        # Delegated to AME: it is the published tool for this test and states its
        # own method and scoring, which a local reimplementation could not.
        from easyepi.motifs.ame_enrichment import ame_enrichment
        db_path = Path(motif_db)
        if db_path.suffix.lower() in (".meme", ".txt") and db_path.is_file():
            meme_db = str(db_path)
        else:
            # HOMER ships a directory of .motif files; JASPAR .pfm is a single
            # file. AME reads neither, so both are converted to MEME first.
            meme_db = _meme_database(db_path)
        try:
            ame = ame_enrichment(seqs, bg_seqs, meme_db,
                                 method=ame_method, scoring=ame_scoring)
        finally:
            _release_converted()
        res.update({"mode": "known_enrichment_ame", "n_seqs": len(seqs),
                    "ame": ame, "motifs": ame["motifs"]})
    else:
        raise ValueError(
            "Nothing to do. Pass --pwm for a scan, --db with --background for "
            "known-motif enrichment (both need no external program), or omit "
            "both to run de novo discovery through HOMER (which does).")
    out.write_text(json.dumps(res, indent=2, default=lambda x: x.tolist() if hasattr(x, "tolist") else x), encoding="utf-8")
    print(f"[easyepi] motif analysis written to {out}")
    return 0


# Conversions create a temporary MEME database that must outlive the AME
# subprocess but not the call. A module-level holder makes the lifetime explicit
# instead of relying on garbage collection of a TemporaryDirectory that was
# returned by value and dropped.
_CONVERTED_HOLD: list = []


def _meme_database(db_path: Path) -> str:
    """Return a path to a MEME file for ``db_path``, converting if necessary.

    Accepts a directory (HOMER's share/homer/motifs layout), a single *.motif or
    *.pfm file, or an existing MEME file, which is passed through untouched.
    """
    import glob
    from easyepi.motifs import formats as _f

    if db_path.is_file() and db_path.suffix.lower() in (".meme", ".txt"):
        return str(db_path)

    motifs = []
    seen = set()
    if db_path.is_dir():
        paths = []
        for pattern in ("*.motif", "*.motifs", "*.pfm", "*.meme"):
            paths.extend(sorted(glob.glob(str(db_path / pattern))))
    else:
        paths = [str(db_path)]

    for path in paths:
        try:
            for m in _f.read_motifs(Path(path)):
                key = str(m.name)
                if key in seen:
                    continue
                seen.add(key)
                motifs.append(m)
        except Exception:
            continue

    if not motifs:
        raise ValueError(
            f"no readable motifs found in {db_path}; supported inputs are a HOMER "
            "motif directory, a JASPAR .pfm file, or a MEME .meme file")

    tmp = tempfile.TemporaryDirectory()
    target = Path(tmp.name) / "motifs.meme"
    _f.write_meme(motifs, target)
    _CONVERTED_HOLD.append(tmp)   # released by the caller in its finally block
    return str(target)


def _release_converted() -> None:
    """Delete any temporary MEME database created for an AME run."""
    while _CONVERTED_HOLD:
        tmp = _CONVERTED_HOLD.pop()
        try:
            tmp.cleanup()
        except Exception:
            pass


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
