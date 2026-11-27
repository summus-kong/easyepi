"""Official HOMER delegation plus a self-contained known-motif enrichment.

De novo discovery is delegated to HOMER; there is no bundled substitute. The
known-motif enrichment below needs no external program: it scans PWMs with
numpy and tests hits with a binomial tail, and it is what `tf-annotate` uses.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from easyepi.core.stats import benjamini_hochberg
from easyepi.motifs import formats as F
from easyepi.motifs import scan_sequences

_BASES = "ACGT"
_BASE_TO_I = {b: i for i, b in enumerate(_BASES)}
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def _rc(seq: str) -> str:
    return "".join(_COMP.get(c, "N") for c in reversed(seq))


# ---------------------------------------------------------------------------
# Sequence I/O
# ---------------------------------------------------------------------------
def read_fasta(path: str) -> List[str]:
    """Read a FASTA into a list of sequence strings (uppercase)."""
    seqs: List[str] = []
    buf: List[str] = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if buf:
                seqs.append("".join(buf).upper())
            buf = []
        else:
            buf.append(line.upper())
    if buf:
        seqs.append("".join(buf).upper())
    return seqs


# ---------------------------------------------------------------------------
# Per-motif background model (the key to HOMER-matching enrichment)
# ---------------------------------------------------------------------------
def _per_motif_model(pwm: np.ndarray, bg_seqs: Sequence[str],
                     threshold_frac: float = 0.99) -> float:
    """Compute a motif-specific score threshold from the BACKGROUND sequences.

    Uses the background-score distribution's high percentile (e.g. 0.99) as the
    "hit" threshold, an exploratory empirical threshold, not HOMER calibration. This is much
    more sensitive than a global z-score for long / GC-rich PWMs.
    """
    bg_scores = np.asarray(scan_sequences(pwm, bg_seqs)["scores"], dtype=np.float64)
    bg_scores = bg_scores[np.isfinite(bg_scores)]
    if bg_scores.size == 0:
        return float("inf")
    return float(np.percentile(bg_scores, 100 * threshold_frac))


def _log_p_binomial(hits: int, n: int, p0: float) -> float:
    """-log10 binomial upper-tail P(X >= hits) given hit rate p0."""
    from scipy.stats import binom
    if n <= 0:
        return 0.0
    if p0 <= 0:
        return float(-np.log10(1.0 / max(n, 1))) if hits > 0 else 0.0
    if p0 >= 1.0:
        return float(-np.log10(max(binom.sf(hits - 1, n, p0), 1e-300)))
    return float(-np.log10(max(binom.sf(hits - 1, n, p0), 1e-300)))


def _gc_fraction(seq: str) -> float:
    """GC content fraction of a sequence (ignores N)."""
    gc = sum(1 for c in seq if c in "GC")
    n = sum(1 for c in seq if c in "ACGT")
    return gc / n if n else 0.5


def _gc_bin(seq: str, n_bins: int = 10) -> int:
    """Assign a sequence to a GC-content bin index (0..n_bins-1)."""
    return min(n_bins - 1, int(_gc_fraction(seq) * n_bins))


def known_motif_enrichment_homer(target_seqs: Sequence[str], bg_seqs: Sequence[str],
                                 motifs: Sequence[F.Motif],
                                 bg_model_threshold: float = 0.99,
                                 mscore: float = 3.0,
                                 gc_norm: bool = False,
                                 gc_bins: int = 10) -> List[Dict]:
    """Legacy exploratory enrichment; results are not equivalent to HOMER.

    For each known motif:
      1. compute its log-odds score in target and background (ZOOPS: best hit);
      2. a sequence "has" the motif if its best score exceeds HOMER's
         ``-mscore`` threshold (default 3.0); when ``bg_model_threshold`` < 1,
         a background-derived percentile is used as a floor so long/GC-rich
         motifs need a genuine enrichment (not just background popularity);
      3. binomial enrichment -> log-P;
      4. BH-corrected q-value.

    By default ``gc_norm`` is False: as in HOMER the background sequences are
    already GC/composition-matched to the peaks (via ``-bg``/``-norm`` GC-matched
    background), so a *separate* GC-stratification tends to over-correct and
    inflate low-complexity repeats. Set ``gc_norm=True`` only if the provided
    background is NOT composition-matched.
    """
    target = [s.upper() for s in target_seqs if s]
    bg = [s.upper() for s in bg_seqs if s]
    n_t = len(target) or 1
    n_b = len(bg) or 1

    t_bins = np.array([_gc_bin(s, gc_bins) for s in target], dtype=int) if gc_norm else None
    b_bins = np.array([_gc_bin(s, gc_bins) for s in bg], dtype=int) if gc_norm else None

    rows: List[Dict] = []
    for m in motifs:
        pwm = m.to_log_odds()
        # Background-derived floor threshold (per-motif null), if requested.
        floor = _per_motif_model(pwm, bg, bg_model_threshold) if bg_model_threshold < 1.0 else -np.inf
        thr = max(floor, mscore)
        t_scores = np.asarray(scan_sequences(pwm, target)["scores"], dtype=np.float64)
        b_scores = np.asarray(scan_sequences(pwm, bg)["scores"], dtype=np.float64)
        t_hit = t_scores >= thr
        b_hit = b_scores >= thr

        if gc_norm:
            t_like = b_like = t_denom = 0
            for b in range(gc_bins):
                tm = (t_bins == b); bm = (b_bins == b)
                if not tm.any() or not bm.any():
                    continue
                bin_b_hits = int(b_hit[bm].sum()); bin_b_n = int(bm.sum())
                bin_t_hits = int(t_hit[tm].sum()); bin_t_n = int(tm.sum())
                p0 = bin_b_hits / bin_b_n if bin_b_n else 0.0
                t_like += bin_t_hits; t_denom += bin_t_n
                b_like += int(round(p0 * bin_t_n))
            p0 = (b_like / t_denom) if t_denom else 0.0
            logp = _log_p_binomial(t_like, t_denom, p0)
            t_frac = t_like / t_denom if t_denom else 0.0
            b_frac = b_like / t_denom if t_denom else 0.0
        else:
            t_like = int(t_hit.sum()); b_like = int(b_hit.sum())
            p0 = (b_like / n_b) if n_b else 0.0
            logp = _log_p_binomial(t_like, n_t, p0)
            t_frac = t_like / n_t; b_frac = b_like / n_b

        if t_frac <= 0:
            lo = -10.0
        elif b_frac <= 0:
            lo = 10.0
        else:
            lo = float(np.log2(t_frac / b_frac))
        rows.append({
            "name": m.name, "id": m.id, "tf": m.tf, "consensus": m.consensus_seq(),
            "n_target_with_motif": int(t_like),
            "pct_target": round(100 * t_frac, 2), "pct_background": round(100 * b_frac, 2),
            "enrichment_fold": round(2 ** lo, 2) if lo != 10.0 else float("inf"),
            "pvalue": float(10 ** -logp) if logp < 300 else 0.0,
            "log_p": round(logp, 3), "log_odds": round(lo, 3),
            "score_threshold": round(float(thr), 3),
        })
    pvals = np.array([max(r["pvalue"], 1e-300) for r in rows], dtype=np.float64)
    qvals = benjamini_hochberg(pvals)
    for r, q in zip(rows, qvals):
        r["qvalue"] = round(float(q), 4)
    rows.sort(key=lambda r: (-r["log_p"], -r["log_odds"]))
    return rows


# ---------------------------------------------------------------------------
# Multiple PFM database merge (HOMER + optional JASPAR/HOCOMOCO/Cis-BP)
# ---------------------------------------------------------------------------
def load_known_motifs(homer_known: Optional[str] = None,
                      extra_pfm_dirs: Optional[Sequence[str]] = None,
                      extra_pfm_files: Optional[Sequence[str]] = None,
                      dedup_threshold: float = 0.9) -> List[F.Motif]:
    """Load HOMER's known motifs plus any optional external TF PFM libraries.

    Args:
        homer_known: path to HOMER ``known.motifs`` (or ``all.motifs``).
        extra_pfm_dirs: directories scanned recursively for ``.pfm``/``.jaspar``/
            ``.meme``/``.motif`` files (JASPAR / HOCOMOCO / Cis-BP exports).
        extra_pfm_files: explicit PFM files.
        dedup_threshold: reverse-complement-aware similarity to drop near-dup.
    """
    motifs: List[F.Motif] = []
    if homer_known and Path(homer_known).exists():
        motifs.extend(F.read_motifs(homer_known))
    sources = []
    for d in (extra_pfm_dirs or []):
        for ext in (".pfm", ".jaspar", ".meme", ".motif", ".motifs"):
            sources.extend(sorted(Path(d).rglob(f"*{ext}")))
    for f in (extra_pfm_files or []):
        if Path(f).exists():
            sources.append(Path(f))
    for src in sources:
        try:
            motifs.extend(F.read_motifs(src))
        except Exception as exc:  # skip unreadable libs
            print(f"[easyepi] skip PFM {src}: {exc}")
    return F.dedup_merge(motifs, threshold=dedup_threshold)


# ---------------------------------------------------------------------------
# Writer for HOMER-format results
# ---------------------------------------------------------------------------
def write_known_results(rows: List[Dict], path: str) -> None:
    """Write ``knownResults.txt``-style tab-delimited table."""
    header = ["Motif Name", "Consensus", "P-value", "Log P-value",
              "q-value (Benjamini)", "# Target with Motif", "% Target",
              "# Background with Motif", "% Background"]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join([
                r["name"], r["consensus"], "%.1e" % r["pvalue"],
                "%.3f" % (-r["log_p"]) if r["log_p"] else "0", "%.4f" % r.get("qvalue", 0),
                "-", "%.2f%%" % r["pct_target"], "-", "%.2f%%" % r["pct_background"],
            ]) + "\n")


def write_homer_results(rows: List[Dict], path: str) -> None:
    """Write ``homerResults.txt``-style table for de novo motifs."""
    header = ["Motif Name", "Consensus", "P-value", "Log P-value", "q-value",
              "% Target", "% Background", "IC (bits)"]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for i, r in enumerate(rows):
            fh.write("\t".join([
                f"motif_{i+1}", r["consensus"], "%.1e" % r["pvalue"],
                "%.3f" % (-r["log_p"]) if r["log_p"] else "0", "-",
                "%.2f%%" % r["pct_target"], "%.2f%%" % r["pct_background"],
                "%.2f" % r["ic_bits"],
            ]) + "\n")


# ---------------------------------------------------------------------------
# Top-level workflow
# ---------------------------------------------------------------------------
def run_homer_workflow(peak_fasta: str, bg_fasta: Optional[str], out_dir: str,
                       known_motif_files: Optional[Sequence[str]] = None,
                       extra_pfm_dirs: Optional[Sequence[str]] = None,
                       k: int = 8, n_denovo: int = 10, n_known: int = 25,
                       bg_model_threshold: float = 0.99,
                       mscore: float = 3.0,
                       gc_norm: bool = False, gc_bins: int = 10,
                       homer_known: Optional[str] = None,
                       backend: str = "auto",
                       homer_dir: Optional[str] = None) -> Dict:
    """De novo discovery through HOMER, plus an optional built-in known enrichment.

    De novo is HOMER's job and only HOMER's: the package ships no substitute, so
    there is no mode in which a caller receives heuristic numbers believing them
    to be HOMER's. ``backend`` is kept for API compatibility:

        ``auto`` / ``homer``  run HOMER; absent or failing raises.
        ``native``            not supported for de novo any more; the known
                              enrichment still works and is reached with
                              ``--db``/``--pwm`` through ``easyepi motifs``.

    This makes HOMER a requirement for this function, not for the package: every
    other command, and `easyepi motifs --pwm`/`--db`, runs without it.
    """
    from easyepi.motifs import external_tools as ET

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    if backend == "native":
        raise ValueError(
            "de novo motif discovery has no built-in engine any more: the "
            "homer_like implementation was removed because it did not reproduce "
            "HOMER (3.90 bp mean shared consensus, and its EM variant claimed "
            "every sequence as a site). Install HOMER and use --backend homer, "
            "or scan known motifs with --pwm/--db, which needs no external "
            "program.")
    if backend not in {"auto", "homer"}:
        raise ValueError(f"Unknown motif backend: {backend}")

    if not bg_fasta:
        raise ValueError("HOMER FASTA mode requires an explicit background FASTA")

    homer_result = ET.run_homer_motifs(
        peak_fasta, bg_fasta, str(out), k=k, homer_dir=homer_dir)
    if homer_result.get("returncode") != 0:
        raise RuntimeError(
            "HOMER findMotifs.pl failed (exit "
            f"{homer_result.get('returncode')}); stderr: "
            f"{homer_result.get('stderr', '')[-800:]}")

    homer_result["backend"] = "homer"
    homer_result["delegated"] = True
    homer_result["validated_equivalent_to"] = "HOMER"
    homer_result["experimental"] = False
    homer_result["out_dir"] = str(out)
    return homer_result
