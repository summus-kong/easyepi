"""Annotate peaks with transcription factors by scanning known PWMs in situ.

What this does, stated precisely
--------------------------------
For every peak, each supplied known motif is scanned along both strands of the
peak sequence. A motif is reported for a peak when the best site in that peak
scores above a threshold derived from the **exact** score distribution of the
PWM under a uniform background. The output is therefore "this motif has a
statistically compatible binding site inside this peak", which is the same
question HOMER's ``annotatePeaks.pl`` and FIMO answer. It is *not* evidence that
the named protein binds there in the assayed cells: motif identity is shared
between paralogues, chromatin context is ignored, and no footprint is computed.

Threshold method (no invented statistic)
----------------------------------------
A PWM log-odds score is a sum of independent per-column contributions. Under a
uniform base composition the score distribution is therefore the exact
convolution of the four per-column outcomes, which is computed here by dynamic
programming. The reported threshold is the smallest score whose upper-tail
probability is <= ``fpr``. This is the same construction FIMO uses (Grant et al.
2011, Bioinformatics 27:1017) and it needs no fitted null and no background
sample. ``fpr`` is a per-motif, per-sequence false-positive rate under that
uniform model; it is not a genome-wide FDR.

When a background FASTA is supplied, each TF additionally gets a Fisher exact
test of "fraction of peaks with a hit" vs "fraction of background sequences with
a hit" at the same threshold, with Benjamini-Hochberg correction across TFs.

References
----------
* Grant CE, Bailey TL, Noble WS (2011) FIMO. Bioinformatics 27:1017-1018.
* Heinz S et al. (2010) HOMER. Mol Cell 38:576-589.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from easyepi.motifs.formats import Motif, read_motifs

UNIFORM_BG = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)
_BASES = "ACGT"
_BASE_BYTES = [ord(c) for c in _BASES]


# ---------------------------------------------------------------------------
# Exact score distribution of one PWM
# ---------------------------------------------------------------------------
def pwm_score_distribution(pwm: np.ndarray,
                           bg: np.ndarray | None = None,
                           max_outcomes: int = 20_000) -> tuple[np.ndarray, np.ndarray]:
    """Score distribution of a log-odds PWM under a uniform sequence model.

    Each column contributes one of four scores with probability ``bg[b]``; the
    total is their convolution. The reachable state count grows up to 4x per
    column, and a real 19 bp motif has 4**19 ~ 2.7e11 distinct sums, so tracking
    them exactly is impossible. The score axis is therefore binned into at most
    ``max_outcomes`` equal-width bins spanning the motif's full score range, and
    the convolution is carried out over bin indices.

    Binning the score axis is exactly what FIMO does (Grant et al. 2011), and the
    returned threshold is conservative to within one bin width. Returns
    ``(scores, probs)`` ascending with probabilities summing to 1.
    """
    pwm = np.asarray(pwm, dtype=np.float64)
    if pwm.ndim != 2 or pwm.shape[0] != 4 or pwm.shape[1] == 0:
        raise ValueError("PWM must be a non-empty (4, length) matrix")
    if not np.isfinite(pwm).all():
        raise ValueError("PWM must be finite")
    bg = UNIFORM_BG if bg is None else np.asarray(bg, dtype=np.float64)
    if bg.shape != (4,) or bg.sum() <= 0:
        raise ValueError("background must be 4 base probabilities")
    bg = bg / bg.sum()
    if max_outcomes < 256:
        raise ValueError("max_outcomes must be >= 64")

    # Bin each column independently around that column's own minimum, then carry
    # the exact per-column offset alongside the bin index. Binning against a single
    # global range would shift every column by a different amount and bias the
    # whole distribution, which silently misplaces the threshold.
    span = float(np.ptp(pwm, axis=0).sum())
    if span <= 0:
        return np.array([float(pwm.sum())]), np.array([1.0])

    width = span / max_outcomes
    col_min = pwm.min(axis=0)                     # (L,)
    offsets = np.round(col_min / width).astype(np.int64)   # integer bin offsets
    bins = np.rint((pwm - col_min[None, :]) / width).astype(np.int64)
    bins = np.clip(bins, 0, max_outcomes)

    states = bins[:, 0]
    probs = bg.copy()
    for j in range(1, pwm.shape[1]):
        col = bins[:, j]
        new_states = (states[:, None] + col[None, :]).ravel()
        new_probs = (probs[:, None] * bg[None, :]).ravel()
        uniq, inverse = np.unique(new_states, return_inverse=True)
        states = uniq
        probs = np.bincount(inverse, weights=new_probs, minlength=uniq.size)

    total = probs.sum()
    if total <= 0:
        raise ValueError("degenerate PWM score distribution")
    # map bin index back to a score: add the accumulated per-column offsets
    scores = (states + int(offsets.sum())) * width
    return scores, probs / total


def pwm_threshold(pwm: np.ndarray, fpr: float = 1e-4,
                  bg: np.ndarray | None = None) -> tuple[float, float]:
    """Smallest score whose upper-tail probability is <= ``fpr``.

    Returns ``(threshold, achieved_fpr)``. Raises when ``fpr`` is not in (0, 1)
    or when the PWM is so short that no site can reach the requested rate.
    """
    if not (0.0 < fpr < 1.0):
        raise ValueError("fpr must be strictly between 0 and 1")
    scores, probs = pwm_score_distribution(pwm, bg)
    tail = np.cumsum(probs[::-1])[::-1]  # P(score >= s)
    ok = np.nonzero(tail <= fpr)[0]
    if ok.size == 0:
        raise ValueError(
            f"no PWM score reaches the requested fpr={fpr:g}; "
            f"the best achievable rate is {probs[-1]:g}")
    idx = int(ok[0])
    return float(scores[idx]), float(tail[idx])


def site_pvalue(pwm: np.ndarray, score: float,
                bg: np.ndarray | None = None,
                table: tuple[np.ndarray, np.ndarray] | None = None) -> float:
    """Upper-tail probability of observing at least ``score`` under the PWM null.

    Pass ``table`` (the ``(scores, tail)`` pair from :func:`tail_table`) when
    scoring many sites of the same motif: recomputing the convolution per hit is
    O(hits x L x 4) work, which dominates a genome-scale run.
    """
    if table is None:
        table = tail_table(pwm, bg)
    scores, tail = table
    idx = np.searchsorted(scores, score, side="left")
    if idx >= scores.size:
        return 0.0
    return float(tail[idx])


def tail_table(pwm: np.ndarray, bg: np.ndarray | None = None
               ) -> tuple[np.ndarray, np.ndarray]:
    """Precompute ``(scores, P(X >= score))`` once per motif."""
    scores, probs = pwm_score_distribution(pwm, bg)
    tail = np.cumsum(probs[::-1])[::-1]
    return scores, tail


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def read_peak_intervals(path: str | Path) -> list[dict[str, Any]]:
    """Read BED/narrowPeak/broadPeak as peak records with stable identifiers."""
    records: list[dict[str, Any]] = []
    opener = __import__("gzip").open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end <= start:
                continue
            name = fields[3] if len(fields) > 3 and fields[3] else f"peak_{lineno + 1}"
            records.append({"chrom": fields[0], "start": start, "end": end,
                            "name": name, "row": lineno + 1})
    if not records:
        raise ValueError(f"no intervals read from {path}")
    return records


def fetch_sequences(records: Sequence[dict[str, Any]], fasta_path: str | Path,
                    max_chroms: int | None = None) -> list[str]:
    """Pull one sequence per interval from an indexed FASTA.

    Requires ``pysam`` (the ``bam`` extra) and a ``.fai`` beside the FASTA.
    Missing contigs raise rather than silently returning empty sequences, because
    a silent empty sequence would be scored as "no motif" and hide a naming
    mismatch between the peak file and the genome build.
    """
    from easyepi.core.io import _require
    pysam = _require("pysam", "pip install 'easyepi[bam]'")
    fasta = pysam.FastaFile(str(fasta_path))
    try:
        available = set(fasta.references)
        missing = sorted({r["chrom"] for r in records} - available)
        if missing:
            raise ValueError(
                "peak chromosomes absent from the FASTA: " + ", ".join(missing[:10]) +
                (" ..." if len(missing) > 10 else "") +
                f" (FASTA has e.g. {sorted(available)[:5]})")
        seqs = []
        for rec in records:
            seq = fasta.fetch(rec["chrom"], rec["start"], rec["end"])
            seqs.append(seq.upper())
    finally:
        fasta.close()
    return seqs


def _score_all(seqs: Sequence[str], pwm: np.ndarray
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Best score and position of one PWM across every sequence, both strands.

    The scan is a correlation: with ``x[b, i]`` the one-hot indicator of base ``b``
    at position ``i``, the window starting at ``w`` scores
    ``sum_b sum_j pwm[b, j] * x[b, w + j]``. Correlating each base channel with its
    own PWM row is O(4 * n) in a C-level loop, instead of O(W * L) per strand with
    a materialised (W, L) temporary.

    The minus strand is the *reverse complement*: base ``c`` of the reversed
    sequence is the complement of the forward base, so channel ``c`` is
    ``onehot[inv_comp[c]]`` reversed in position. Reversing positions without
    complementing would score the reverse string against the forward PWM and
    silently report a wrong offset.

    Returns ``(score, forward_offset, strand)`` per sequence, with the offset
    always expressed on the forward sequence.
    """
    pwm = np.asarray(pwm, dtype=np.float64)
    L = pwm.shape[1]
    n_seq = len(seqs)
    best = np.full(n_seq, -np.inf)
    best_pos = np.full(n_seq, -1, dtype=np.int64)
    best_strand = np.full(n_seq, "+", dtype="<U1")
    # A->T, C->G, G->C, T->A as base indices (complement is an involution)
    inv_comp = [3, 2, 1, 0]

    for i, seq in enumerate(seqs):
        n = len(seq)
        if n < L:
            continue
        code = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
        onehot = np.zeros((4, n), dtype=np.float64)
        known = np.zeros(n, dtype=bool)
        for b in range(4):
            hit = code == _BASE_BYTES[b]
            onehot[b] = hit
            known |= hit

        def best_window(channels, known_mask):
            """Best offset/score for a set of base channels over this sequence."""
            bad_prefix = np.concatenate(([0], np.cumsum(~known_mask)))
            bad = (bad_prefix[L:] - bad_prefix[:-L]) > 0
            sc = np.zeros(n - L + 1)
            for c in range(4):
                sc += np.correlate(channels[c], pwm[c], mode="valid")
            sc = np.where(bad, -np.inf, sc)
            if not np.isfinite(sc).any():
                return None
            pos = int(np.argmax(sc))
            return pos, float(sc[pos])

        fwd = best_window(onehot, known)
        rev_channels = [onehot[inv_comp[c]][::-1] for c in range(4)]
        rev = best_window(rev_channels, known[::-1])

        for strand, res in (("+", fwd), ("-", rev)):
            if res is None:
                continue
            pos, score = res
            offset = pos if strand == "+" else n - pos - L
            if score > best[i]:
                best[i], best_pos[i], best_strand[i] = score, offset, strand
    return best, best_pos, best_strand


def _read_fasta_sequences(path: str | Path) -> list[str]:
    """Minimal FASTA reader for background sequence sets."""
    seqs: list[str] = []
    buf: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf).upper())
                    buf = []
            else:
                buf.append(line)
    if buf:
        seqs.append("".join(buf).upper())
    return seqs


def load_motif_database(spec: str | Path) -> list[Motif]:
    """Load one motive database; a comma-separated list is merged in order."""
    motifs: list[Motif] = []
    for item in str(spec).split(","):
        item = item.strip()
        if not item:
            continue
        motifs.extend(read_motifs(Path(item)))
    if not motifs:
        raise ValueError(f"no motifs read from {spec!r}")
    return motifs


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------
def annotate_sequences(seqs: Sequence[str], motifs: Sequence[Motif],
                       fpr: float = 1e-4, top: int = 5) -> dict[str, Any]:
    """Best above-threshold site of every motif in every sequence.

    Returns ``{"per_peak": [...], "summary": [...], "thresholds": {...}}``.
    ``per_peak`` is a long table (one row per peak x motif hit, best ``top``
    motifs per peak by score). Coordinates in ``start``/``end`` are offsets
    inside the supplied sequence, always on the forward strand.
    """

    seqs = [str(s).upper() for s in seqs]
    if not seqs:
        raise ValueError("no sequences to annotate")
    if top < 1:
        raise ValueError("top must be >= 1")

    ann_motifs = len(motifs)

    thresholds: dict[str, dict[str, float]] = {}
    hits_by_peak: list[list[dict[str, Any]]] = [[] for _ in seqs]
    skipped: dict[str, dict[str, Any]] = {}

    for motif in motifs:
        pwm = motif.to_log_odds()
        try:
            threshold, achieved = pwm_threshold(pwm, fpr)
        except ValueError as exc:
            # A short motif cannot reach a stringent per-sequence rate. Record it
            # and carry on: one unresolvable motif must not abort the whole run,
            # but it must never be reported as "no site" either.
            scores, probs = pwm_score_distribution(pwm)
            skipped[motif.name] = {"motif_id": motif.id, "length": motif.length,
                                   "best_achievable_fpr": float(probs[-1]),
                                   "reason": str(exc)}
            continue
        thresholds[motif.name] = {"motif_id": motif.id, "threshold": threshold,
                                  "fpr": fpr, "achieved_fpr": achieved,
                                  "length": motif.length,
                                  "consensus": motif.consensus_seq()}
        scores, positions, strands = _score_all(seqs, pwm)
        table = tail_table(pwm)
        for i, score in enumerate(scores):
            if not np.isfinite(score) or score < threshold:
                continue
            pos = int(positions[i])
            if pos < 0:
                continue
            hits_by_peak[i].append({
                "tf": motif.name,
                "motif_id": motif.id,
                "consensus": motif.consensus_seq(),
                "score": round(float(score), 4),
                "pvalue": site_pvalue(pwm, float(score), None, table=table),
                "strand": str(strands[i]),
                "offset_start": pos,
                "offset_end": pos + motif.length,
            })

    per_peak: list[dict[str, Any]] = []
    for i, hits in enumerate(hits_by_peak):
        hits.sort(key=lambda h: (-h["score"], h["tf"]))
        for rank, hit in enumerate(hits[:top], start=1):
            per_peak.append({"seq_index": i, "rank": rank, **hit})

    summary: list[dict[str, Any]] = []
    for motif in motifs:
        if motif.name in skipped:
            continue
        col = thresholds[motif.name]
        rows = [h for peak in hits_by_peak for h in peak if h["tf"] == motif.name]
        n_with = sum(1 for peak in hits_by_peak if any(h["tf"] == motif.name for h in peak))
        best = max((h["score"] for h in rows), default=float("nan"))
        summary.append({
            "tf": motif.name,
            "motif_id": motif.id,
            "consensus": col["consensus"],
            "motif_length": col["length"],
            "score_threshold": round(col["threshold"], 4),
            "n_peaks_with_site": n_with,
            "pct_peaks_with_site": round(100.0 * n_with / len(seqs), 2),
            "n_sites": len(rows),
            "best_score": None if not rows else round(best, 4),
        })
    summary.sort(key=lambda r: (-r["n_peaks_with_site"], -r["best_score"] if r["best_score"] is not None else 0.0))
    return {"per_peak": per_peak, "summary": summary, "thresholds": thresholds,
            "skipped_motifs": skipped,
            "n_sequences": len(seqs), "n_motifs": ann_motifs}


def enrichment_test(annotation: dict[str, Any], n_peaks: int,
                    bg_present: dict[str, int], n_bg: int) -> list[dict[str, Any]]:
    """Fisher exact enrichment of per-TF peak hit rate vs a background set.

    ``bg_present`` maps motif name to the number of background sequences that
    contain an above-threshold site. p-values are BH-corrected across the motifs
    actually tested.
    """
    from scipy.stats import fisher_exact
    from easyepi.core.stats import benjamini_hochberg

    present = {row["tf"]: row["n_peaks_with_site"] for row in annotation["summary"]}
    rows = []
    for row in annotation["summary"]:
        name = row["tf"]
        a = present[name]
        b = n_peaks - a
        c = bg_present.get(name, 0)
        d = n_bg - c
        if a + c == 0:
            continue
        odds, p = fisher_exact([[a, b], [c, d]])
        rows.append({"tf": name, "motif_id": row["motif_id"],
                     "n_peaks_with_site": a,
                     "pct_peaks_with_site": row["pct_peaks_with_site"],
                     "n_bg_with_site": c,
                     "pct_bg_with_site": round(100.0 * c / n_bg, 2) if n_bg else None,
                     "odds_ratio": None if not np.isfinite(odds) else round(float(odds), 4),
                     "pvalue": float(p)})
    if rows:
        padj = benjamini_hochberg([r["pvalue"] for r in rows])
        for row, q in zip(rows, padj):
            row["padj"] = None if not np.isfinite(q) else float(q)
        rows.sort(key=lambda r: (r["pvalue"], -(r["odds_ratio"] or 0)))
    return rows


# ---------------------------------------------------------------------------
# CLI driver
# ---------------------------------------------------------------------------
def run_tf_annotate(args) -> int:
    """CLI driver: peaks -> sequences -> per-motif exact thresholds -> tables."""
    from easyepi._plots import selected_formats
    from easyepi.footprint import fig_tf_annotation

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    records = read_peak_intervals(args.peaks)
    motifs = load_motif_database(args.motifs)
    seqs = fetch_sequences(records, args.genome_fasta)
    usable = [(rec, seq) for rec, seq in zip(records, seqs) if len(seq) > 0]
    if not usable:
        raise ValueError("every peak sequence is empty; check the genome build")
    records = [rec for rec, _ in usable]
    seqs = [seq for _, seq in usable]

    ann = annotate_sequences(seqs, motifs, fpr=float(args.fpr), top=int(args.top))

    per_peak_path = out / "tf_per_peak.tsv"
    with per_peak_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["peak_id", "chrom", "start", "end", "rank", "tf", "motif_id",
                         "consensus", "score", "site_pvalue", "strand",
                         "site_start", "site_end"])
        for row in ann["per_peak"]:
            rec = records[row["seq_index"]]
            writer.writerow([rec["name"], rec["chrom"], rec["start"], rec["end"],
                             row["rank"], row["tf"], row["motif_id"],
                             row["consensus"], row["score"], row["pvalue"],
                             row["strand"],
                             rec["start"] + row["offset_start"],
                             rec["start"] + row["offset_end"]])

    summary_path = out / "tf_summary.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["tf", "motif_id", "consensus", "motif_length",
                         "score_threshold", "n_peaks_with_site",
                         "pct_peaks_with_site", "n_sites", "best_score"])
        for row in ann["summary"]:
            writer.writerow([row["tf"], row["motif_id"], row["consensus"],
                             row["motif_length"], row["score_threshold"],
                             row["n_peaks_with_site"], row["pct_peaks_with_site"],
                             row["n_sites"],
                             "" if row["best_score"] is None else row["best_score"]])

    enrichment_rows: list[dict[str, Any]] = []
    result: dict[str, Any] = {
        "peaks_file": str(args.peaks),
        "genome_fasta": str(args.genome_fasta),
        "motif_sources": [s.strip() for s in str(args.motifs).split(",") if s.strip()],
        "fpr": float(args.fpr),
        "top": int(args.top),
        "n_peaks_input": len(read_peak_intervals(args.peaks)),
        "n_peaks_used": len(seqs),
        "n_motifs": ann["n_motifs"],
        "score_model": "exact PWM convolution under uniform background",
        "thresholds": ann["thresholds"],
        "skipped_motifs": ann["skipped_motifs"],
        "summary": ann["summary"],
    }
    if ann["skipped_motifs"]:
        print(f"[easyepi] {len(ann['skipped_motifs'])} motif(s) skipped: no site can "
              f"reach fpr={args.fpr:g} (see skipped_motifs in tf_annotate.json)")

    if args.background:
        bg_seqs = _read_fasta_sequences(args.background)
        bg_seqs = [s for s in bg_seqs if s]
        if not bg_seqs:
            raise ValueError(f"no sequences read from background FASTA {args.background}")
        bg_ann = annotate_sequences(bg_seqs, motifs, fpr=float(args.fpr), top=int(args.top))
        bg_present = {row["tf"]: row["n_peaks_with_site"] for row in bg_ann["summary"]}
        enrichment_rows = enrichment_test(ann, len(seqs), bg_present, len(bg_seqs))
        result["background"] = {"fasta": str(args.background), "n_sequences": len(bg_seqs)}
        result["enrichment"] = enrichment_rows
        enrich_path = out / "tf_enrichment.tsv"
        with enrich_path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["tf", "motif_id", "n_peaks_with_site", "pct_peaks_with_site",
                             "n_bg_with_site", "pct_bg_with_site", "odds_ratio",
                             "pvalue", "padj"])
            for row in enrichment_rows:
                writer.writerow([row["tf"], row["motif_id"], row["n_peaks_with_site"],
                                 row["pct_peaks_with_site"], row["n_bg_with_site"],
                                 row["pct_bg_with_site"],
                                 "" if row["odds_ratio"] is None else row["odds_ratio"],
                                 row["pvalue"],
                                 "" if row.get("padj") is None else row["padj"]])

    (out / "tf_annotate.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    formats = selected_formats(args)
    if formats:
        fig_tf_annotation(ann["summary"], n_peaks=len(seqs),
                          title=f"Known-motif annotation of {len(seqs)} peaks",
                          save=str(out / "tf_annotation"), formats=list(formats))

    n_hit_peaks = len({row["seq_index"] for row in ann["per_peak"]})
    print(f"[easyepi] tf-annotate: {ann['n_motifs']} motifs x {len(seqs)} peaks; "
          f"{n_hit_peaks} peaks carry at least one site at fpr<={args.fpr:g}")
    for row in ann["summary"][:5]:
        print(f"    {row['tf']:12s} {row['pct_peaks_with_site']:5.1f}% of peaks  "
              f"consensus={row['consensus']}")
    if enrichment_rows:
        top = enrichment_rows[0]
        print(f"[easyepi] most enriched: {top['tf']} "
              f"(padj={top.get('padj')}, odds={top['odds_ratio']})")
    print(f"[easyepi] wrote {per_peak_path.name}, {summary_path.name}, tf_annotate.json")
    return 0
