"""Composition-matched background regions.

Why not just shift the target sites
-----------------------------------
A uniform shift keeps a locus on the same chromosome and roughly the same context,
but it does not preserve GC content, and GC is the strongest single predictor of
both accessibility signal and CpG-island proximity. On a mark like H3K4me3 the
shifted control lands in AT-rich gene deserts about as often as in GC-rich
promoters, so a target-vs-shifted contrast understates the real difference and its
p-value is meaningless.

What this module does
---------------------
1. **Bin the genome** into fixed-width windows (default 1000 bp).
2. **Measure GC** on each window from the genome FASTA.
3. **Match** every target region to a candidate window with a GC fraction inside a
   tolerance window, a minimum distance from any target (so a control cannot
   overlap the thing it controls for), and the same chromosome by default.
4. **Sample** without replacement, so one candidate cannot serve two targets.

This is the standard composition-matching construction used for enrichment
backgrounds (see e.g. the GC-matched shuffle in ChIPseeker/DeepTools workflows).
It is still a *sequence-composition* match: it does not match mappability, repeat
content, chromatin state or expression, and it cannot create GC-matched controls
for targets in regions absent from the binned genome.

References
----------
* Yu G, Wang LG, He QY (2015) ChIPseeker. Bioinformatics 31:2382 (GC-matched
  background framing and its stated limits).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterator, Sequence

import numpy as np


def gc_fraction(seq: str) -> float:
    """GC fraction of a sequence, ignoring anything that is not A/C/G/T."""
    if not seq:
        return 0.0
    upper = seq.upper()
    gc = upper.count("G") + upper.count("C")
    at = upper.count("A") + upper.count("T")
    total = gc + at
    return float(gc / total) if total else 0.0


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(name, sequence)`` from a FASTA, uppercased."""
    name: str | None = None
    buf: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf).upper()


def genome_windows(fasta_path: str | Path, window: int = 1000,
                   chromosomes: Sequence[str] | None = None,
                   min_gc_known: float = 0.5
                   ) -> dict[str, list[tuple[int, int, float]]]:
    """GC content of fixed-width windows across the genome.

    Windows where more than ``1 - min_gc_known`` of bases are N are dropped: their
    GC cannot be measured, so they cannot be matched on GC.
    """
    if window < 1:
        raise ValueError("window must be >= 1")
    keep = set(chromosomes) if chromosomes else None
    out: dict[str, list[tuple[int, int, float]]] = {}
    for name, seq in iter_fasta(fasta_path):
        if keep is not None and name not in keep:
            continue
        rows: list[tuple[int, int, float]] = []
        n = len(seq)
        for start in range(0, n, window):
            chunk = seq[start:start + window]
            if not chunk:
                continue
            acgt = sum(chunk.count(b) for b in "ACGT")
            if acgt < min_gc_known * len(chunk):
                continue
            rows.append((start, min(n, start + len(chunk)), gc_fraction(chunk)))
        if rows:
            out[name] = rows
    if not out:
        raise ValueError(f"no usable windows found in {fasta_path}")
    return out


def match_background(targets: Sequence[tuple[str, int, int]],
                     windows: dict[str, list[tuple[int, int, float]]],
                     tolerance: float = 0.05, min_distance: int = 2000,
                     same_chrom: bool = True, seed: int = 0,
                     exclude: Sequence[tuple[str, int, int]] | None = None,
                     max_tries_per_target: int = 200
                     ) -> tuple[list[tuple[str, int, int]], dict[str, Any]]:
    """Pick one GC- and width-matched control window per target region.

    ``tolerance`` is the allowed absolute difference in GC fraction. Candidates
    must lie at least ``min_distance`` bp from every target (and every entry of
    ``exclude``) so a control cannot silently overlap the signal it controls for.
    Each candidate is used at most once.

    Returns ``(controls, summary)``. Unmatched targets are reported rather than
    silently dropped, because a control set that quietly loses half the hard cases
    biases the comparison.
    """
    if not (0.0 < tolerance < 1.0):
        raise ValueError("tolerance must be in (0, 1)")
    if min_distance < 0:
        raise ValueError("min_distance must be non-negative")

    rng = np.random.default_rng(seed)
    target_gc: list[float] = []
    by_chrom_cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chrom, rows in windows.items():
        arr = np.asarray(rows, dtype=np.float64)
        by_chrom_cache[chrom] = (arr[:, 0], arr[:, 1], arr[:, 2])

    blocked: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e in list(targets) + list(exclude or []):
        blocked.setdefault(chrom, []).append((int(s), int(e)))
    blocked_arr = {c: np.asarray(sorted(v), dtype=np.int64) for c, v in blocked.items()}

    # GC of every target is fixed, so resolve it once and reuse it in the retry
    # loop; recomputing per attempt turned a 200-try search into 200 lookups.
    target_gc_map: dict[tuple[str, int, int], float] = {}
    for chrom, s0, e0 in targets:
        target_gc_map[(chrom, int(s0), int(e0))] = _target_gc(
            chrom, (int(s0) + int(e0)) // 2, windows)

    used: set[tuple[str, int, int]] = set()
    # ``exclude`` is a hard exclusion: those windows can never be controls, whatever
    # min_distance is. Treating them as merely "blocked by distance" made the
    # parameter silently do nothing when min_distance was 0.
    for chrom_x, s_x, e_x in (exclude or []):
        used.add((str(chrom_x), int(s_x), int(e_x)))
    controls: list[tuple[str, int, int]] = []
    unmatched: list[dict[str, Any]] = []

    for chrom, s, e in targets:
        # A control must never be the target itself, even when min_distance is 0
        # (the distance test is on interval starts and would pass for distance 0).
        self_key = (chrom, int(s), int(e))
        used.add(self_key)
        if chrom not in by_chrom_cache:
            unmatched.append({"chrom": chrom, "start": s, "end": e,
                              "reason": "chromosome not in genome windows"})
            continue
        starts, ends, gcs = by_chrom_cache[chrom]
        width = int(e - s)
        want_gc = target_gc_map[(chrom, int(s), int(e))]
        # candidate windows of comparable width: use the window whose span is
        # closest to the target width, then match GC within tolerance of it
        # Enumerate the eligible candidates instead of sampling blindly: with a
        # small candidate pool, random retries can miss the only valid window and
        # report a false "no match". Sampling then picks among real options.
        blocked_chrom = blocked_arr.get(chrom)
        cand_starts, cand_ends, cand_gcs = starts, ends, gcs
        width_slack = max(int(0.5 * width), 1)
        ok = np.abs((cand_ends - cand_starts) - width) <= width_slack
        ok &= np.abs(cand_gcs - want_gc) <= tolerance
        if blocked_chrom is not None and blocked_chrom.size:
            gaps = np.maximum(
                np.maximum(blocked_chrom[:, None, 0] - cand_ends[None, :],
                           cand_starts[None, :] - blocked_chrom[:, None, 1]), 0)
            ok &= ~(gaps < min_distance).any(axis=0)
        pick = None
        if ok.any():
            idxs = np.nonzero(ok)[0]
            free = [int(i) for i in idxs
                    if (chrom, int(cand_starts[i]), int(cand_ends[i])) not in used]
            if free:
                # one drawn index selects one whole candidate; drawing twice would
                # mix the start of one window with the end of another
                chosen = free[int(rng.integers(0, len(free)))]
                pick = (chrom, int(cand_starts[chosen]), int(cand_ends[chosen]))
        if pick is None:
            unmatched.append({"chrom": chrom, "start": s, "end": e,
                              "reason": "no GC- and distance-matched window found"})
            continue
        used.add(pick)
        controls.append(pick)

    summary = {
        "n_targets": len(targets),
        "n_matched": len(controls),
        "n_unmatched": len(unmatched),
        "frac_matched": round(len(controls) / len(targets), 4) if targets else 0.0,
        "tolerance": tolerance,
        "min_distance": min_distance,
        "unmatched": unmatched[:50],
    }
    return controls, summary


def _target_gc(chrom: str, centre: int,
               windows: dict[str, list[tuple[str, int, int, float]]]) -> float:
    """GC fraction of the window containing ``centre`` (binary search).

    A linear scan over genome windows is O(n) per target and becomes the dominant
    cost on a real genome; the starts are sorted, so a search is enough.
    """
    rows = windows.get(chrom) or []
    if not rows:
        return 0.0
    starts = np.fromiter((r[0] for r in rows), dtype=np.int64, count=len(rows))
    idx = int(np.searchsorted(starts, centre, side="right")) - 1
    if idx < 0:
        return 0.0
    _s, e, gc = rows[idx]
    return float(gc) if centre < e else 0.0


def write_bed(regions: Sequence[tuple[str, int, int]], path: str | Path,
              names: Sequence[str] | None = None) -> Path:
    """Write regions as BED3/BED4."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        for i, (chrom, s, e) in enumerate(regions):
            name = names[i] if names and i < len(names) else f"bg_{i+1}"
            fh.write(f"{chrom}\t{int(s)}\t{int(e)}\t{name}\n")
    return out
