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

#: GC is bucketed at 1/200 resolution for the candidate index. A tolerance below
#: half a bucket still sees its neighbours, because the index is queried over the
#: whole tolerance range rather than at one point.
_GC_RESOLUTION = 0.005
_GC_BUCKETS = int(round(1.0 / _GC_RESOLUTION)) + 1


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

    ``window`` sets the sampling grid, *not* the maximum control width: a target
    wider than one window is matched to a control spanning several adjacent
    windows, so choose a window no larger than the narrowest target. Passing the
    smallest target width keeps the grid fine enough that length matching has
    candidates to choose from.
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


def _blocked_within(cand_starts: np.ndarray, cand_ends: np.ndarray,
                    blocked_starts: np.ndarray, blocked_ends: np.ndarray,
                    min_distance: int) -> np.ndarray:
    """True where a blocked interval is closer than ``min_distance``.

    Exactly equivalent to the outer product it replaces::

        gaps = np.maximum(np.maximum(bs[None, :] - ce[:, None],
                                     cs[:, None] - be[None, :]), 0)
        return (gaps < min_distance).any(axis=1)

    Both comparisons are strict, and the original clamps the gap at zero, so
    ``min_distance = 0`` excludes nothing at all -- including a control that
    overlaps its own target. That is the existing contract; it is preserved here
    rather than quietly improved, and the zero case is short-circuited so the
    intent is visible instead of resting on a `maximum` with zero.

    For positive ``min_distance`` the blocked intervals are sorted by start, so
    the ones satisfying the first branch form a prefix found by one search, and
    the second branch is answered by a running maximum of ends over that prefix.
    """
    cand_starts = np.asarray(cand_starts, dtype=np.int64)
    cand_ends = np.asarray(cand_ends, dtype=np.int64)
    if cand_starts.size == 0:
        return np.zeros(0, dtype=bool)
    if min_distance <= 0 or blocked_starts.size == 0:
        return np.zeros(cand_starts.size, dtype=bool)

    # Branch 1: bs - ce < md  <=>  bs < ce + md, strict on the left.
    upper = np.searchsorted(blocked_starts, cand_ends + min_distance, side="left")
    # Branch 2 within that prefix: cs - be < md  <=>  be > cs - md, strict.
    run_max_end = np.maximum.accumulate(blocked_ends)
    idx = np.maximum(upper - 1, 0)
    reach = np.where(upper > 0, run_max_end[idx], -1)
    return (upper > 0) & (reach > (cand_starts - min_distance))


def _interval_gc(prefix: np.ndarray, starts: np.ndarray, ends: np.ndarray,
                 gcs: np.ndarray, s0: int, e0: int) -> float | None:
    """GC fraction of ``[s0, e0)`` from a window-grid prefix sum.

    Exact for any interval, including partial windows at each end, so the GC that
    selection tests is the GC a caller would measure on the returned sequence.
    """
    if e0 <= s0:
        return None
    lo = int(np.searchsorted(starts, s0, side="right")) - 1
    if lo < 0:
        lo = 0
    hi = int(np.searchsorted(ends, e0, side="left"))
    if hi <= lo:
        hi = min(lo + 1, starts.size)
    total = 0.0
    weighted = 0.0
    for i in range(lo, min(hi + 1, starts.size)):
        a = max(int(starts[i]), s0)
        b = min(int(ends[i]), e0)
        if b > a:
            total += b - a
            weighted += float(gcs[i]) * (b - a)
    return (weighted / total) if total else None


def read_accessible_bed(path: str | Path
                        ) -> dict[str, np.ndarray]:
    """Read an accessible-region BED into per-chromosome sorted intervals.

    Accepts BED3+ and skips track/browser/comment lines. Overlapping intervals are
    merged, because a window inside two overlapping accessible regions is simply
    accessible, and leaving them separate would make the containment test depend
    on which interval happened to be found first.
    """
    raw: dict[str, list[tuple[int, int]]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end > start:
                raw.setdefault(fields[0], []).append((start, end))
    if not raw:
        raise ValueError(f"no accessible regions read from {path}")

    merged: dict[str, np.ndarray] = {}
    for chrom, pairs in raw.items():
        pairs.sort()
        out: list[tuple[int, int]] = []
        for start, end in pairs:
            if out and start <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], end))
            else:
                out.append((start, end))
        merged[chrom] = np.asarray(out, dtype=np.int64)
    return merged


def _contained(chrom: str, start: int, end: int,
               accessible: dict[str, np.ndarray]) -> bool:
    """True when ``[start, end)`` lies wholly inside one accessible interval."""
    intervals = accessible.get(chrom)
    if intervals is None or intervals.size == 0:
        return False
    idx = int(np.searchsorted(intervals[:, 0], start, side="right")) - 1
    if idx < 0:
        return False
    return int(intervals[idx, 1]) >= end


def accessible_fraction(chrom: str, regions: Sequence[tuple[str, int, int]],
                        accessible: dict[str, np.ndarray]) -> float:
    """Fraction of ``regions`` that lie wholly inside the accessible universe.

    Reported so a reader can see how restrictive the accessible set was: a
    universe that contains 15% of the genome and 95% of the targets is not an
    accessibility match, it is a restatement of the targets.
    """
    subset = [(s, e) for c, s, e in regions if c == chrom]
    if not subset:
        return 0.0
    hits = sum(1 for s, e in subset if _contained(chrom, int(s), int(e), accessible))
    return hits / len(subset)


def match_background(targets: Sequence[tuple[str, int, int]],
                     windows: dict[str, list[tuple[int, int, float]]],
                     tolerance: float = 0.05, min_distance: int = 2000,
                     same_chrom: bool = True, seed: int = 0,
                     exclude: Sequence[tuple[str, int, int]] | None = None,
                     accessible: dict[str, np.ndarray] | None = None,
                     max_tries_per_target: int = 200
                     ) -> tuple[list[tuple[str, int, int]], dict[str, Any]]:
    """Pick one GC- and width-matched control interval per target region.

    ``tolerance`` is the allowed absolute difference in GC fraction. Candidates
    must lie at least ``min_distance`` bp from every target (and every entry of
    ``exclude``) so a control cannot silently overlap the signal it controls for.
    Each candidate is used at most once.

    ``accessible`` restricts controls to a supplied accessible universe (see
    :func:`read_accessible_bed`). This is the difference between a
    composition-matched null and an accessibility-matched one: for ATAC data,
    sampling from the whole genome compares peaks against mostly-inaccessible
    sequence, which no plausible motif enrichment can be read against. When
    supplied, every control lies wholly inside an accessible interval.

    Returns ``(controls, summary)``. Unmatched targets are reported rather than
    silently dropped, because a control set that quietly loses half the hard cases
    biases the comparison. ``summary["background_design"]`` names the construction
    so a downstream report cannot confuse the two.
    """
    if not (0.0 < tolerance < 1.0):
        raise ValueError("tolerance must be in (0, 1)")
    if min_distance < 0:
        raise ValueError("min_distance must be non-negative")

    # A background drawn from the targets it controls for is circular, and a
    # p-value computed against it is meaningless. The check is on coverage rather
    # than on file identity, because an accessible set that merely contains every
    # peak is equally circular whatever it was named.
    if accessible is not None:
        covered = [accessible_fraction(c, targets, accessible) for c in
                   {t[0] for t in targets}]
        if covered and min(covered) >= 0.99:
            raise ValueError(
                "the accessible universe contains every target region, so a "
                "control drawn from it would come from the peaks being tested. "
                "Supply an accessibility call that is independent of these "
                "targets, or omit --accessible to use a composition-matched "
                "background and report it as such.")

    rng = np.random.default_rng(seed)
    by_chrom_cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    # GC-bucket index per chromosome: bucket b holds the window indices whose GC
    # lies in [b/res, (b+1)/res). Without it, every target scans every window.
    gc_buckets: dict[str, dict[int, np.ndarray]] = {}
    chrom_end_cache: dict[str, int] = {}
    prefix: dict[str, np.ndarray] = {}
    for chrom, rows in windows.items():
        arr = np.asarray(rows, dtype=np.float64)
        by_chrom_cache[chrom] = (arr[:, 0], arr[:, 1], arr[:, 2])
        # The last window's end is the only bound a control needs, and it does not
        # change between targets. Recomputing it per target cost 373M element
        # visits on a chr1-scale grid.
        chrom_end_cache[chrom] = int(arr[:, 1].max())
        # prefix[i] = sum of gc * width over windows 0..i-1, so the GC of any run
        # of whole windows is one subtraction.
        widths = arr[:, 1] - arr[:, 0]
        prefix[chrom] = np.concatenate(([0.0], np.cumsum(arr[:, 2] * widths)))
        buckets: dict[int, list[int]] = {}
        bucket_ids = np.clip(arr[:, 2] / _GC_RESOLUTION, 0,
                             _GC_BUCKETS - 1).astype(np.int64)
        for i, b in enumerate(bucket_ids):
            # When an accessible universe is given, only windows starting inside
            # it enter the index at all. Full containment of the control interval
            # is still verified per draw: a window can start in an accessible
            # region and the control can still run past its end.
            if accessible is not None and not _contained(
                    chrom, int(arr[i, 0]), int(arr[i, 0]) + 1, accessible):
                continue
            buckets.setdefault(int(b), []).append(i)
        gc_buckets[chrom] = {b: np.asarray(v, dtype=np.int64)
                             for b, v in buckets.items()}

    blocked: dict[str, list[tuple[int, int]]] = {}
    for chrom, s, e in list(targets) + list(exclude or []):
        blocked.setdefault(chrom, []).append((int(s), int(e)))

    # Sorted by start, so "is any blocked interval within min_distance of this
    # candidate" becomes an interval search instead of an outer product.
    blocked_idx: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, pairs in blocked.items():
        arr = np.asarray(sorted(pairs), dtype=np.int64)
        blocked_idx[chrom] = (arr[:, 0], arr[:, 1])

    # GC of every target is fixed, so resolve it once and reuse it in the retry
    # loop; recomputing per attempt turned a 200-try search into 200 lookups.
    # The GC a target must be matched *on* is the GC of its own interval, not the
    # GC of the single window containing its midpoint. For a target wider than a
    # window the two differ, and comparing a control's interval GC to a window's
    # GC left 4 of 197 controls outside the tolerance their report claims.
    target_gc_map: dict[tuple[str, int, int], float] = {}
    for chrom, s0, e0 in targets:
        key = (chrom, int(s0), int(e0))
        _st, _en, _gc = by_chrom_cache.get(chrom, (None, None, None))
        if _st is None:
            target_gc_map[key] = 0.0
            continue
        exact = _interval_gc(prefix[chrom], _st, _en, _gc, int(s0), int(e0))
        target_gc_map[key] = (_target_gc(chrom, (int(s0) + int(e0)) // 2, windows)
                              if exact is None else exact)

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
        # Only the GC buckets overlapping [want_gc - tolerance, want_gc + tolerance]
        # can contain a match, so gather those indices first. Everything below
        # then runs on a few hundred candidates instead of every window.
        lo_b = max(0, int((want_gc - tolerance) // _GC_RESOLUTION))
        hi_b = min(_GC_BUCKETS - 1, int((want_gc + tolerance) // _GC_RESOLUTION))
        per_chrom = gc_buckets.get(chrom, {})
        cand_idx = [per_chrom[b] for b in range(lo_b, hi_b + 1) if b in per_chrom]
        if not cand_idx:
            unmatched.append({"chrom": chrom, "start": s, "end": e,
                              "reason": "no window in the requested GC band"})
            continue
        sub = np.concatenate(cand_idx)
        # The control is exactly the target's width, starting at a window
        # boundary, and its GC is the GC of that exact interval. Snapping the width
        # to whole windows changed 119 of 197 control widths, and testing the
        # window-weighted average let 144 of 197 controls outside the tolerance the
        # report then claims -- both measured against the true interval GC.
        pref = prefix.get(chrom)
        n_windows = starts.size
        if pref is None or n_windows == 0:
            unmatched.append({"chrom": chrom, "start": s, "end": e,
                              "reason": "genome grid unusable for this chromosome"})
            continue
        blocked_chrom = blocked_idx.get(chrom)
        starts_b, ends_b = blocked_chrom if blocked_chrom is not None else (
            np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64))

        # Rejection sampling over the GC-band buckets: draw a window start, compute
        # the exact GC of [start, start + width), accept on GC, room and distance.
        # O(1) per attempt, instead of building a boolean over every window.
        pick = None
        for _ in range(max_tries_per_target):
            chosen = int(sub[int(rng.integers(0, sub.size))])
            first = int(np.searchsorted(starts, starts[chosen]))
            cs_c = int(starts[first])
            ce_c = cs_c + width
            if ce_c > chrom_end_cache.get(chrom, ce_c):
                continue
            if accessible is not None and not _contained(chrom, cs_c, ce_c,
                                                         accessible):
                continue
            span_gc = _interval_gc(pref, starts, ends, gcs, cs_c, ce_c)
            if span_gc is None or abs(span_gc - want_gc) > tolerance:
                continue
            cand = (chrom, cs_c, ce_c)
            if cand in used:
                continue
            if starts_b.size and _blocked_within(
                    np.array([cs_c]), np.array([ce_c]),
                    starts_b, ends_b, min_distance)[0]:
                continue
            pick = cand
            break
        if pick is None:
            unmatched.append({"chrom": chrom, "start": s, "end": e,
                              "reason": "no GC-matched window of sufficient width "
                                        "and distance was available"})
            continue
        # Report the control at the target's width so GC is measured over a
        # comparable number of bases. When the target is wider than one window,
        # extend across adjacent windows on the same chromosome, which is what
        # makes a 7.7 kb peak matchable at all.
        # The pick is already exactly the target's width.
        _pc, cs, ce = pick
        chrom_end = chrom_end_cache.get(chrom, ce)
        if ce > chrom_end:
            ce = chrom_end
            cs = max(0, ce - width)
        control = (chrom, int(cs), int(ce))
        if control in used or control == (chrom, int(s), int(e)):
            # the resized interval can collide with an earlier pick or the target
            unmatched.append({"chrom": chrom, "start": s, "end": e,
                              "reason": "resized control collided with an "
                                        "existing region"})
            continue
        used.add(control)
        used.add(pick)
        controls.append(control)

    summary = {
        "background_design": ("accessible, GC- and length-matched" if accessible
                              is not None else "GC- and length-matched"),
        "accessible_universe_supplied": accessible is not None,
        "n_targets": len(targets),
        "n_matched": len(controls),
        "n_unmatched": len(unmatched),
        "frac_matched": round(len(controls) / len(targets), 4) if targets else 0.0,
        "tolerance": tolerance,
        "min_distance": min_distance,
        "unmatched": unmatched[:50],
    }
    return controls, summary


_GC_START_CACHE: dict[int, tuple[np.ndarray, list]] = {}


def _target_gc(chrom: str, centre: int,
               windows: dict[str, list[tuple[str, int, int, float]]]) -> float:
    """GC fraction of the window containing ``centre`` (binary search).

    A linear scan over genome windows is O(n) per target and dominates on a real
    genome, so the starts are searched. The rebuilt start array is cached per
    window list: a 200 bp grid on chr1 is 1.25M rows, and rebuilding it once per
    target turned a fast lookup into 300 full array constructions.
    """
    rows = windows.get(chrom) or []
    if not rows:
        return 0.0
    key = id(rows)
    cached = _GC_START_CACHE.get(key)
    if cached is None or cached[1] is not rows:
        starts = np.fromiter((r[0] for r in rows), dtype=np.int64, count=len(rows))
        _GC_START_CACHE[key] = (starts, rows)
        cached = (starts, rows)
    starts = cached[0]
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
