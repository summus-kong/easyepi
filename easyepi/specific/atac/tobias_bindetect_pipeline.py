"""1:1 reimplementation of TOBIAS BINDetect's orchestration layer.

The statistics and normalisation are already ported and verified. What remains is
the plumbing that produces their inputs, transcribed from ``scan_and_score`` and
``process_tfbs``:

**Motif scanning** (`utils/motifs.py`)
    PSSM built from the motif with pseudocount 0.01
    threshold = MOODS.tools.threshold_from_p(pssm, bg, motif_pvalue, 4)
    MOODS scanner with core 7, run on both strands

**Site scoring** (`bindetect_functions.py:scan_and_score`)
    the score attributed to a TFBS is the footprint value at its MIDPOINT:
        pos = TFBS.start - region.start + int(motif_length / 2)

**Background sampling** (same function)
    `random.seed(region_length)` then
    `random.sample(range(reglen), max(1, int(reglen / 200)))`
    so the sampled positions are a deterministic function of region length alone.

**Per-TF output** (`bindetect_functions.py:process_tfbs`)
    normalise -> clip negatives to 0 -> round to 5 -> bound if score > threshold
    log2fc = round(log2((c1 + pseudo) / (c2 + pseudo)), 5)
    observed per-TF statistic is averaged per peak before the differential test

This module deliberately keeps the MOODS dependency: the p-value thresholding is
part of BINDetect's definition, and easyepi's own best-hit scanner answers a
different question (where is the best match, not which sites pass a p-value).
"""
from __future__ import annotations

import random

import numpy as np

try:                                              # pragma: no cover - import guard
    import MOODS.scan
    import MOODS.tools
except Exception:                                 # noqa: BLE001
    MOODS = None

#: Pseudocount used when converting a motif to a PSSM (`get_pssm(ps=0.01)`).
PSSM_PSEUDOCOUNT = 0.01

#: MOODS core model size used by the reference scanner.
MOODS_CORE = 7

#: Default motif p-value for scanning.
DEFAULT_MOTIF_PVALUE = 0.0001

#: One background sample per this many bp of region length.
RANDOM_WINDOW = 200

#: Uniform background assumed by the reference.
UNIFORM_BG = (0.25, 0.25, 0.25, 0.25)


def pwm_to_pssm(pwm, pseudocount: float = PSSM_PSEUDOCOUNT, bg=UNIFORM_BG):
    """Log-odds matrix, laid out (4, L) with rows A, C, G, T.

    Mirrors `OneMotif.get_pssm` in TOBIAS, including two details that are easy to
    miss: the pseudocount is added as ``ps * bg`` and the result is *renormalised*
    per column, and the log is the natural log rather than log2. Using log2 inflates
    every entry by ln(2), which shifts MOODS thresholds as well.
    """
    pwm = np.asarray(pwm, dtype=float)
    if pwm.ndim != 2 or pwm.shape[0] != 4:
        raise ValueError("PWM must have 4 rows (A, C, G, T)")

    bg_col = np.asarray(bg, dtype=float).reshape((-1, 1))
    pfm = pwm / pwm.sum(axis=0, keepdims=True)
    pseudo_vector = pseudocount * bg_col
    pfm_pseudocount = (pfm + pseudo_vector) / (pfm + pseudo_vector).sum(
        axis=0, keepdims=True)

    pssm = np.log(pfm_pseudocount) - np.log(bg_col)
    return np.nan_to_num(pssm)


def motif_threshold(pssm, bg=UNIFORM_BG, pvalue: float = DEFAULT_MOTIF_PVALUE):
    """``MOODS.tools.threshold_from_p`` for one PSSM."""
    if MOODS is None:
        raise ImportError("MOODS is required to reproduce BINDetect's thresholds")
    pssm_tuple = tuple(tuple(float(v) for v in row) for row in np.asarray(pssm))
    return MOODS.tools.threshold_from_p(pssm_tuple, tuple(bg), pvalue, 4)


def scan_region(sequence: str, pssms, thresholds, bg=UNIFORM_BG):
    """All motif hits in one sequence, both strands.

    Returns ``(motif_index, start, end, strand, score)``.

    Each motif is scanned with its own scanner. MOODS returns
    ``tuple[strand] -> tuple[match, ...]`` where a ``match`` carries only ``pos`` and
    ``score`` -- no motif index -- so a combined scan cannot attribute hits back to
    motifs. Scanning per motif keeps that information explicit.
    """
    if MOODS is None:
        raise ImportError("MOODS is required to reproduce BINDetect's scanning")

    seq = sequence.upper()
    hits = []
    for motif_idx, (pssm, thr) in enumerate(zip(pssms, thresholds)):
        pssm_tuple = tuple(tuple(float(v) for v in row) for row in np.asarray(pssm))
        scanner = MOODS.scan.Scanner(MOODS_CORE)
        scanner.set_motifs([pssm_tuple], tuple(bg), [float(thr)])
        width = len(pssm_tuple[0])

        scanned = scanner.scan(seq)
        for strand_idx, matches in enumerate(scanned):
            strand = "+" if strand_idx == 0 else "-"
            for match in matches:
                hits.append((motif_idx, int(match.pos), int(match.pos) + width,
                             strand, float(match.score)))
    return hits


def midpoint_offset(tfbs_start: int, region_start: int, motif_length: int) -> int:
    """Index used to read the footprint score for a TFBS.

    Reference: ``pos = TFBS.start - region.start + int(motif_length / 2)``.
    """
    return int(tfbs_start - region_start + int(motif_length / 2))


def background_positions(region_length: int):
    """Deterministic background sample for one region.

    The reference seeds the RNG with the region length, so the positions depend on
    nothing else -- two runs over the same region give identical samples.
    """
    random.seed(region_length)
    n = max(1, int(region_length / RANDOM_WINDOW))
    return random.sample(range(region_length), n)


def score_tfbs(hits, motif_lengths, signals: dict):
    """Attach per-condition scores to each hit, read at the TFBS midpoint.

    ``hits`` carry *region-relative* coordinates (``start``/``end`` already had
    ``region.start`` subtracted), and ``signals`` maps condition name -> array of
    footprint values over the same region, so the midpoint indexes the array
    directly:

        pos = TFBS.start - region.start + int(motif_length / 2)

    with the subtraction already applied when the hits were built.
    """
    rows = []
    for motif_idx, start, end, strand, _score in hits:
        mid = midpoint_offset(start, 0, motif_lengths[motif_idx])
        row = {"motif": motif_idx, "start": start, "end": end, "strand": strand,
               "mid": mid}
        for cond, values in signals.items():
            idx = min(max(mid, 0), len(values) - 1)
            row[cond] = float(values[idx])
        rows.append(row)
    return rows


def finalize_tfbs(rows, norm_objects: dict, thresholds: dict, pseudo: float,
                  comparisons):
    """Normalise, clip, classify bound/unbound, and add log2 fold changes.

    Transcribed from ``process_tfbs``: the score is normalised then clipped at 0,
    rounded to 5 decimals, and compared strictly greater than the threshold.
    """
    for row in rows:
        for cond, norm in norm_objects.items():
            original = row[cond]
            value = float(norm.normalize(original))
            value = value if value > 0 else 0.0
            row[cond + "_score"] = float("{0:.5f}".format(value))
            row[cond + "_bound"] = 1 if row[cond + "_score"] > thresholds[cond] else 0

        for cond1, cond2 in comparisons:
            base = "{0}_{1}".format(cond1, cond2)
            fc = float(np.log2((row[cond1 + "_score"] + pseudo)
                               / (row[cond2 + "_score"] + pseudo)))
            row[base + "_log2fc"] = float("{0:.5f}".format(fc))
    return rows


def mean_log2fc_per_peak(rows, peak_ids, base: str):
    """Observed per-TF statistic: the log2fc averaged within each peak."""
    ids = np.asarray(peak_ids)
    vals = np.asarray([r[base + "_log2fc"] for r in rows], dtype=float)
    if ids.size == 0:
        return np.array([])
    order = np.argsort(ids, kind="stable")
    ids_sorted, vals_sorted = ids[order], vals[order]
    out, start = [], 0
    for i in range(1, len(ids_sorted) + 1):
        if i == len(ids_sorted) or ids_sorted[i] != ids_sorted[start]:
            out.append(float(np.mean(vals_sorted[start:i])))
            start = i
    return np.asarray(out)


def sample_signal(bw, chrom: str, start: int, end: int):
    """Footprint values over ``[start, end)`` as a numpy array.

    Mirrors ``region.get_signal``: missing values come back as ``NaN`` so that the
    downstream ``nanmean``-style statistics skip them rather than treating them as
    zero.
    """
    if bw is None:
        raise ValueError("a bigWig handle is required")
    if end <= start:
        return np.empty(0, dtype=float)
    values = bw.values(chrom, start, end, numpy=True)
    return np.asarray(values, dtype=float)


def background_signal(signals: dict, region_length: int):
    """Background values sampled at the reference's deterministic positions."""
    positions = background_positions(region_length)
    return {cond: np.asarray(values, dtype=float)[positions]
            for cond, values in signals.items()}


def tfbs_bed_row(row, cond_names, motif_name: str, peak_id: str = "."):
    """One row of the per-site BED, in the reference's column order.

    ``chrom start end name score strand`` from the peak, then the per-condition
    normalised scores, the bound flags, and the log2 fold changes -- matching the
    order `process_tfbs` appends them in.
    """
    fields = [row.get("chrom", "."), str(row.get("start", 0)), str(row.get("end", 0)),
              "{0}_{1}".format(motif_name, peak_id), ".", row.get("strand", "+")]
    for cond in cond_names:
        fields.append("{0:.3f}".format(float(row.get(cond + "_score", 0.0))))
    for cond in cond_names:
        fields.append(str(int(row.get(cond + "_bound", 0))))
    for key in sorted(k for k in row if k.endswith("_log2fc")):
        fields.append("{0:.3f}".format(float(row[key])))
    return "\t".join(fields)


def observed_statistics(rows, cond_names, comparison):
    """Per-TF summary numbers used by the permutation test.

    The reference reduces the per-site table to the mean score per condition and the
    mean log2 fold change; those three numbers drive `log2fc_pvalue`. Returns plain
    floats so the result can be written to the overview table as-is.
    """
    cond1, cond2 = comparison
    base = "{0}_{1}".format(cond1, cond2)
    stats = {}
    for cond in cond_names:
        vals = np.asarray([r.get(cond + "_score", np.nan) for r in rows], dtype=float)
        finite = vals[np.isfinite(vals)]
        stats[cond + "_mean_score"] = float(finite.mean()) if finite.size else float("nan")
    fcs = np.asarray([r.get(base + "_log2fc", np.nan) for r in rows], dtype=float)
    finite_fc = fcs[np.isfinite(fcs)]
    stats[base + "_mean_log2fc"] = (float(finite_fc.mean())
                                    if finite_fc.size else float("nan"))
    stats["n_sites"] = int(len(rows))
    return stats
