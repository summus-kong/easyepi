"""Replicate concordance: do the samples that should agree actually agree?

Why this is its own module
--------------------------
A differential result is only as good as its replicates. Two samples assigned to the
same group must resemble each other more than they resemble the other group, or the
comparison is measuring batch, library quality or mislabelling rather than biology.
That check is cheap and it is not part of any single-sample metric, so it lives here.

Two measures, because they fail differently
-------------------------------------------
* **Signal correlation** (Spearman by default) across features. Sensitive to global
  quality differences: a low-input library correlates poorly with everything.
* **Region-set overlap (Jaccard)** on called regions. Sensitive to threshold and depth:
  two good replicates of unequal depth can still share few exact regions.

They are reported separately rather than folded into one score, because a sample that
fails one and passes the other points at a specific cause, and averaging that away
loses the diagnosis.

What is deliberately not claimed
--------------------------------
There is no universal pass threshold. What counts as concordant depends on the assay,
the depth and the question, so this module reports the numbers and a within-group
versus between-group summary, and leaves the judgement to the caller. A threshold
table for a specific project is a project decision.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

__all__ = ["pairwise_concordance", "replicate_concordance",
           "read_region_sets", "json_safe"]


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman correlation, ranking with ties averaged (no scipy dependency)."""
    if a.size < 2:
        return float("nan")
    ra = _rank(a)
    rb = _rank(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = float(np.sqrt((ra ** 2).sum() * (rb ** 2).sum()))
    if denom == 0:
        return float("nan")
    return float((ra * rb).sum() / denom)


def _rank(x: np.ndarray) -> np.ndarray:
    """Average ranks, so ties do not bias the correlation."""
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1, dtype=float)
    # average tied ranks
    sorted_x = x[order]
    i = 0
    while i < len(sorted_x):
        j = i
        while j + 1 < len(sorted_x) and sorted_x[j + 1] == sorted_x[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return ranks


def _jaccard(set_a: set, set_b: set) -> float:
    if not set_a and not set_b:
        return float("nan")
    union = len(set_a | set_b)
    return float(len(set_a & set_b) / union) if union else float("nan")


def _comparison(groups: Optional[Dict[str, str]], a: str, b: str) -> str:
    """Label a pair as within or between group, or ``unlabelled`` when there are none.

    Without groups every pair would otherwise be labelled ``between_group``, and the
    summary would then report "0 within-group pairs" as though replicates had been
    compared and failed -- rather than as though no grouping was supplied.
    """
    if not groups:
        return "unlabelled"
    ga, gb = groups.get(a), groups.get(b)
    if ga is None or gb is None:
        return "unlabelled"
    return "within_group" if ga == gb else "between_group"


def _read_bed_regions(path) -> set:
    """Read a BED/narrowPeak/broadPeak file as a set of ``(chrom, start, end)``.

    Only the first three columns are kept, so BED3, BED6, narrowPeak and
    broadPeak all compare on the same terms -- a Jaccard between a BED6 and a
    narrowPeak is otherwise defined by whichever extra columns each happens to
    have. Intervals with non-positive width are dropped: samtools and MACS3 both
    emit them when a peak collapses to a point, and they would count as regions
    that can never overlap anything.
    """
    import gzip

    text = str(path)
    opener = gzip.open if text.endswith(".gz") else open
    regions = set()
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start, end = int(parts[1]), int(parts[2])
                except ValueError:
                    continue
                if end > start:
                    regions.add((parts[0], start, end))
    except OSError as exc:
        raise ValueError(f"cannot read peak file {path}: {exc}") from exc
    return regions


def read_region_sets(specs: Sequence[str], *, require_exists: bool = True
                     ) -> Dict[str, set]:
    """Build ``{sample_id: region set}`` from ``SAMPLE=PATH`` or bare paths.

    ``SAMPLE=PATH`` names the sample explicitly. A bare path takes its sample id
    from the filename with the track extension removed, which is the same rule
    ``qc_inputs`` uses, so a folder of ``A1.narrowPeak`` files needs no mapping.

    The explicit form exists because a peak file's name is not always its sample
    name -- ``A1_peaks.narrowPeak`` would otherwise become sample ``A1_peaks``,
    and the design table would then fail to match with a confusing message about
    a sample nobody wrote down.
    """
    from easyepi.formats import PEAK_EXTS

    # Accept repeated flags and comma lists: --peaks A=a.bed,B=b.bed A2=c.bed
    #
    # argparse with ``action="append", nargs="+"`` hands over a *list of lists*
    # (one inner list per flag occurrence), while a direct API caller passes a
    # flat list of strings. Both arrive here, so flatten by type rather than by
    # depth: calling str() on an inner list produced the single "path"
    # "['/data/A1.narrowPeak'", which then failed a not-found check with a
    # message quoting a Python list at the user.
    flat = []

    def _collect(value):
        if isinstance(value, (list, tuple)):
            for item in value:
                _collect(item)
            return
        for piece in str(value).split(","):
            piece = piece.strip()
            if piece:
                flat.append(piece)

    for entry in specs or ():
        _collect(entry)
    if not flat:
        raise ValueError("no --peaks entries were given")

    out: Dict[str, set] = {}
    for entry in flat:
        if "=" in entry:
            sample, path = entry.split("=", 1)
            sample, path = sample.strip(), path.strip()
            if not sample or not path:
                raise ValueError(
                    f"--peaks entry {entry!r} must be SAMPLE=PATH with both parts nonempty")
        else:
            sample, path = "", entry
        if require_exists:
            from pathlib import Path as _Path

            if not _Path(path).exists():
                raise FileNotFoundError(f"--peaks file not found: {path}")
        if not sample:
            sample = _stem_for_peaks(path, PEAK_EXTS)
        if sample in out:
            raise ValueError(
                f"--peaks gives sample {sample!r} more than once; a sample can only "
                "have one peak set")
        out[sample] = _read_bed_regions(path)
    return out


def _stem_for_peaks(path, exts) -> str:
    """Sample id from a filename, stripping the longest known peak extension.

    ``a.bed.gz`` must reduce to ``a``, not ``a.bed``.
    """
    from pathlib import Path as _Path

    name = _Path(path).name
    for ext in sorted(exts, key=len, reverse=True):
        if name.endswith(ext):
            return name[: -len(ext)]
    return _Path(path).stem


def json_safe(value):
    """Replace NaN/inf with ``None`` so the output is strict JSON.

    ``json.dumps`` writes bare ``NaN``, which is not JSON: a strict parser (and
    ``jq``, and most other languages' loaders) rejects the whole document. The
    concordance code returns NaN deliberately -- an empty region set against an
    empty region set has no Jaccard, and fewer than three finite features give no
    correlation -- so the conversion belongs here, at the boundary, where "no
    value" becomes ``null`` rather than being silently rewritten to 0.
    """
    import math

    if isinstance(value, float):
        return None if not math.isfinite(value) else value
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def pairwise_concordance(
    sample_ids: Sequence[str],
    matrix: np.ndarray,
    *,
    groups: Optional[Dict[str, str]] = None,
    region_sets: Optional[Dict[str, set]] = None,
    method: str = "spearman",
) -> List[Dict[str, Any]]:
    """Every pairwise comparison, each tagged with whether it is within a group.

    ``matrix`` is features x samples, matching the count-matrix layout used elsewhere.
    ``region_sets`` maps a sample id to a set of called regions; when absent the
    Jaccard column is ``None`` rather than 0, because "not computed" and "no overlap"
    mean different things.
    """
    if method not in ("spearman", "pearson"):
        raise ValueError(f"method must be 'spearman' or 'pearson', got {method!r}")
    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2:
        raise ValueError("matrix must be two-dimensional (features x samples)")
    if values.shape[1] != len(sample_ids):
        raise ValueError(
            f"matrix has {values.shape[1]} columns but {len(sample_ids)} sample ids")

    out: List[Dict[str, Any]] = []
    for i, j in combinations(range(len(sample_ids)), 2):
        a, b = values[:, i], values[:, j]
        keep = np.isfinite(a) & np.isfinite(b)
        if keep.sum() < 3:
            rho = float("nan")
        elif method == "spearman":
            rho = _spearman(a[keep], b[keep])
        else:
            aa, bb = a[keep] - a[keep].mean(), b[keep] - b[keep].mean()
            denom = float(np.sqrt((aa ** 2).sum() * (bb ** 2).sum()))
            rho = float((aa * bb).sum() / denom) if denom else float("nan")

        entry: Dict[str, Any] = {
            "sample_a": sample_ids[i],
            "sample_b": sample_ids[j],
            "correlation": rho,
            "method": method,
            "comparison": _comparison(groups, sample_ids[i], sample_ids[j]),
        }
        if groups:
            entry["group_a"] = groups.get(sample_ids[i])
            entry["group_b"] = groups.get(sample_ids[j])
        if region_sets is not None:
            entry["jaccard"] = _jaccard(region_sets.get(sample_ids[i], set()),
                                        region_sets.get(sample_ids[j], set()))
        else:
            entry["jaccard"] = None
        out.append(entry)
    return out


def replicate_concordance(
    sample_ids: Sequence[str],
    matrix: np.ndarray,
    *,
    groups: Optional[Dict[str, str]] = None,
    region_sets: Optional[Dict[str, set]] = None,
    method: str = "spearman",
) -> Dict[str, Any]:
    """Pairwise concordance plus a within- versus between-group summary.

    The summary answers the question the numbers are for: are replicates closer to each
    other than to the other group? A negative answer is the finding -- it means the
    grouping is not reflected in the data, and a differential test on it would be
    reporting something else.
    """
    pairs = pairwise_concordance(sample_ids, matrix, groups=groups,
                                 region_sets=region_sets, method=method)
    within = [p for p in pairs if p["comparison"] == "within_group"]
    between = [p for p in pairs if p["comparison"] == "between_group"]
    unlabelled = [p for p in pairs if p["comparison"] == "unlabelled"]

    def _mean(values):
        vals = [v for v in values if v is not None and np.isfinite(v)]
        return float(np.mean(vals)) if vals else None

    within_rho = _mean(p["correlation"] for p in within)
    between_rho = _mean(p["correlation"] for p in between)
    within_j = _mean(p["jaccard"] for p in within if p["jaccard"] is not None)
    between_j = _mean(p["jaccard"] for p in between if p["jaccard"] is not None)

    summary: Dict[str, Any] = {
        "pairs": len(pairs),
        "within_group_pairs": len(within),
        "between_group_pairs": len(between),
        "within_group_mean_correlation": within_rho,
        "between_group_mean_correlation": between_rho,
        "grouping_supported": None,
    }
    if within_rho is not None and between_rho is not None:
        # Strictly greater, not "not worse": equal means the grouping carries no signal.
        summary["grouping_supported"] = bool(within_rho > between_rho)
    if within_j is not None or between_j is not None:
        summary["within_group_mean_jaccard"] = within_j
        summary["between_group_mean_jaccard"] = between_j
    if unlabelled and not groups:
        summary["note"] = (
            "no groups supplied, so every pair is an all-pairs comparison; pass a "
            "design with a group column to judge replicate concordance")
    elif not within:
        # One sample per group: concordance within a group cannot be assessed at all,
        # and saying so is more useful than returning an empty average.
        summary["note"] = (
            "no within-group pairs; replicate concordance needs at least two samples "
            "per group")
    return {"summary": summary, "pairs": pairs}
