"""Time-course specificity and trend analysis.

Two questions are answered here, and they are deliberately kept separate:

  * how *specific* is a feature to one time point? (TAU)
  * does it trend monotonically over time? (Spearman)

TAU is a descriptive statistic, not a test. It says nothing about whether the
shape is reproducible, only how concentrated the signal is on a single time
point. Treating it as a p-value is the most common way this number gets
misread, so no p-value is attached to it anywhere in this module.

Definition (Yanai et al. 2005, the standard formulation)::

    x_hat_i = x_i / max_j(x_j)
    tau     = sum_i (1 - x_hat_i) / (n - 1)

with ``tau`` in [0, 1]. ``tau == 1`` means every time point but one is zero,
``tau == 0`` means a perfectly flat profile.

Three preconditions are enforced rather than assumed, because each one silently
produces a meaningless number if ignored:

  * expression must be non-negative -- an absolute value would manufacture
    specificity that is not there;
  * at least three time points -- with two, tau is a deterministic function of
    the ratio and cannot rank features;
  * an all-zero profile has no maximum, so its tau is undefined, not 0.0.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

MIN_TIME_POINTS = 3


class TimeCourseError(ValueError):
    """Raised when the design cannot support the requested analysis."""


def read_time_design(path, sample_ids):
    """Read a design TSV and return time per sample plus the ordered time levels.

    The design must carry ``sample_id`` and ``time``. ``time`` may be numeric
    (``0, 2, 6, 24``) or ordered labels (``t0, t2``); numeric values are kept as
    floats so they can be regressed on, labels are used only for grouping.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise TimeCourseError("design TSV is empty")
    if "sample_id" not in rows[0]:
        raise TimeCourseError("design TSV must contain a sample_id column")
    if "time" not in rows[0]:
        raise TimeCourseError(
            "design TSV must contain a time column for time-course analysis")

    if len(sample_ids) != len(set(sample_ids)):
        raise TimeCourseError("duplicate sample IDs in the count matrix")
    by_id = {r["sample_id"]: r for r in rows}
    if len(by_id) != len(rows):
        raise TimeCourseError("duplicate sample IDs in the design")

    missing = [s for s in sample_ids if s not in by_id]
    extra = [s for s in by_id if s not in sample_ids]
    if missing or extra:
        raise TimeCourseError(
            "design/count sample mismatch; missing=%s, extra=%s" % (missing, extra))

    raw = [str(by_id[s]["time"]).strip() for s in sample_ids]
    if any(not v for v in raw):
        raise TimeCourseError("every design row needs a non-empty time value")

    numeric = True
    values = []
    for v in raw:
        try:
            values.append(float(v))
        except ValueError:
            numeric = False
            break
    if numeric:
        time_values = np.array(values, dtype=float)
    else:
        time_values = None

    levels = list(dict.fromkeys(raw))
    if len(levels) < MIN_TIME_POINTS:
        raise TimeCourseError(
            "time-course analysis needs at least %d distinct time points, got %d (%s)"
            % (MIN_TIME_POINTS, len(levels), ", ".join(levels)))
    return raw, levels, time_values


def tau_specificity(per_time):
    """TAU for one profile given as a vector over time points.

    Returns ``None`` when the profile is all zero: there is no maximum to
    normalise by, and reporting 0.0 would read as "no specificity" when the
    truth is "not measurable".
    """
    x = np.asarray(per_time, dtype=float)
    if x.size < MIN_TIME_POINTS:
        raise TimeCourseError(
            "tau needs at least %d time points" % MIN_TIME_POINTS)
    if np.any(x < 0):
        raise TimeCourseError(
            "tau is undefined for negative expression; fix the input rather than "
            "taking absolute values")
    peak = x.max()
    if peak <= 0:
        return None
    x_hat = x / peak
    return float((1.0 - x_hat).sum() / (x.size - 1))


def stage_specificity(per_time):
    """One specificity score per stage, for a profile ordered by stage.

    Returns ``None`` for an all-zero profile: there is no scale to normalise
    against, and 0.0 would be misread as "measured, no specificity".

    ``max(result)`` equals :func:`tau_specificity` for the same profile. That
    is deliberate -- see the module docstring -- so the per-stage table and the
    global tau can never tell contradictory stories.
    """
    x = np.asarray(per_time, dtype=float)
    if x.size < MIN_TIME_POINTS:
        raise TimeCourseError(
            "stage specificity needs at least %d stages" % MIN_TIME_POINTS)
    if np.any(x < 0):
        raise TimeCourseError(
            "specificity is undefined for negative expression; fix the input "
            "rather than taking absolute values")
    if x.max() <= 0:
        return None
    n = x.size
    total_except = x.sum() - x
    # x_s == 0 means the stage has nothing to be specific about.
    with np.errstate(divide="ignore", invalid="ignore"):
        score = 1.0 - total_except / ((n - 1) * x)
    score = np.where(x > 0, score, 0.0)
    return [float(min(1.0, max(0.0, v))) for v in score]


def replicate_agreement(mat, assign, n_levels):
    """Per-feature spread of tau across replicates, and stage-call agreement.

    ``mat`` is features x samples, ``assign`` maps each sample to a stage index.
    Returns ``(tau_min, tau_max, top_stage_consistent, max_stage_cv)`` per
    feature. A stage with a single replicate contributes no spread, so its CV is
    ``None`` rather than 0 -- reporting 0 would claim agreement that was never
    measured.
    """
    n_feat = mat.shape[0]
    tau_min = np.full(n_feat, np.nan)
    tau_max = np.full(n_feat, np.nan)
    consistent = np.zeros(n_feat, dtype=bool)
    stage_cv = np.full(n_feat, np.nan)

    # Which samples belong to each stage, so a per-replicate profile can be
    # assembled by taking one column per stage.
    members = [np.flatnonzero(assign == s) for s in range(n_levels)]
    n_rep = min(len(m) for m in members) if members else 0
    if n_rep == 0:
        return tau_min, tau_max, consistent, stage_cv

    for feat in range(n_feat):
        taus = []
        tops = []
        for r in range(n_rep):
            profile = np.array([mat[feat, members[s][r]] for s in range(n_levels)],
                               dtype=float)
            tau = tau_specificity(profile)
            # An all-zero replicate has no specificity to spread, so it must not
            # enter min()/max(). Keeping the None would also have to mean "the
            # spread is undefined", which is what leaving these cells NaN says.
            if tau is None:
                continue
            taus.append(tau)
            tops.append(int(np.argmax(profile)))
        if taus:
            tau_min[feat] = min(taus)
            tau_max[feat] = max(taus)
        # Agreement is only meaningful when at least one replicate had a
        # measurable peak; with none, the stage call is unknown, not agreed.
        consistent[feat] = len(tops) > 0 and len(set(tops)) == 1

        # Spread of the winning stage, measured on the mean profile that tau
        # actually uses.
        mean_profile = np.array(
            [mat[feat, members[s]].mean() for s in range(n_levels)], dtype=float)
        winner = int(np.argmax(mean_profile))
        reps = mat[feat, members[winner]].astype(float)
        if reps.size > 1 and reps.mean() > 0:
            stage_cv[feat] = float(reps.std(ddof=1) / reps.mean())
    return tau_min, tau_max, consistent, stage_cv


def timecourse_summary(mat, region_ids, sample_ids, raw_time, levels, time_values):
    """Per-feature TAU, peak time and (optionally) a monotonic trend test.

    ``mat`` is features x samples. Samples are collapsed to per-time means
    before TAU so that uneven replicate counts do not weight a time point more
    heavily than its measurement warrants.
    """
    if mat.size == 0:
        return []
    if np.any(mat < 0):
        raise TimeCourseError(
            "counts contain negative values (%g minimum); tau requires "
            "non-negative expression" % float(mat.min()))

    index = {lvl: i for i, lvl in enumerate(levels)}
    # `raw_time` arrives in sample order (see read_time_design), so the level
    # index lines up positionally; looking it up by sample id would need a
    # dict and silently mangle the assignment.
    assign = np.array([index[t] for t in raw_time])
    per_time = np.empty((mat.shape[0], len(levels)), dtype=float)
    for lvl, col in index.items():
        per_time[:, col] = mat[:, assign == col].mean(axis=1)

    # Diagnostics that a mean profile cannot carry. Computed once for the whole
    # matrix; the per-feature values are read off inside the loop below.
    rep_tau_min, rep_tau_max, rep_consistent, rep_stage_cv = replicate_agreement(
        mat, assign, len(levels))

    rows = []
    for feat, rid in enumerate(region_ids):
        profile = per_time[feat]
        tau = tau_specificity(profile)
        stage_scores = stage_specificity(profile)
        # Replicate spread travels with the row so the ranked table can be read
        # without a second lookup.
        rep_min = rep_tau_min[feat]
        rep_max = rep_tau_max[feat]
        rep_cv = rep_stage_cv[feat]
        if tau is None:
            rows.append({
                "feature": rid,
                "tau": None,
                "stage_specificity": None,
                "top_stage": "",
                "top_stage_score": None,
                "peak_time": "",
                "second_time": "",
                "peak_second_ratio": None,
                "rho": None,
                "p": None,
            })
            continue

        order = np.argsort(profile)[::-1]
        top, second = order[0], order[1] if len(order) > 1 else order[0]
        top_value = float(profile[top])
        second_value = float(profile[second])
        ratio = (top_value / second_value) if second_value > 0 else None

        rho = p_value = None
        if time_values is not None:
            # `time_values` is per *sample*; `profile` is per *time point*.
            # Collapse the time axis the same way the profile was collapsed,
            # otherwise the two vectors have different lengths.
            axis = np.array(
                [float(np.mean(time_values[assign == col]))
                 for col in range(len(levels))])
            rho, p_value = _spearman(axis, profile)

        rows.append({
            "feature": rid,
            "tau": tau,
            # Spread of tau across replicates and agreement of the stage call.
            # Both are None when a stage has fewer than two replicates, which is
            # the honest answer: nothing was measured to spread.
            "tau_rep_min": rep_min,
            "tau_rep_max": rep_max,
            "top_stage_consistent": bool(rep_consistent[feat]),
            "max_stage_cv": rep_cv,
            "stage_specificity": stage_scores,
            "top_stage": levels[int(np.argmax(stage_scores))],
            "top_stage_score": max(stage_scores),
            "peak_time": levels[top],
            "second_time": levels[second],
            "peak_second_ratio": ratio,
            "rho": rho,
            "p": p_value,
        })
    return rows


def _spearman(x, y):
    """Spearman rho and its two-sided p-value, without pulling in scipy.

    Tie handling follows the usual mid-rank convention. Returns (None, None)
    when either vector is constant, because rho is undefined there rather than
    equal to zero.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < MIN_TIME_POINTS:
        raise TimeCourseError(
            "a trend test needs at least %d time points" % MIN_TIME_POINTS)
    rx, ry = _rank(x), _rank(y)
    if rx.std() == 0 or ry.std() == 0:
        return None, None
    rho = float(np.corrcoef(rx, ry)[0, 1])
    n = x.size
    if n < 3:
        return rho, None
    # t approximation, exact enough for n >= 3 and conservative at small n.
    if abs(rho) >= 1.0:
        return rho, 0.0
    t = rho * np.sqrt((n - 2) / (1.0 - rho * rho))
    p = _two_sided_t(t, n - 2)
    return rho, p


def _rank(v):
    """Mid-ranks, so tied time points do not bias the correlation."""
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(v.size, dtype=float)
    i = 0
    while i < v.size:
        j = i
        while j + 1 < v.size and v[order[j + 1]] == v[order[i]]:
            j += 1
        mid = (i + j) / 2.0 + 1.0
        ranks[order[i:j + 1]] = mid
        i = j + 1
    return ranks


def _two_sided_t(t, df):
    """Two-sided p for Student's t, via the regularised incomplete beta."""
    if df <= 0:
        return None
    x = df / (df + t * t)
    return float(_betainc(df / 2.0, 0.5, x))


def _betainc(a, b, x):
    """Regularised incomplete beta I_x(a, b) by continued fraction.

    Small and dependency-free; the time-course test only needs a handful of
    tenths of a percent of accuracy and the series converges quickly for the
    parameter range that arises here.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = _lgamma(a + b) - _lgamma(a) - _lgamma(b)
    front = np.exp(lbeta + a * np.log(x) + b * np.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def _betacf(a, b, x, max_iter=200, eps=3e-14):
    """Lentz's algorithm for the incomplete beta continued fraction."""
    tiny = 1e-30
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _lgamma(x):
    """log Gamma via the Lanczos approximation (no scipy dependency)."""
    coef = [76.18009172947146, -86.50532032941677, 24.01409824083091,
            -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
    y = x
    tmp = x + 5.5
    tmp -= (x + 0.5) * np.log(tmp)
    ser = 1.000000000190015
    for c in coef:
        y += 1.0
        ser += c / y
    return -tmp + np.log(2.5066282746310005 * ser / x)


def benjamini_hochberg(pvalues):
    """BH-adjusted p-values, preserving order and leaving NAs in place."""
    idx = [i for i, p in enumerate(pvalues) if p is not None]
    if not idx:
        return list(pvalues)
    ps = np.array([pvalues[i] for i in idx], dtype=float)
    n = ps.size
    order = np.argsort(ps)
    ranked = ps[order]
    adj = ranked * n / (np.arange(n) + 1.0)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = list(pvalues)
    for pos, original in enumerate(order):
        out[idx[original]] = float(adj[pos])
    return out


def write_timecourse(rows, out_path, levels, time_values):
    """Write the per-feature table, sorted by tau so specificity reads first."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ordered = sorted(
        rows,
        key=lambda r: (r["tau"] is None, -(r["tau"] or 0.0), r["feature"]))

    padj = benjamini_hochberg([r["p"] for r in ordered])

    columns = (["feature", "tau"] + ["spec_%s" % lvl for lvl in levels] +
               ["top_stage", "top_stage_score", "padj", "peak_time",
                "second_time", "peak_second_ratio", "rho", "p"])
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row, q in zip(ordered, padj):
            scores = row.get("stage_specificity")
            writer.writerow(
                [row["feature"], _fmt(row["tau"])] +
                [_fmt(None if scores is None else scores[i])
                 for i in range(len(levels))] +
                [row.get("top_stage", ""), _fmt(row.get("top_stage_score")),
                 _fmt(q),
                 row["peak_time"],
                 row["second_time"],
                 _fmt(row["peak_second_ratio"]),
                 _fmt(row["rho"]),
                 _fmt(row["p"])])
    _write_long_form(ordered, out_path, levels)
    return ordered, padj


def _write_long_form(ordered, out_path, levels):
    """One row per feature per stage, for people who prefer tidy data.

    A cell is blank when the feature has no measurable specificity at all, so
    "no signal" never masquerades as "score 0".
    """
    path = Path(str(out_path)[:-4] + ".by_stage.tsv") if str(out_path).endswith(
        ".tsv") else Path(str(out_path) + ".by_stage.tsv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature", "stage", "specificity", "is_top_stage"])
        for row in ordered:
            scores = row.get("stage_specificity")
            for i, lvl in enumerate(levels):
                writer.writerow([
                    row["feature"],
                    lvl,
                    _fmt(None if scores is None else scores[i]),
                    "" if scores is None else int(i == int(np.argmax(scores))),
                ])


def _fmt(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return "%.6g" % value
    return str(value)
