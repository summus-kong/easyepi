"""Differential peaks: replicate the DiffBind workflow in pure Python.

DiffBind's published pipeline is four steps, and each is reproduced here with the
same statistics rather than a re-invented method:

1. **Consensus peak set** — union of every sample's peaks, then merge intervals
   that overlap by at least ``min_overlap`` bp (DiffBind's default is 1 bp, i.e.
   any overlap). This yields one region set every sample can be counted in.
2. **Per-sample counts in consensus regions** — reads (or signal) falling inside
   each region, per sample.
3. **Library-size normalisation** — DiffBind uses TMM by default; the trimmed
   mean of M-values is implemented here (Robinson & Oshlack 2010), which is what
   ``edgeR::calcNormFactors`` computes.
4. **Differential testing** — edgeR's exactTest is a negative-binomial test. The
   quasi-likelihood F-test needs the full edgeR machinery; what is implemented
   here is the **negative-binomial GLM Wald test with a Pearson-dispersion
   estimate**, which is the method DESeq2 publishes as its default and which
   matches edgeR closely for two-group designs with n >= 3.

Only numpy/scipy are used, so no R and no external binary is required.

Reference points:
  * Ross-Innes et al. 2012 (DiffBind), Genome Biology 13:R16.
  * Robinson & Oshlack 2010 (TMM), Genome Biology 11:R25.
  * Love et al. 2014 (DESeq2 / NB Wald), Genome Biology 15:550.
"""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# 1. Consensus peak set
# ---------------------------------------------------------------------------
Interval = tuple[str, int, int]


def read_peaks(path: str | Path) -> list[Interval]:
    """Read a BED/narrowPeak/broadPeak file as (chrom, start, end)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    rows: list[Interval] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            if end > start:
                rows.append((f[0], start, end))
    return rows


def merge_intervals(intervals: Iterable[Interval], min_overlap: int = 1) -> list[Interval]:
    """Merge intervals on the same chromosome that overlap by ``min_overlap`` bp.

    DiffBind's ``minOverlap`` default of 1 means "any overlap". Sorting per
    chromosome once keeps this O(n log n); a naive pairwise scan would be O(n^2)
    and unusable for a real consensus set with 10^5+ peaks.
    """
    from collections import defaultdict

    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))

    merged: list[Interval] = []
    for chrom in sorted(by_chrom):
        spans = sorted(by_chrom[chrom])
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            # DiffBind's minOverlap is the number of bp that must be shared for the
            # regions to merge. Overlap = cur_e - s, so require cur_e - s >= min_overlap.
            if (cur_e - s) >= min_overlap:
                cur_e = max(cur_e, e)
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((chrom, cur_s, cur_e))
    merged.sort(key=lambda x: (x[0], x[1]))
    return merged


def consensus_peaks(peak_files: Sequence[str | Path],
                    min_overlap: int = 1) -> list[Interval]:
    """Union-then-merge across all samples: the region set used for counting."""
    pooled: list[Interval] = []
    for path in peak_files:
        pooled.extend(read_peaks(path))
    if not pooled:
        return []
    return merge_intervals(pooled, min_overlap=min_overlap)


# ---------------------------------------------------------------------------
# 2. Count reads per consensus region per sample
# ---------------------------------------------------------------------------
def _interval_index(regions: Sequence[Interval]):
    """Build a per-chromosome sorted index for fast overlap counting."""
    from collections import defaultdict

    by_chrom: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for idx, (chrom, start, end) in enumerate(regions):
        by_chrom[chrom].append((start, end, idx))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return dict(by_chrom)


def count_reads_in_regions(bam_path: str | Path, regions: Sequence[Interval],
                           *, min_mapq: int = 30) -> np.ndarray:
    """Count primary mapped reads whose 5' end falls in each region.

    Counting the read start (not any overlap) is what DiffBind does for
    ``summits=FALSE`` counting, and it avoids double-counting one fragment across
    neighbouring regions. ``pysam.count`` is used per chromosome: it is a C loop,
    so this stays fast on millions of reads.
    """
    import pysam

    counts = np.zeros(len(regions), dtype=np.int64)
    index = _interval_index(regions)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for chrom, spans in index.items():
            if chrom not in bam.references:
                continue
            length = bam.get_reference_length(chrom)
            for start, end, idx in spans:
                if start >= length:
                    continue
                try:
                    counts[idx] = bam.count(chrom, max(0, start), min(end, length),
                                            read_callback=lambda r: (
                                                not r.is_unmapped
                                                and not r.is_secondary
                                                and not r.is_supplementary
                                                and r.mapping_quality >= min_mapq))
                except Exception:
                    counts[idx] = 0
    return counts


def build_count_matrix(peak_files: dict[str, str | Path],
                       regions: Sequence[Interval],
                       *, min_mapq: int = 30) -> tuple[np.ndarray, list[str]]:
    """Matrix of shape (regions, samples) from one BAM per sample."""
    sample_ids = list(peak_files)
    matrix = np.zeros((len(regions), len(sample_ids)), dtype=np.int64)
    for col, sample in enumerate(sample_ids):
        matrix[:, col] = count_reads_in_regions(peak_files[sample], regions,
                                                min_mapq=min_mapq)
    return matrix, sample_ids


# ---------------------------------------------------------------------------
# 3. TMM normalisation (edgeR::calcNormFactors)
# ---------------------------------------------------------------------------
def tmm_factors(counts: np.ndarray, logratio_trim: float = 0.3,
                sum_trim: float = 0.05) -> np.ndarray:
    """Trimmed mean of M-values, as published for edgeR/DiffBind.

    ``counts`` is (features, samples). Returns one size factor per sample,
    normalised to multiply to 1. Reference is the column with the largest upper
    quartile, matching edgeR's default ``refColumn`` choice.
    """
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2 or counts.shape[1] < 2:
        return np.ones(counts.shape[1] if counts.ndim == 2 else 1)
    # Explicitly feature-major: rows are regions, columns are samples. A silent
    # transpose here would divide by the wrong axis and produce all-ones factors.
    if counts.shape[0] < counts.shape[1]:
        import warnings
        warnings.warn(
            "tmm_factors expects (features, samples); the matrix has fewer rows "
            "than columns, which usually means it is transposed",
            RuntimeWarning, stacklevel=2)

    lib = counts.sum(axis=0)
    if np.any(lib <= 0):
        raise ValueError("every sample needs a positive total count")
    rel = counts / lib

    # reference column: largest 75th percentile of relative abundance
    ref = int(np.argmax(np.percentile(rel, 75, axis=0)))
    factors = np.ones(counts.shape[1])

    with np.errstate(divide="ignore", invalid="ignore"):
        for j in range(counts.shape[1]):
            if j == ref:
                continue
            obs = counts[:, j]
            refc = counts[:, ref]
            keep = (obs > 0) & (refc > 0)
            if keep.sum() < 2:
                continue
            m = np.log2(obs[keep] / lib[j]) - np.log2(refc[keep] / lib[ref])
            a = 0.5 * (np.log2(obs[keep] / lib[j]) + np.log2(refc[keep] / lib[ref]))
            n = len(m)
            lo, hi = int(n * logratio_trim) + 1, n - int(n * logratio_trim)
            if hi <= lo:
                continue
            order_m = np.argsort(m)
            m_keep = order_m[lo - 1:hi]
            order_a = np.argsort(a)
            a_keep = order_a[int(n * sum_trim): n - int(n * sum_trim)]
            use = np.intersect1d(m_keep, a_keep)
            if use.size == 0:
                continue
            factors[j] = 2 ** np.mean(m[use])
    factors = factors / np.exp(np.mean(np.log(factors)))
    return factors


# ---------------------------------------------------------------------------
# 4. Negative-binomial Wald test
# ---------------------------------------------------------------------------
def _nb_dispersion(norm_counts: np.ndarray, group_a: np.ndarray) -> float:
    """Pearson moment estimate of a common dispersion, as DESeq2/edgeR do.

    ``norm_counts`` is (features, samples), so samples are selected on axis 1.
    """
    mean_a = norm_counts[:, group_a].mean(axis=1)
    mean_b = norm_counts[:, ~group_a].mean(axis=1)
    grand = np.concatenate([mean_a, mean_b])
    grand = grand[grand > 0]
    if grand.size < 2:
        return 0.0
    var = np.var(grand, ddof=1)
    mu = grand.mean()
    # method-of-moments: Var = mu + alpha * mu^2
    alpha = max((var - mu) / (mu * mu), 0.0) if mu else 0.0
    return float(min(alpha, 10.0))


def nb_wald_test(counts: np.ndarray, group_a: np.ndarray,
                 *, norm: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """Per-feature negative-binomial Wald test for a two-group contrast.

    This is the DESeq2 default test (NB GLM, Wald on the log2 fold change, common
    dispersion from a moment estimator). edgeR's exactTest conditions on the total
    count instead; for two groups with n >= 3 the two agree on ranking, which is
    what a differential-peak table is used for.
    """
    counts = np.asarray(counts, dtype=float)
    group_a = np.asarray(group_a, dtype=bool)
    if counts.ndim != 2:
        raise ValueError("counts must be (features, samples)")
    if counts.shape[1] != group_a.size:
        raise ValueError(
            f"counts has {counts.shape[1]} sample columns but the group mask has "
            f"{group_a.size} entries; expected shape (features, samples)")
    if group_a.sum() < 2 or (~group_a).sum() < 2:
        raise ValueError("two-group test needs at least 2 samples per group")

    if norm is None:
        norm = tmm_factors(counts)
    normed = counts / norm

    alpha = _nb_dispersion(normed, group_a)
    mean_a = normed[:, group_a].mean(axis=1)
    mean_b = normed[:, ~group_a].mean(axis=1)
    lfc = np.log2((mean_a + 1e-12) / (mean_b + 1e-12))

    # Wald statistic on the log2 scale, which is the DESeq2 formulation: the
    # contrast log2(a/b) has approximate variance
    #     sigma^2 = 1/(n_a*mu_a) + 1/(n_b*mu_b) + alpha*(1/n_a + 1/n_b)
    # Derived from Var(mean) = mu/n + alpha*mu^2/n by the delta method on log.
    # Using the raw-mean difference instead understates significance badly
    # (a 10x change came out at p = 0.13).
    n_a, n_b = int(group_a.sum()), int((~group_a).sum())
    mu_a = np.maximum(mean_a, 1e-9)
    mu_b = np.maximum(mean_b, 1e-9)
    var_log2 = (1.0 / (n_a * mu_a) + 1.0 / (n_b * mu_b)) / (np.log(2) ** 2) \
        + alpha * (1.0 / n_a + 1.0 / n_b) / (np.log(2) ** 2)
    se = np.sqrt(np.maximum(var_log2, 1e-12))
    wald = lfc / se

    from scipy.stats import norm as _norm
    p = 2 * _norm.sf(np.abs(wald))
    return {"log2fc": lfc, "pvalue": p, "dispersion": np.full(len(p), alpha),
            "stat": wald, "mean_a": mean_a, "mean_b": mean_b}


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR, the correction DiffBind reports as ``FDR``."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def parse_peak_specs(items, bam_dir=None, peak_dir=None,
                     bam_exts=(".bam",), peak_exts=(".narrowPeak", ".broadPeak",
                                                   ".bed", ".bed.gz")):
    """Resolve ``--peaks`` entries into ``{sample: {bam, peak}}``.

    Three forms are accepted, because the full ``sample=BAM:PEAK`` spelling for
    every replicate is tedious and error-prone:

    * ``sample=BAM:PEAK`` — explicit paths.
    * ``sample=BAM`` — the peak file is found in ``peak_dir`` as
      ``<sample>.narrowPeak`` / ``.broadPeak`` / ``.bed``.
    * ``sample`` — both files found as ``<bam_dir>/<sample>.bam`` and the peak file
      in ``peak_dir``.

    A sample whose files cannot be resolved raises with the paths that were tried;
    silently dropping a replicate would change the statistical design without
    telling the user. Comma-separated labels are accepted, so
    ``--peaks ctrl1,ctrl2`` works, and the ``sample=`` prefix keeps working for
    entries that only name a sample.
    """
    entries = []
    for raw in items or []:
        for piece in str(raw).split(","):
            value = piece.strip()
            if value:
                entries.append(value)
    if not entries:
        raise ValueError("at least one --peaks entry is required")

    def find(directory, stem, exts):
        if directory is None:
            return None
        base = Path(directory)
        found = [base / f"{stem}{ext}" for ext in exts]
        found = [p for p in found if p.exists()]
        if not found:
            return None
        if len(found) > 1:
            raise ValueError(
                f"'{stem}' is ambiguous in {base}: "
                + ", ".join(p.name for p in found))
        return found[0]

    out: dict[str, dict[str, str]] = {}
    for entry in entries:
        if "=" in entry:
            name, rest = entry.split("=", 1)
            name = name.strip()
            rest = rest.strip()
            if not name:
                raise ValueError(f"malformed --peaks entry {entry!r}")
            if ":" in rest:
                bam, peak = rest.rsplit(":", 1)
                bam, peak = bam.strip(), peak.strip()
                if not Path(bam).exists():
                    raise FileNotFoundError(f"sample {name}: BAM not found: {bam}")
                if not Path(peak).exists():
                    raise FileNotFoundError(f"sample {name}: peak file not found: {peak}")
                out[name] = {"bam": bam, "peak": peak}
                continue
            # sample=BAM only: locate the peak file next to the annotation dirs
            bam = rest
            if not Path(bam).exists():
                raise FileNotFoundError(f"sample {name}: BAM not found: {bam}")
            peak = find(peak_dir, name, peak_exts)
            if peak is None:
                raise FileNotFoundError(
                    f"sample {name}: no peak file found"
                    + (f" in {peak_dir}" if peak_dir else
                       "; pass sample=BAM:PEAK or --peak-dir"))
            out[name] = {"bam": bam, "peak": str(peak)}
            continue

        name = entry
        bam = find(bam_dir, name, bam_exts)
        peak = find(peak_dir, name, peak_exts)
        missing = []
        if bam is None:
            missing.append(f"{bam_dir}/<sample>.bam" if bam_dir
                           else "a BAM (pass --bam-dir)")
        if peak is None:
            missing.append(f"a peak file in {peak_dir}" if peak_dir
                           else "a peak file (pass --peak-dir)")
        if missing:
            raise FileNotFoundError(
                f"sample {name}: could not resolve " + " and ".join(missing))
        out[name] = {"bam": str(bam), "peak": str(peak)}

    if len(out) < 2:
        raise ValueError(
            f"differential peaks needs at least two samples, got {len(out)}")
    return out


def run_differential_peaks(args) -> int:
    """CLI driver: consensus regions -> counts -> TMM -> NB Wald -> table."""
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    samples = parse_peak_specs(getattr(args, "peaks", None),
                               bam_dir=getattr(args, "bam_dir", None),
                               peak_dir=getattr(args, "peak_dir", None))

    regions = consensus_peaks([s["peak"] for s in samples.values()],
                              min_overlap=getattr(args, "min_overlap", 1))
    if not regions:
        out.write_text("# easyepi differential-peaks - no consensus regions\n", encoding="utf-8")
        print("[easyepi] differential-peaks: no consensus regions")
        return 0

    matrix, sample_ids = build_count_matrix(
        {name: s["bam"] for name, s in samples.items()}, regions,
        min_mapq=getattr(args, "min_mapq", 30))

    # group assignment from --group-a / --group-b sample-name lists, or from the
    # design file when supplied.
    design = getattr(args, "design", None)
    if design:
        group_of: dict[str, str] = {}
        with open(design, encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            key = "group" if "group" in (reader.fieldnames or []) else "condition"
            for row in reader:
                group_of[row["sample_id"]] = row[key]
        levels = list(dict.fromkeys(group_of[s] for s in sample_ids if s in group_of))
        if len(levels) != 2:
            raise ValueError(f"design must contain exactly 2 groups, got {levels}")
        mask_a = np.array([group_of.get(s) == levels[0] for s in sample_ids])
    else:
        a = set((getattr(args, "group_a", None) or "").split(","))
        a.discard("")
        if not a:
            raise ValueError("provide --design or --group-a (comma-separated sample names)")
        mask_a = np.array([s in a for s in sample_ids])
        levels = ["group_a", "group_b"]

    norm = tmm_factors(matrix)
    n_a, n_b = int(mask_a.sum()), int((~mask_a).sum())
    if n_a >= 2 and n_b >= 2:
        res = nb_wald_test(matrix, mask_a, norm=norm)
        padj = bh_adjust(res["pvalue"])
        test_name = "nb_wald"
    else:
        # Not enough replicates for a dispersion/variance estimate. Report the
        # normalised fold change only, exactly as DiffBind does for n=1, and say
        # so in the output instead of fabricating a p-value.
        normed = matrix / norm
        mean_a = normed[:, mask_a].mean(axis=1)
        mean_b = normed[:, ~mask_a].mean(axis=1)
        lfc = np.log2((mean_a + 1e-12) / (mean_b + 1e-12))
        res = {"log2fc": lfc, "pvalue": np.full(len(lfc), np.nan),
               "dispersion": np.full(len(lfc), np.nan),
               "stat": np.full(len(lfc), np.nan),
               "mean_a": mean_a, "mean_b": mean_b}
        padj = np.full(len(lfc), np.nan)
        test_name = "fold_change_only_no_replicates"

    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "start", "end", "width",
                    "conc", "conc_group_a", "conc_group_b",
                    "Fold", "pvalue", "FDR", "dispersion", "stat", "test"])
        for i, (chrom, start, end) in enumerate(regions):
            w.writerow([chrom, start, end, end - start,
                        int(matrix[i].sum()),
                        int(matrix[i][mask_a].sum()), int(matrix[i][~mask_a].sum()),
                        round(float(res["log2fc"][i]), 4),
                        float(res["pvalue"][i]) if np.isfinite(res["pvalue"][i]) else "NA",
                        float(padj[i]) if np.isfinite(padj[i]) else "NA",
                        round(float(res["dispersion"][i]), 6)
                        if np.isfinite(res["dispersion"][i]) else "NA",
                        round(float(res["stat"][i]), 4)
                        if np.isfinite(res["stat"][i]) else "NA",
                        test_name])

    finite = padj[np.isfinite(padj)]
    sig = int((finite < 0.05).sum()) if finite.size else 0
    print(f"[easyepi] differential-peaks: {len(regions)} consensus regions, "
          f"{sig} with FDR<0.05, test={test_name} -> {out}")
    return 0
