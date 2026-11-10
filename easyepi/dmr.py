"""Differential methylation regions between two groups of bisulfite BAMs.

Scope, stated up front
----------------------
Two groups, any number of BAMs per group, one comparison. That is what the request
asked for and it is also the comparison almost every study makes. Multi-factor
designs, covariates and paired analyses are out of scope and are not silently
approximated.

Why per-cytosine, not per-CpG-dyad
----------------------------------
A CpG dyad carries a methylated C on each strand, and the two are separate
observations in the alignments. Merging them hides strand-specific effects (hemi-
methylation, asymmetric maintenance) and depends on the aligner's coordinate
convention for the reverse C. Sites are therefore kept strand-specific throughout and
the strand is a column in every output.

Test used, and its exact limitation
-----------------------------------
Per site this computes a **Wald test on the log odds ratio** of methylation
proportions, with the Haldane-Anscombe correction (0.5 added to each cell) so that
zero counts do not produce an infinite statistic:

    p1 = (mA + 0.5) / (nA + 1)          p2 = (mB + 0.5) / (nB + 1)
    LOR = logit(p1) - logit(p2)
    SE  = sqrt(1/(mA+0.5) + 1/(nA-mA+0.5) + 1/(mB+0.5) + 1/(nB-mB+0.5))

This is **not** a beta-binomial model. Real methylation data is overdispersed
relative to binomial, so the Wald test is anti-conservative at high coverage and
unstable at low coverage. Both are mitigated by requiring a minimum coverage per
group (``--min-coverage``, default 5) and a minimum absolute difference
(``--min-diff``), and by merging into regions where the region-level test is
recomputed from summed counts. It is not a substitute for DSS/bsseq/metilene on a
publication-grade analysis; it is a fast, transparent screen, and the docstring in
``docs/RELIABILITY.md`` says so.

Region calling
--------------
Significant sites on the same chromosome and strand are merged when separated by at
most ``--max-gap`` bp. Each region's own test is recomputed from the summed counts,
so a region is reported with its own effect size and p-value rather than with the
minimum p-value of its members.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


def _require(name: str, hint: str):
    from easyepi.core.io import _require as _r
    return _r(name, hint)


# ---------------------------------------------------------------------------
# Per-site counts
# ---------------------------------------------------------------------------
def site_counts(bam_path: str | Path, *, chroms: Sequence[str] | None = None,
                contexts: Iterable[str] = ("CpG",), min_mapq: int = 10,
                skip_duplicates: bool = True,
                max_reads: int | None = None):
    """Methylated / total counts per cytosine, one row per (chrom, pos, strand).

    Chromosomes are flushed one at a time so peak memory is bounded by the largest
    contig rather than by the genome: a whole-genome dict of 28 M sites would cost
    several GB, a single chromosome is a few tens of MB.

    Returns ``(frame, stats)``; ``stats`` names the tag encoding and the read and
    site counts so a suspiciously small number is visible rather than silent.
    """
    pysam = _require("pysam", "pip install 'easyepi[bam]'")
    import pandas as pd

    from easyepi.methylation_bam import (_calls_from_mm, _calls_from_xm,
                                         has_methylation_tags)

    probe = has_methylation_tags(bam_path)
    if not (probe["has_mm"] or probe["has_xm"]):
        raise ValueError(
            f"{bam_path} carries neither an MM/ML tag nor a Bismark XM tag, so it "
            "contains no methylation calls. A bisulfite alignment writes one of "
            "them; a BAM without either is not usable for methylation analysis.")

    want = {c.strip() for c in contexts}
    unknown = want - {"CpG", "CHG", "CHH"}
    if unknown:
        raise ValueError(
            f"unknown methylation context(s) {sorted(unknown)}; use CpG, CHG or CHH")

    rows: list[tuple[str, int, str, int, int]] = []
    reads = used = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        targets = list(chroms) if chroms else list(bam.references)
        for chrom in targets:
            if chrom not in bam.references:
                raise ValueError(
                    f"chrom {chrom!r} is not in {bam_path} "
                    f"(has e.g. {list(bam.references)[:5]})")
            # (pos, strand) -> [methylated, total]
            tally: dict[tuple[int, str], list[int]] = {}
            for rec in bam.fetch(chrom):
                reads += 1
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                if skip_duplicates and rec.is_duplicate:
                    continue
                if rec.mapping_quality < min_mapq:
                    continue
                used += 1
                gen = (_calls_from_mm(rec, chrom) if rec.has_tag("MM")
                       else _calls_from_xm(rec, chrom))
                for call in gen:
                    if call.context not in want:
                        continue
                    key = (call.pos, call.strand)
                    cell = tally.get(key)
                    if cell is None:
                        tally[key] = [1 if call.methylated else 0, 1]
                    else:
                        cell[1] += 1
                        if call.methylated:
                            cell[0] += 1
                if max_reads is not None and used >= max_reads:
                    break
            for (pos, strand), (nmeth, ntot) in tally.items():
                rows.append((chrom, int(pos), strand, int(nmeth), int(ntot)))
            if max_reads is not None and used >= max_reads:
                break

    frame = pd.DataFrame(rows, columns=["chrom", "pos", "strand", "n_meth", "n_total"])
    if not frame.empty:
        frame = frame.sort_values(["chrom", "pos", "strand"], kind="stable")
    stats = {"bam": str(bam_path),
             "encoding": "MM/ML" if probe["has_mm"] else "XM",
             "reads_seen": reads, "reads_used": used,
             "sites": int(len(frame)), "contexts": sorted(want)}
    return frame.reset_index(drop=True), stats


def group_counts(samples: Sequence[str | Path], *, label: str = "group",
                 contexts: Iterable[str] = ("CpG",), min_mapq: int = 10,
                 chroms: Sequence[str] | None = None,
                 max_reads: int | None = None):
    """Sum per-site counts across the samples of one group.

    Samples are summed, not averaged: with unequal coverage an average would weight
    a low-coverage replicate as heavily as a deep one, which changes the result.
    """
    import pandas as pd

    if not samples:
        raise ValueError(f"{label} has no BAMs")

    per_sample = []
    sample_stats = []
    for bam in samples:
        frame, stats = site_counts(bam, chroms=chroms, contexts=contexts,
                                   min_mapq=min_mapq, max_reads=max_reads)
        per_sample.append(frame.set_index(["chrom", "pos", "strand"]))
        sample_stats.append(stats)

    if not per_sample:
        empty = pd.DataFrame(columns=["n_meth", "n_total"])
        return empty, {"samples": [], "sites": 0}

    stacked = pd.concat(per_sample, axis=1, keys=range(len(per_sample)))
    meth_cols = [c for c in stacked.columns if c[1] == "n_meth"]
    total_cols = [c for c in stacked.columns if c[1] == "n_total"]
    summed = pd.DataFrame({
        "n_meth": stacked[meth_cols].fillna(0).sum(axis=1),
        "n_total": stacked[total_cols].fillna(0).sum(axis=1),
    })
    summed = summed[summed["n_total"] > 0]
    return summed, {"samples": [s["bam"] for s in sample_stats],
                    "sites": int(len(summed)),
                    "per_sample": sample_stats}


# ---------------------------------------------------------------------------
# Per-site test
# ---------------------------------------------------------------------------
def _wald_log_odds(m_a, n_a, m_b, n_b):
    """Vectorised Wald test on the log odds ratio with a 0.5 continuity correction.

    Returns ``(diff, statistic, pvalue)`` where ``diff`` is the methylation-level
    difference ``pB - pA``. The correction keeps zero cells finite; without it a site
    methylated in every read of one group and none of the other gives an infinite
    statistic and a zero p-value, which is not a real result at that depth.
    """
    from scipy.stats import norm

    m_a = np.asarray(m_a, dtype=float)
    n_a = np.asarray(n_a, dtype=float)
    m_b = np.asarray(m_b, dtype=float)
    n_b = np.asarray(n_b, dtype=float)

    p_a = (m_a + 0.5) / (n_a + 1.0)
    p_b = (m_b + 0.5) / (n_b + 1.0)
    diff = p_b - p_a

    with np.errstate(divide="ignore", invalid="ignore"):
        lor = np.log(p_a / (1.0 - p_a)) - np.log(p_b / (1.0 - p_b))
        var = (1.0 / (m_a + 0.5) + 1.0 / (n_a - m_a + 0.5)
               + 1.0 / (m_b + 0.5) + 1.0 / (n_b - m_b + 0.5))
        se = np.sqrt(var)
        stat = lor / se
    pvalue = 2.0 * norm.sf(np.abs(stat))
    bad = ~np.isfinite(stat)
    stat = np.where(bad, np.nan, stat)
    pvalue = np.where(bad, np.nan, pvalue)
    return diff, stat, pvalue


def differential_sites(a_counts, b_counts, *, min_coverage: int = 5,
                       min_diff: float = 0.0, alpha: float = 0.05):
    """Compare two groups site by site and return the tested table.

    ``min_coverage`` is applied **per group**, not to the sum: a site covered 1x in
    one group and 50x in the other cannot estimate that group's level, and allowing
    it through would let a single read drive a "significant" call.

    ``padj`` is Benjamini-Hochberg across the tested sites. Sites failing a filter
    are dropped before correction, so the FDR refers to the sites that were actually
    testable and ``stats`` reports how many were dropped and why.
    """
    import pandas as pd

    from easyepi.core.stats import benjamini_hochberg

    joined = a_counts.join(b_counts, how="outer", lsuffix="_a", rsuffix="_b")
    joined = joined.fillna(0.0)
    joined = joined.rename(columns={"n_meth_a": "m_a", "n_total_a": "n_a",
                                    "n_meth_b": "m_b", "n_total_b": "n_b"})
    for col in ("m_a", "n_a", "m_b", "n_b"):
        joined[col] = joined[col].astype(np.int64)

    n_joined = len(joined)
    keep = ((joined["n_a"] >= min_coverage) & (joined["n_b"] >= min_coverage))
    tested = joined[keep].copy()
    stats = {"sites_union": int(n_joined), "sites_tested": int(len(tested)),
             "dropped_low_coverage": int(n_joined - len(tested)),
             "min_coverage": int(min_coverage)}

    if tested.empty:
        tested["diff"] = []
        tested["statistic"] = []
        tested["pvalue"] = []
        tested["padj"] = []
        stats.update({"significant": 0, "min_diff": float(min_diff),
                      "alpha": float(alpha)})
        return tested.reset_index(), stats

    tested["level_a"] = tested["m_a"] / tested["n_a"]
    tested["level_b"] = tested["m_b"] / tested["n_b"]
    diff, stat, pvalue = _wald_log_odds(tested["m_a"], tested["n_a"],
                                        tested["m_b"], tested["n_b"])
    tested["diff"] = diff
    tested["statistic"] = stat
    tested["pvalue"] = pvalue

    valid = tested["pvalue"].notna()
    padj = np.full(len(tested), np.nan)
    if valid.any():
        padj[valid.to_numpy()] = benjamini_hochberg(
            tested.loc[valid, "pvalue"].to_numpy())
    tested["padj"] = padj

    tested["significant"] = ((tested["padj"] < alpha)
                             & (tested["diff"].abs() >= min_diff))
    stats.update({"significant": int(tested["significant"].sum()),
                  "min_diff": float(min_diff), "alpha": float(alpha),
                  "sites_untestable": int((~valid).sum())})
    return tested.reset_index(), stats


# ---------------------------------------------------------------------------
# Region merging
# ---------------------------------------------------------------------------
def merge_dmr_regions(tested, *, max_gap: int = 100, min_sites: int = 1,
                      alpha: float = 0.05, min_diff: float = 0.0):
    """Merge adjacent significant sites into regions and retest each region.

    Sites are grouped by (chrom, strand) and merged when consecutive positions are
    at most ``max_gap`` apart, so a run of significant CpGs in a CpG island becomes
    one region while two islands separated by a few hundred bp stay separate.

    The region's p-value comes from the **summed** counts of its members, not from
    the smallest member p-value: summing independent observations increases the
    evidence, and quoting the minimum would multiply-test within the region.
    """
    import pandas as pd

    from easyepi.core.stats import benjamini_hochberg

    if tested.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "n_sites",
                                     "m_a", "n_a", "m_b", "n_b", "level_a",
                                     "level_b", "diff", "pvalue", "padj"]), {}

    sig = tested[tested["padj"].notna() & (tested["padj"] < alpha)
                 & (tested["diff"].abs() >= min_diff)].copy()
    if sig.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "n_sites",
                                     "m_a", "n_a", "m_b", "n_b", "level_a",
                                     "level_b", "diff", "pvalue", "padj"]), {
            "regions": 0, "sites_used": 0}

    regions: list[dict[str, Any]] = []
    for (chrom, strand), block in sig.groupby(["chrom", "strand"], sort=True):
        block = block.sort_values("pos")
        pos = block["pos"].to_numpy()
        if pos.size == 0:
            continue
        breaks = np.flatnonzero(np.diff(pos) > max_gap) + 1
        for chunk in np.split(np.arange(pos.size), breaks):
            if chunk.size == 0:
                continue
            sub = block.iloc[chunk]
            if len(sub) < min_sites:
                continue
            regions.append({
                "chrom": chrom,
                "start": int(sub["pos"].min()),
                "end": int(sub["pos"].max()) + 1,
                "strand": strand,
                "n_sites": int(len(sub)),
                "m_a": int(sub["m_a"].sum()), "n_a": int(sub["n_a"].sum()),
                "m_b": int(sub["m_b"].sum()), "n_b": int(sub["n_b"].sum()),
            })

    if not regions:
        return pd.DataFrame(columns=["chrom", "start", "end", "strand", "n_sites",
                                     "m_a", "n_a", "m_b", "n_b", "level_a",
                                     "level_b", "diff", "pvalue", "padj"]), {
            "regions": 0, "sites_used": int(len(sig))}

    frame = pd.DataFrame(regions)
    frame["level_a"] = frame["m_a"] / frame["n_a"]
    frame["level_b"] = frame["m_b"] / frame["n_b"]
    diff, stat, pvalue = _wald_log_odds(frame["m_a"], frame["n_a"],
                                        frame["m_b"], frame["n_b"])
    frame["diff"] = diff
    frame["statistic"] = stat
    frame["pvalue"] = pvalue

    valid = frame["pvalue"].notna()
    padj = np.full(len(frame), np.nan)
    if valid.any():
        padj[valid.to_numpy()] = benjamini_hochberg(
            frame.loc[valid, "pvalue"].to_numpy())
    frame["padj"] = padj

    frame = frame.sort_values(["chrom", "start", "strand"], kind="stable")
    # Direction comes from the region's own recomputed effect size. Using the
    # members' mean diff (as an earlier version did) lets a signed average disagree
    # with the summed-count statistic, so a region could be labelled "hyper" while
    # its own test says otherwise.
    frame["direction"] = np.where(frame["diff"] > 0, "hyper", "hypo")
    stats = {"regions": int(len(frame)), "sites_used": int(len(sig)),
             "max_gap": int(max_gap), "min_sites": int(min_sites),
             "significant_regions": int((frame["padj"] < alpha).sum())}
    cols = ["chrom", "start", "end", "strand", "n_sites", "m_a", "n_a", "m_b",
            "n_b", "level_a", "level_b", "diff", "direction", "statistic",
            "pvalue", "padj"]
    return frame[cols].reset_index(drop=True), stats


# ---------------------------------------------------------------------------
# CLI driver
# ---------------------------------------------------------------------------
def _split_bams(values: Sequence[str] | None, label: str) -> list[str]:
    """Flatten repeated and comma-separated BAM paths for one group."""
    out: list[str] = []
    for raw in values or []:
        for piece in str(raw).split(","):
            path = piece.strip()
            if not path:
                continue
            if not Path(path).exists():
                raise FileNotFoundError(f"{label}: BAM not found: {path}")
            out.append(path)
    if not out:
        raise ValueError(f"{label} is empty; pass at least one BAM")
    dupes = sorted({p for p in out if out.count(p) > 1})
    if dupes:
        raise ValueError(f"{label} lists the same BAM twice: {', '.join(dupes)}")
    return out


def run_dmr(args) -> int:
    """CLI driver: two groups of bisulfite BAMs -> site table + DMR table."""
    import pandas as pd

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    group_a = _split_bams(getattr(args, "group_a", None), "--group-a")
    group_b = _split_bams(getattr(args, "group_b", None), "--group-b")
    overlap = set(group_a) & set(group_b)
    if overlap:
        raise ValueError(
            "the same BAM appears in both groups: " + ", ".join(sorted(overlap)))

    chroms = None
    if getattr(args, "chrom", None):
        chroms = [c.strip() for c in str(args.chrom).split(",") if c.strip()]

    print(f"[easyepi] dmr: {len(group_a)} vs {len(group_b)} sample(s)")
    for label, group in ((getattr(args, "label_a", None) or "A", group_a),
                         (getattr(args, "label_b", None) or "B", group_b)):
        print(f"  group {label}: " + ", ".join(Path(b).name for b in group))

    counts_a, stats_a = group_counts(group_a, label="--group-a",
                                     contexts=(getattr(args, "context", None) or "CpG",),
                                     min_mapq=int(args.min_mapq), chroms=chroms)
    counts_b, stats_b = group_counts(group_b, label="--group-b",
                                     contexts=(getattr(args, "context", None) or "CpG",),
                                     min_mapq=int(args.min_mapq), chroms=chroms)
    print(f"  sites: A={stats_a['sites']:,}  B={stats_b['sites']:,}")

    tested, test_stats = differential_sites(
        counts_a, counts_b, min_coverage=int(args.min_coverage),
        min_diff=float(args.min_diff), alpha=float(args.alpha))
    print(f"  tested {test_stats['sites_tested']:,} site(s); "
          f"{test_stats['significant']:,} significant "
          f"(dropped {test_stats['dropped_low_coverage']:,} below "
          f"{args.min_coverage}x in a group)")

    site_path = out.with_name(out.stem + ".sites.tsv")
    tested.to_csv(site_path, sep="\t", index=False)

    regions, region_stats = merge_dmr_regions(
        tested, max_gap=int(args.max_gap), min_sites=int(args.min_sites),
        alpha=float(args.alpha), min_diff=float(args.min_diff))
    regions.to_csv(out, sep="\t", index=False)

    summary = {
        "groups": {"a": {"bams": group_a, "label": getattr(args, "label_a", None) or "A"},
                   "b": {"bams": group_b, "label": getattr(args, "label_b", None) or "B"}},
        "parameters": {"context": getattr(args, "context", None) or "CpG",
                       "min_coverage": int(args.min_coverage),
                       "min_mapq": int(args.min_mapq),
                       "min_diff": float(args.min_diff),
                       "alpha": float(args.alpha),
                       "max_gap": int(args.max_gap),
                       "min_sites": int(args.min_sites)},
        "counts": {"group_a": stats_a, "group_b": stats_b},
        "site_test": test_stats,
        "regions": region_stats,
        "method": "per-site Wald log-odds (Haldane-Anscombe), BH FDR, "
                  "gap-merged regions retested on summed counts",
        "limitation": "binomial-based, not beta-binomial; overdispersion is not "
                      "modelled, so p-values are anti-conservative at high coverage",
        "outputs": {"sites": str(site_path), "regions": str(out)},
    }
    summary_path = out.with_name(out.stem + ".json")
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False),
                            encoding="utf-8")

    print(f"  regions: {region_stats.get('regions', 0)} "
          f"({region_stats.get('significant_regions', 0)} at padj < {args.alpha})")
    if not regions.empty:
        top = regions.reindex(regions["padj"].abs().sort_values().index).head(3)
        for _, row in top.iterrows():
            print(f"    {row['chrom']}:{int(row['start'])}-{int(row['end'])} "
                  f"{row['strand']} {row['direction']:5s} "
                  f"diff={row['diff']:+.3f} padj={row['padj']:.2e} "
                  f"sites={int(row['n_sites'])}")
    print(f"[easyepi] wrote {out.name}, {site_path.name}, {summary_path.name}")
    return 0
