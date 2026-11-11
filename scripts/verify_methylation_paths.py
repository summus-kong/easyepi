#!/usr/bin/env python3
"""Are the two methylation QC paths equivalent?

Two code paths survive from the merge of the ``*-qc`` aliases:

* ``qc --assay methylation`` reads a **bisulfite BAM** (MM/ML or Bismark XM tags) and
  summarises the beta values with :func:`easyepi.assay_metrics.beta_distribution`.
* ``easyepi.methylation.summarize_methylation`` reads a **table** (beta or
  methylated/unmethylated counts) and reports its own field set.

Both claim to describe the same thing, so they are run on the *same* underlying
methylation levels: a synthetic BAM is written with known per-site states, its beta
values are extracted, and the identical vector is fed to the table summariser. Any
difference is then a real inconsistency between the two paths rather than a
difference in input data.

What is compared
----------------
``n_sites``, ``mean_beta``, ``median_beta`` and the hypomethylated / hypermethylated
fractions. The two paths name things differently (``fraction_*`` vs ``pct_*``), so the
comparison maps fields explicitly instead of assuming identical keys.

What is not
-----------
The two paths do not compute identical *sets* of statistics: the BAM path adds a
histogram and a bimodality flag, the table path adds coverage summaries. That is
acceptable. What must not differ is a statistic computed from the same beta vector.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUT = Path("/data1/GaoLab/houoweibo/easyepi_work/validation_meth")


def build_bam(path: Path, states: list[tuple[int, bool]]) -> Path:
    """One single-base read per recorded site, XM-coded, on chr1."""
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 100_000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, (pos, methylated) in enumerate(sorted(states)):
            rec = pysam.AlignedSegment()
            rec.query_name = f"r{i}"
            rec.query_sequence = "C"
            rec.flag = 0
            rec.reference_id = 0
            rec.reference_start = int(pos)
            rec.mapping_quality = 60
            rec.cigartuples = [(0, 1)]
            rec.query_qualities = pysam.qualitystring_to_array("I")
            rec.set_tag("XM", "Z" if methylated else "z")
            out.write(rec)
    pysam.index(str(path))
    return path


def main() -> int:
    from easyepi.assay_metrics import beta_distribution
    from easyepi.methylation import read_methylation_table, summarize_methylation
    from easyepi.methylation_bam import (cpg_level_methylation, read_methylation_calls,
                                         calls_to_dataframe)

    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(3)

    # A realistic bimodal methylation profile: most CpGs near 0 or 1, some between.
    n_sites = 4000
    levels = np.where(rng.random(n_sites) < 0.45,
                      rng.beta(0.6, 6.0, n_sites),      # low
                      rng.beta(6.0, 0.6, n_sites))      # high
    levels = np.clip(levels, 0.0, 1.0)

    # Each site becomes 10 reads whose methylation count reproduces its level.
    states: list[tuple[int, bool]] = []
    pos = 1000
    truth: dict[int, float] = {}
    for level in levels:
        n_meth = int(round(level * 10))
        truth[pos] = n_meth / 10.0
        for k in range(10):
            states.append((pos, k < n_meth))
        pos += 3
    bam = build_bam(OUT / "meth.bam", states)
    print(f"wrote {bam.name}: {len(states)} reads over {len(truth)} sites")

    # --- path 1: BAM -------------------------------------------------------
    calls, stats = read_methylation_calls(bam, chrom="chr1", contexts=("CpG",),
                                          max_reads=10_000_000)
    frame = calls_to_dataframe(calls, {"chr1": 100_000})
    sites, beta_stats = cpg_level_methylation(frame, min_coverage=1)
    bam_beta = sites["beta"].to_numpy(dtype=float)
    bam_summary = beta_distribution(bam_beta)
    print(f"BAM path  : {len(bam_beta)} sites, mean {bam_summary['mean_beta']:.4f}")

    # --- path 2: table -----------------------------------------------------
    table_path = OUT / "meth.tsv"
    sites[["chrom", "start", "end", "beta"]].to_csv(table_path, sep="\t", index=False)
    table = read_methylation_table(table_path)
    table_summary = summarize_methylation(table, min_coverage=0)
    print(f"table path: {len(table)} rows, mean {table_summary.get('mean_beta'):.4f}")

    # Sanity: the beta vector the BAM path produced must match the planted levels.
    planted = np.array([truth[p] for p in sites["start"]], dtype=float)
    reconstruction = float(np.max(np.abs(bam_beta - planted)))
    print(f"BAM beta vs planted levels: max abs difference {reconstruction:.2e}")

    def pct(d, key):
        v = d.get(key)
        return None if v is None else float(v)

    comparison = {
        "n_sites": {"bam": bam_summary["n_sites"], "table": table_summary.get("n_sites")},
        "mean_beta": {"bam": round(bam_summary["mean_beta"], 6),
                      "table": round(float(table_summary.get("mean_beta", float("nan"))), 6)},
        "median_beta": {"bam": round(bam_summary["median_beta"], 6),
                        "table": round(float(table_summary.get("median_beta", float("nan"))), 6)},
        "fraction_hypomethylated": {
            "bam": round(pct(bam_summary, "fraction_hypomethylated"), 6),
            "table": round(float(table_summary.get("pct_hypomethylated", float("nan"))) / 100.0, 6)
            if table_summary.get("pct_hypomethylated") is not None else None,
        },
        "fraction_hypermethylated": {
            "bam": round(pct(bam_summary, "fraction_hypermethylated"), 6),
            "table": round(float(table_summary.get("pct_hypermethylated", float("nan"))) / 100.0, 6)
            if table_summary.get("pct_hypermethylated") is not None else None,
        },
    }
    print("\nfield comparison (bam vs table):")
    for key, pair in comparison.items():
        print(f"  {key:26s} {pair['bam']}  vs  {pair['table']}")

    # --- gates -------------------------------------------------------------
    gates: dict[str, bool] = {}
    gates["beta_reconstructs_planted_levels"] = reconstruction < 1e-9
    gates["same_site_count"] = comparison["n_sites"]["bam"] == comparison["n_sites"]["table"]
    for key in ("mean_beta", "median_beta", "fraction_hypomethylated",
                "fraction_hypermethylated"):
        b, t = comparison[key]["bam"], comparison[key]["table"]
        gates[f"{key}_agrees"] = (b is not None and t is not None
                                  and abs(b - t) < 1e-6)

    report = {"bam_stats": stats, "sites": int(len(bam_beta)),
              "comparison": comparison, "gates": gates,
              "note": "both paths consumed the identical beta vector; the BAM path's "
                      "own read-to-beta step is checked separately against the "
                      "planted levels"}
    (OUT / "equivalence.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("\ngates:")
    for name, ok in gates.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    print(f"\nreport: {OUT / 'equivalence.json'}")
    return 0 if all(gates.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
