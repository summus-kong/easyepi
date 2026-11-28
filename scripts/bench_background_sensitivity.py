"""Run the background sensitivity analysis end to end on real peak sequences.

This is the evidence for the roadmap's M3: does the choice of background change a
motif enrichment ranking, and by how much.

Three configurations, all against the same target sequences and the same motifs:

  genomic            random genome windows, composition-unmatched
  gc_length_matched  `easyepi background` output (GC- and length-matched)
  accessible_matched same, restricted to a supplied accessible universe

Every configuration records its own GC and length distributions so the comparison
can be read alongside the question of whether the nulls actually differed.

Usage (on the analysis host, with the easyepi environment):
    python scripts/bench_background_sensitivity.py \
        --peaks peaks.bed --genome-fasta genome.fa \
        --motifs share/homer/motifs --accessible accessible.bed \
        --chromosomes chr1 --out bg_sensitivity/ --seed 0
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np


def _read_bed(path: str | Path) -> List[tuple[str, int, int]]:
    regions: List[tuple[str, int, int]] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            try:
                regions.append((fields[0], int(fields[1]), int(fields[2])))
            except ValueError:
                continue
    if not regions:
        raise ValueError(f"no regions in {path}")
    return regions


def _fetch(genome_fasta: str, regions: Sequence[tuple[str, int, int]]) -> List[str]:
    """Sequences for the regions; empty strings where the FASTA lacks a contig."""
    try:
        import pysam
    except Exception:
        pysam = None

    if pysam is not None:
        out: List[str] = []
        with pysam.FastaFile(genome_fasta) as fa:
            for chrom, start, end in regions:
                try:
                    seq = fa.fetch(chrom, int(start), int(end))
                except Exception:
                    seq = ""
                out.append((seq or "").upper())
        return out

    from easyepi.gc_background import iter_fasta
    wanted: Dict[str, List[int]] = {}
    for i, (chrom, _s, _e) in enumerate(regions):
        wanted.setdefault(chrom, []).append(i)
    out = [""] * len(regions)
    for name, seq in iter_fasta(genome_fasta):
        for i in wanted.get(name, []):
            _c, s, e = regions[i]
            out[i] = seq[int(s):int(e)]
    return out


def _random_genomic(targets: Sequence[tuple[str, int, int]],
                    genome_fasta: str,
                    chromosomes: Sequence[str],
                    rng: random.Random) -> List[tuple[str, int, int]]:
    """Random windows of the targets' widths, anywhere on the chosen chromosomes.

    Deliberately unmatched: this is the null a naive analysis uses, and the point
    of the comparison is to show what it costs.
    """
    try:
        import pysam
        with pysam.FastaFile(genome_fasta) as fa:
            lengths = {c: fa.get_reference_length(c) for c in chromosomes
                       if c in set(fa.references)}
    except Exception:
        from easyepi.gc_background import iter_fasta
        lengths = {}
        for name, seq in iter_fasta(genome_fasta):
            if name in set(chromosomes):
                lengths[name] = len(seq)

    usable = [c for c in chromosomes if lengths.get(c, 0) > 1000]
    if not usable:
        raise ValueError("no usable chromosome for random background")

    out: List[tuple[str, int, int]] = []
    for _chrom, s, e in targets:
        width = int(e - s)
        for _ in range(200):  # bounded retries per target
            c = rng.choice(usable)
            limit = lengths[c] - width
            if limit <= 0:
                continue
            start = rng.randint(0, limit)
            out.append((c, start, start + width))
            break
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--peaks", required=True, help="target regions (BED)")
    parser.add_argument("--genome-fasta", required=True)
    parser.add_argument("--motifs", required=True,
                        help="motif file or directory of HOMER .motif files")
    parser.add_argument("--accessible", default=None,
                        help="accessible BED; enables the third configuration")
    parser.add_argument("--chromosomes", default=None,
                        help="comma-separated chromosomes to consider")
    parser.add_argument("--window", type=int, default=1000)
    parser.add_argument("--tolerance", type=float, default=0.05)
    parser.add_argument("--min-distance", type=int, default=2000)
    parser.add_argument("--max-motifs", type=int, default=200,
                        help="cap the motif library so the run stays interactive")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    from easyepi.background_report import compare_distributions
    from easyepi.gc_background import (gc_fraction, genome_windows,
                                       match_background, read_accessible_bed)
    from easyepi.motifs import formats as F
    from easyepi.motifs.background_sensitivity import compare_rankings
    from easyepi.motifs.homer_workflow import known_motif_enrichment_homer

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    chromosomes = ([c.strip() for c in args.chromosomes.split(",") if c.strip()]
                   if args.chromosomes else None)
    targets = _read_bed(args.peaks)
    if chromosomes:
        targets = [t for t in targets if t[0] in set(chromosomes)]
    if not targets:
        raise SystemExit("no target regions left after chromosome filtering")

    # --- motifs -----------------------------------------------------------
    motif_path = Path(args.motifs)
    if motif_path.is_dir():
        files = sorted(motif_path.glob("*.motif"))[: args.max_motifs]
        motifs = []
        for f in files:
            try:
                motifs.extend(F.read_homer_motifs(str(f)))
            except Exception:
                continue
    else:
        motifs = F.read_motifs(motif_path)[: args.max_motifs]
    if not motifs:
        raise SystemExit(f"no motifs parsed from {args.motifs}")
    print(f"[sensitivity] {len(motifs)} motif(s) from {args.motifs}")

    target_seqs = _fetch(args.genome_fasta, targets)
    target_seqs = [s for s in target_seqs if s]
    print(f"[sensitivity] {len(target_seqs)} target sequence(s)")

    # --- background configurations ---------------------------------------
    backgrounds: Dict[str, List[tuple[str, int, int]]] = {}

    rng = random.Random(args.seed)
    chromosome_pool = chromosomes or sorted({t[0] for t in targets})
    backgrounds["genomic"] = _random_genomic(targets, args.genome_fasta,
                                             chromosome_pool, rng)

    windows = genome_windows(args.genome_fasta, window=int(args.window),
                             chromosomes=chromosomes)
    ctrl, summary = match_background(targets, windows,
                                     tolerance=float(args.tolerance),
                                     min_distance=int(args.min_distance),
                                     seed=int(args.seed))
    backgrounds["gc_length_matched"] = ctrl
    matched_summary = summary

    if args.accessible:
        accessible = read_accessible_bed(args.accessible)
        ctrl_acc, summary_acc = match_background(
            targets, windows, tolerance=float(args.tolerance),
            min_distance=int(args.min_distance), seed=int(args.seed),
            accessible=accessible)
        backgrounds["accessible_matched"] = ctrl_acc
        matched_summary_acc = summary_acc
    else:
        matched_summary_acc = None

    # --- enrichment per configuration ------------------------------------
    rankings: Dict[str, List[Dict]] = {}
    composition: Dict[str, Any] = {}
    limits: Dict[str, Any] = {}
    for label, regions in backgrounds.items():
        seqs = [s for s in _fetch(args.genome_fasta, regions) if s]
        if not seqs:
            print(f"[sensitivity] {label}: no background sequences, skipped")
            continue
        limits[label] = {
            "n_regions": len(regions),
            "n_sequences_fetched": len(seqs),
            "median_length": float(np.median([len(s) for s in seqs])),
        }
        composition[label] = compare_distributions(
            [gc_fraction(s) for s in target_seqs],
            [gc_fraction(s) for s in seqs],
            [len(s) for s in target_seqs],
            [len(s) for s in seqs],
        )
        rows = known_motif_enrichment_homer(target_seqs, seqs, motifs)
        rankings[label] = rows
        print(f"[sensitivity] {label}: {len(seqs)} background sequences, "
              f"{len(rows)} motif rows")

    if len(rankings) < 2:
        raise SystemExit("fewer than two usable background configurations")

    # The base is the *unmatched* null: the comparison question is "what changes
    # when the null is made appropriate", so the naive one has to be the reference.
    ordered = {k: rankings[k] for k in
               ("genomic", "gc_length_matched", "accessible_matched")
               if k in rankings}
    comparison = compare_rankings(ordered, top_n=10)

    payload = {
        "params": {
            "peaks": str(args.peaks),
            "genome_fasta": str(args.genome_fasta),
            "motifs": str(args.motifs),
            "accessible": str(args.accessible) if args.accessible else None,
            "chromosomes": chromosomes,
            "window": int(args.window),
            "tolerance": float(args.tolerance),
            "min_distance": int(args.min_distance),
            "seed": int(args.seed),
            "n_target_sequences": len(target_seqs),
            "n_motifs": len(motifs),
        },
        "background_limits": limits,
        "composition": composition,
        "matching": {
            "gc_length_matched": {
                "n_matched": matched_summary.get("n_matched"),
                "frac_matched": matched_summary.get("frac_matched"),
                "design": matched_summary.get("background_design"),
            },
            **({"accessible_matched": {
                "n_matched": matched_summary_acc.get("n_matched"),
                "frac_matched": matched_summary_acc.get("frac_matched"),
                "design": matched_summary_acc.get("background_design"),
            }} if matched_summary_acc else {}),
        },
        "ranking_comparison": comparison,
        "rankings": {k: [{"name": r.get("name"), "log_p": r.get("log_p"),
                          "enrichment_fold": r.get("enrichment_fold"),
                          "qvalue": r.get("qvalue")} for r in v[:50]]
                     for k, v in ordered.items()},
    }
    out_file = out_dir / "background_sensitivity.json"
    out_file.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    # --- report ----------------------------------------------------------
    print()
    for label, cmp in comparison["comparisons"].items():
        rho = cmp["spearman_rho"]
        rho_txt = "undefined" if rho is None else f"{rho:.3f}"
        print(f"[sensitivity] genomic vs {label}: "
              f"spearman rho = {rho_txt}, "
              f"top-10 overlap = {cmp['top_overlap']['n']}/10, "
              f"{cmp['n_shared_motifs']} shared motifs")
        for move in cmp["largest_rank_moves"][:3]:
            print(f"    {move['name'][:38]:38s} rank {move['rank_in_base']:>3} -> "
                  f"{move['rank_here']:>3}  (delta {move['delta']:+d})")
    print()
    print(f"[sensitivity] wrote {out_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
