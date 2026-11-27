"""`easyepi background`: build a matched control set and report whether it matched.

The command exists because "GC-matched background" is a claim, and a claim needs
an artefact. This writes the control BED plus a JSON report that a reader can
check: the two GC and length distributions, a KS test on each, the fraction of
targets that received a control, and the seed.

The final BED is not enough on its own: a run that matched 40% of targets and
dropped the rest produces a perfectly ordinary-looking control file.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

import numpy as np

from easyepi.background_report import background_report
from easyepi.gc_background import (accessible_fraction, genome_windows,
                                   match_background, read_accessible_bed,
                                   write_bed)


def _read_bed(path: str | Path) -> list[tuple[str, int, int]]:
    """Read BED3+ into (chrom, start, end), skipping track/comment lines."""
    regions: list[tuple[str, int, int]] = []
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
        raise ValueError(f"no regions read from {path}")
    return regions


def _fetch(genome_fasta: str, regions: Sequence[tuple[str, int, int]]) -> list[str]:
    """Sequences for the regions, or empty strings where the FASTA lacks them.

    pysam is optional, so the fallback reads the FASTA once into memory. Either
    way a missing contig yields "" rather than aborting: the report then shows a
    GC distribution computed from fewer sequences, which is visible, instead of
    failing after a long window scan.
    """
    try:
        import pysam
    except Exception:
        pysam = None

    if pysam is not None:
        out: list[str] = []
        with pysam.FastaFile(genome_fasta) as fa:
            for chrom, start, end in regions:
                try:
                    seq = fa.fetch(chrom, int(start), int(end))
                except Exception:
                    seq = ""
                out.append((seq or "").upper())
        return out

    from easyepi.gc_background import iter_fasta

    wanted: dict[str, list[int]] = {}
    for i, (chrom, _s, _e) in enumerate(regions):
        wanted.setdefault(chrom, []).append(i)
    out = [""] * len(regions)
    for name, seq in iter_fasta(genome_fasta):
        for i in wanted.get(name, []):
            _c, s, e = regions[i]
            out[i] = seq[int(s):int(e)]
    return out


def run_background(args) -> int:
    """CLI driver: targets -> matched controls -> BED + report."""
    from easyepi.logging_utils import RunRecorder, setup_logging

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger = setup_logging("background", out_path.parent,
                           verbose=getattr(args, "verbose", False),
                           quiet=getattr(args, "quiet", False),
                           log_file=getattr(args, "log_file", None))

    targets = _read_bed(args.peaks)
    exclude = _read_bed(args.exclude) if getattr(args, "exclude", None) else None
    accessible = (read_accessible_bed(args.accessible)
                  if getattr(args, "accessible", None) else None)
    chromosomes = None
    if getattr(args, "chromosomes", None):
        chromosomes = [c.strip() for c in args.chromosomes.split(",") if c.strip()]

    with RunRecorder("background", out_path.parent, logger, render_figures=False) as run:
        logger.info("reading genome windows from %s", args.genome_fasta)
        windows = genome_windows(args.genome_fasta, window=int(args.window),
                                 chromosomes=chromosomes)

        controls, match_summary = match_background(
            targets, windows,
            tolerance=float(args.tolerance),
            min_distance=int(args.min_distance),
            same_chrom=not getattr(args, "any_chrom", False),
            seed=int(args.seed),
            exclude=exclude,
            accessible=accessible,
        )

        # Only measure GC on sequences when something was matched; otherwise the
        # report would compare a real distribution against an empty one.
        target_seqs = control_seqs = None
        if controls:
            target_seqs = _fetch(args.genome_fasta, targets)
            control_seqs = _fetch(args.genome_fasta, controls)

        report = background_report(
            targets, controls,
            target_seqs=target_seqs, control_seqs=control_seqs,
            match_summary=match_summary,
        )
        # The design has to travel with the result: a reader cannot otherwise tell
        # an accessibility-matched background from a composition-matched one, and
        # the two support different conclusions.
        report["background_design"] = match_summary["background_design"]
        if accessible is not None:
            # Per chromosome, because a multi-chromosome run has a different
            # coverage on each and a single number would hide the empty ones.
            chroms = sorted({t[0] for t in targets})
            coverage = {c: round(accessible_fraction(c, targets, accessible), 4)
                        for c in chroms}
            report["accessible"] = {
                "bed": str(args.accessible),
                "intervals_note": "overlapping accessible intervals merged on load",
                # This is the fraction of TARGETS that lie wholly inside the
                # universe, which answers "how much of the signal is in the
                # accessible set", not "did matching work". It can be near zero
                # while every control is accessible, when the universe's intervals
                # are shorter than the peaks -- a control only has to fit inside an
                # interval, a target of the same width does not.
                "target_fraction_wholly_inside_universe_per_chrom": coverage,
                "target_fraction_wholly_inside_universe_overall": round(
                    sum(accessible_fraction(c, targets, accessible) *
                        sum(1 for t in targets if t[0] == c) for c in chroms)
                    / len(targets), 4) if targets else 0.0,
                "what_this_number_means":
                    "fraction of target REGIONS wholly contained in the accessible "
                    "universe; it measures how much signal lies in accessible "
                    "sequence, not whether controls were drawn from it. Controls "
                    "are always drawn from the universe when this field is present.",
            }
        report["params"] = {
            "peaks": str(args.peaks),
            "genome_fasta": str(args.genome_fasta),
            "window": int(args.window),
            "tolerance": float(args.tolerance),
            "min_distance": int(args.min_distance),
            "same_chrom": not getattr(args, "any_chrom", False),
            "seed": int(args.seed),
            "exclude": str(args.exclude) if getattr(args, "exclude", None) else None,
            "accessible": str(args.accessible) if getattr(args, "accessible", None) else None,
        }

        write_bed(controls, out_path)
        report_path = out_path.with_suffix(out_path.suffix + ".json")
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

        # The numbers a reader needs, printed rather than buried in the JSON.
        gc = report["distributions"]["gc"]
        if gc.get("target"):
            print("[easyepi] GC   target mean %.3f / background mean %.3f  "
                  "KS p=%.3g" % (
                      gc["target"]["mean"], gc["background"]["mean"],
                      (gc["ks"]["pvalue"] if gc["ks"]["pvalue"] is not None else float("nan"))))
        else:
            print("[easyepi] GC not measured: %s" % gc.get("status", "unavailable"))
        ln = report["distributions"]["length"]
        print("[easyepi] len  target median %s / background median %s" % (
            ln["target"]["median"], ln["background"]["median"]))
        print("[easyepi] background design: %s" % match_summary["background_design"])
        if accessible is not None:
            print("[easyepi] accessible universe: %d interval(s) on %d chromosome(s); "
                  "every control lies wholly inside it"
                  % (sum(len(v) for v in accessible.values()), len(accessible)))
        print("[easyepi] matched %d/%d targets (%.0f%%)" % (
            match_summary["n_matched"], match_summary["n_targets"],
            100 * match_summary["frac_matched"]))
        if report.get("warning"):
            print("[easyepi] WARNING: %s" % report["warning"])

        gc_block = report["distributions"]["gc"]
        run.set_metrics(
            n_targets=match_summary["n_targets"],
            n_matched=match_summary["n_matched"],
            frac_matched=match_summary["frac_matched"],
            gc_ks_pvalue=(gc_block.get("ks", {}).get("pvalue")
                          if gc_block.get("target") else None),
        )
        logger.info("wrote %s and %s", out_path, report_path)

    print(f"[easyepi] background regions written to {out_path}")
    print(f"[easyepi] match report written to {report_path}")
    return 0
