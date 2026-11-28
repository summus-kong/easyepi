"""Easyepi command line interface.

Built with the standard-library ``argparse`` so the CLI needs no extra
dependency beyond the core package. Every analysis module is exposed as a
subcommand; running ``easyepi <module> --help`` shows module options.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from easyepi import __version__
from easyepi.assays import accepted_spellings

#: Assay labels accepted by ``bam2bw``/``qc``. Kept next to the version import so
#: both parsers read one definition; a second copy would drift.
#: One working command per subcommand, shown when its required arguments are
#: missing. argparse prints a usage wall that a first-time user cannot act on;
#: an example they can paste is what actually unblocks them.
MINIMAL_EXAMPLE: dict[str, str] = {
    "bam2bw": "easyepi bam2bw sample.bam -o tracks/ --assay atac",
    "qc": "easyepi qc --assay atac --bam sample.bam --out qc.json",
    "peaks": "easyepi peaks sample.bam -o peaks/ --assay atac",
    "annotate": "easyepi annotate --peaks peaks/sample/sample_peaks.narrowPeak "
                "--gtf genes.gtf --out annot/",
    "tf-annotate": "easyepi tf-annotate --peaks p.narrowPeak --genome-fasta g.fa "
                   "--motifs jaspar.pfm --out motifs/",
    "motifs": "easyepi motifs --peaks p.narrowPeak --genome-fasta g.fa "
              "--motifs jaspar.pfm --out motifs/",
    "footprint": "easyepi footprint --bigwig cut.bw --target sites.bed "
                 "--chrom chr1 --start 0 --end 1000000 --out fp/",
    "nucleosome": "easyepi nucleosome --bigwig mnase.bw "
                  "--chrom chr1 --start 0 --end 1000000 --out nuc/",
    "nucleosome-positioning": "easyepi nucleosome-positioning --bam paired.bam "
                              "--chrom chr1 --start 0 --end 1000000 --out np/",
    "gene-panel": "easyepi gene-panel --bigwig a.bw,b.bw --gtf genes.gtf "
                  "--gene AGRN,TP73 --out panel/ --figure-format pdf",
    "differential-peaks": "easyepi differential-peaks --peaks ctrl1,ctrl2 "
                          "--bam-dir bams --peak-dir peaks --group-a ctrl1 "
                          "--out diff.tsv",
    "dmr": "easyepi dmr --group-a ctrl.bam --group-b trt.bam --out dmr/regions.tsv",
    "histone-states": "easyepi histone-states --signal H3K4me3=a.bg,H3K27ac=b.bg "
                      "--chrom chr1 --start 0 --end 5000000 --out states/",
    "integrate": "easyepi integrate --gene-bed genes.bed --track S=a.bw --out int/",
    "normalize": "easyepi normalize --matrix counts.tsv --factors f.tsv --out norm.tsv",
    "diff": "easyepi diff --counts counts.tsv --design design.tsv --method deseq2 --out diff/deseq2.tsv",
    "concordance": "easyepi concordance --peaks A1=peaks/A1.bed,A2=peaks/A2.bed "
                   "--design design.tsv --out conc/",
    "viz": "easyepi viz --regions regions.bed --out viz/ --figure-format pdf",
    "timecourse": "easyepi timecourse --counts counts.tsv --design design.tsv --out tau.tsv",
    "washu": "easyepi washu --manifest tracks.tsv --out washu/",
    "init": "easyepi init /data/bams --out project.toml",
    "validate-project": "easyepi validate-project --project project.toml --out v/",
    "make-qc-report": "easyepi make-qc-report --project project.toml --out batch/",
    "merge-qc-reports": "easyepi merge-qc-reports \"batch_*/qc_report.json\" --out study/",
}


class EasyepiParser(argparse.ArgumentParser):
    """ArgumentParser that ends a usage error with a runnable example."""

    def error(self, message: str):
        program = self.prog.split()[-1] if self.prog else ""
        example = MINIMAL_EXAMPLE.get(program)
        text = f"{self.prog}: {message}"
        if example:
            text += ("\n  a working command looks like:\n    " + example)
        text += ("\n  run with --help for the full option list")
        self.print_usage(sys.stderr)
        self.exit(2, text + "\n")


#: Every spelling the suite accepts, generated from :mod:`easyepi.assays` so a
#: new alias reaches every command at once. Canonical labels come first, which is
#: what `--help` shows most prominently.
ASSAY_CHOICES_FOR_CLI = accepted_spellings()


def _cmd_merge_qc_reports(args):
    from easyepi.qc_report_merge import collect_reports, merge_reports

    paths, unresolved = collect_reports(args.reports)
    if unresolved:
        import sys as _sys
        print(f"[easyepi] no qc_report.json matched: {', '.join(unresolved)}",
              file=_sys.stderr)
    if not paths:
        import sys as _sys
        print("[easyepi] nothing to merge", file=_sys.stderr)
        return 1
    merged = merge_reports(paths)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(merged, indent=2, ensure_ascii=False), encoding="utf-8")
    summary = merged["verdict_summary"]["samples"]
    print(f"[easyepi] merged {merged['n_reports']} report(s), "
          f"{merged['n_samples']} sample(s) -> {out}")
    if not getattr(args, "no_html", False):
        from easyepi.qc_html import render_merged_html

        # Batch verdicts are a table, not a set of figures: this report carries
        # every sample from every input report so one page answers "which batches
        # and which samples warned".
        html_path = out.with_suffix(".html")
        try:
            render_merged_html(merged, html_path)
        except Exception as exc:  # noqa: BLE001 - report only, never fail the merge
            print(f"[easyepi] HTML report skipped: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
        else:
            print(f"[easyepi] batch QC report written to {html_path}")
    print(f"[easyepi] verdicts: {summary}")
    if merged["problems"]:
        print(f"[easyepi] {len(merged['problems'])} report(s) skipped", file=sys.stderr)
    return 0


def _cmd_diffpeaks(args):
    from easyepi.peaks.diffpeaks import run_differential_peaks
    return run_differential_peaks(args)


def _cmd_histone_states(args):
    from easyepi.histone_states_cli import run_histone_states
    return run_histone_states(args)


def _cmd_gene_panel(args):
    from easyepi.viz.panels import run_gene_panel
    return run_gene_panel(args)


def _cmd_nucleosome_positioning(args):
    from easyepi.nucleosome_positioning import run_positioning
    return run_positioning(args)


def _cmd_tf_annotate(args):
    from easyepi.tf_annotate import run_tf_annotate
    return run_tf_annotate(args)






def _cmd_dmr(args):
    from easyepi.dmr import run_dmr
    return run_dmr(args)


def _cmd_bam2bw(args):
    from easyepi.bam2bw_cli import run_bam2bw
    return run_bam2bw(args)


def _cmd_init(args):
    from easyepi.quickflow import run_init
    return run_init(args)


def _cmd_qc(args):
    from easyepi.qc import run_qc

    # No `_render_cli_figures` here: `run_qc` renders its own figures after
    # writing the JSON (the plot helpers read the written payload), and it knows
    # the report needs a raster format. Calling the shared renderer as well ran a
    # second pass with the bare default -- the pass that wrote PDF-only files and
    # left the HTML report with no images.
    return run_qc(args)


def _cmd_normalize(args):
    from easyepi.normalize import run_normalize
    result = run_normalize(args)
    _render_cli_figures("normalize", args)
    return result


def _cmd_peaks_v2(args):
    from easyepi.peaks_cli import run_peaks_command
    return run_peaks_command(args)


def _cmd_peaks(args):
    from easyepi.peaks import run_peaks
    result = run_peaks(args)
    _render_cli_figures("peaks", args)
    return result


def _cmd_annotate(args):
    from easyepi.annotate import run_annotate
    result = run_annotate(args)
    _render_cli_figures("annotation", args)
    return result


def _cmd_motifs(args):
    from easyepi.motifs import run_motifs
    result = run_motifs(args)
    _render_cli_figures("motifs", args)
    return result


def _cmd_footprint(args):
    """TF footprint: target vs background accessibility around motif sites."""
    from easyepi.footprint import (extract_profiles, fig_footprint,
                                   footprint_profile, footprint_significance,
                                   local_signal, read_site_bed)
    from easyepi.footprint_engines import tobias_port_score, tobias_binary_score
    from easyepi._plots import selected_formats
    import os, json

    t_sites = read_site_bed(args.target)
    b_sites = read_site_bed(args.background) if args.background else None

    chrom = args.chrom
    start, end = args.start, args.end
    signal = local_signal(args.bigwig, chrom, start, end)
    offset = -start  # to convert genome coords -> signal-array index
    N = len(signal)
    target_idx = [pos + offset for chrom_, pos in t_sites if chrom_ == chrom]
    target_idx = [i for i in target_idx if 0 <= i < N]
    bg_idx = None
    if b_sites:
        bg_idx = [pos + offset for chrom_, pos in b_sites if chrom_ == chrom]
        bg_idx = [i for i in bg_idx if 0 <= i < N]
    if not target_idx:
        raise ValueError(
            f"no target site falls inside {chrom}:{start}-{end}; "
            "widen --start/--end or check the chromosome name")

    prof = footprint_profile(signal, target_idx,
                             flank=args.flank, background_sites=bg_idx)

    # The TOBIAS-style score is reported next to the native one so the two are
    # comparable from a single JSON file. `native` computes it too (it is cheap)
    # but only `port`/`tobias` make it the headline result.
    engine = getattr(args, "engine", "native")
    port_score = None
    if engine in ("port", "tobias"):
        try:
            if engine == "port":
                port_score = tobias_port_score(signal, target_idx, args.flank)
            else:
                port_score = tobias_binary_score(args.bigwig, args.target,
                                                 args.chrom, args.start,
                                                 args.end)
        except (ImportError, OSError, ValueError) as exc:
            raise SystemExit(
                f"[easyepi] footprint --engine {engine} unavailable: {exc}\n"
                "          install TOBIAS for --engine tobias, or use "
                "--engine native") from exc

    # Per-locus shoulder-vs-centre contrast, tested against the background set.
    target_mat = extract_profiles(signal, target_idx, 0, args.flank)
    bg_mat = extract_profiles(signal, bg_idx, 0, args.flank) if bg_idx else None
    sig = footprint_significance(target_mat, bg_mat)

    os.makedirs(args.out, exist_ok=True)
    res_path = os.path.join(args.out, "footprint.json")
    with open(res_path, "w", encoding="utf-8") as fh:
        payload = {"score": prof["score"], "n_sites": prof["n_sites"],
                   "chrom": chrom, "start": start, "end": end, "flank": args.flank,
                   "engine": engine,
                   "target_profile": prof["target_profile"].tolist(),
                   "background_profile": prof["background_profile"].tolist(),
                   "significance": {k: v for k, v in sig.items() if k != "target_dips"}}
        if port_score is not None:
            payload["tobias_engine"] = port_score
        json.dump(payload, fh, indent=2)

    p_mw = sig.get("p_mannwhitney")
    p_w = sig.get("p_wilcoxon")
    title = args.title
    if p_mw is not None and p_mw == p_mw:
        title = f"{title}\nMann-Whitney p={p_mw:.2e}"
    formats = selected_formats(args)
    fig_footprint(prof, title=title,
                  save=os.path.join(args.out, "footprint") if formats else None,
                  formats=list(formats))
    print(f"[easyepi] footprint score={prof['score']:.3f} n_sites={prof['n_sites']} "
          f"chrom={chrom}:{start}-{end}")
    if p_mw is not None and p_mw == p_mw:
        print(f"[easyepi] target-vs-background centre contrast: "
              f"mean={sig['target_mean_dip']:.4f} Mann-Whitney p={p_mw:.3e}")
    else:
        print("[easyepi] no background set supplied: no target-vs-background test")
    if p_w is not None and p_w == p_w:
        print(f"[easyepi] centre-vs-shoulder (greater) Wilcoxon p={p_w:.3e} "
              f"(a dip is positive; negative scores mean the centre is a peak)")
    if port_score is not None:
        print(f"[easyepi] TOBIAS-style footprint score={port_score['score']:.4f} "
              f"({port_score['n_sites']} sites, engine={engine})")
        if port_score.get("degenerate"):
            print("[easyepi] note: signal has no negative values, so the TOBIAS "
                  "score equals the flank level by definition")
    print(f"[easyepi] wrote {res_path}")
    return 0


def _cmd_nucleosome(args):
    from easyepi.nucleosome import run_nucleosome
    result = run_nucleosome(args)
    _render_cli_figures("nucleosome", args)
    return result


def _cmd_timecourse(args):
    from easyepi.timecourse_cli import run_timecourse
    return run_timecourse(args)


def _cmd_diff(args):
    from easyepi.differential import run_diff
    result = run_diff(args)
    _render_cli_figures("differential", args)
    return result


def _cmd_concordance(args):
    from easyepi.concordance_cli import run_concordance

    result = run_concordance(args)
    _render_cli_figures("concordance", args)
    return result


def _cmd_viz(args):
    from easyepi.viz import run_viz
    return run_viz(args)


def _render_cli_figures(module, args):
    """Reuse shared plotting functions; record plotting failures separately from analysis failures."""
    try:
        from easyepi._plots import render_module_figures, selected_formats
        render_module_figures(module, args.out, selected_formats(args))
    except Exception as exc:
        print(f"[easyepi] figure rendering skipped: {exc}", file=sys.stderr)


def _cmd_meta(args):
    print(f"Easyepi version: {__version__}")
    print("python:", sys.version.split()[0])


def _cmd_validate_project(args):
    from easyepi.project_commands import validate_project_command
    return validate_project_command(args)


def _cmd_make_qc_report(args):
    from easyepi.project_commands import make_qc_report_command
    return make_qc_report_command(args)


def _cmd_assay_qc(args):
    from easyepi.assay_commands import run_assay_qc
    return run_assay_qc(args, args.assay_name)


def _add_tss_args(parser) -> None:
    """Register the TSS options on a QC parser.

    Included for the same reason as ``--mito-pattern``: the five assay aliases build
    their own parser, so a flag added to ``qc`` alone is absent there and the caller's
    ``getattr(..., None)`` hides it.
    """
    parser.add_argument(
        "--tss", default=None,
        help="TSV of transcription start sites (chrom, tss[, strand]) for TSS "
             "enrichment; requires --bam and an ATAC assay")
    parser.add_argument(
        "--tss-flank", type=int, default=2000,
        help="half-width of the TSS window in bp (default 2000)")


def _add_mito_pattern_arg(parser) -> None:
    """Register ``--mito-pattern`` on a QC parser.

    Shared by ``qc`` and the five assay aliases: the aliases build their own parser,
    so adding the flag in one place only is how it came to be missing on the alias
    path while `getattr(args, "mito_pattern", None)` made it look merely unset.
    """
    parser.add_argument(
        "--mito-pattern", default=None,
        help="regular expression matching mitochondrial contig names, e.g. "
             "'^(chr)?M(T)?$'. Required for mitochondrial_fraction: the contig is "
             "spelled chrM, MT, chrMT or an accession depending on the reference, "
             "and no single spelling is assumed")


def _add_runtime_args(parser):
    """Add shared project, logging and reproducibility parameters."""
    parser.add_argument("--project", default=None, help="Project TOML or TSV sample manifest")
    parser.add_argument("--out", required=True, help="directory to write the run outputs into; created if absent, and a status.json recording the resolved inputs lands at its root")
    parser.add_argument("--genome", default=None, help="reference genome FASTA")
    parser.add_argument("--annotation", default=None, help="gene annotation (GTF/GFF) used to restrict and label the TSS and feature summaries; without it only whole-genome metrics are reported")
    parser.add_argument("--blacklist", default=None,
                        help="ENCODE-style blacklist BED; enables blacklist_fraction QC")
    parser.add_argument("--marker", default=None,
                        help="histone mark label (e.g. H3K27me3); selects marker-aware "
                             "QC thresholds")
    parser.add_argument("--qc-thresholds", default=None, dest="qc_thresholds",
                        help="JSON file overriding QC warn thresholds, e.g. {\"frip\": {\"op\": \"lt\", \"value\": 0.3}}")
    parser.add_argument("--genome-build", default=None, help="genome build, e.g. hg38")
    parser.add_argument("--bam", action="append", default=None, help="BAM; repeat for command-line samples")
    parser.add_argument("--bigwig", action="append", default=None, help="BigWig; repeat for command-line samples")
    parser.add_argument("--peak", action="append", default=None, help="BED/narrowPeak; repeat for command-line samples")
    parser.add_argument("--assay", default=None, help="assay label for direct command-line input")
    parser.add_argument("--no-check-files", action="store_true", help="only validate schema, not file existence")
    parser.add_argument("--verbose", action="store_true", help="write detailed logs/run.log")
    parser.add_argument("--log-file", default=None, help="also write a detailed log to this path; omitted, only status.json records the run")
    parser.add_argument("--quiet", action="store_true", help="disable panel log output")
    parser.add_argument("--threads", type=int, default=1, help="threads recorded for this module")
    parser.add_argument("--seed", type=int, default=1, help="random seed recorded for reproducibility")
    parser.add_argument("--dry-run", action="store_true", help="validate and resolve inputs without analysis")
    parser.add_argument("--print-config", action="store_true", help="print resolved configuration after parsing")
    from easyepi._plots import add_figure_args
    add_figure_args(parser)

def _cmd_integrate(args):
    """Integrate gene-centric BigWig, peak and RNA result tables."""
    from easyepi.integration import run_integration
    from easyepi.viz.panels import expand_csv_values

    tracks = {}
    for item in expand_csv_values(args.track):
        if "=" not in item:
            raise ValueError("--track must be formatted as name=path.bw")
        name, path = item.split("=", 1)
        if not name or not path or name in tracks:
            raise ValueError("--track names must be non-empty and unique")
        tracks[name] = path
    result = run_integration(args.gene_bed, tracks, args.out, clusters=args.clusters, seed=args.seed,
                             peaks=args.peaks, rna_results=args.rna_results,
                             rna_gene_column=args.rna_gene_column, max_link_distance=args.max_link_distance)
    print(f"[easyepi] integrated {result['genes_clustered']} genes across {len(tracks)} track(s)")
    print(f"[easyepi] manifest: {result['outputs']['manifest']}")
    return 0


def _cmd_washu(args):
    from easyepi.washu import build_washu_hub, write_washu_readme
    from easyepi.washu import read_manifest
    samples = read_manifest(args.manifest)
    res = build_washu_hub(samples, args.out, genome=args.genome, base_url=args.base_url)
    readme = write_washu_readme(args.out, samples, genome=args.genome)
    print(f"[easyepi] WashU hub: {res['hub_path']}")
    print(f"[easyepi] browser URL: {res['browser_url']}")
    print(f"[easyepi] README: {readme}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the argparse tree for every analysis module."""
    parser = EasyepiParser(
        prog="easyepi",
        description="Easyepi 鈥?unified bulk epigenomics downstream analysis suite.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="module", required=True,
                                parser_class=EasyepiParser)



    # Guided entry points come first in --help: they are what a first-time
    # user needs, and burying them under 25 subcommands is the main
    # usability complaint.
    # The chain starts here: peaks, QC and tracks all read what this writes.

    p_b2 = sub.add_parser(
        "bam2bw",
        help="convert BAM(s) to BigWig coverage; the first step of any analysis")
    p_b2.add_argument("bam", nargs="+",
                      help="coordinate-sorted BAM file(s)")
    p_b2.add_argument("-o", "--out", required=True,
                      help="directory receiving one signal track per sample plus status.json")
    p_b2.add_argument("--output-file", default=None, dest="output_file",
                      help="explicit BigWig path; only valid with a single BAM")
    p_b2.add_argument("--assay", default=None, choices=list(ASSAY_CHOICES_FOR_CLI),
                      help="sets sensible defaults for shift/bin-size/normalize")
    p_b2.add_argument("--shift", default=None, choices=["none", "tn5"],
                      help="tn5 applies the +4/-5 cut-site shift (ATAC/CUT&Tag)")
    p_b2.add_argument("--bin-size", type=int, default=None, dest="bin_size",
                      help="bin width in bp")
    p_b2.add_argument("--normalize", default=None,
                      choices=["raw", "cpm", "rpkm"],
                      help="value scaling; cpm/rpkm make samples comparable")
    p_b2.add_argument("--chrom", default=None,
                      help="restrict to one chromosome (default: all)")
    p_b2.add_argument("--min-mapq", type=int, default=30, dest="min_mapq",
                      help="minimum mapping quality for a read to be counted (default 30)")
    p_b2.add_argument("--max-depth", type=float, default=None, dest="max_depth",
                      help="cap per-base depth before binning")
    p_b2.add_argument("--name", action="append", default=None,
                      metavar="BAM=NAME",
                      help="override a sample name (default: BAM file stem)")
    p_b2.add_argument("--jobs", "-j", type=int, default=None,
                      help="worker processes for coverage accumulation; "
                           "default: every core the scheduler allows, 1 = "
                           "in-process")
    p_b2.add_argument("--manifest", default=None,
                      help="write a JSON manifest of the produced tracks")
    p_b2.set_defaults(func=_cmd_bam2bw)



    p_i = sub.add_parser(
        "init",
        help="write a project file from a folder, detecting samples and assays")
    p_i.add_argument("source", help="folder holding BAM/BigWig/peak files")
    p_i.add_argument("--out", required=True, help="path of the project TOML to create; refuses to overwrite an existing file")
    p_i.add_argument("--project-id", default=None, dest="project_id",
                     help="identifier written to project.json; defaults to the output directory name")
    p_i.add_argument("--genome-build", default=None, dest="genome_build",
                     help="hg38 / hg19 / mm10 / mm39; inferred from paths if omitted")
    p_i.add_argument("--annotation", default=None, help="GTF/GFF path to record")
    p_i.add_argument("--motif-database", default=None, dest="motif_database",
                     help="motif file to record (JASPAR/HOMER/MEME)")
    p_i.set_defaults(func=_cmd_init)

    p_qc = sub.add_parser("qc", help="QC and statistical summary")
    p_qc.add_argument("--bam", action="append", default=None,
                      help="BAM file(s); repeatable and/or comma-separated")
    p_qc.add_argument("--bigwig", action="append", default=None,
                      help="BigWig signal(s); repeatable and/or comma-separated. "
                           "Enables assay-level metrics")
    p_qc.add_argument("--peak", action="append", default=None,
                      help="BED/narrowPeak(s); repeatable and/or comma-separated. "
                           "Enables FRiP and peak metrics")
    p_qc.add_argument("--input-dir", default=None, dest="input_dir",
                      help="folder scanned for *.bam / *.bw / *.bigWig and peak "
                           "files; each sample becomes <stem>. Saves listing every "
                           "file by hand")
    p_qc.add_argument("--out", required=True, help="path of the JSON report to write; the parent directory must already exist")
    # Assay is optional: when omitted the command stays a generic BAM summary.
    # The choices are listed so `--help` shows them and argparse rejects a typo,
    # which is what used to produce a report with an empty assay_qc.
    # TSS enrichment: opt-in, because it needs a TSS table and a different
    # reference than the other metrics. A gene-level table and a transcript-level
    # table give materially different scores, so the caller supplies the file
    # rather than having one guessed from a GTF.
    _add_mito_pattern_arg(p_qc)
    _add_tss_args(p_qc)
    p_qc.add_argument("--no-html", action="store_true", dest="no_html",
                      help="do not write the self-contained HTML report beside "
                           "the JSON")
    p_qc.add_argument("--assay", default=None, choices=list(ASSAY_CHOICES_FOR_CLI),
                      help="assay label: atac, chip, histone, mnase or "
                           "methylation; CUT&Tag/CUT&RUN and common synonyms are "
                           "accepted. Enables the same verdict block as the "
                           "assay-specific commands")
    p_qc.add_argument("--project", default=None, help="Project TOML or TSV sample manifest")
    p_qc.add_argument("--annotation", default=None, help="GTF/GFF annotation for TSS enrichment")
    p_qc.add_argument("--blacklist", default=None, help="ENCODE-style blacklist BED")
    p_qc.add_argument("--qc-thresholds", default=None, dest="qc_thresholds",
                      help="JSON file overriding QC warn thresholds")
    p_qc.add_argument("--max-reads", type=int, default=200000,
                      help="cap on reads inspected per metric (default: 200000)")
    p_qc.add_argument("--verbose", action="store_true",
                       help="write a detailed log/run.log alongside status.json")
    p_qc.add_argument("--quiet", action="store_true",
                       help="suppress the progress panel; errors are still reported")
    from easyepi._plots import add_figure_args
    add_figure_args(p_qc)
    p_qc.set_defaults(func=_cmd_qc)

    p_dm = sub.add_parser(
        "dmr",
        help="differential methylation regions between two groups of "
             "bisulfite BAMs (per-site Wald test, BH FDR, gap-merged regions)")
    p_dm.add_argument("--group-a", action="append", required=True,
                      dest="group_a", metavar="BAM[,BAM...]",
                      help="repeatable or comma-separated bisulfite BAMs for group A")
    p_dm.add_argument("--group-b", action="append", required=True,
                      dest="group_b", metavar="BAM[,BAM...]",
                      help="repeatable or comma-separated bisulfite BAMs for group B")
    p_dm.add_argument("--label-a", default=None, dest="label_a",
                      help="display name for group A (default A)")
    p_dm.add_argument("--label-b", default=None, dest="label_b",
                      help="display name for group B (default B)")
    p_dm.add_argument("--context", default="CpG",
                      help="methylation context to test: CpG, CHG or CHH")
    p_dm.add_argument("--chrom", default=None,
                      help="restrict to these chromosomes (comma-separated)")
    p_dm.add_argument("--min-coverage", type=int, default=5,
                      dest="min_coverage",
                      help="minimum reads per group at a site (default 5)")
    p_dm.add_argument("--min-mapq", type=int, default=10, dest="min_mapq",
                      help="minimum mapping quality for a read to be counted (default 10, lower than the peak commands because bisulfite aligners assign lower MAPQ)")
    p_dm.add_argument("--min-diff", type=float, default=0.0,
                      dest="min_diff",
                      help="minimum absolute methylation difference (0-1)")
    p_dm.add_argument("--alpha", type=float, default=0.05,
                      help="FDR threshold for significant sites and regions")
    p_dm.add_argument("--max-gap", type=int, default=100, dest="max_gap",
                      help="merge significant sites at most this many bp apart")
    p_dm.add_argument("--pooled", action="store_true",
                      help="force the pooled-binomial site test even when both "
                           "groups have >=2 replicates. It ignores "
                           "between-replicate variance, so p-values are "
                           "anti-conservative with heterogeneous replicates; "
                           "use only to reproduce an older run.")
    p_dm.add_argument("--no-smoothing", action="store_true", dest="no_smoothing",
                      help="use DSS's unsmoothed site test (DMLtest smoothing=FALSE): "
                           "the level and variance come from each site alone, not "
                           "from a 500 bp window. It has its own variance floor "
                           "(1e-05, one-sided), so p-values are not comparable with "
                           "the default smoothed run. Smoothed is the default because "
                           "pooling neighbouring CpGs is what makes single-site "
                           "coverage usable.")
    p_dm.add_argument("--min-sites", type=int, default=1, dest="min_sites",
                      help="minimum significant sites per reported region")
    p_dm.add_argument("--out", required=True,
                      help="path of the per-site DMR table; the run summary JSON and the region table are written beside it")
    p_dm.set_defaults(func=_cmd_dmr)

    p_norm = sub.add_parser("normalize", help="signal / count normalization")
    p_norm.add_argument("--matrix", required=True, help="TSV: region ID then sample columns")
    p_norm.add_argument("--factors", required=True, help="TSV: sample_id and value; measured library sizes or explicit multipliers")
    p_norm.add_argument("--method", choices=["cpm","spikein","rpgc"], default="cpm",
                        help="normalisation: cpm (library size), spikein (exogenous scale factor) or rpgc (reads per genomic content); default cpm")
    p_norm.add_argument("--effective-genome-size", type=int, default=0,
                       help="mappable genome length in bp used by rpgc; 0 means derive it from the BigWig chrom sizes, so set it explicitly when the track covers only part of the genome")
    p_norm.add_argument("--fragment-length", type=float,
                        help="average fragment length in bp; with rpgc this scales by read length / fragment length. Omit to use read length alone")
    p_norm.add_argument("--out", required=True,
                     help="directory for the normalised tables plus status.json")
    add_figure_args(p_norm)
    p_norm.set_defaults(func=_cmd_normalize)

    p_pk = sub.add_parser(
        "peaks",
        help="call peaks with MACS3 (imported as a library)")
    p_pk.add_argument("bam", nargs="+", help="treatment BAM(s), coordinate-sorted")
    p_pk.add_argument("-c", "--control", action="append", default=None,
                      help="control BAM(s); repeatable or comma-separated")
    p_pk.add_argument("-o", "--out", required=True, help="directory to write the run outputs into; created if absent, and a status.json recording the resolved inputs lands at its root")
    from easyepi.peaks_cli import ASSAY_CHOICES as PEAK_ASSAY_CHOICES

    # MACS3 has no bisulfite mode, so `methylation` is not offered here; the list
    # comes from the shared vocabulary restricted to the supported assays.
    p_pk.add_argument("--assay", default=None,
                      choices=list(accepted_spellings(PEAK_ASSAY_CHOICES)),
                      help="sets format and broad/narrow defaults")
    p_pk.add_argument("--format", default=None,
                      help="MACS3 format: BAM, BAMPE, BED, BEDPE (default from --assay)")
    p_pk.add_argument("--genome", default="hs",
                      help="MACS3 genome size: hs, mm, dm or a number (default hs)")
    p_pk.add_argument("--qvalue", type=float, default=0.01,
                      help="q-value cutoff (default 0.01)")
    p_pk.add_argument("--extsize", type=int, default=None,
                       help="extend each read to this length in bp (MACS3 --extsize); omit to use the MACS3 model")
    p_pk.add_argument("--shift-size", type=int, default=None, dest="shift_size",
                      help="shift read starts by this many bp (MACS3 --shift, commonly -100); omit to use the MACS3 model")
    p_pk.add_argument("--broad", action="store_true", default=None,
                      help="broad peaks (default from --assay)")
    p_pk.add_argument("--macs-arg", action="append", default=None,
                      dest="macs_arg", metavar="'--opt VALUE'",
                      help="extra MACS3 options passed through verbatim, e.g. "
                           "--macs-arg=\"--nomodel --extsize 200\"; repeatable")
    p_pk.add_argument("--keep-dup-all", action="store_true", dest="keep_dup_all",
                      help="keep all duplicate reads (MACS3 default is auto)")
    p_pk.set_defaults(func=_cmd_peaks_v2)

    p_an = sub.add_parser("annotate", help="peak annotation and enrichment")
    p_an.add_argument("--peaks", required=True,
                       help="peak file (narrowPeak / broadPeak / BED) whose intervals are annotated")
    p_an.add_argument("--gtf", default=None, help="GTF/GFF gene annotation")
    p_an.add_argument("--backend", default="native", choices=["native", "pybedtools", "pyranges"],
                       help="interval engine; pybedtools and pyranges must be installed, native needs neither (default native)")
    p_an.add_argument("--mode", default="midpoint", choices=["midpoint", "overlap"],
                       help="attribute a peak by its midpoint or by any overlap with the feature (default midpoint)")
    p_an.add_argument("--min-overlap", type=float, default=0.0,
                       help="fraction of the peak that must overlap the feature in overlap mode (default 0.0); ignored in midpoint mode")
    # The promoter window was hardcoded to TSS +/-3000, which is ChIPseeker's
    # default and therefore the right default here; exposing it lets a user match
    # a different ChIPseeker tssRegion, or ask about the core promoter only.
    p_an.add_argument("--promoter-upstream", type=int, default=3000,
                      help="bases upstream of the TSS counted as promoter "
                           "(default 3000, as ChIPseeker)")
    p_an.add_argument("--promoter-downstream", type=int, default=3000,
                      help="bases downstream of the TSS counted as promoter "
                           "(default 3000, as ChIPseeker)")
    # Mirrors getOption("ChIPseeker.downstreamDistance"), default 300.
    p_an.add_argument("--downstream-distance", type=int, default=300,
                      dest="downstream_distance",
                      help="distance to the downstream gene end still reported "
                           "as Downstream (default 300, as ChIPseeker)")
    p_an.add_argument("--no-bins-summary", action="store_false", dest="bins_summary",
                      help="do not write annotation_stat.tsv")
    # Output names are configurable because the feature table is often joined
    # with other results by filename; a fixed name forces renaming after the run.
    p_an.add_argument("--annotated-name", default="annotated.tsv",
                      dest="annotated_name",
                      help="filename for the per-peak table "
                           "(default annotated.tsv)")
    p_an.add_argument("--stat-name", default="annotation_stat.tsv",
                      dest="stat_name",
                      help="filename for the ChIPseeker feature distribution "
                           "(default annotation_stat.tsv)")
    p_an.add_argument("--add-flank-info", action="store_true",
                      dest="add_flank_info",
                      help="append flank_txIds/flank_geneIds/flank_gene_distances "
                           "for features within --flank-distance "
                           "(as ChIPseeker addFlankGeneInfo=TRUE)")
    p_an.add_argument("--flank-distance", type=int, default=5000,
                      help="window in bp used by --add-flank-info "
                           "(default 5000, as ChIPseeker flankDistance)")
    p_an.add_argument("--same-strand", action="store_true",
                      dest="same_strand",
                      help="only consider features on the peak's own strand, "
                           "as annotatePeak(sameStrand=TRUE)")
    p_an.add_argument("--level", default="transcript",
                      choices=("transcript", "gene"),
                      help="annotation unit, as annotatePeak(level=): "
                           "transcript (default) or gene")
    p_an.add_argument("--out", required=True,
                   help="directory for the annotation tables; one TSV per input plus status.json at its root")
    add_figure_args(p_an)
    p_an.set_defaults(func=_cmd_annotate)

    p_mo = sub.add_parser("motifs", help="motif discovery and footprint")
    p_mo.add_argument(
        "--peaks", default=None,
        help="peak file the run is about. Required when motifs are discovered de "
             "novo (the regions come from it); not needed for --pwm scanning or "
             "for --db enrichment, which take their sequences from --fasta and "
             "record --peaks only as provenance")
    p_mo.add_argument("--pwm", default=None, help="JASPAR .pfm PWM file")
    p_mo.add_argument("--fasta", default=None, help="peak sequences FASTA (target)")
    p_mo.add_argument("--background", default=None, help="background sequences FASTA")
    p_mo.add_argument("--db", default=None,
                      help="known motif database (HOMER .motif dir, JASPAR .pfm or "
                           "MEME .meme). Enrichment runs through AME; non-MEME input "
                           "is converted first")
    p_mo.add_argument("--k", type=int, default=8, help="de novo k-mer length (default 8)")
    p_mo.add_argument(
        "--backend", default="auto", choices=["auto", "homer"],
        help="de novo discovery runs HOMER; there is no built-in engine, so a "
             "missing HOMER is reported instead of substituted. Scan a known PWM "
             "with --pwm or test known enrichment with --db to avoid needing "
             "any external program")
    p_mo.add_argument("--homer-dir", default=None, help="path to HOMER share/homer dir (migrated install)")
    p_mo.add_argument(
        "--ame-method", default="fisher",
        choices=["fisher", "3dmhg", "4dmhg", "ranksum", "pearson", "spearman"],
        help="AME statistical test for known-motif enrichment (default fisher)")
    p_mo.add_argument(
        "--ame-scoring", default="avg",
        choices=["avg", "max", "sum", "totalhits"],
        help="AME sequence scoring for known-motif enrichment (default avg)")
    p_mo.add_argument("--out", required=True,
                   help="directory for the motif tables and figures plus status.json")
    add_figure_args(p_mo)
    p_mo.set_defaults(func=_cmd_motifs)

    p_bg = sub.add_parser(
        "background",
        help="GC- and length-matched control regions for an enrichment test")
    p_bg.add_argument("--peaks", required=True,
                      help="target regions (BED/narrowPeak/broadPeak)")
    p_bg.add_argument("--genome-fasta", required=True,
                      help="genome FASTA the regions are taken from")
    p_bg.add_argument("--window", type=int, default=1000,
                      help="width of the genome windows to sample controls from "
                           "(default 1000)")
    p_bg.add_argument("--tolerance", type=float, default=0.05,
                      help="allowed absolute GC-fraction difference (default 0.05)")
    p_bg.add_argument("--min-distance", type=int, default=2000,
                      help="minimum bp between a control and any target or "
                           "excluded region (default 2000)")
    p_bg.add_argument("--exclude", default=None,
                      help="BED of regions that can never be controls, e.g. an "
                           "ENCODE blacklist")
    p_bg.add_argument("--accessible", default=None,
                      help="BED of accessible regions (e.g. a union peak set) to "
                           "restrict controls to. For ATAC data a genome-wide "
                           "background compares peaks against mostly-"
                           "inaccessible sequence; this makes the null an "
                           "accessibility-matched one. Every control is required "
                           "to lie wholly inside an accessible interval, and a "
                           "universe that already contains every target is "
                           "refused as circular")
    p_bg.add_argument("--any-chrom", action="store_true",
                      help="allow controls from a different chromosome; the "
                           "default keeps them on the same chromosome")
    p_bg.add_argument("--chromosomes", default=None,
                      help="comma-separated chromosomes to consider")
    p_bg.add_argument("--seed", type=int, default=0,
                      help="sampling seed, recorded so a background is "
                           "reproducible (default 0)")
    p_bg.add_argument("--out", required=True,
                      help="output BED for the matched control regions")
    add_figure_args(p_bg)

    p_fp = sub.add_parser("footprint", help="TF footprint from accessibility signal")
    p_fp.add_argument("--bigwig", required=True, help="per-base accessibility/cut bigWig")
    p_fp.add_argument("--target", required=True, help="BED of target motif sites (chrom start end)")
    p_fp.add_argument("--background", default=None, help="BED of background/control sites")
    p_fp.add_argument("--chrom", required=True, help="chromosome for the probe region")
    p_fp.add_argument("--start", type=int, required=True, help="probe region start (bp)")
    p_fp.add_argument("--end", type=int, required=True, help="probe region end (bp)")
    p_fp.add_argument("--flank", type=int, default=100, help="flank around each site (bp, default 100)")
    p_fp.add_argument("--title", default="TF footprint",
                       help="figure title (default \"TF footprint\"); set it per comparison so published panels are distinguishable")
    p_fp.add_argument("--out", required=True, help="directory to write the run outputs into; created if absent, and a status.json recording the resolved inputs lands at its root")
    p_fp.add_argument("--engine", choices=("native", "port", "tobias"), default="native",
                      help="footprint scorer: native (default, proportion contrast), port (bundled TOBIAS reimplementation, no TOBIAS needed), or tobias (external binary)")
    add_figure_args(p_fp)
    p_fp.set_defaults(func=_cmd_footprint)

    p_tf = sub.add_parser(
        "tf-annotate",
        help="annotate peaks with transcription factors: scan known PWMs in the "
             "peak sequences with exact per-motif score thresholds")
    p_tf.add_argument("--peaks", required=True,
                      help="BED/narrowPeak/broadPeak to annotate")
    p_tf.add_argument("--genome-fasta", required=True, dest="genome_fasta",
                      help="indexed genome FASTA (.fai required); peak sequences "
                           "are fetched from it")
    p_tf.add_argument("--motifs", required=True,
                      help="motif database: JASPAR .pfm/.jaspar, HOMER .motif or "
                           "MEME .meme; comma-separate to merge several")
    p_tf.add_argument("--fpr", type=float, default=1e-4,
                      help="per-motif per-sequence false-positive rate used for the "
                           "exact threshold (default 1e-4)")
    p_tf.add_argument("--background", default=None,
                      help="background-sequences FASTA; enables the Fisher "
                           "enrichment test and BH correction")
    p_tf.add_argument("--top", type=int, default=5,
                      help="hits reported per peak (default 5)")
    p_tf.add_argument("--out", required=True,
                   help="directory for the TF annotation tables and status.json")
    add_figure_args(p_tf)
    p_tf.set_defaults(func=_cmd_tf_annotate)

    p_nuc = sub.add_parser(
        "nucleosome",
        help="nucleosome occupancy / nucleosome-free fraction from a signal track")
    p_nuc.add_argument("--bigwig", required=True,
                       help="occupancy track (bigWig); for dyad calls from "
                            "paired-end BAM use nucleosome-positioning")
    p_nuc.add_argument("--chrom", required=True,
                       help="contig name as written in the bigWig, e.g. chr1")
    p_nuc.add_argument("--start", type=int, default=0,
                       help="0-based, inclusive start (default 0)")
    p_nuc.add_argument("--end", type=int, required=True,
                       help="0-based, exclusive end of the region to summarise")
    p_nuc.add_argument("--bin-size", type=int, default=10, dest="bin_size",
                       help="bin width in bp (default 10)")
    p_nuc.add_argument("--out", required=True,
                       help="directory for the occupancy metrics, the repeat estimate and status.json")
    add_figure_args(p_nuc)
    p_nuc.set_defaults(func=_cmd_nucleosome)

    p_np = sub.add_parser(
        "nucleosome-positioning",
        help="call nucleosome positions from paired-end fragment dyad density")
    p_np.add_argument("--bam", required=True,
                      help="coordinate-sorted paired-end BAM (MNase-seq or ATAC-seq). "
                           "Single-end input is refused: one cut site gives no fragment "
                           "span, so the dyad cannot be located. Use `easyepi nucleosome` "
                           "for a measured occupancy signal instead")
    p_np.add_argument("--chrom", required=True,
                      help="contig name as written in the BAM header, e.g. chr1")
    p_np.add_argument("--start", type=int, required=True,
                      help="0-based, inclusive start of the region to call")
    p_np.add_argument("--end", type=int, required=True,
                      help="0-based, exclusive end of the region to call")
    p_np.add_argument("--bin-size", type=int, default=10, dest="bin_size",
                      help="dyad density bin width in bp (default 10); spacing and "
                           "fuzziness are reported in bp, not bins")
    p_np.add_argument("--min-length", type=int, default=120, dest="min_length",
                      help="minimum fragment length kept (default 120 bp)")
    p_np.add_argument("--max-length", type=int, default=200, dest="max_length",
                      help="maximum fragment length kept (default 200 bp)")
    p_np.add_argument("--smooth-bins", type=int, default=3, dest="smooth_bins",
                      help="running-mean window in bins applied to the dyad density "
                           "before peak calling (default 3)")
    p_np.add_argument("--mad-k", type=float, default=2.0, dest="mad_k",
                      help="call threshold in scaled MADs above the median")
    p_np.add_argument("--min-spacing", type=int, default=100, dest="min_spacing",
                      help="minimum bp between accepted calls (default 100)")
    p_np.add_argument("--min-height", type=float, default=0.0, dest="min_height",
                      help="absolute floor on the smoothed dyad density a call must "
                           "reach, in density units (default 0 = no absolute floor; "
                           "the relative floor below still applies)")
    p_np.add_argument("--min-relative-height", type=float, default=0.05,
                      dest="min_relative_height",
                      help="fraction of the strongest peak a call must reach "
                           "(default 0.05); guards against noise when MAD is 0")
    p_np.add_argument("--fuzziness-fraction", type=float, default=0.5,
                      dest="fuzziness_fraction",
                      help="fraction of a dyad cluster's signal the reported span must "
                           "cover, measured outward from the peak (default 0.5)")
    p_np.add_argument("--fuzziness-radius", type=int, default=150,
                      dest="fuzziness_radius",
                      help="how far either side of a peak the span search may reach, "
                           "in bp (default 150, about one nucleosome plus linker). "
                           "Values near the neighbour spacing push the span to the cap")
    p_np.add_argument("--max-fragments", type=int, default=2000000,
                      dest="max_fragments",
                      help="cap on fragments read from the region before density is "
                           "built (default 2000000); raise it for deep sequencing")
    p_np.add_argument("--min-mapq", type=int, default=30, dest="min_mapq",
                      help="minimum mapping quality for a fragment to be used "
                           "(default 30)")
    p_np.add_argument("--out", required=True,
                   help="directory for nucleosome_positions.bed, nucleosome_positions.json and status.json")
    add_figure_args(p_np)
    p_np.set_defaults(func=_cmd_nucleosome_positioning)

    p_gp = sub.add_parser(
        "gene-panel",
        help="multi-locus x multi-mark track panel: one column per gene, "
             "one row per BigWig, with the gene model underneath")
    p_gp.add_argument("--bigwig", action="append", required=True,
                      metavar="PATH[,PATH...] or NAME=PATH",
                      help="BigWig path(s), absolute or relative. A single path "
                           "is labelled from its filename stem; NAME=PATH "
                           "overrides it. Comma lists accepted, e.g. a.bw,b.bw")
    p_gp.add_argument("--gtf", default=None,
                      help="GTF/GFF gene annotation used for gene bodies and exons")
    p_gp.add_argument("--bed", default=None,
                      help="alternative to --gtf: chrom/start/end/name[/score/strand]")
    p_gp.add_argument("--gene", action="append", default=None,
                      help="gene symbol(s) to draw; repeatable and/or comma-separated, "
                           "e.g. --gene AGRN,TP73. Omit to take the first --max-genes genes")
    p_gp.add_argument("--max-genes", type=int, default=6, dest="max_genes",
                       help="number of genes drawn when --gene is omitted (default 6)")
    p_gp.add_argument("--flank", type=int, default=2000,
                      help="bp padding each side of the gene body (default 2000)")
    p_gp.add_argument("--bin-size", type=int, default=25, dest="bin_size",
                       help="bin width in bp for the signal profile (default 25)")
    p_gp.add_argument("--independent-ylim", action="store_true",
                      dest="independent_ylim",
                      help="scale each panel separately instead of one y range "
                           "per track row")
    p_gp.add_argument("--title", default=None,
                       help="figure title; omitted when not given")
    p_gp.add_argument("--out", required=True,
                   help="directory for the panel tables and figures plus status.json")
    add_figure_args(p_gp)
    p_gp.set_defaults(func=_cmd_gene_panel)



    p_dp = sub.add_parser(
        "differential-peaks",
        help="DiffBind-style differential peaks in pure Python "
             "(consensus peaks -> TMM -> negative-binomial Wald)")
    p_dp.add_argument("--peaks", action="append", required=True,
                      metavar="SAMPLE[=BAM[:PEAK]]",
                      help="repeatable and/or comma-separated. Three forms: "
                           "SAMPLE=BAM:PEAK (explicit), SAMPLE=BAM (peak file "
                           "found via --peak-dir), or SAMPLE alone (both found "
                           "via --bam-dir and --peak-dir), e.g. --peaks ctrl1,ctrl2")
    p_dp.add_argument("--bam-dir", default=None, dest="bam_dir",
                      help="folder holding <sample>.bam for bare-label --peaks entries")
    p_dp.add_argument("--peak-dir", default=None, dest="peak_dir",
                      help="folder holding <sample>.narrowPeak / .broadPeak / .bed")
    p_dp.add_argument("--group-a", default=None,
                      help="comma-separated sample names forming group A")
    p_dp.add_argument("--design", default=None,
                      help="design TSV (sample_id + group) instead of --group-a")
    p_dp.add_argument("--min-overlap", type=int, default=1, dest="min_overlap",
                      help="bp overlap needed to merge peaks into one region (default 1)")
    p_dp.add_argument("--min-mapq", type=int, default=30, dest="min_mapq",
                       help="minimum mapping quality for a count to enter the model (default 30)")
    p_dp.add_argument("--out", required=True, help="path of the TSV to write; the parent directory must already exist")
    p_dp.set_defaults(func=_cmd_diffpeaks)

    p_hs = sub.add_parser(
        "histone-states",
        help="Chromatin states (ChromHMM) in pure Python: binarise -> "
             "multivariate-Bernoulli HMM -> Viterbi segmentation -> annotation")
    p_hs.add_argument("--signal", action="append", required=True,
                      metavar="MARK=PATH",
                      help="repeatable; bedGraph/bigWig-derived signal per mark")
    p_hs.add_argument("--chrom", required=True,
                       help="contig name as written in the BAM header, e.g. chr1")
    p_hs.add_argument("--start", type=int, default=0,
                       help="0-based, inclusive start (default 0)")
    p_hs.add_argument("--end", type=int, required=True,
                       help="0-based, exclusive end of the region to segment")
    p_hs.add_argument("--resolution", type=int, default=200,
                      help="bin size in bp (ChromHMM default 200)")
    p_hs.add_argument("--states", type=int, default=6,
                       help="number of chromatin states to learn (default 6)")
    p_hs.add_argument("--binarize", choices=["poisson", "peaks"], default="poisson",
                      help="how marks become presence calls: poisson (thresholded counts) or peaks (use the peak file); default poisson")
    p_hs.add_argument("--poisson-cutoff", type=float, nargs="+", default=[0.99],
                      help="one Poisson tail probability per mark, or a single value reused for all marks (default 0.99)",
                      dest="poisson_cutoff")
    p_hs.add_argument("--max-iter", type=int, default=200, dest="max_iter",
                      help="EM iteration cap for the state model (default 200)")
    p_hs.add_argument("--n-init", type=int, default=8, dest="n_init",
                      help="random HMM restarts; the best log-likelihood wins "
                           "(default 8)")
    p_hs.add_argument("--seed", type=int, default=0,
                       help="random seed for state initialisation (default 0, so reruns agree)")
    p_hs.add_argument("--out", required=True,
                   help="directory for the state tables, the segment BED and status.json")
    p_hs.set_defaults(func=_cmd_histone_states)

    p_df = sub.add_parser("diff", help="design-aware differential analysis")
    p_df.add_argument("--counts", required=True, help="count matrix TSV")
    p_df.add_argument("--design", required=True, help="design TSV with sample_id and group/condition")
    p_df.add_argument("--method", default="mannwhitney",
                      choices=["mannwhitney", "ttest", "wilcoxon", "deseq2"],
                      help="independent-group test (default: mannwhitney); "
                           "wilcoxon requires --pair; deseq2 fits a negative-binomial "
                           "count model and needs raw integer counts plus >=2 "
                           "replicates per group")
    p_df.add_argument("--pair", default=None, help="design column containing paired subject IDs")
    p_df.add_argument("--equal-var", action="store_true", help="assume equal variance for t-test")
    # deseq2-only options. They are accepted for every method so a project file can
    # carry them, but the adapter is what enforces their meaning.
    p_df.add_argument("--design-factors", default=None,
                      help="comma-separated design columns to add as additive terms "
                           "(e.g. batch,sex); deseq2 only")
    p_df.add_argument("--ref-level", default=None,
                      help="reference level for the group factor; defaults to the "
                           "second level alphabetically; deseq2 only")
    p_df.add_argument("--lfc-shrink", action="store_true",
                      help="apply apeGLM-style LFC shrinkage; deseq2 only")
    p_df.add_argument("--no-independent-filter", action="store_true",
                      help="test every feature, skipping DESeq2 independent "
                           "filtering; deseq2 only")
    p_df.add_argument("--out", required=True,
                   help="path of the comparison table to write; the parent directory must already exist")
    add_figure_args(p_df)
    p_df.set_defaults(func=_cmd_diff)

    p_tc = sub.add_parser(
        "timecourse",
        help="time-course specificity (TAU) and monotonic trend per feature",
        description=(
            "Score how concentrated each feature is on a single time point "
            "(TAU) and test for a monotonic trend over time (Spearman). TAU is "
            "descriptive: it carries no p-value and says nothing about whether "
            "the shape reproduces between replicates."),
        epilog=(
            "Design TSV columns: sample_id and time, both required. time may be "
            "numeric (0, 2, 6, 24) or ordered labels; at least 3 distinct points "
            "are required. Replicates at the same time are averaged before TAU, "
            "and one column is written per time point."),
    )
    p_tc.add_argument("--counts", required=True,
                      help="count matrix TSV, features x samples, first column is the feature id")
    p_tc.add_argument("--design", required=True,
                      help="sample sheet TSV with sample_id and time; time points "
                           "are never inferred from column order")
    p_tc.add_argument("--out", required=True,
                      help="path of the per-feature table; a sibling <out>.json "
                           "summary and <out>.by_stage.tsv are written next to it")
    p_tc.set_defaults(func=_cmd_timecourse)

    p_cc = sub.add_parser(
        "concordance",
        help="replicate concordance: correlation and region overlap across samples")
    # Two mutually exclusive input layers rather than two commands: peak files
    # give correlation and Jaccard, a count matrix gives correlation only, and
    # the output shape is identical, so one command with one output contract is
    # less to learn than two.
    p_cc.add_argument("--peaks", action="append", default=None, nargs="+",
                      help="SAMPLE=peakfile (or a bare path, taking the sample id "
                           "from the filename); repeatable and/or comma-separated. "
                           "Enables Jaccard as well as correlation")
    p_cc.add_argument("--matrix", default=None,
                      help="features x samples TSV (first column: feature id), the "
                           "same layout as `diff --counts`; correlation only")
    p_cc.add_argument("--design", required=True,
                      help="design TSV with sample_id and group/condition")
    p_cc.add_argument("--method", default="spearman",
                      choices=["spearman", "pearson"],
                      help="correlation method (default: spearman)")
    p_cc.add_argument("--min-samples", type=int, default=2, dest="min_samples",
                      help="refuse to run with fewer samples than this (default 2)")
    p_cc.add_argument("--out", required=True,
                      help="directory for concordance.json, the pairwise TSV and the figures; created if absent")
    p_cc.add_argument("--verbose", action="store_true",
                      help="write detailed logs/run.log")
    p_cc.add_argument("--log-file", default=None,
                      help="also write a detailed log to this path; omitted, only status.json records the run")
    p_cc.add_argument("--quiet", action="store_true",
                      help="disable panel log output")
    add_figure_args(p_cc)
    p_cc.set_defaults(func=_cmd_concordance)

    p_vz = sub.add_parser("viz", help="visualization")
    p_vz.add_argument("--regions", required=True,
                       help="region file (chrom/start/end[/name]) plotted as rows")
    p_vz.add_argument("--out", required=True,
                   help="directory for the figures plus status.json")
    add_figure_args(p_vz)
    p_vz.set_defaults(func=_cmd_viz)

    p_meta = sub.add_parser("meta", help="show package info")
    p_meta.set_defaults(func=_cmd_meta)

    p_project = sub.add_parser("validate-project", help="validate project.toml/manifest and write resolved inputs")
    _add_runtime_args(p_project)
    p_project.set_defaults(func=_cmd_validate_project)

    p_merge = sub.add_parser(
        "merge-qc-reports",
        help="merge several qc_report.json files into one batch-level verdict summary")
    p_merge.add_argument("reports", nargs="+",
                         help="qc_report.json paths, directories, or globs "
                              "(e.g. 'batch_*/qc_report.json')")
    p_merge.add_argument("--out", required=True, help="path of the merged report JSON to write")
    p_merge.add_argument("--no-html", action="store_true", dest="no_html",
                         help="do not write the batch HTML report beside the JSON")
    p_merge.set_defaults(func=_cmd_merge_qc_reports)

    p_qcr = sub.add_parser("make-qc-report", help="summarize input QC only; no motif/differential/RNA conclusions")
    _add_runtime_args(p_qcr)
    # Same knobs as the assay commands so all three QC paths accept identical flags.
    p_qcr.add_argument("--max-reads", type=int, default=200000,
                       help="cap on reads inspected per metric (default: 200000)")
    p_qcr.set_defaults(func=_cmd_make_qc_report)

    # The five assay QC commands were five names for one implementation, so they
    # are registered as aliases of `qc --assay` and hidden from --help. Each one
    # behaves exactly as before, including its exit code and output files; the
    # only change is that `qc --assay atac` is now the documented spelling.
    # The five assay QC commands were five names for one implementation, so they
    # are kept as thin aliases of `qc --assay`. Each still parses the same flags
    # and writes the same outputs; the documented spelling is `qc --assay`.
    for command, assay, help_text in (
        ("atac-qc", "atac", "alias of `qc --assay atac`"),
        ("chip-qc", "chip", "alias of `qc --assay chip`"),
        ("histone-qc", "histone", "alias of `qc --assay histone`"),
        ("methylation-qc", "methylation", "alias of `qc --assay methylation`"),
        ("mnase-qc", "mnase", "alias of `qc --assay mnase`"),
    ):
        p_assay = sub.add_parser(command, help=help_text)
        _add_runtime_args(p_assay)
        p_assay.add_argument(
            "--max-reads", type=int, default=200000,
            help="cap on reads inspected for FRiP/fragment metrics (default: 200000)")
        _add_tss_args(p_assay)
        _add_mito_pattern_arg(p_assay)
        # _add_runtime_args already defines --assay; the alias only fixes its
        # default so the user cannot accidentally run atac-qc with chip thresholds.
        for action in p_assay._actions:
            if action.dest == "assay":
                action.default = assay
        p_assay.set_defaults(func=_cmd_assay_qc, assay_name=assay)

    p_int = sub.add_parser("integrate", help="gene-centric multi-track integration for QCed BigWig, peak and RNA results")
    p_int.add_argument("--gene-bed", required=True, help="BED4+ gene intervals; BED6 strand is used when available")
    p_int.add_argument("--track", action="append", required=True, metavar="NAME=BIGWIG",
                       help="named BigWig track; repeat for each assay/sample summary")
    p_int.add_argument("--out", required=True, help="directory to write the run outputs into; created if absent, and a status.json recording the resolved inputs lands at its root")
    p_int.add_argument("--clusters", type=int, default=4, help="K-means cluster count (default: 4)")
    p_int.add_argument("--seed", type=int, default=1, help="random seed recorded in manifest")
    p_int.add_argument("--peaks", default=None, help="optional contrast-specific BED to link to nearest TSS")
    p_int.add_argument("--max-link-distance", type=int, default=100000, help="maximum nearest-TSS distance in bp")
    p_int.add_argument("--rna-results", default=None, help="optional RNA result TSV; merged by gene ID")
    p_int.add_argument("--rna-gene-column", default="gene_id", help="RNA result gene-ID column")
    add_figure_args(p_int)
    p_int.set_defaults(func=_cmd_integrate)

    p_washu = sub.add_parser("washu", help="export WashU Epigenome Browser tracks & URL")
    p_washu.add_argument("--manifest", required=True, help="TSV: name<tab>bigwig<tab>group[<tab>peak_file]")
    p_washu.add_argument("--out", required=True, help="directory receiving the track hub JSON and its README")
    p_washu.add_argument("--genome", default="hg38", help="WashU genome name (default hg38)")
    p_washu.add_argument("--base-url", default=None, help="http(s) base URL serving the track files")
    p_washu.set_defaults(func=_cmd_washu)
    for command in sub.choices.values():
        if not any(action.dest == "verbose" for action in command._actions):
            command.add_argument("--verbose", action="store_true", help="Write detailed logs/run.log")
            command.add_argument("--log-file", help="Detailed file log path")
            command.add_argument("--quiet", action="store_true", help="Suppress panel logging")
    return parser


def main(argv=None) -> int:
    """Parse arguments, open the run record and dispatch to the selected module."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.module == "background":
            from easyepi.background_cli import run_background
            return run_background(args)
        if args.module in {"qc", "normalize", "annotate", "motifs", "footprint",
                           "nucleosome", "diff", "viz", "integrate", "washu",
                           "tf-annotate", "nucleosome-positioning", "gene-panel"}:
            from easyepi.logging_utils import RunRecorder, setup_logging
            from pathlib import Path
            output = Path(args.out)
            directory = output.parent if output.suffix.lower() in {".tsv", ".csv", ".json", ".txt"} else output
            logger = setup_logging(args.module, directory, verbose=getattr(args,"verbose",False),
                                   quiet=getattr(args,"quiet",False), log_file=getattr(args,"log_file",None))
            with RunRecorder(args.module, directory, logger, render_figures=False) as run:
                result = args.func(args)
                # Modules that build a summary return a dict; scalar modules
                # still return an exit code. Record the dict instead of letting
                # metrics.json end up empty.
                if isinstance(result, dict):
                    run.set_metrics(**result)
                    return 0
                return int(result or 0)
        return int(args.func(args) or 0)
    except ImportError as exc:
        print(f"[easyepi] Missing optional dependency: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover - top-level catch
        print(f"[easyepi] Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
