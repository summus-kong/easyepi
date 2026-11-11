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
    "nucleosome": "easyepi nucleosome --bigwig mnase.bw --out nuc/",
    "nucleosome-positioning": "easyepi nucleosome-positioning --bam paired.bam "
                              "--chrom chr1 --start 0 --end 1000000 --out np/",
    "gene-panel": "easyepi gene-panel --bigwig a.bw,b.bw --gtf genes.gtf "
                  "--gene AGRN,TP73 --out panel/ --figure-format png",
    "differential-peaks": "easyepi differential-peaks --peaks ctrl1,ctrl2 "
                          "--bam-dir bams --peak-dir peaks --group-a ctrl1 "
                          "--out diff.tsv",
    "dmr": "easyepi dmr --group-a ctrl.bam --group-b trt.bam --out dmr/regions.tsv",
    "histone-states": "easyepi histone-states --signal H3K4me3=a.bg,H3K27ac=b.bg "
                      "--chrom chr1 --start 0 --end 5000000 --out states/",
    "integrate": "easyepi integrate --gene-bed genes.bed --track S=a.bw --out int/",
    "normalize": "easyepi normalize --matrix counts.tsv --factors f.tsv --out norm.tsv",
    "diff": "easyepi diff --project project.toml --out diff/",
    "viz": "easyepi viz --regions regions.bed --out viz/ --figure-format png",
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


#: Canonical labels plus the shell-safe spellings. `&` is a shell operator, so
#: `--assay cut&tag` must be quoted; accepting `cuttag` and `cutrun` unquoted is
#: what the assay pages tell a user to type. Both forms resolve to one table key
#: (see ``bam2bw_cli._ASSAY_SPELLINGS``).
ASSAY_CHOICES_FOR_CLI = ("atac", "cut&tag", "cut&run", "cuttag", "cutrun",
                         "chip", "histone", "mnase", "methylation")


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
    result = run_qc(args)
    _render_cli_figures("qc", args)
    return result


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

    # Per-locus shoulder-vs-centre contrast, tested against the background set.
    target_mat = extract_profiles(signal, target_idx, 0, args.flank)
    bg_mat = extract_profiles(signal, bg_idx, 0, args.flank) if bg_idx else None
    sig = footprint_significance(target_mat, bg_mat)

    os.makedirs(args.out, exist_ok=True)
    res_path = os.path.join(args.out, "footprint.json")
    with open(res_path, "w", encoding="utf-8") as fh:
        json.dump({"score": prof["score"], "n_sites": prof["n_sites"],
                   "chrom": chrom, "start": start, "end": end, "flank": args.flank,
                   "target_profile": prof["target_profile"].tolist(),
                   "background_profile": prof["background_profile"].tolist(),
                   "significance": {k: v for k, v in sig.items() if k != "target_dips"}},
                  fh, indent=2)

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
    print(f"[easyepi] wrote {res_path}")
    return 0


def _cmd_nucleosome(args):
    from easyepi.nucleosome import run_nucleosome
    result = run_nucleosome(args)
    _render_cli_figures("nucleosome", args)
    return result


def _cmd_diff(args):
    from easyepi.differential import run_diff
    result = run_diff(args)
    _render_cli_figures("differential", args)
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


def _add_runtime_args(parser):
    """Add shared project, logging and reproducibility parameters."""
    parser.add_argument("--project", default=None, help="Project TOML or TSV sample manifest")
    parser.add_argument("--out", required=True, help="output directory")
    parser.add_argument("--genome", default=None, help="reference genome FASTA")
    parser.add_argument("--annotation", default=None, help="GTF/GFF annotation")
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
    parser.add_argument("--log-file", default=None, help="optional detailed log file path")
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
                      help="output directory for one track per sample")
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
    p_b2.add_argument("--min-mapq", type=int, default=30, dest="min_mapq")
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
    p_i.add_argument("--out", required=True, help="project TOML to write")
    p_i.add_argument("--project-id", default=None, dest="project_id")
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
    p_qc.add_argument("--out", required=True, help="output report path")
    # Assay is optional: when omitted the command stays a generic BAM summary.
    p_qc.add_argument("--assay", default=None,
                      help="assay label (ATAC-seq, ChIP-seq, histone, MNase-seq, "
                           "methylation); enables the same verdict block as the "
                           "assay-specific commands")
    p_qc.add_argument("--project", default=None, help="Project TOML or TSV sample manifest")
    p_qc.add_argument("--annotation", default=None, help="GTF/GFF annotation for TSS enrichment")
    p_qc.add_argument("--blacklist", default=None, help="ENCODE-style blacklist BED")
    p_qc.add_argument("--qc-thresholds", default=None, dest="qc_thresholds",
                      help="JSON file overriding QC warn thresholds")
    p_qc.add_argument("--max-reads", type=int, default=200000,
                      help="cap on reads inspected per metric (default: 200000)")
    p_qc.add_argument("--verbose", action="store_true")
    p_qc.add_argument("--quiet", action="store_true")
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
    p_dm.add_argument("--min-mapq", type=int, default=10, dest="min_mapq")
    p_dm.add_argument("--min-diff", type=float, default=0.0,
                      dest="min_diff",
                      help="minimum absolute methylation difference (0-1)")
    p_dm.add_argument("--alpha", type=float, default=0.05,
                      help="FDR threshold for significant sites and regions")
    p_dm.add_argument("--max-gap", type=int, default=100, dest="max_gap",
                      help="merge significant sites at most this many bp apart")
    p_dm.add_argument("--min-sites", type=int, default=1, dest="min_sites",
                      help="minimum significant sites per reported region")
    p_dm.add_argument("--out", required=True,
                      help="DMR table (a .sites.tsv and .json are written beside it)")
    p_dm.set_defaults(func=_cmd_dmr)

    p_norm = sub.add_parser("normalize", help="signal / count normalization")
    p_norm.add_argument("--matrix", required=True, help="TSV: region ID then sample columns")
    p_norm.add_argument("--factors", required=True, help="TSV: sample_id and value; measured library sizes or explicit multipliers")
    p_norm.add_argument("--method", choices=["cpm","spikein","rpgc"], default="cpm")
    p_norm.add_argument("--effective-genome-size", type=int, default=0)
    p_norm.add_argument("--fragment-length", type=float)
    p_norm.add_argument("--out", required=True)
    add_figure_args(p_norm)
    p_norm.set_defaults(func=_cmd_normalize)

    p_pk = sub.add_parser(
        "peaks",
        help="call peaks with MACS3 (imported as a library)")
    p_pk.add_argument("bam", nargs="+", help="treatment BAM(s), coordinate-sorted")
    p_pk.add_argument("-c", "--control", action="append", default=None,
                      help="control BAM(s); repeatable or comma-separated")
    p_pk.add_argument("-o", "--out", required=True, help="output directory")
    # `methylation` is dropped: ASSAY_CHOICES_FOR_CLI is shared with `bam2bw`, but
    # MACS3 has no bisulfite mode, so offering it here promised a mode that
    # the command rejects. CUT&Tag keeps the bam2bw spelling; peaks_cli
    # normalises it to its own table key.
    p_pk.add_argument("--assay", default=None,
                      choices=[a for a in ASSAY_CHOICES_FOR_CLI
                               if a != "methylation"],
                      help="sets format and broad/narrow defaults")
    p_pk.add_argument("--format", default=None,
                      help="MACS3 format: BAM, BAMPE, BED, BEDPE (default from --assay)")
    p_pk.add_argument("--genome", default="hs",
                      help="MACS3 genome size: hs, mm, dm or a number (default hs)")
    p_pk.add_argument("--qvalue", type=float, default=0.01,
                      help="q-value cutoff (default 0.01)")
    p_pk.add_argument("--extsize", type=int, default=None)
    p_pk.add_argument("--shift-size", type=int, default=None, dest="shift_size")
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
    p_an.add_argument("--peaks", required=True)
    p_an.add_argument("--gtf", default=None, help="GTF/GFF gene annotation")
    p_an.add_argument("--backend", default="native", choices=["native", "pybedtools", "pyranges"])
    p_an.add_argument("--mode", default="midpoint", choices=["midpoint", "overlap"])
    p_an.add_argument("--min-overlap", type=float, default=0.0)
    p_an.add_argument("--out", required=True)
    add_figure_args(p_an)
    p_an.set_defaults(func=_cmd_annotate)

    p_mo = sub.add_parser("motifs", help="motif discovery and footprint")
    p_mo.add_argument("--peaks", required=True)
    p_mo.add_argument("--pwm", default=None, help="JASPAR .pfm PWM file")
    p_mo.add_argument("--fasta", default=None, help="peak sequences FASTA (target)")
    p_mo.add_argument("--background", default=None, help="background sequences FASTA")
    p_mo.add_argument("--db", default=None, help="known motif database (HOMER .motif / JASPAR .pfm / MEME .meme)")
    p_mo.add_argument("--k", type=int, default=8, help="de novo k-mer length (default 8)")
    p_mo.add_argument("--backend", default="auto", choices=["auto", "native", "homer"],
                      help="auto/homer require HOMER; native is a legacy exploratory method")
    p_mo.add_argument("--homer-dir", default=None, help="path to HOMER share/homer dir (migrated install)")
    p_mo.add_argument("--out", required=True)
    add_figure_args(p_mo)
    p_mo.set_defaults(func=_cmd_motifs)

    p_fp = sub.add_parser("footprint", help="TF footprint from accessibility signal")
    p_fp.add_argument("--bigwig", required=True, help="per-base accessibility/cut bigWig")
    p_fp.add_argument("--target", required=True, help="BED of target motif sites (chrom start end)")
    p_fp.add_argument("--background", default=None, help="BED of background/control sites")
    p_fp.add_argument("--chrom", required=True, help="chromosome for the probe region")
    p_fp.add_argument("--start", type=int, required=True, help="probe region start (bp)")
    p_fp.add_argument("--end", type=int, required=True, help="probe region end (bp)")
    p_fp.add_argument("--flank", type=int, default=100, help="flank around each site (bp, default 100)")
    p_fp.add_argument("--title", default="TF footprint", help="plot title")
    p_fp.add_argument("--out", required=True, help="output directory")
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
    p_tf.add_argument("--out", required=True)
    add_figure_args(p_tf)
    p_tf.set_defaults(func=_cmd_tf_annotate)

    p_nuc = sub.add_parser("nucleosome", help="nucleosome occupancy / positioning")
    p_nuc.add_argument("--bam", action="append")
    p_nuc.add_argument("--out", required=True)
    add_figure_args(p_nuc)
    p_nuc.set_defaults(func=_cmd_nucleosome)

    p_np = sub.add_parser(
        "nucleosome-positioning",
        help="call nucleosome positions from paired-end fragment dyad density")
    p_np.add_argument("--bam", required=True,
                      help="coordinate-sorted paired-end BAM (MNase-seq or ATAC-seq)")
    p_np.add_argument("--chrom", required=True)
    p_np.add_argument("--start", type=int, required=True)
    p_np.add_argument("--end", type=int, required=True)
    p_np.add_argument("--bin-size", type=int, default=10, dest="bin_size")
    p_np.add_argument("--min-length", type=int, default=120, dest="min_length",
                      help="minimum fragment length kept (default 120 bp)")
    p_np.add_argument("--max-length", type=int, default=200, dest="max_length",
                      help="maximum fragment length kept (default 200 bp)")
    p_np.add_argument("--smooth-bins", type=int, default=3, dest="smooth_bins")
    p_np.add_argument("--mad-k", type=float, default=2.0, dest="mad_k",
                      help="call threshold in scaled MADs above the median")
    p_np.add_argument("--min-spacing", type=int, default=100, dest="min_spacing",
                      help="minimum bp between accepted calls (default 100)")
    p_np.add_argument("--min-height", type=float, default=0.0, dest="min_height")
    p_np.add_argument("--min-relative-height", type=float, default=0.05,
                      dest="min_relative_height",
                      help="fraction of the strongest peak a call must reach "
                           "(default 0.05); guards against noise when MAD is 0")
    p_np.add_argument("--fuzziness-fraction", type=float, default=0.5,
                      dest="fuzziness_fraction")
    p_np.add_argument("--fuzziness-radius", type=int, default=500,
                      dest="fuzziness_radius")
    p_np.add_argument("--max-fragments", type=int, default=2000000,
                      dest="max_fragments")
    p_np.add_argument("--min-mapq", type=int, default=30, dest="min_mapq")
    p_np.add_argument("--out", required=True)
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
    p_gp.add_argument("--max-genes", type=int, default=6, dest="max_genes")
    p_gp.add_argument("--flank", type=int, default=2000,
                      help="bp padding each side of the gene body (default 2000)")
    p_gp.add_argument("--bin-size", type=int, default=25, dest="bin_size")
    p_gp.add_argument("--independent-ylim", action="store_true",
                      dest="independent_ylim",
                      help="scale each panel separately instead of one y range "
                           "per track row")
    p_gp.add_argument("--title", default=None)
    p_gp.add_argument("--out", required=True)
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
    p_dp.add_argument("--min-mapq", type=int, default=30, dest="min_mapq")
    p_dp.add_argument("--out", required=True, help="output TSV")
    p_dp.set_defaults(func=_cmd_diffpeaks)

    p_hs = sub.add_parser(
        "histone-states",
        help="Chromatin states (ChromHMM) in pure Python: binarise -> "
             "multivariate-Bernoulli HMM -> Viterbi segmentation -> annotation")
    p_hs.add_argument("--signal", action="append", required=True,
                      metavar="MARK=PATH",
                      help="repeatable; bedGraph/bigWig-derived signal per mark")
    p_hs.add_argument("--chrom", required=True)
    p_hs.add_argument("--start", type=int, default=0)
    p_hs.add_argument("--end", type=int, required=True)
    p_hs.add_argument("--resolution", type=int, default=200,
                      help="bin size in bp (ChromHMM default 200)")
    p_hs.add_argument("--states", type=int, default=6)
    p_hs.add_argument("--binarize", choices=["poisson", "peaks"], default="poisson")
    p_hs.add_argument("--poisson-cutoff", type=float, nargs="+", default=[0.99],
                      dest="poisson_cutoff")
    p_hs.add_argument("--max-iter", type=int, default=200, dest="max_iter")
    p_hs.add_argument("--n-init", type=int, default=8, dest="n_init",
                      help="random HMM restarts; the best log-likelihood wins "
                           "(default 8)")
    p_hs.add_argument("--seed", type=int, default=0)
    p_hs.add_argument("--out", required=True)
    p_hs.set_defaults(func=_cmd_histone_states)

    p_df = sub.add_parser("diff", help="design-aware differential analysis")
    p_df.add_argument("--counts", required=True, help="count matrix TSV")
    p_df.add_argument("--design", required=True, help="design TSV with sample_id and group/condition")
    p_df.add_argument("--method", default="mannwhitney", choices=["mannwhitney", "ttest", "wilcoxon"],
                      help="independent-group test (default: mannwhitney); wilcoxon requires --pair")
    p_df.add_argument("--pair", default=None, help="design column containing paired subject IDs")
    p_df.add_argument("--equal-var", action="store_true", help="assume equal variance for t-test")
    p_df.add_argument("--out", required=True)
    add_figure_args(p_df)
    p_df.set_defaults(func=_cmd_diff)

    p_vz = sub.add_parser("viz", help="visualization")
    p_vz.add_argument("--regions", required=True)
    p_vz.add_argument("--out", required=True)
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
    p_merge.add_argument("--out", required=True, help="output merged JSON path")
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
    p_int.add_argument("--out", required=True, help="output directory")
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
    p_washu.add_argument("--out", required=True, help="output dir for hub JSON + README")
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
        if args.module in {"qc", "normalize", "annotate", "motifs", "footprint",
                           "nucleosome", "diff", "viz", "integrate", "washu",
                           "tf-annotate", "nucleosome-positioning", "gene-panel"}:
            from easyepi.logging_utils import RunRecorder, setup_logging
            from pathlib import Path
            output = Path(args.out)
            directory = output.parent if output.suffix.lower() in {".tsv", ".csv", ".json", ".txt"} else output
            logger = setup_logging(args.module, directory, verbose=getattr(args,"verbose",False),
                                   quiet=getattr(args,"quiet",False), log_file=getattr(args,"log_file",None))
            with RunRecorder(args.module, directory, logger, render_figures=False):
                return int(args.func(args) or 0)
        return int(args.func(args) or 0)
    except ImportError as exc:
        print(f"[easyepi] Missing optional dependency: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover - top-level catch
        print(f"[easyepi] Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
