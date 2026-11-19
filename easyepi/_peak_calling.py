"""Assay-aware command construction for the official MACS3 executable."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from easyepi._runtime import add_project_args, prepare_run
from easyepi.logging_utils import run_external


def _signal_input(sample):
    """Return the track to call peaks on: BAM first, else bigWig/bedGraph signal.

    Processed histone experiments frequently ship a normalised bigWig rather than
    raw reads. MACS3 accepts those as BED tracks, so refusing them would make the
    command unusable for exactly the data most people have.
    """
    for attr in ("bam", "bigwig"):
        value = getattr(sample, attr, None)
        if value:
            return str(value)
    raise ValueError(
        f"sample {getattr(sample, 'sample_id', '?')} needs a BAM or a bigWig/bedGraph track")


def _is_paired_end_bam(path) -> bool | None:
    """Detect paired-end alignments from the BAM header/flags.

    Returns True/False, or None when the file is not a readable BAM. Processed
    histone experiments on ENCODE are frequently single-end, and passing BAMPE for
    them makes MACS3 reinterpret reads as improper pairs, so the mode must be
    derived from the data rather than assumed.
    """
    import pysam

    if Path(str(path)).suffix.lower() != ".bam":
        return None
    try:
        with pysam.AlignmentFile(str(path), "rb") as bam:
            # The first primary alignment decides: a BAM is either paired or not,
            # it cannot be mixed. Avoiding a full scan keeps this cheap on 1 GB+ files.
            for rec in bam.fetch(until_eof=True):
                if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                    continue
                return bool(rec.is_paired)
        return False
    except Exception:
        return None


def _macs3_format(path, *, paired_end, qvalue):
    """Map an input path to the MACS3 ``-f`` value it actually is."""
    name = Path(str(path)).name.lower()
    # Compressed tracks carry the meaningful suffix one level in
    # (``x.bedGraph.gz``), so strip a compression extension before deciding.
    for comp in (".gz", ".bz2", ".zst", ".xz"):
        if name.endswith(comp):
            name = name[: -len(comp)]
            break
    suffix = Path(name).suffix
    if suffix in {".bw", ".bigwig"}:
        # A bigWig is not a MACS3 format. Callers must convert first; failing
        # loudly here beats letting MACS3 emit a confusing binary-parse error.
        raise ValueError(
            f"MACS3 cannot read bigWig directly ({path}). Convert it to bedGraph "
            "first (bigWigToBedGraph) and pass the .bedGraph file as --bigwig.")
    if suffix in {".bedgraph", ".bdg", ".bed"}:
        return "BED"
    if suffix in {".bam", ".sam", ".cram"}:
        return "BAMPE" if paired_end else "BAM"
    # Unknown extension: treat as BED rather than silently claiming BAMPE, which
    # produced a confusing MACS3 parse failure on a compressed bedGraph.
    return "BED"


def macs3_command(executable, sample, output, *, genome_size, assay, broad=False,
                  paired_end: bool | None = None, qvalue=0.01, control=None):
    """Build MACS3 arguments for BAM, BED or bedGraph input.

    ``paired_end=None`` (the default) asks the BAM itself; pass True/False to
    override for a deliberately mislabelled file.
    """
    if not 0 < qvalue < 1:
        raise ValueError("qvalue must be between zero and one")
    track = _signal_input(sample)
    if paired_end is None:
        detected = _is_paired_end_bam(track)
        paired_end = True if detected is None else detected
    fmt = _macs3_format(track, paired_end=paired_end, qvalue=qvalue)
    command = [executable, "callpeak", "-t", track, "-f", fmt, "-g", str(genome_size),
               "-n", sample.sample_id, "--outdir", str(output), "-q", str(qvalue)]
    if control:
        command += ["-c", control]
    if broad:
        command += ["--broad"]
    if fmt == "BED":
        # A coverage track has no fragment-length model: tell MACS3 to skip it and
        # use a fixed extension, otherwise peak widths are arbitrary.
        command += ["--nomodel", "--extsize", "147"]
    elif assay == "ATAC" and not paired_end:
        command += ["--nomodel", "--shift", "-100", "--extsize", "200"]
    return command


def main(argv=None, *, assay="ATAC", broad=False):
    """Shared MACS3 peak-calling driver parameterised by assay and peak mode."""
    parser = argparse.ArgumentParser(description=f"{assay} peak calling with MACS3")
    add_project_args(parser)
    parser.add_argument("--genome-size", required=True,
                        help="MACS3 effective genome size: e.g. hs, mm, or a numeric value")
    parser.add_argument("--single-end", action="store_true",
                        help="Force single-end mode. By default the pairing is "
                             "auto-detected from the BAM header.")
    parser.add_argument("--paired-end", action="store_true",
                        help="Force paired-end (BAMPE) mode regardless of the BAM flags.")
    parser.add_argument("--qvalue", type=float, default=0.01)
    parser.add_argument("--macs3", default="macs3", help="MACS3 executable name or path")
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, f"{assay}-PEAKS")
    with recorder as run:
        executable = shutil.which(args.macs3)
        if executable is None:
            raise FileNotFoundError("MACS3 is required for peak calling; install macs3 separately")
        samples = config.samples
        missing = [s.sample_id for s in samples
                   if not (getattr(s, "bam", None) or getattr(s, "bigwig", None))]
        if missing:
            raise ValueError(
                "Every selected sample requires a BAM or a bigWig/bedGraph track; "
                f"missing for: {', '.join(missing)}")
        manifest = []
        for sample in samples:
            sample_dir = Path(args.out) / sample.sample_id
            command = macs3_command(executable, sample, sample_dir,
                                    genome_size=args.genome_size, assay=assay, broad=broad,
                                    paired_end=(False if args.single_end
                                                else True if args.paired_end
                                                else None), qvalue=args.qvalue,
                                    control=sample.control_bam)
            run_external(command, run.logger)
            extension = "broadPeak" if broad else "narrowPeak"
            peak = sample_dir / f"{sample.sample_id}_peaks.{extension}"
            if not peak.is_file():
                raise FileNotFoundError(f"MACS3 did not produce the expected output: {peak}")
            manifest.append({"sample_id": sample.sample_id, "peak": str(peak), "command": command})
        path = Path(args.out) / f"{assay.lower()}_peaks_manifest.json"
        path.write_text(json.dumps({"backend": "MACS3", "samples": manifest}, indent=2), encoding="utf-8")
        run.set_metrics(samples=len(manifest), backend="MACS3", output=str(path))
    return 0
