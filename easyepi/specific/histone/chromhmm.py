"""Prepare inputs and run the official ChromHMM BinarizeBed/LearnModel commands."""
from __future__ import annotations
import argparse
import json
import shutil
from pathlib import Path
from easyepi.logging_utils import RunRecorder, run_external, setup_logging
from easyepi.project import load_project, validate_project, write_inputs_tsv
from easyepi._plots import add_figure_args, selected_formats


def build_commands(args, cellmark, input_dir, binary, model):
    """Assemble the ChromHMM command sequence for this run mode."""
    commands = []
    java = shutil.which("java") or "java"
    prefix = [java, "-mx4000M", "-jar", str(args.chromhmm_jar)]
    if args.run_binarize:
        if not args.genome:
            raise ValueError("BinarizeBed requires --genome chromosome lengths")
        commands.append(prefix + ["BinarizeBed", "-peaks", "-b", str(args.resolution),
                                  str(Path(args.genome).resolve()), str(input_dir), str(cellmark), str(binary)])
    if args.run_learn:
        if not args.genome_build:
            raise ValueError("LearnModel requires --genome-build")
        commands.append(prefix + ["LearnModel", "-b", str(args.resolution), "-p", str(args.threads), "-s", str(args.seed),
                                  str(binary), str(model), str(args.states), args.genome_build])
    return commands


def main(argv=None):
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Prepare or execute official ChromHMM stages")
    parser.add_argument("--project", help="Project TOML/TSV; required for preparation/binarization")
    parser.add_argument("--out", required=True)
    parser.add_argument("--genome-build")
    parser.add_argument("--chromhmm-jar")
    parser.add_argument("--genome", help="Chromosome length TSV, not a FASTA")
    parser.add_argument("--binary-dir", help="Existing binarized directory for standalone learning")
    parser.add_argument("--states", type=int, default=15)
    parser.add_argument("--resolution", type=int, default=200)
    parser.add_argument("--run-binarize", action="store_true")
    parser.add_argument("--run-learn", action="store_true")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--no-check-files", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Write command plan without executing Java")
    add_figure_args(parser)
    args = parser.parse_args(argv)
    if min(args.states,args.resolution,args.threads) < 1:
        raise ValueError("states, resolution and threads must be positive")
    out = Path(args.out).resolve()
    logger = setup_logging("CHROMHMM",out,verbose=args.verbose,quiet=args.quiet,log_file=args.log_file)
    with RunRecorder("CHROMHMM",out,logger,figure_formats=selected_formats(args)) as run:
        cellmark, input_dir = out / "cellmarkfiletable.tsv", out / "input_beds"
        binary = Path(args.binary_dir).resolve() if args.binary_dir else out / "binarized"
        rows = []
        needs_inputs = args.run_binarize or not args.run_learn
        if needs_inputs:
            if not args.project:
                raise ValueError("--project is required for preparation/binarization")
            config = load_project(args.project)
            validate_project(config,check_files=not args.no_check_files)
            args.genome_build = args.genome_build or config.genome_build
            config.write_resolved(out)
            write_inputs_tsv(config,out)
            input_dir.mkdir(exist_ok=True)
            seen = set()
            for index,sample in enumerate(config.samples):
                cell = str(sample.metadata.get("cell") or sample.sample_id)
                mark = sample.metadata.get("mark") or sample.metadata.get("target")
                if not sample.peak or not mark:
                    raise ValueError("Each ChromHMM sample requires peak and explicit mark (or target); use cell to combine marks")
                if any("\t" in value or "\n" in value for value in (cell,str(mark))):
                    raise ValueError("ChromHMM cell and mark labels must not contain tabs or newlines")
                if (cell,str(mark)) in seen:
                    raise ValueError("Duplicate cell/mark pair; merge replicates explicitly before ChromHMM")
                seen.add((cell,str(mark)))
                filename = f"track_{index+1}.bed"
                rows.append({"cell":cell,"mark":str(mark),"bed":filename,"source":sample.peak})
                if not args.dry_run:
                    shutil.copyfile(sample.peak,input_dir/filename)
            with cellmark.open("w",encoding="utf-8") as handle:
                for row in rows:
                    handle.write(f"{row['cell']}\t{row['mark']}\t{row['bed']}\n")
        if args.run_learn and not args.run_binarize and not args.binary_dir:
            raise ValueError("Standalone learning requires --binary-dir")
        if (args.run_learn or args.run_binarize) and not args.chromhmm_jar:
            raise ValueError("Execution requires --chromhmm-jar")
        commands = build_commands(args,cellmark,input_dir,binary,out/"learned_model")
        (out/"chromhmm_manifest.json").write_text(json.dumps({"tracks":rows,"commands":commands,"cellmark":str(cellmark)},indent=2),encoding="utf-8")
        if not args.dry_run:
            for command in commands:
                run_external(command,logger)
        run.set_metrics(tracks=len(rows),commands=commands,dry_run=args.dry_run)
        logger.info("Wrote ChromHMM manifest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
