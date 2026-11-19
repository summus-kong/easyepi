"""Run established external motif discovery or scanning without heuristic fallback."""
import argparse
import json
from pathlib import Path
from easyepi.logging_utils import RunRecorder, setup_logging
from easyepi.motifs.external_tools import run_meme_discovery, run_fimo_scan, run_homer_motifs


def main(argv=None):
    """Command-line entry point for this module."""
    parser=argparse.ArgumentParser(description="Independent MEME discovery, FIMO scanning or HOMER analysis")
    parser.add_argument("--backend",choices=["meme","fimo","homer"],required=True)
    parser.add_argument("--fasta",required=True)
    parser.add_argument("--background",help="Matched background FASTA, required for HOMER")
    parser.add_argument("--database",help="MEME motif database, required for FIMO")
    parser.add_argument("--motifs",type=int,default=5)
    parser.add_argument("--min-width",type=int,default=8)
    parser.add_argument("--max-width",type=int,default=20)
    parser.add_argument("--out",required=True)
    parser.add_argument("--verbose",action="store_true")
    parser.add_argument("--log-file")
    parser.add_argument("--quiet",action="store_true")
    args=parser.parse_args(argv)
    logger=setup_logging("MOTIF-ANALYSIS",args.out,verbose=args.verbose,quiet=args.quiet,log_file=args.log_file)
    with RunRecorder("MOTIF-ANALYSIS",args.out,logger,render_figures=False) as run:
        if args.backend=="meme":
            result=run_meme_discovery(args.fasta,str(Path(args.out)/"meme"),args.motifs,args.min_width,args.max_width)
        elif args.backend=="fimo":
            if not args.database:
                raise ValueError("FIMO requires --database")
            result=run_fimo_scan(args.database,args.fasta,str(Path(args.out)/"fimo"))
        else:
            if not args.background:
                raise ValueError("HOMER FASTA analysis requires explicit --background")
            result=run_homer_motifs(args.fasta,args.background,str(Path(args.out)/"homer"))
        (Path(args.out)/"motif_backend.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
        run.set_metrics(**result)
        logger.info("Completed %s",args.backend)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
