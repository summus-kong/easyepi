"""Count filtered alignment overlaps or integrate BigWig signal over fixed regions."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from easyepi._runtime import add_project_args, prepare_run
from easyepi.differential_matrix import bigwig_region_sum


def main(argv=None):
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Region matrix: read alignment counts or integrated signal")
    add_project_args(parser)
    parser.add_argument("--regions", required=True, help="Non-overlapping BED3/4 regions")
    parser.add_argument("--source", choices=["bam", "bigwig"], default="bam")
    parser.add_argument("--min-mapq", type=int, default=30)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "REGION-MATRIX")
    with recorder as run:
        regions = []
        with open(args.regions, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError("Invalid BED coordinates")
                regions.append((chrom, start, end, fields[3] if len(fields)>3 else f"region_{len(regions)+1}"))
        if not regions or len({r[3] for r in regions}) != len(regions):
            raise ValueError("Regions must be non-empty with unique IDs")
        ordered = sorted(regions)
        if any(a[0] == b[0] and a[2] > b[1] for a,b in zip(ordered,ordered[1:])):
            raise ValueError("Regions must not overlap")
        result = pd.DataFrame(index=pd.Index([r[3] for r in regions], name="region"))
        for sample in config.samples:
            source = getattr(sample, args.source)
            if not source:
                raise ValueError(f"{sample.sample_id} requires {args.source}")
            if args.source == "bam":
                import pysam
                with pysam.AlignmentFile(source, "rb") as bam:
                    def accept(read):
                        """Validate and normalise one command-line argument value."""
                        return (not read.is_unmapped and not read.is_secondary and not read.is_supplementary
                                and not read.is_duplicate and not read.is_qcfail and read.mapping_quality >= args.min_mapq)
                    result[sample.sample_id] = [bam.count(c,s,e,read_callback=accept) for c,s,e,_ in regions]
            else:
                import pyBigWig
                with pyBigWig.open(source) as bw:
                    result[sample.sample_id] = [bigwig_region_sum(bw,c,s,e) for c,s,e,_ in regions]
        stem = "alignment_counts" if args.source == "bam" else "integrated_signal"
        output = Path(args.out) / f"{stem}.tsv"
        result.to_csv(output, sep="\t")
        pd.DataFrame(regions,columns=["chrom","start","end","region"]).to_csv(Path(args.out)/"regions.tsv",sep="\t",index=False)
        run.set_metrics(regions=len(regions),samples=len(config.samples),source=args.source,
                        units="read alignment overlaps" if args.source == "bam" else "integrated BigWig units",
                        output=str(output))
        run.logger.info("Wrote %s",output.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
