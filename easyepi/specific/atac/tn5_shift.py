"""Independent tn5 shift analysis entry."""
from __future__ import annotations

import argparse
from pathlib import Path

from easyepi._runtime import add_project_args, prepare_run


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi ATAC Tn5 shift")
    add_project_args(parser)
    parser.add_argument("--min-mapq", type=int, default=30)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "ATAC-TN5-SHIFT")
    with recorder as run:
        import pysam
        n = 0
        for sample in config.samples:
            if not sample.bam:
                continue
            output = Path(args.out) / f"{sample.sample_id}.tn5.bed"
            with pysam.AlignmentFile(sample.bam, "rb") as bam, output.open("w", encoding="utf-8") as out:
                for read in bam.fetch(until_eof=True):
                    if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                            or read.is_qcfail or read.mapping_quality < args.min_mapq or read.reference_start < 0):
                        continue
                    pos = (read.reference_end - 5) if read.is_reverse else (read.reference_start + 4)
                    chrom = bam.get_reference_name(read.reference_id)
                    if 0 <= pos < bam.get_reference_length(chrom):
                        out.write(f"{chrom}\t{pos}\t{pos + 1}\t{sample.sample_id}_{n}\t1\t{ '-' if read.is_reverse else '+'}\n")
                        n += 1
            run.logger.info("Wrote %s", output.name)
        run.set_metrics(cut_sites=n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
