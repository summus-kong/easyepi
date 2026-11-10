"""Annotate methylation intervals with candidate nearest TSS links."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from easyepi.methylation import read_methylation_table
from easyepi._runtime import add_project_args, prepare_run


def _read_tss(gtf: str) -> dict[str, list[tuple[int, str]]]:
    result: dict[str, list[tuple[int, str]]] = {}
    with open(gtf, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            try:
                start, end = int(p[3]) - 1, int(p[4])
            except ValueError:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', p[8]))
            gene = attrs.get("gene_id", attrs.get("transcript_id", "unknown"))
            result.setdefault(p[0], []).append((end - 1 if p[6] == "-" else start, gene))
    return result


def main(argv=None) -> int:
    """Command-line entry point for this module."""
    parser = argparse.ArgumentParser(description="Easyepi methylation nearest-TSS annotation")
    add_project_args(parser)
    parser.add_argument("--annotation", required=True, help="GTF with gene features and gene_id attributes")
    parser.add_argument("--max-distance", type=int, default=100000)
    args = parser.parse_args(argv)
    config, recorder = prepare_run(args, "METHYLATION-ANNOTATE")
    with recorder as run:
        tss = _read_tss(args.annotation)
        out = Path(args.out)
        n = 0
        for sample in config.samples:
            if not sample.methylation:
                continue
            frame = read_methylation_table(sample.methylation).copy()
            genes, distances = [], []
            for chrom, pos in zip(frame["chrom"], frame["start"]):
                options = tss.get(chrom, [])
                if not options:
                    genes.append(None); distances.append(None); continue
                nearest = min(options, key=lambda x: abs(x[0] - int(pos)))
                dist = int(pos) - nearest[0]
                genes.append(nearest[1] if abs(dist) <= args.max_distance else None)
                distances.append(dist if abs(dist) <= args.max_distance else None)
            frame["nearest_gene_id"] = genes
            frame["signed_distance_to_tss"] = distances
            frame.to_csv(out / f"{sample.sample_id}.methylation_annotated.tsv", sep="\t", index=False)
            n += 1
        run.set_metrics(samples=n, max_distance=args.max_distance)
        run.logger.info("Wrote %d annotated methylation tables", n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

