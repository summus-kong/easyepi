# easyepi dmr

```bash
easyepi dmr --group-a ctrl.bam --group-b trt.bam --out dmr/regions.tsv
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--group-a` | required |  | repeatable or comma-separated bisulfite BAMs for group A |
| `--group-b` | required |  | repeatable or comma-separated bisulfite BAMs for group B |
| `--label-a` | optional |  | display name for group A (default A) |
| `--label-b` | optional |  | display name for group B (default B) |
| `--context` | optional | `CpG` | methylation context to test: CpG, CHG or CHH |
| `--chrom` | optional |  | restrict to these chromosomes (comma-separated) |
| `--min-coverage` | optional | `5` | minimum reads per group at a site (default 5) |
| `--min-mapq` | optional | `10` | minimum mapping quality for a read to be counted (default 10, lower than the peak commands because bisulfite aligners assign lower MAPQ) |
| `--min-diff` | optional |  | minimum absolute methylation difference (0-1) |
| `--alpha` | optional | `0.05` | FDR threshold for significant sites and regions |
| `--max-gap` | optional | `100` | merge significant sites at most this many bp apart |
| `--pooled` | optional |  | force the pooled-binomial site test even when both groups have >=2 replicates. It ignores between-replicate variance, so p-values are anti-conservative with heterogeneous replicates; use only to reproduce an older run. |
| `--no-smoothing` | optional |  | use DSS's unsmoothed site test (DMLtest smoothing=FALSE): the level and variance come from each site alone, not from a 500 bp window. It has its own variance floor (1e-05, one-sided), so p-values are not comparable with the default smoothed run. Smoothed is the default because pooling neighbouring CpGs is what makes single-site coverage usable. |
| `--min-sites` | optional | `1` | minimum significant sites per reported region |
| `--out` | required |  | output DMR table file path (.sites.tsv and .json are written beside it) |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
