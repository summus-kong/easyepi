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
| `--min-mapq` | optional | `10` |  |
| `--min-diff` | optional |  | minimum absolute methylation difference (0-1) |
| `--alpha` | optional | `0.05` | FDR threshold for significant sites and regions |
| `--max-gap` | optional | `100` | merge significant sites at most this many bp apart |
| `--min-sites` | optional | `1` | minimum significant sites per reported region |
| `--out` | required |  | DMR table (a .sites.tsv and .json are written beside it) |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
