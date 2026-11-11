# easyepi peaks

```bash
easyepi peaks tracks/sample.bw -o peaks/ --assay atac
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--bam` | optional |  | BAM input(s); repeatable or comma-separated |
| `--bigwig` | optional |  | BigWig input(s); repeatable or comma-separated |
| `--assay` | optional |  | selects calling defaults (shift, broad output) |
| `--broad` | optional |  | write a broadPeak file (default comes from --assay) |
| `--backend` | optional | `auto` | caller to use; auto prefers MACS3 when installed |
| `--chrom` | optional |  | restrict to these chromosomes (comma-separated) |
| `--bin-size` | optional | `200` |  |
| `--threshold` | optional | `5.0` | fold-change threshold for the built-in caller |
| `--max-qvalue` | optional | `0.01` |  |
| `--out` | required |  |  |
| `--figure-format` | optional |  | Optional export format; repeat to export multiple formats |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
