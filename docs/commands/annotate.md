# easyepi annotate

```bash
easyepi annotate --peaks peaks/sample/sample_peaks.narrowPeak --gtf genes.gtf --out annot/
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--peaks` | required |  |  |
| `--gtf` | optional |  | GTF/GFF gene annotation |
| `--backend` | optional | `native` |  |
| `--mode` | optional | `midpoint` |  |
| `--min-overlap` | optional |  |  |
| `--out` | required |  |  |
| `--figure-format` | optional |  | Optional export format; repeat to export multiple formats. Defaults to pdf (vector) when a figure is produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
