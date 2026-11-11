# easyepi viz

```bash
easyepi viz --regions regions.bed --out viz/ --figure-format pdf
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--regions` | required |  |  |
| `--out` | required |  |  |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
