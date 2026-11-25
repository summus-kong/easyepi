# easyepi timecourse

```bash
easyepi timecourse --counts counts.tsv --design design.tsv --out tau.tsv
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--counts` | required |  | count matrix TSV, features x samples, first column is the feature id |
| `--design` | required |  | sample sheet TSV with sample_id and time; time points are never inferred from column order |
| `--out` | required |  | path of the per-feature table; a sibling <out>.json summary and <out>.by_stage.tsv are written next to it |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

Generated from the live parser by `scripts/gen_docs.py`.
