# easyepi concordance

```bash
easyepi concordance --peaks A1=peaks/A1.bed,A2=peaks/A2.bed --design design.tsv --out conc/
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--peaks` | optional |  | SAMPLE=peakfile (or a bare path, taking the sample id from the filename); repeatable and/or comma-separated. Enables Jaccard as well as correlation |
| `--matrix` | optional |  | features x samples TSV (first column: feature id), the same layout as `diff --counts`; correlation only |
| `--design` | required |  | design TSV with sample_id and group/condition |
| `--method` | optional | `spearman` | correlation method (default: spearman) |
| `--min-samples` | optional | `2` | refuse to run with fewer samples than this (default 2) |
| `--out` | required |  | output directory for concordance.json, concordance_pairs.tsv and figures |
| `--verbose` | optional |  | write detailed logs/run.log |
| `--log-file` | optional |  | optional detailed log file path |
| `--quiet` | optional |  | disable panel log output |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |

## Two inputs, two resolutions

`--peaks` takes one peak file per sample and gives correlation **and**
region-set Jaccard. `--matrix` takes a features x samples table (the same
layout as `diff --counts`) and gives correlation only, because a count
matrix carries no genomic coordinates. Passing both is refused.

`--design` is required and must have `sample_id` plus `group` or
`condition`. Groups are never inferred from column order or filenames.

## What is deliberately not produced

No `qc_verdict` and no pass/fail threshold. `grouping_supported` is a
strict comparison of the mean within-group correlation against the mean
between-group correlation, and it is a *group-level* statement: measured on
synthetic designs with real group structure, swapping one label in six
still reports `True`. It does not identify a mislabelled sample, and the
command prints that caveat on every run rather than only on odd ones.

A value that could not be computed is `null` in the JSON and an empty
field in the TSV -- never `0`. Two empty peak sets share no regions, but
that is not a Jaccard of zero; it is no Jaccard at all.

## Outputs

* `concordance.json` -- summary, per-pair records, the correlation matrix
  (design-ordered) and any caveats that applied.
* `concordance_pairs.tsv` -- the same pairs long-form, for pandas/R/awk.
* `figures/correlation_heatmap.pdf` -- ordered by group with separators
  drawn where the group changes. Not clustered: clustering reorders samples
  to maximise block structure, which destroys the evidence the figure is
  for.

Generated from the live parser by `scripts/gen_docs.py`.
