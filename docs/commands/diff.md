# easyepi diff

```bash
easyepi diff --counts counts.tsv --design design.tsv --method deseq2 --out diff/deseq2.tsv
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--counts` | required |  | count matrix TSV |
| `--design` | required |  | design TSV with sample_id and group/condition |
| `--method` | optional | `mannwhitney` | independent-group test (default: mannwhitney); wilcoxon requires --pair; deseq2 fits a negative-binomial count model and needs raw integer counts plus >=2 replicates per group |
| `--pair` | optional |  | design column containing paired subject IDs |
| `--equal-var` | optional |  | assume equal variance for t-test |
| `--design-factors` | optional |  | comma-separated design columns to add as additive terms (e.g. batch,sex); deseq2 only |
| `--ref-level` | optional |  | reference level for the group factor; defaults to the second level alphabetically; deseq2 only |
| `--lfc-shrink` | optional |  | apply apeGLM-style LFC shrinkage; deseq2 only |
| `--no-independent-filter` | optional |  | test every feature, skipping DESeq2 independent filtering; deseq2 only |
| `--out` | required |  | output TSV for the comparison table |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

## Choosing a method

`mannwhitney` (default), `ttest` and `wilcoxon` compare the values in
the count matrix directly. They are signal comparisons: no count model,
no covariates, no dispersion estimate. Use them to look at data, not to
conclude from it.

`deseq2` fits a negative-binomial GLM through PyDESeq2, which is the
inference path: it models counts, absorbs covariates supplied with
`--design-factors`, and reports shrunken fold changes on request.

### What `--method deseq2` requires

Raw integer counts. The matrix is rejected, with the reason, if it holds
non-integer values (which usually means normalised signal), negative
counts, NaN or inf, fewer than four samples, or fewer than two
replicates per group. A `--design-factors` column with a single level is
rejected too, because it cannot be estimated.

### Direction of `log2FoldChange`

Reported as `log2(group_a / group_b)`: positive means higher in the
first group. PyDESeq2 itself returns the reciprocal for a given
contrast operand order, so the adapter passes the operands reversed
deliberately -- do not "fix" that without re-running the sign test in
`tests/test_deseq2_adapter.py`.

`<out>.deseq2.json` records the model that was fitted: PyDESeq2 version,
design formula, contrast, replicate counts, dispersion fit and how many
features were filtered out before testing.

Generated from the live parser by `scripts/gen_docs.py`.
