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
| `--promoter-upstream` | optional | `3000` | bases upstream of the TSS counted as promoter (default 3000, as ChIPseeker) |
| `--promoter-downstream` | optional | `3000` | bases downstream of the TSS counted as promoter (default 3000, as ChIPseeker) |
| `--no-bins-summary` | optional | `True` | do not write annotation_bins.tsv |
| `--out` | required |  | output directory for annotation tables |
| `--figure-format` | optional |  | Export format for figures; repeat to export several. Default: pdf (vector) |
| `--no-figures` | optional |  | Do not write figures; only tables and JSON are produced |
| `--verbose` | optional |  | Write detailed logs/run.log |
| `--log-file` | optional |  | Detailed file log path |
| `--quiet` | optional |  | Suppress panel logging |

## Genomic-region bins

Each peak is assigned to the first bin, in this order, whose interval
contains the peak midpoint:

    Promoter > 5'UTR > 3'UTR > Exon > Intron > Intergenic

`Promoter` is the window around the TSS set by `--promoter-upstream` and
`--promoter-downstream` (3000 each by default, the ChIPseeker `tssRegion`
default). A UTR is the exonic, non-coding piece of an exon: between the
transcript end and the first coding base (5'), or between the last coding
base and the polyadenylation site (3'). An exon further along the transcript
than the last CDS block is reported as `Exon`, not UTR.

The UTRs precede `Exon` on purpose: an exon interval always contains its own
UTRs, so with `Exon` first the UTR bins could never be filled.

`annotation_bins.tsv` is written beside the per-peak table with the bin
counts and fractions, which are the numbers usually quoted.
`--no-bins-summary` suppresses it.

Generated from the live parser by `scripts/gen_docs.py`.
