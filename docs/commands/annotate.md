# easyepi annotate

```bash
easyepi annotate --peaks peaks/sample/sample_peaks.narrowPeak --gtf genes.gtf --out annot/
```

## Options

| flag | required | default | description |
| --- | --- | --- | --- |
| `--peaks` | required |  | peak file (narrowPeak / broadPeak / BED) whose intervals are annotated |
| `--gtf` | optional |  | GTF/GFF gene annotation |
| `--backend` | optional | `native` | interval engine; pybedtools and pyranges must be installed, native needs neither (default native) |
| `--mode` | optional | `midpoint` | attribute a peak by its midpoint or by any overlap with the feature (default midpoint) |
| `--min-overlap` | optional |  | fraction of the peak that must overlap the feature in overlap mode (default 0.0); ignored in midpoint mode |
| `--promoter-upstream` | optional | `3000` | bases upstream of the TSS counted as promoter (default 3000, as ChIPseeker) |
| `--promoter-downstream` | optional | `3000` | bases downstream of the TSS counted as promoter (default 3000, as ChIPseeker) |
| `--downstream-distance` | optional | `300` | distance to the downstream gene end still reported as Downstream (default 300, as ChIPseeker) |
| `--no-bins-summary` | optional | `True` | do not write annotation_stat.tsv |
| `--annotated-name` | optional | `annotated.tsv` | filename for the per-peak table (default annotated.tsv) |
| `--stat-name` | optional | `annotation_stat.tsv` | filename for the ChIPseeker feature distribution (default annotation_stat.tsv) |
| `--add-flank-info` | optional |  | append flank_txIds/flank_geneIds/flank_gene_distances for features within --flank-distance (as ChIPseeker addFlankGeneInfo=TRUE) |
| `--flank-distance` | optional | `5000` | window in bp used by --add-flank-info (default 5000, as ChIPseeker flankDistance) |
| `--same-strand` | optional |  | only consider features on the peak's own strand, as annotatePeak(sameStrand=TRUE) |
| `--level` | optional | `transcript` | annotation unit, as annotatePeak(level=): transcript (default) or gene |
| `--out` | required |  | directory for the annotation tables; one TSV per input plus status.json at its root |
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

`annotation_stat.tsv` is written beside the per-peak table with the
feature distribution: percentages per category, following ChIPseeker's
`getAnnoStat`. Promoter sub-categories stay separate rows, exons and
introns are split into `1st`/`Other`, and `Downstream` is its own row;
the category count follows `--promoter-upstream`/`--promoter-downstream`
and the data itself. `--no-bins-summary` suppresses the file.

Both table names are configurable via `--annotated-name` and `--stat-name`.

Generated from the live parser by `scripts/gen_docs.py`.
