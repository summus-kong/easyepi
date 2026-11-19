# ATAC-seq

ATAC-seq measures accessibility through Tn5 insertion. The signal of interest is the
insertion site, not the read body, so `bam2bw --assay atac` applies the standard
+4/−5 shift before accumulating.

## Commands

```bash
easyepi bam2bw atac.bam -o tracks/ --assay atac --jobs 8
easyepi peaks atac.bam -o peaks/ --assay atac
easyepi qc --assay atac --bam atac.bam --bigwig tracks/atac.bw \
    --peak peaks/atac/atac_peaks.narrowPeak --out qc/atac.json
easyepi footprint --bigwig tracks/atac.bw --target peaks/atac/atac_peaks.narrowPeak \
    --chrom chr1 --start 0 --end 5000000 --flank 100 --out fp/
```

## Metrics and what they mean

| Metric | Good range | If out of range |
| --- | --- | --- |
| FRiP (genome-wide) | 0.2–0.5 | low: over-digestion, or peaks called over a sub-region |
| Fragment periodicity | visible 10 bp ladder | lost: nuclei prep or digestion problem |
| Nucleosome ratio | mononucleosome dominant | dinucleosome-heavy: under-digestion |
| TSS enrichment | > 5 | low: dead cells or failed transposition |

`qc` reports FRiP as `interpretable: false` when the peak set spans only a small
fraction of the mapped genome — a regional numerator over a global denominator is not
a FRiP, and the message says so.

## Assumptions

- the BAM is coordinate-sorted and indexed
- reads are paired-end for fragment-length metrics; single-end still yields FRiP and
  signal enrichment
- `--assay atac` implies Tn5 shift, 10 bp bins and CPM; override with `--shift none`
  if the BAM was already shifted
