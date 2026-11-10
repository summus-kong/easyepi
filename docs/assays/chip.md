# ChIP-seq

Two shapes of signal, handled by the same commands with different flags.

## Transcription factor (narrow peaks)

```bash
easyepi bam2bw chip_tf.bam -o tracks/ --assay chip --jobs 8
easyepi peaks tracks/chip_tf.bw -o peaks/ --assay chip
easyepi qc --assay chip --bam chip_tf.bam --bigwig tracks/chip_tf.bw \
    --peak peaks/chip_tf.narrowPeak --out qc/chip_tf.json
easyepi tf-annotate --peaks peaks/chip_tf.narrowPeak --genome-fasta genome.fa \
    --motifs jaspar.pfm --fpr 1e-4 --out motifs/
```

## Histone mark (broad peaks and states)

```bash
easyepi bam2bw h3k27me3.bam -o tracks/ --assay histone --jobs 8
easyepi peaks tracks/h3k27me3.bw -o peaks/ --assay histone --broad
easyepi qc --assay histone --bam h3k27me3.bam --bigwig tracks/h3k27me3.bw \
    --peak peaks/h3k27me3.broadPeak --marker H3K27me3 --out qc/h3k27me3.json
```

`--marker` selects the threshold table for that mark: H3K27me3 and H3K9me3 are
expected to be broad and low, H3K27ac narrow and high, so one FRiP cut-off cannot
apply to both.

## Chromatin states across marks

```bash
easyepi histone-states \
    --signal H3K4me3=a.bedGraph,H3K27ac=b.bedGraph,H3K4me1=c.bedGraph \
    --chrom chr1 --start 0 --end 5000000 --resolution 200 --states 6 \
    --n-init 8 --out states/
```

Reads bedGraph rather than BigWig. `states.json` records the restart count and final
log-likelihood so you can judge whether EM converged; a low `n_init` with a poor
likelihood usually means more restarts are needed.

## Assumptions

- input is ChIP over input/IP control; no control track is required
- `--assay chip` implies no shift, 25 bp bins and CPM
- for CUT&Tag and CUT&RUN use `--assay cut&tag` instead: the signal is the cut site
