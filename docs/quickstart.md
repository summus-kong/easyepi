# Quick start

Five commands, each independent. Run only the ones you need.

## 1. BAM to BigWig

```bash
easyepi bam2bw sample.bam -o tracks/ --assay atac
```

`--assay` picks the sensible defaults: for ATAC it applies the Tn5 +4/−5 shift, 10 bp
bins and CPM normalisation. Override any of them with `--shift`, `--bin-size`,
`--normalize`.

Multiple BAMs are accepted, and `--jobs` controls worker processes:

```bash
easyepi bam2bw rep1.bam rep2.bam -o tracks/ --assay chip --jobs 8
```

## 2. Peaks

```bash
easyepi peaks tracks/sample.bw -o peaks/ --assay atac
```

Writes `sample.narrowPeak` (or `.broadPeak` with `--broad`).

## 3. Quality control

```bash
easyepi qc --assay atac --bam sample.bam --bigwig tracks/sample.bw \
    --peak peaks/sample.narrowPeak --out qc/sample.json
```

`qc` is the single QC entry point; `--assay` selects which metrics apply —
`atac`, `chip`, `histone`, `mnase` or `methylation`. Metrics whose inputs are absent
are skipped and reported as `not_evaluated`, never as passing.

## 4. Annotation

```bash
easyepi annotate peaks/sample.narrowPeak --gtf genes.gtf -o annot/
```

## 5. Figures

```bash
easyepi gene-panel --bigwig tracks/sample.bw --gtf genes.gtf \
    --gene AGRN,TP73 --out panel/ --figure-format png
```

Figures are only written when `--figure-format` is given.

## Bisulfite data

A bisulfite BAM is the only input; methylation calls are read from its `MM`/`ML` or
Bismark `XM` tags.

```bash
easyepi qc --assay methylation --bam ctrl1.bam --out qc/ctrl1.json
easyepi dmr --group-a ctrl1.bam,ctrl2.bam --group-b trt1.bam,trt2.bam \
    --out dmr/regions.tsv
```

## When something goes wrong

| Message | Meaning |
| --- | --- |
| `input not found: ...` | a path is wrong; the message shows the resolved absolute path and the nearest existing filename |
| `pysam is required for ...` | optional dependency missing; the message prints the exact install command |
| `no methylation calls found in this BAM` | the BAM has neither `MM`/`ML` nor `XM` |
| `a working command looks like: ...` | a required argument is missing; the example is copy-pasteable |

Next: [End-to-end examples](../EXAMPLES.md) for a complete sequence per assay.
