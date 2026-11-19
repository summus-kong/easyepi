# End-to-end examples

Every example is a copy-pasteable sequence for one assay, from an aligned BAM to the
final table or figure. They assume:

- `easyepi` is on `PATH` (otherwise prefix with `python -m easyepi.cli`)
- a genome build and a matching annotation are available
- figures are requested explicitly with `--figure-format`

Nothing here requires a project file; `--project` is an optional convenience.

---

## ATAC-seq

BAM -> BigWig (Tn5 cut sites) -> peaks -> annotation -> footprint -> figure.

```bash
# 1. BAM -> BigWig, Tn5-shifted, CPM-normalised, 8 workers
easyepi bam2bw atac_rep1.bam -o tracks/ --assay atac --jobs 8

# 2. Peaks from the BAM (MACS3 reads alignments, not tracks)
easyepi peaks atac_rep1.bam -o peaks/ --assay atac

# 3. Library QC: fragment periodicity, FRiP, signal enrichment
easyepi qc --assay atac \
    --bam atac_rep1.bam \
    --bigwig tracks/atac_rep1.bw \
    --peak peaks/atac_rep1/atac_rep1_peaks.narrowPeak \
    --out qc/atac_rep1.json

# 4. Peaks to genes
easyepi annotate --peaks peaks/atac_rep1/atac_rep1_peaks.narrowPeak \
    --gtf genes.gtf --out annot/

# 5. Accessibility profile around the peaks (nucleosome-depleted region)
easyepi footprint --bigwig tracks/atac_rep1.bw \
    --target peaks/atac_rep1/atac_rep1_peaks.narrowPeak \
    --chrom chr1 --start 0 --end 5000000 \
    --flank 100 --out fp/ --figure-format pdf
```

**What to expect.** `qc` FRiP for a good ATAC library is 0.2-0.5 genome-wide. The
footprint step reports a Mann-Whitney p for the centre-vs-background contrast; a
value near 0.05 means the shifted background is too similar to the target, usually
because the peak set and the signal cover different regions.

**Common failure.** `footprint` with `--start 0 --end 5000000` but peaks spread over
all chromosomes finds no target sites in the window and stops with a message naming
the region. Widen the window or match it to the peaks.

---

## ChIP-seq (transcription factor, narrow peaks)

```bash
# 1. BAM -> BigWig, read-body depth (no Tn5 shift), CPM
easyepi bam2bw chip_tf.bam -o tracks/ --assay chip --jobs 8

# 2. Narrow peaks with MACS3, via the imported library
easyepi peaks chip_tf.bam -o peaks/ --assay chip

# 3. QC: FRiP and enrichment at the peaks
easyepi qc --assay chip \
    --bam chip_tf.bam --bigwig tracks/chip_tf.bw \
    --peak peaks/chip_tf/chip_tf_peaks.narrowPeak --out qc/chip_tf.json

# 4. Which TFs can bind these peaks
easyepi tf-annotate \
    --peaks peaks/chip_tf/chip_tf_peaks.narrowPeak \
    --genome-fasta genome.fa \
    --motifs jaspar_subset.pfm \
    --fpr 1e-4 --top 5 --out motifs/ --figure-format pdf
```

**What to expect.** `tf-annotate` needs a `.fai` beside the FASTA. A 4 bp motif
cannot reach `--fpr 1e-4` and is reported in `skipped_motifs` rather than as
"no site"; raise `--fpr` if a short motif matters.

---

## ChIP-seq (histone marks) and chromatin states

```bash
# Per-mark peaks: broad marks need --broad
easyepi bam2bw h3k27me3.bam -o tracks/ --assay histone --jobs 8
easyepi peaks h3k27me3.bam -o peaks/ --assay histone --broad

# Chromatin states across several marks
easyepi bam2bw h3k4me3.bam h3k27ac.bam h3k4me1.bam -o tracks/ --assay histone --jobs 8

easyepi histone-states \
    --signal H3K4me3=tracks/h3k4me3.bg,H3K27ac=tracks/h3k27ac.bg,H3K4me1=tracks/h3k4me1.bg \
    --chrom chr1 --start 0 --end 5000000 --resolution 200 --states 6 \
    --n-init 8 --out states/

# Multi-locus figure: one column per gene, one row per mark
easyepi gene-panel \
    --bigwig tracks/h3k4me3.bw,tracks/h3k27ac.bw,tracks/h3k4me1.bw \
    --gtf genes.gtf --gene AGRN,PLCH2,SKI,TP73 \
    --flank 3000 --bin-size 50 --out panel/ --figure-format pdf
```

**Note on `histone-states`.** It reads bedGraph, not BigWig; convert first
(`bigWigToBedGraph`) or point `--signal` at bedGraph you already have. The command
records the restart count and final log-likelihood in `states.json`, which is how
you tell whether the EM converged.

---

## CUT&Tag / CUT&RUN

Same path as ATAC: the signal is the Tn5 cut site.

```bash
easyepi bam2bw cuttag.bam -o tracks/ --assay cut\&tag --jobs 8
easyepi peaks cuttag.bam -o peaks/ --assay cut\&tag
easyepi qc --assay chip --bam cuttag.bam \
    --bigwig tracks/cuttag.bw --peak peaks/cuttag/cuttag_peaks.narrowPeak \
    --out qc/cuttag.json
```

`--assay cut&tag` sets shift=tn5 and bin=10; the quoting is needed because `&` is a
shell operator. Use `--assay cuttag` or quote the value.

---

## MNase-seq

Occupancy and periodicity from a signal track; dyad calls need paired-end fragments.

```bash
# Occupancy and periodicity from the signal
easyepi bam2bw mnase.bam -o tracks/ --assay mnase --jobs 8
easyepi nucleosome --bam mnase.bam --out nuc/ --figure-format pdf

# Dyad calling from paired-end fragments (single-end is refused)
easyepi nucleosome-positioning \
    --bam mnase_paired.bam \
    --chrom chr1 --start 0 --end 2000000 \
    --bin-size 10 --min-length 120 --max-length 200 \
    --out positioning/ --figure-format pdf
```

**Why single-end is refused.** A single-end read reports one cut site and no
fragment span, so the dyad cannot be located. The error says so and points to
`easyepi nucleosome`, which works from an occupancy signal instead.

---

## WGBS / bisulfite

Methylation is read from the alignments' `MM`/`ML` or Bismark `XM` tags, so the BAM
is the only input.

```bash
# 1. Per-sample methylation QC (beta distribution, coverage)
easyepi qc --assay methylation --bam ctrl1.bam --out qc/ctrl1.json

# 2. Two-group differential methylation
easyepi dmr \
    --group-a ctrl1.bam,ctrl2.bam \
    --group-b trt1.bam,trt2.bam \
    --label-a control --label-b treated \
    --min-coverage 10 --min-diff 0.2 --alpha 0.05 \
    --out dmr/regions.tsv
```

Outputs `regions.tsv`, `regions.sites.tsv`, `regions.json`.

**What to expect.** With `--min-coverage 10 --min-diff 0.2` the region count is in
the hundreds to low thousands for a real contrast. A `regions.tsv` with zero rows
usually means the coverage filter removed everything; `regions.json` reports
`dropped_low_coverage` so you can tell that apart from "no difference".

**Limitation.** The per-site test is binomial, not beta-binomial, so p-values are
anti-conservative at high coverage. For publication-grade DMR calling, DSS, bsseq or
metilene remain the reference; this is a transparent fast screen.

---

## Two-condition comparison (ATAC or ChIP)

```bash
# Per-sample tracks and peaks
for s in ctrl1 ctrl2 trt1 trt2; do
    easyepi bam2bw ${s}.bam -o tracks/ --assay atac --jobs 4
    easyepi peaks ${s}.bam -o peaks/ --assay atac
done

# Differential accessibility with a negative-binomial model
easyepi differential-peaks \
    --peaks ctrl1,ctrl2,trt1,trt2 \
    --bam-dir bams --peak-dir peaks \
    --group-a ctrl1,ctrl2 \
    --out diff/differential.tsv

# Plot the result
easyepi viz --regions diff/differential.tsv --out figures/ --figure-format pdf
```

`--peaks` also accepts the explicit `SAMPLE=BAM:PEAK` form and `SAMPLE=BAM` (peak
file found via `--peak-dir`). Samples that cannot be resolved raise rather than
silently changing the replicate count.

---

## Region-set similarity and intersections

```python
from easyepi.viz.genomeviews import jaccard_bp, plot_peak_overlap
from easyepi.viz.panels import plot_upset
from easyepi.peaks.diffpeaks import read_peaks

sets = {name: read_peaks(path) for name, path in
        [("H3K4me3", "peaks/h3k4me3/h3k4me3_peaks.narrowPeak"),
         ("H3K27ac", "peaks/h3k27ac/h3k27ac_peaks.narrowPeak"),
         ("H3K27me3", "peaks/h3k27me3/h3k27me3_peaks.broadPeak")]}

# base-pair Jaccard is width-robust; interval counts overstate broad domains
print(jaccard_bp(sets["H3K4me3"], sets["H3K27ac"]))

plot_peak_overlap(list(sets.values()), save="figures/overlap", formats=["png"])
plot_upset(sets, save="figures/upset", formats=["png"])
```

---

## A project file, when you have many samples

```bash
# Detect samples from a folder and write a project file
easyepi init /data/bams --out project.toml --genome-build hg38 --annotation genes.gtf

# Edit the generated file, then every command can take --project instead of paths
easyepi qc --assay atac --project project.toml --out qc/
easyepi validate-project --project project.toml --out validation/
```

`init` records the assay it inferred for each sample and prints it, so a wrong guess
is visible before any analysis runs.

---

## Common errors and what they mean

| Message | Cause | Fix |
| --- | --- | --- |
| `input not found: ...` | path typo, or a relative path from the wrong directory | the message prints the resolved absolute path and the nearest existing filename |
| `pysam is required for reading BAM alignments` | optional dependency absent | the message prints the exact `pip install` for the running interpreter |
| `no methylation calls found in this BAM` | the BAM has neither `MM`/`ML` nor `XM` | confirm it is the aligner output, not a derivative that stripped the tag |
| `nucleosome positioning needs paired-end fragments` | single-end BAM | use paired-end data, or `easyepi nucleosome` for occupancy |
| `no PWM score reaches the requested fpr` | motif too short | raise `--fpr`, or accept that a 4 bp motif cannot reach `1e-4` |
| `peaks span only X% of mapped reference` | peak set and BAM cover different regions | FRiP is not comparable; re-call peaks over the same universe |
