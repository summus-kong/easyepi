# Complete usage guide

## 1. Choose a step

Every entry is an independent Python module:

```bash
python -m easyepi.universal.MODULE --help
python -m easyepi.specific.ASSAY.MODULE --help
```

Do not run all modules automatically. FASTQ alignment, read trimming and reference downloading are outside the present scope. Start with aligned BAM, normalized BigWig, peak BED, methylation tables or previously computed RNA results. Choose one assay per project manifest for assay-specific steps.

## 2. Linux installation and dependencies

Linux is the supported runtime. Use Python 3.11 in an isolated venv or the provided conda environment; no native Windows implementation is planned.

```bash
conda env create -f environment.yml
conda activate easyepi
python -m pip install -e ".[dev]"
```

| Feature | Python extra | External software |
| --- | --- | --- |
| Tables, statistical summaries, plotting | Base: NumPy, SciPy, pandas, Matplotlib | None |
| BAM QC, overlap counts, Tn5 cuts | `bam`: pysam | None beyond pysam/HTSlib |
| BigWig reading/integration | `bigwig`: pyBigWig | None |
| Known motif enrichment | `motif`: Biopython | None |
| Multi-track clustering | `stats`: scikit-learn, statsmodels | None |
| TF/regulator ULM estimates | `activity`: decoupler >=1.9.2,<2 | None; numba/llvmlite are transitive dependencies |
| Peak calling | `bam` for adjacent QC/count steps | MACS3 |
| Chromatin states | None additional | Java and official ChromHMM.jar |
| De novo motifs / FIMO | None additional | MEME Suite |
| Optional HOMER | None additional | HOMER and Perl |
| Optional TOBIAS footprint | Separate tool installation | TOBIAS, with bias-corrected input |

`python -m pip install -e ".[all]"` installs the Python extras; it does not install external executables, genomes or motif databases. The base package can be imported without pysam, pyBigWig, decoupler or scikit-learn. Missing optional dependencies fail explicitly when the corresponding analysis is requested.

## 3. Unified sample input

TOML is recommended for stable sample metadata and reference paths. CLI flags express per-run settings. Paths in TOML/TSV resolve relative to the manifest. Use ASCII sample IDs starting with a letter/digit, followed by letters, digits, underscores, periods or hyphens.

```toml
[project]
project_id = "atac-time-course"
genome_build = "hg38"

[[samples]]
sample_id = "T0_R1"
assay = "ATAC-seq"
bam = "bams/T0_R1.bam"
bigwig = "tracks/T0_R1.bw"
peak = "peaks/T0_R1.narrowPeak"
time = 0

[[samples]]
sample_id = "T24_R1"
assay = "ATAC-seq"
bam = "bams/T24_R1.bam"
bigwig = "tracks/T24_R1.bw"
peak = "peaks/T24_R1.narrowPeak"
time = 24
```

Remove unavailable file fields rather than inventing paths. Optional fields include `group`, `batch`, `subject_id`, `control_bam`, `spikein_bam` and `methylation`. Genome and annotation references belong in `[references]`; do not declare references you do not use. ChromHMM additionally requires explicit `mark` and an optional `cell` identifier.

A tab-separated manifest is also accepted:

```text
sample_id	assay	bam	group
A1	ATAC-seq	bams/A1.bam	A
A2	ATAC-seq	bams/A2.bam	A
B1	ATAC-seq	bams/B1.bam	B
B2	ATAC-seq	bams/B2.bam	B
```

The displayed separators represent actual tabs. No experiment/control labels are required for QC, counts, annotation, motifs or descriptive profiles. `time` is metadata: no longitudinal or batch-adjusted statistical model is implemented yet. A two-group comparison cannot replace a repeated-measures time-course model.

## 4. Validation, logs and plotting

```bash
python -m easyepi.universal.validate_project --project project.toml --out results/validate --verbose
python -m easyepi.specific.atac.qc --project project.toml --out results/qc --log-file results/qc/details.log
python -m easyepi.universal.make_qc_report --project project.toml --out results/qc_report
```

Reports summarize input QC only. They do not run downstream analyses or produce biological conclusions. Available metrics are input-level summaries; ATAC TSS enrichment, FRiP, library complexity and replicate concordance are not fully implemented here.

Runs write `status.json`, `metrics.json` and `command.txt`; project-aware steps also write `resolved_config.json` and `inputs.tsv`. Status records contain Python and selected package versions. `--verbose` adds `logs/run.log`; `--log-file` chooses its path; `--quiet` suppresses panel logs. External stdout/stderr is captured in detailed logs. Command provenance uses process arguments; when invoking `main(argv)` programmatically, also record the passed argument list in your notebook or workflow.

`--no-check-files` checks manifest schema without guaranteeing that an analysis can run. Shared project modules accept `--dry-run` to validate inputs and write metadata without analysis. ChromHMM dry runs additionally produce a Java command plan. Dry runs create an output directory and metadata files.

No image is exported by default. Where supported:

```bash
python -m easyepi.specific.atac.qc --project project.toml --out results/qc --figure-format pdf --figure-format svg
```

Use `easyepi.plot` functions to obtain Figures and customize labels, colors and export formats. Automatic generic plots are summaries, not validated assay-specific diagnostic figures. Standard volcano, heatmap, bar, profile and distribution functions are shared; specialized motif logos and profile plots remain in their domain modules.

## 5. ATAC and ChIP/CUT&Tag/CUT&RUN

Coordinate-sorted indexed BAM files must share the same genome assembly and contig naming. Resolve multimapping, mitochondrial reads, blacklists, library complexity and replicate QC before confirmatory analyses.

```bash
python -m easyepi.specific.atac.tn5_shift --project atac.toml --out results/cuts --min-mapq 30
python -m easyepi.specific.atac.peaks --project atac.toml --genome-size hs --out results/atac_peaks
python -m easyepi.specific.chip.narrow_peaks --project chip.toml --genome-size hs --out results/chip_peaks
python -m easyepi.specific.chip.broad_peaks --project histone.toml --genome-size hs --out results/broad_peaks
```

Peak calling requires MACS3 on PATH (or `--macs3 /path/to/macs3`). The default `BAMPE` mode uses observed fragment boundaries and does not apply an extra Tn5 shift. Use `--single-end` only for single-end BAM: ATAC then uses MACS3 `--nomodel --shift -100 --extsize 200`, while ChIP uses MACS3 model defaults. The standalone Tn5 step produces +4/-5 cut sites for signal inspection, not input required by the BAMPE caller. No automatic CUT&RUN/CUT&Tag-specific SEACR selection is implemented.

`control_bam` is forwarded per sample when supplied. A missing control is allowed by MACS3, but suitability depends on the assay and experimental design. `--genome-size` is MACS3 effective genome size (such as `hs`, `mm` or a number); it is not a FASTA file or HOMER genome identifier.

Update peak paths in a subsequent manifest with the produced `*_peaks.narrowPeak`/`broadPeak` paths. No hidden pipeline mutates the input manifest.

```bash
python -m easyepi.specific.atac.consensus --project peaks.toml --min-samples 2 --out results/consensus
python -m easyepi.specific.atac.count --project atac.toml --regions results/consensus/atac_consensus.bed --source bam --out results/counts
python -m easyepi.specific.chip.count --project chip.toml --regions consensus.bed --source bigwig --out results/signal
```

Consensus selects genomic positions covered by peaks from at least N distinct samples and merges adjacent retained segments. Its BED score reports minimum support over each merged interval. It is an interval-support summary, not IDR.

BAM matrices count primary, nonsupplementary, nonduplicate, QC-passing alignments with MAPQ >=30 by default. **These are read-alignment overlap counts, not paired fragments**: both mates can count, and an alignment may overlap multiple separate regions. Output is `alignment_counts.tsv`, with region ID followed by sample columns; `regions.tsv` retains coordinates. Use `--min-mapq` to set a different threshold. BigWig mode writes `integrated_signal.tsv` with signal units, not read counts. Unknown BigWig contigs raise errors; uncovered bases on existing contigs contribute zero.

## 6. Normalize measured matrices

A factor TSV has `sample_id` and `value`. For CPM, value is measured library size; for spike-in, value is an explicit positive multiplier. Sample IDs are matched, not inferred from row order.

```bash
python -m easyepi.universal.normalize_signal --matrix alignment_counts.tsv --factors library_sizes.tsv --method cpm --out results/normalized
```

RPGC is supported only for per-base coverage input and requires measured fragment counts in the factors file, `--effective-genome-size` and `--fragment-length`. Do not apply this coverage formula to arbitrary peak counts. No missing library size or spike-in quantity is replaced with a fabricated default.

## 7. Differential comparisons

The existing SciPy tests are exploratory, replicated two-group comparisons. They do not estimate count dispersion, fit batch/time covariates or substitute for PyDESeq2/edgeR/DESeq2. Use biologically independent replicates (minimum two per group; larger studies are preferable).

```text
sample_id	group	subject
A1	A	1
A2	A	2
B1	B	1
B2	B	2
```

```bash
easyepi diff --counts normalized_matrix.tsv --design design.tsv --method mannwhitney --out results/differential.tsv
easyepi diff --counts normalized_matrix.tsv --design design.tsv --method wilcoxon --pair subject --out results/paired.tsv
python -m easyepi.specific.atac.differential --project replicated_signals.toml --method ttest --out results/signal_comparison
```

The general CLI aligns design IDs to matrix columns and reorders paired samples by subject. Welch's t-test is the default t-test variant unless equal variance is explicitly requested. Missing observations remain missing; paired tests remove whole pairs. BH correction excludes untested features. For count/signal comparisons `log2FC = log2((mean(group_a)+1)/(mean(group_b)+1))`; group order follows matrix/manifest first appearance. The assay-specific signal comparison records direction in `comparison.json`.

Assay-specific differential scripts integrate normalized BigWig signal over the union of peaks. All regions are retained unless `--max-regions` is explicitly provided. Explicit subsampling is seeded and must be reported. Paired testing is available through `easyepi diff` with an explicit design; assay-specific signal modules accept independent tests only.

## 8. Peak and methylation annotation

```bash
python -m easyepi.universal.peak_annotation --peaks peaks.bed --gtf genes.gtf --out results/annotation.tsv
python -m easyepi.specific.methylation.annotate --project methylation.toml --annotation genes.gtf --out results/methylation_annotation
```

Nearest-TSS annotations are candidate gene links, not experimentally confirmed regulatory interactions. The methylation annotation step requires GTF gene features with `gene_id`, and supports `--max-distance`. Motif-to-TF database identities are separate from genomic peak annotation; motif enrichment does not annotate a peak's nearest gene automatically.

## 3b. Scanning a folder instead of listing files

```bash
easyepi qc --input-dir /data/tracks --assay histone --out results/qc.json
```

`--input-dir` groups `*.bam`, `*.bw` / `*.bigWig`, peak files
(`*.narrowPeak`/`*.broadPeak`/`*.bed`/`*.bed.gz`) and methylation tables
(`*.cov`/`*.bedMethyl`) into samples by filename stem. `a.bed.gz` becomes sample
`a`, and `*.bai`/`*.csi` indexes are ignored.

Tracks are attached to samples **by name, not by list position**, so a sample
without a BigWig no longer shifts every later pairing. Giving one sample two files
of the same kind is an error rather than an arbitrary pick. Explicit `--bam` /
`--bigwig` / `--peak` values override the scan for the same sample and kind, and
all three accept comma-separated lists.

If nothing is recognised, the command fails; it does not write an empty report.

### Differential peaks with short sample labels

```bash
easyepi differential-peaks \
  --peaks ctrl1,ctrl2,treat1,treat2 \
  --bam-dir bams --peak-dir peaks \
  --group-a ctrl1,ctrl2 --out results/diff.tsv
```

Three `--peaks` forms are accepted: `SAMPLE=BAM:PEAK` (explicit), `SAMPLE=BAM`
(peak file located via `--peak-dir`), and bare `SAMPLE` (both located via
`--bam-dir` and `--peak-dir`). Samples that cannot be resolved, or that match more
than one file, raise instead of quietly changing the replicate count.

## 9. Motif scan, enrichment and discovery

```bash
python -m easyepi.universal.motif_scan --fasta peaks.fa --pwm TF.pfm --out results/scan
python -m easyepi.universal.motif_enrichment --target peaks.fa --background matched.fa --pwm TF1.pfm --pwm TF2.pfm --fpr 0.0001 --out results/enrichment
python -m easyepi.universal.motif_analysis --backend meme --fasta peaks.fa --out results/denovo
python -m easyepi.universal.motif_analysis --backend fimo --fasta peaks.fa --database motifs.meme --out results/fimo
```

Scanning accepts a single JASPAR motif (optional header and A/C/G/T labels) or a plain four-row PFM. It returns the best score, strand and zero-based forward-sequence start per sequence. The low-level scanner does not report all sites; use FIMO for complete, calibrated site calls and genomic coordinate conversion.

Known enrichment uses Biopython score distributions under a uniform background to choose a per-strand score threshold, one presence/absence call per sequence, Fisher's exact test and BH across supplied motifs. The `--fpr` value is a model threshold, not a study-wide binding error rate. Target/background sequences must be distinct and matched for GC, length, mappability and tested genomic universe. Matching is user-supplied; automatic background generation is not implemented. Closely related TFs often share motifs, so motif identity does not establish the active family member.

MEME/FIMO require MEME Suite. Optional HOMER uses `--backend homer --background matched.fa`, requires its installation and never silently falls back to a native heuristic. Historical low-level native discovery helpers remain explicitly exploratory and are excluded from recommended workflows.

## 9b. TF annotation of peaks (peak -> TF)

```bash
easyepi tf-annotate --peaks H3K4me3_peaks.narrowPeak \
    --genome-fasta /ref/GRCh38/genome.fa \
    --motifs jaspar_subset.pfm --fpr 1e-4 --out results/tf
```

Peak sequences are fetched from the indexed FASTA (a `.fai` must exist), each
supplied PWM is scanned on both strands, and a motif counts as present when its
best site exceeds a threshold taken from the **exact PWM score distribution**
under a uniform background, binned on the score axis. This is the FIMO
construction (Grant et al. 2011); `--fpr` is a per-motif, per-sequence
false-positive rate, not a genome-wide FDR.

Outputs:

| File | Contents |
| --- | --- |
| `tf_per_peak.tsv` | best `--top` motifs per peak with genomic site coordinates |
| `tf_summary.tsv` | per-TF peak coverage, threshold, best score |
| `tf_annotate.json` | thresholds, skipped motifs, run parameters |

A motif too short to reach `--fpr` (a 4 bp motif tops out at 0.25^4 = 0.0039) is
listed in `skipped_motifs` instead of being reported as having no site. Supply
`--background matched.fa` to add a Fisher exact enrichment test with BH
correction.

**Interpretation limits.** A motif hit means the peak contains a sequence
compatible with that PWM. It is not evidence that the named protein occupies the
site in the assayed cells: paralogues share motifs, chromatin context is ignored,
and no footprint is computed here.

## 11b. Nucleosome positioning (paired-end)

```bash
easyepi nucleosome-positioning --bam mnase.paired.bam \
    --chrom chr1 --start 0 --end 200000 \
    --bin-size 10 --min-length 120 --max-length 200 \
    --out results/positioning --figure-format pdf
```

Positions are called from the **fragment dyad density**: each fragment contributes
its midpoint (the midpoint of the two MNase cut sites), the midpoints are binned,
and local maxima of the smoothed curve above `median + mad_k * MAD` are accepted,
spaced by at least `--min-spacing` bp. A relative floor (`--min-relative-height`,
default 0.05 of the strongest peak) prevents single-fragment noise from being called
when MAD is 0, which is normal for a sparse dyad track.

| Output | Contents |
| --- | --- |
| `nucleosome_positions.bed` | one interval per call, ±one bin around the dyad |
| `nucleosome_positions.json` | calls, fuzziness, spacing, periodicity, parameters |

Single-end BAMs are refused: a single-end read gives one cut site and no fragment
span, so the dyad position is not identifiable. The error says so and points to
`easyepi nucleosome`, which works from a measured occupancy signal instead.

This is not a template-based or kinetic model of digestion, it does not resolve
overlapping nucleosomes, and it must not be applied to ATAC-seq as if ATAC reported
protection.

## 12b. Multi-locus panels, UpSet and GC-matched background

### Several genes x several marks on one figure

```bash
easyepi gene-panel \
  --bigwig H3K4me3=/data/bigwigs/H3K4me3.bw,H3K27ac=/data/bigwigs/H3K27ac.bw,H3K4me1=/data/bigwigs/H3K4me1.bw,H3K27me3=/data/bigwigs/H3K27me3.bw \
  --gtf gencode.annotation.gtf.gz \
  --gene AGRN,PLCH2,SKI,TP73 \
  --flank 3000 --bin-size 50 \
  --out results/panel --figure-format pdf
```

Both `--bigwig` and `--gene` accept comma-separated values, and the flags may still be repeated.

`--bigwig-dir` lets `--bigwig` take bare labels: `H3K4me3` resolves to `<dir>/H3K4me3.bw` or `<dir>/H3K4me3.bigWig`. A label that matches nothing, or matches more than one file, is an error rather than a silently missing track. `NAME=PATH` remains available for files outside one folder.

One column per gene, one row per mark, gene model (exons plus a strand arrow)
underneath, and genomic coordinates in kb on the bottom row. `--gtf` supplies gene
bodies and exons; `--bed` accepts `chrom/start/end/name` instead. A requested gene
symbol that is absent from the annotation raises, so a silently short figure cannot
pass as a complete one.

y-limits are computed per panel from a robust quantile, so a saturated locus does
not flatten the others. `--independent-ylim` scales each panel on its own.

### UpSet and sample correlation

```python
from easyepi.viz.panels import plot_upset, plot_sample_correlation
plot_upset({"H3K4me3": peaks_a, "H3K27ac": peaks_b, "H3K4me1": peaks_c},
           save="results/upset", formats=["png"])
plot_sample_correlation(counts_frame, method="pearson",
                        save="results/corr", formats=["png"])
```

UpSet counts **elementary union segments**, not whole intervals, so a partially
overlapping pair is split into its unique and shared parts. `plot_sample_correlation`
drops constant features and constant samples, then optionally draws a dendrogram.

### GC-matched background

```python
from easyepi.gc_background import genome_windows, match_background, write_bed
windows = genome_windows("genome.fa", window=1000, chromosomes=["chr1"])
controls, summary = match_background(peaks, windows, tolerance=0.05,
                                     min_distance=2000, seed=0)
write_bed(controls, "results/background.bed")
print(summary["n_unmatched"], "targets had no GC-matched control")
```

A uniform shift does not preserve GC, and GC is the strongest single predictor of
accessibility, so a shifted control understates the difference at GC-rich sites.
This matcher requires the same GC fraction within `tolerance`, a comparable window
width, an interval distance of at least `min_distance` from every target and from
every `exclude` region, and uses each window at most once. Unmatched targets are
listed in the summary instead of being dropped, because silently losing the hard
cases biases the comparison. It matches sequence composition only; mappability,
repeat content and chromatin state are not matched.

## 10. TF activity estimates

```bash
python -m easyepi.universal.tf_activity --matrix normalized_accessibility.tsv --network tf_peak_network.tsv --min-targets 5 --out results/activity
python -m easyepi.universal.tf_activity_differential --activity results/activity/tf_activity.tsv --design design.tsv --out results/activity_comparison
```

The matrix has feature/peak IDs in the first column and sample columns. The network TSV has `source` (TF), `target` (matching feature ID), `weight` (signed or positive numeric value). For a binary motif-peak network, use weight 1 for motif-containing peaks and omit absent edges; ULM treats absent network weights as zero. Include the full prespecified feature universe in the matrix. Build the network from quality-controlled motif matches outside this adapter; its construction is not automatic yet.

Estimation calls the published decoupler ULM implementation, returning per-sample `activity`, `pvalue`, per-sample BH `padj` and `method`. It does not reproduce chromVAR, correct GC/accessibility bias automatically or prove TF binding. Use an independently justified network, normalization and controls. Differential activity currently supports independent Mann–Whitney tests only, requires complete sample matching and >=2 replicates per group, and reports delta = group2 minus group1.

The adapter pins decoupler to the verified 1.9 API. The Linux environment and CI pin compatible numerical dependencies. The local host smoke check does not replace a Linux run; no Windows compatibility workaround is part of the software.

## 11. DNA methylation

Input TSV requires `chrom`, zero-based `start`, optional exclusive `end`, and either `beta` in [0,1] or integer `methylated`/`unmethylated` counts. Coverage is derived from counts or supplied as `coverage`. Percentage input is rejected; convert percentages explicitly. Coverage and beta must agree with counts when both are supplied. Duplicate intervals and invalid coordinates are rejected.

```bash
python -m easyepi.specific.methylation.qc --project methylation.toml --out results/methylation_qc
python -m easyepi.specific.methylation.profile --project methylation.toml --min-coverage 5 --out results/methylation_profile
python -m easyepi.specific.methylation.differential --project replicated_methylation.toml --min-coverage 5 --out results/site_comparison
```

For beta-only tables with no measured coverage, use `--min-coverage 0`; unknown coverage is not fabricated. Differential output is an exploratory site-wise rank comparison, **not a coverage-aware beta-binomial model or DMR caller**. Each group needs >=2 replicates; sites with insufficient measurements remain untested. `delta_beta = group2 - group1`. Bismark one-based coverage formats, bedMethyl percentages and arrays require explicit conversion before input.

## 12. Histone chromatin states

```toml
[[samples]]
sample_id = "C1_H3K27ac"
assay = "histone"
cell = "C1"
mark = "H3K27ac"
peak = "peaks/C1_H3K27ac.bed"
```

Use the same `cell` across different marks for a biological condition. Duplicate cell/mark pairs are rejected; combine replicates using a documented strategy first.

```bash
python -m easyepi.specific.histone.chromhmm_prepare --project histone.toml --out results/chromhmm_prepare
python -m easyepi.specific.histone.chromhmm_binarize --project histone.toml --genome chrom.sizes --chromhmm-jar /tools/ChromHMM.jar --resolution 200 --out results/chromhmm_binary
python -m easyepi.specific.histone.chromhmm_train --binary-dir results/chromhmm_binary/binarized --genome-build hg38 --chromhmm-jar /tools/ChromHMM.jar --resolution 200 --states 15 --threads 8 --seed 1 --out results/chromhmm_model
python -m easyepi.specific.histone.chromhmm_annotate --states segmentation.bed --labels state_labels.tsv --out results/state_summary
```

Preparation writes a headerless CellMark table and stages input BED files. Binarization uses official `BinarizeBed -peaks -b ...`, appropriate for peak interval inputs, followed by optional LearnModel with the same bin size, seed and threads. `--genome` here is a chromosome-length TSV, not a FASTA. LearnModel can produce its own HTML and diagnostic images; those are native ChromHMM outputs.

State annotation summarizes coverage and optionally joins manually reviewed `state,label` definitions. It does not infer biological labels automatically. Inspect emissions and genomic enrichments before calling a state promoter, enhancer or repressed.

## 13. MNase

```bash
python -m easyepi.specific.mnase.qc --project mnase.toml --out results/mnase_qc
python -m easyepi.specific.mnase.profile --signal signal.tsv --out results/mnase_profile
```

`signal.tsv` is a single-column measured numeric signal vector. Occupancy and periodicity are descriptive summaries. `specific.mnase.positioning` now rejects execution because no established positioning backend has been connected. It no longer runs the unrelated ATAC peak caller. Nucleosome positioning and uncertainty estimates remain development work.

## 14. RNA integration

```bash
python -m easyepi.universal.rna_prepare --input rna_results.tsv --gene-column gene_id --out results/rna
python -m easyepi.universal.multiomics_integrate --gene-bed genes.bed --track ATAC=atac.bw --track H3K27ac=h3k27ac.bw --rna-results results/rna/rna_prepared.tsv --clusters 4 --seed 1 --out results/integration
```

Gene BED requires unique IDs and zero-based half-open coordinates; BED6 strand is honored. Integration constructs strand-oriented profiles, standardizes rows, uses scikit-learn KMeans and joins existing RNA result columns. It is descriptive; no RNA differential model, causal network, joint multi-omic latent model or batch correction is fitted. Nearest-TSS links remain candidates. Missing contigs remain missing rather than being silently filled with zero.

## 15. Interpretation checklist

Record genome assembly, annotation/database versions, sample exclusions, independent replicates, control strategy, normalization, test direction, background selection, tool versions and seeds. For publication, benchmark methods on real held-out datasets and compare with the original tool outputs; software tests alone do not establish scientific validity. See [reliability](RELIABILITY.md) and [references](REFERENCES.md).
