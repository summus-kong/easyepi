# easyepi

Bulk epigenomics downstream analysis in pure Python. Every tool takes an aligned
**BAM** and writes ordinary files (BigWig, BED, TSV, JSON, PNG). Peak calling imports
MACS3, so that package must be installed; no other external program is
needed, and no step depends on another having been run first.

```bash
pip install -e ".[bam,viz]"      # pysam, pyBigWig, matplotlib, hmmlearn
```

> **Documentation: <https://summus-kong.github.io/easyepi/>**
>
> [Install](https://summus-kong.github.io/easyepi/install/) ·
> [Quick start](https://summus-kong.github.io/easyepi/quickstart/) ·
> [End-to-end examples](https://summus-kong.github.io/easyepi/EXAMPLES/) ·
> [Command reference](https://summus-kong.github.io/easyepi/commands/) ·
> [Python API](https://summus-kong.github.io/easyepi/PYTHON_API/) ·
> [Verification results](https://summus-kong.github.io/easyepi/VERIFICATION/)

## The chain (each step optional)

```
BAM ──bam2bw──> BigWig ──peaks──> peaks ──annotate──> gene assignment
 │                  │                │
 │                  │                └──tf-annotate──> named TFs
 │                  ├──footprint──> accessibility profile
 │                  ├──gene-panel─> multi-locus figure
 │                  └──differential-peaks──> two-group comparison
 │
 └──qc──> library verdict          qc reads BAM, BigWig, peaks or a bisulfite BAM
```

A minimal run (ChIP/ATAC):

```bash
easyepi bam2bw sample.bam -o tracks/ --assay atac      # 1. BAM -> BigWig
easyepi peaks sample.bam -o peaks/ --assay atac   # 2. call peaks
easyepi annotate --peaks peaks/sample/sample_peaks.narrowPeak \
    --gtf genes.gtf --out annot/
easyepi qc --assay atac --bam sample.bam --out qc.json  # anywhere in the chain
```

Single BAM or many; `--jobs` controls worker processes.

## Assays

| Assay | Status | What it does |
| --- | --- | --- |
| **ATAC-seq** | full | Tn5 shift, FRiP, peak calling, footprint, nucleosome |
| **ChIP-seq** (TF) | full | narrow/broad peaks, FRiP, signal enrichment |
| **ChIP-seq** (histone) | full | broad peaks, chromatin states across marks |
| **CUT&Tag / CUT&RUN** | full | same Tn5 path as ATAC |
| **MNase-seq** | partial | occupancy and periodicity; dyad calling needs paired-end |
| **WGBS / bisulfite** | full | per-cytosine calls and beta values read **from the BAM**; two-group DMR |

Methylation is read from the alignments themselves (`MM`/`ML` tags from samtools
1.10+, or Bismark's `XM` tag), so a bisulfite BAM is the only input -- no separate
`.cov`/`bedMethyl` table has to be produced first. A BAM carrying neither tag is
refused with an explanation rather than treated as unmethylated.

## Commands

### Prepare

| Command | Purpose |
| --- | --- |
| `bam2bw` | BAM -> BigWig. Tn5 shift, CPM/RPKM, blacklist, run-length encoded write, `--jobs` |

### Quality control

| Command | Purpose |
| --- | --- |
| `qc --assay …` | One entry point; `--assay` selects the metric set and thresholds |
| `make-qc-report`, `merge-qc-reports` | Batch and study-level roll-up |

`qc` is the single QC command. `atac-qc`, `chip-qc`, `histone-qc`, `mnase-qc` and
`methylation-qc` remain as aliases of the same implementation.

### Peaks

| Command | Purpose |
| --- | --- |
| `peaks` | MACS3 (imported as a library); narrow or broad; `--macs-arg` passthrough |
| `differential-peaks` | Consensus peaks -> counts -> TMM -> negative binomial -> BH FDR |
| `diff` | Exploratory two-group comparison of a normalised matrix |

### Annotation

| Command | Purpose |
| --- | --- |
| `annotate` | Promoter / exon / intron / intergenic assignment from a GTF |
| `gene-panel` | Multi-locus x multi-mark track figure with gene models |

### Motifs

| Command | Purpose |
| --- | --- |
| `motifs` | Known or de-novo motif discovery in peak sequences |
| `tf-annotate` | Name the TFs a peak can bind, with exact PWM thresholds |

### Structure

| Command | Purpose |
| --- | --- |
| `footprint` | Aggregate accessibility profile around motif sites, with significance |
| `nucleosome` | Occupancy and periodicity from a signal track |
| `nucleosome-positioning` | Dyad calls from paired-end fragments |

### States and integration

| Command | Purpose |
| --- | --- |
| `histone-states` | Chromatin states across marks (ChromHMM equivalent) |
| `dmr` | Differential methylation regions between two groups of bisulfite BAMs |
| `integrate` | Join assay layers over shared regions |
| `normalize` | CPM / spike-in / RPGC scaling |

### Figures

| Command | Purpose |
| --- | --- |
| `viz` | Heatmaps, metagene profiles, sample correlation |
| `washu` | Track manifest for the WashU Epigenome Browser |

### Project

| Command | Purpose |
| --- | --- |
| `init` | Write a project file from a folder of BAMs |
| `validate-project` | Validate a project file and its paths |

## Figures

Analysis never writes a figure unless asked. Functions return Matplotlib Figures:

```python
import pandas as pd
from easyepi.plot import plot_volcano

fig = plot_volcano(pd.read_csv("differential.tsv", sep="\t"))
fig.savefig("volcano.pdf", bbox_inches="tight")
```

Commands accept repeated `--figure-format` (`png`, `pdf`, `svg`, `tiff`, `jpg`).

## Verification

- `pytest tests -q` — 437 tests pass on Linux / Python 3.11
- `maintenance/verify_linux.py` — CLI surface and module imports
- `scripts/bench_coverage.py` — asserts the parallel coverage path produces
  byte-identical tracks to the serial one
- `scripts/verify_macs3_concordance.py` — the MACS3 wrapper produces byte-identical
  peak files to a direct `macs3 callpeak` in three configurations

Measured on a real ENCODE H3K4me3 BAM (chr1, 249 Mb, 25 bp bins):

| Path | Wall clock |
| --- | --- |
| Original per-column pileup | 64 s |
| Chunked, 1 process | 22.6 s |
| Chunked, 48 processes, run-length write | 14.5 s |

The remaining cost is BigWig I/O, not computation (compute 0.6 s, write 4.5 s after
run-length encoding).

## DMR

```bash
easyepi dmr --group-a ctrl1.bam,ctrl2.bam --group-b trt1.bam,trt2.bam \
    --label-a control --label-b treated --out dmr/regions.tsv
```

Two groups, any number of BAMs each. Sites are kept strand-specific; a per-site Wald
test on the log odds ratio with a Haldane-Anscombe correction, BH FDR, then
significant sites within `--max-gap` bp merge into regions whose p-value is
recomputed from the **summed** counts rather than from the smallest member.

Writes `regions.tsv`, `regions.sites.tsv` and `regions.json` (the JSON records every
parameter and the method's limitation).

**Limitation, stated in the output as well:** the test is binomial, not
beta-binomial. Methylation data is overdispersed relative to binomial, so p-values
are anti-conservative at high coverage and unstable at low coverage. `--min-coverage`
(default 5, applied per group) and `--min-diff` mitigate this but do not replace a
beta-binomial caller for publication-grade work.

## Not covered

Single-cell ATAC, Hi-C/TAD, gene-set enrichment (GO/KEGG), and direct comparison
against published callers for the pure-Python algorithms other than MACS3. See
[Reliability and limitations](https://summus-kong.github.io/easyepi/RELIABILITY/) before interpreting results.

## Documentation

The full manual is published at **<https://summus-kong.github.io/easyepi/>**. Every page listed here is built from
the repository by `scripts/gen_docs.py` and deployed on each push to `main`.

### Getting started

| Page | What it covers |
| --- | --- |
| [Home](https://summus-kong.github.io/easyepi/) | What the suite does, the analysis chain, supported assays |
| [Install](https://summus-kong.github.io/easyepi/install/) | Extras, system requirements, what each command needs |
| [Quick start](https://summus-kong.github.io/easyepi/quickstart/) | One BAM to peaks to annotation, end to end |
| [End-to-end examples](https://summus-kong.github.io/easyepi/EXAMPLES/) | Worked runs for ATAC, ChIP, CUT&Tag, MNase, WGBS |

### Per-assay guides

| Page | What it covers |
| --- | --- |
| [ATAC-seq](https://summus-kong.github.io/easyepi/assays/atac/) | Tn5 shift, FRiP, footprinting, nucleosome spacing |
| [ChIP-seq](https://summus-kong.github.io/easyepi/assays/chip/) | TF and histone peaks, narrow vs broad, enrichment |
| [CUT&Tag / CUT&RUN](https://summus-kong.github.io/easyepi/assays/cuttag/) | The Tn5 path shared with ATAC |
| [MNase-seq](https://summus-kong.github.io/easyepi/assays/mnase/) | Occupancy and periodicity |
| [Bisulfite / WGBS](https://summus-kong.github.io/easyepi/assays/wgbs/) | Beta values read from the BAM, two-group DMR |

### Command reference

| Page | What it covers |
| --- | --- |
| [All commands](https://summus-kong.github.io/easyepi/commands/) | Every subcommand, generated from the live CLI |
| [bam2bw](https://summus-kong.github.io/easyepi/commands/bam2bw/) | BAM to BigWig, normalisation, blacklist |
| [qc](https://summus-kong.github.io/easyepi/commands/qc/) | The single QC entry point and its metric sets |
| [peaks](https://summus-kong.github.io/easyepi/commands/peaks/) | MACS3 as a library, narrow and broad |
| [annotate](https://summus-kong.github.io/easyepi/commands/annotate/) | Genomic bins, promoter windows, priority order |
| [tf-annotate](https://summus-kong.github.io/easyepi/commands/tf-annotate/) | Naming the TFs a peak can bind |
| [dmr](https://summus-kong.github.io/easyepi/commands/dmr/) | Differential methylation regions |
| [differential-peaks](https://summus-kong.github.io/easyepi/commands/differential-peaks/) | Consensus peaks to FDR |
| [diff](https://summus-kong.github.io/easyepi/commands/diff/) | Count-model differential analysis (PyDESeq2) |
| [visualisation](https://summus-kong.github.io/easyepi/commands/viz/) | Heatmaps, metagene profiles, correlation |

### Methods and reference

| Page | What it covers |
| --- | --- |
| [Usage and input contracts](https://summus-kong.github.io/easyepi/USAGE/) | What each command accepts and rejects |
| [Python API](https://summus-kong.github.io/easyepi/PYTHON_API/) | Calling the library directly instead of the CLI |
| [Peak calling modes](https://summus-kong.github.io/easyepi/PEAK_MODES/) | Choosing narrow vs broad, and when it matters |
| [QC thresholds](https://summus-kong.github.io/easyepi/QC/) | The shipped cut-offs per assay and how to override them |
| [Motif methods](https://summus-kong.github.io/easyepi/MOTIF_BACKENDS/) | Which motif backends exist and what they assume |
| [Method references](https://summus-kong.github.io/easyepi/REFERENCES/) | The papers each implementation follows |
| [Verification results](https://summus-kong.github.io/easyepi/VERIFICATION/) | Every measured claim, including failures |
| [Reliability and limitations](https://summus-kong.github.io/easyepi/RELIABILITY/) | Where results should not be trusted |
| [Scope and evidence](https://summus-kong.github.io/easyepi/INNOVATION/) | What is novel and what is inherited |

### API reference
Module-level pages, one per public module, generated from docstrings:
[assays](https://summus-kong.github.io/easyepi/api/assays/), [bam2bigwig](https://summus-kong.github.io/easyepi/api/bam2bigwig/),
[coverage_parallel](https://summus-kong.github.io/easyepi/api/coverage_parallel/), [peaks_cli](https://summus-kong.github.io/easyepi/api/peaks_cli/),
[macs_options](https://summus-kong.github.io/easyepi/api/macs_options/), [tool_wrapper](https://summus-kong.github.io/easyepi/api/tool_wrapper/),
[methylation_bam](https://summus-kong.github.io/easyepi/api/methylation_bam/), [dmr](https://summus-kong.github.io/easyepi/api/dmr/),
[tf_annotate](https://summus-kong.github.io/easyepi/api/tf_annotate/), [histone_states](https://summus-kong.github.io/easyepi/api/histone_states/),
[nucleosome_positioning](https://summus-kong.github.io/easyepi/api/nucleosome_positioning/),
[viz.genomeviews](https://summus-kong.github.io/easyepi/api/genomeviews/), [viz.panels](https://summus-kong.github.io/easyepi/api/panels/),
[formats](https://summus-kong.github.io/easyepi/api/formats/), [errors](https://summus-kong.github.io/easyepi/api/errors/).

### Contributing
[Development guide](https://summus-kong.github.io/easyepi/DEVELOPMENT/) · [Releasing](https://summus-kong.github.io/easyepi/PUBLISHING/)

License: BSD-3-Clause.
