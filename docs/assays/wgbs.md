# WGBS / bisulfite

Methylation is read from the alignments themselves: `MM`/`ML` tags (samtools 1.10+) or
Bismark's `XM` tag. A `.cov` or `.bedMethyl` table does not have to be produced first,
so the bisulfite BAM is the only input.

## Per-sample QC

```bash
easyepi qc --assay methylation --bam ctrl1.bam --out qc/ctrl1.json
```

Reports the beta-value distribution, the hypomethylated/hypermethylated fractions and
bimodality. A BAM with neither tag is refused with an explanation rather than treated
as unmethylated.

Betas are per CpG dyad and require `--min-coverage` (5 by default) reads at the site;
sites below it are counted in `methylation_bam.dropped_below_coverage` and left out
rather than averaged from one or two reads. The report records which encoding was
used (`MM/ML` or `XM`) and how many reads and calls it read, so a suspiciously thin
result is visible instead of silent.

A table that was already called (`.cov`, `.bedMethyl`, or any TSV with a `beta`
column) is also accepted when the calls arrive that way: point `--bam` at the BAM and
the project file's `methylation` field at the table. When both are present the table
is used, because it is the input the caller named explicitly.

## Two-group differential methylation

```bash
easyepi dmr \
    --group-a ctrl1.bam,ctrl2.bam \
    --group-b trt1.bam,trt2.bam \
    --label-a control --label-b treated \
    --min-coverage 10 --min-diff 0.2 --alpha 0.05 \
    --out dmr/regions.tsv
```

Writes `regions.tsv`, `regions.sites.tsv` and `regions.json`.

| Parameter | Default | Why it exists |
| --- | --- | --- |
| `--min-coverage` | 5 | applied per group; 1x in one group cannot estimate that group's level |
| `--min-diff` | 0 | a large cohort makes trivial differences significant |
| `--max-gap` | 100 | merges significant sites into regions |
| `--alpha` | 0.05 | BH FDR threshold for sites and regions |

Sites are kept **strand-specific**: the two strands of a CpG are separate
observations, and merging them hides hemi-methylation.

## Limitation

The per-site test is a Wald test on the log odds ratio with a Haldane-Anscombe
correction — binomial, not beta-binomial. Methylation data is overdispersed relative
to binomial, so p-values are anti-conservative at high coverage and unstable at low
coverage. `--min-coverage` and `--min-diff` mitigate this but do not replace DSS,
bsseq or metilene for publication-grade calling. The same sentence is written into
`regions.json` under `limitation`.

## Assumptions

- reads are aligned to the bisulfite-converted genome; the CIGAR is used, not the
  reference sequence
- `--context CpG` is the default; CHG and CHH are supported but sparser
- reading the tags needs `pysam` (`pip install "easyepi[bam]"`); without it the
  command says which install line is missing instead of reporting no methylation
