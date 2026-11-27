# Scale and biological benchmarks

`docs/PUBLISHING.md` step 7 asks for biological *and* performance evidence before
anyone claims method equivalence, runtime improvements or publication readiness.
`docs/RELIABILITY.md` records that biological agreement existed but had never been
written down, and that no performance numbers existed at all. This page is that
record, with the commands needed to reproduce it.

Two rules govern the numbers below, both borrowed from `scripts/bench_coverage.py`:

* a speed claim is accompanied by an equality check, so "faster" cannot quietly
  mean "different";
* every number names the input, the host and the tool versions, because a
  benchmark without those is an anecdote.

## 1. Biological agreement on a real cohort

The comparison is against ChIPseeker, on a shared gene model, so the only
variable is the labelling logic. Data: three donors (CS13-12383, CS14-12408,
CS16-12997) x four histone marks = twelve peak sets, all on chr1.

Basis: verbatim annotation-string equality; `tss_region` = -3000..3000 around the
TSS, `level` = transcript, `downstream_distance` = 300.

| dataset | peaks | label match | distance match |
| --- | --- | --- | --- |
| CS13-12383_H3K27ac | 2974 | 1.000 | 1.000 |
| CS13-12383_H3K27me3 | 709 | 1.000 | 1.000 |
| CS13-12383_H3K4me1 | 763 | 1.000 | 1.000 |
| CS13-12383_H3K4me3 | 1888 | 1.000 | 1.000 |
| CS14-12408_H3K27ac | 3003 | 1.000 | 1.000 |
| CS14-12408_H3K27me3 | 1212 | 1.000 | 1.000 |
| CS14-12408_H3K4me1 | 2201 | 1.000 | 1.000 |
| CS14-12408_H3K4me3 | 2178 | 1.000 | 1.000 |
| CS16-12997_H3K27ac | 3238 | 1.000 | 1.000 |
| CS16-12997_H3K27me3 | 909 | 1.000 | 1.000 |
| CS16-12997_H3K4me1 | 2484 | 1.000 | 1.000 |
| CS16-12997_H3K4me3 | 2087 | 1.000 | 1.000 |

**Interpretation.** With the same gene model, `easyepi` and ChIPseeker agree on
every annotation string and every distance across 23,646 peaks. This is the
strongest agreement evidence the project has, and it is a *like-for-like* result:
the earlier `annotation_agreement.json` figure of 0.725 for H3K4me3 is not a
contradiction. That run compared against a different gene model (UCSC), and
`heart_agreement.json` (0.917) plus `matrix.json` show the disagreement is
concentrated in Promoter-vs-Intron/Exon calls that follow from the TSS set, not
from the labelling rule. `heart_promoter_diagnosis.json` records 1,772 of 2,071
disagreements as `gtf_tss_within_3kb`, i.e. `easyepi` matched the GTF while
ChIPseeker matched a different transcript.

That distinction matters for how this result should be quoted: it establishes
**implementation equivalence given the same annotation**, not that either tool is
biologically correct. Choosing the TSS set remains the user's decision.

## 2. Runtime and memory on real data

Command:

```
python scripts/bench_scale.py \
    --bam <bench_a10>/subset_bam/H3K4me3.bam \
    --peaks <encode_h3k4me3>/peaks_full/H3K4me3/H3K4me3/H3K4me3_peaks.narrowPeak \
    --chrom chr1 --length 20000000 --workers 1,2,4,8 --repeat 2 \
    --json scale_bench.json
```

Input: `H3K4me3.bam` (795,329 reads on chr1), window chr1:0-20,000,000, bin 25 bp,
MAPQ >= 30. Host: 208 cores, Python 3.11.16. Times are best-of-2.

| workers | seconds | speedup | output | work RSS |
| --- | --- | --- | --- | --- |
| 1 | 0.902 | 1.00x | baseline | +38 MB |
| 2 | 0.483 | 1.87x | identical | +37 MB |
| 4 | 0.234 | 3.85x | identical | +38 MB |
| 8 | 0.213 | 4.24x | identical | +35 MB |

Annotation parse: 43,465 peaks in 0.031 s, +38 MB.

**How to read the memory column.** `work_rss` is the rise from the moment the
process has finished importing (about 130 MB) to its high-water mark during the
work. The absolute peak is reported as well in the JSON, but it is dominated by
import cost and comes back identical across worker counts, so it says nothing
about the workload. An earlier version of this script reported only the absolute
peak; the constant 167.9 MB it printed for every configuration was an artefact,
not a measurement.

The honest summary of the memory result is that **binned depth over 20 Mb costs
about 38 MB regardless of worker count**. That is expected: the chunked design
keeps one bin array plus per-chunk buffers, and the chunks are bounded by
`--chunk-mb`, not by how many workers exist.

**What this does not show.** These are scaling numbers on one chromosome window
of one sample. `RELIABILITY.md` item 9 asks for 100k+ peaks and many samples; that
remains open. No external tool was timed against, so this page makes no claim
that `easyepi` is faster than MACS3, deepTools or any R package.

## 3. Tool versions

Recorded alongside the comparison so it can be reproduced or refuted:

```
easyepi        0.6.0 at measurement time (the comparison predates 0.7.0)
python         3.11.16
platform       Linux-4.18.0-513.5.1.el8_9.x86_64-glibc2.28
deepTools      bamCoverage 4.0.0
R              4.5.3
ChIPseeker     1.46.1
csaw           1.44.0
DiffBind       3.20.0
DSS            2.58.0
macs3          unavailable at measurement time (installed later, 3.0.4)
```

The `easyepi` figure is what the comparison actually ran against. Re-running the
biological comparison under 0.7.0 has not been done and is listed as open work.

## Sources

The raw artefacts live outside the repository, under the analysis host's
`bench_a10/` directory:

```
bench_a10/results/verbatim_100.json              twelve-dataset agreement
bench_a10/results/heart_agreement.json           per-dataset label recall
bench_a10/results/heart_promoter_diagnosis.json  why Promoter calls disagreed
bench_a10/results/matrix.json                    agreement under several gene models
bench_a10/results/annotation_agreement.json      ENCODE H3K4me3, UCSC gene model
bench_a10/tool_versions.json                     tool versions above
```

They are not committed because they are regenerated artefacts of a dataset that
does not live in the repository. The commands in sections 1 and 2 are the
reproducible part; the numbers above are what those commands produced.
