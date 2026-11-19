# A10 -- 1:1 reproduction of ChIPseeker's annotation

Question: does EasyEpi's peak annotation reproduce ChIPseeker, on data that is not the
slice the port was developed against?

**Result: 100.0000% label agreement and 100.0000% distance-to-TSS agreement --
23,646 / 23,646 peaks over 12 samples -- comparing the label strings verbatim.**

"Verbatim" is the whole point. ChIPseeker does not emit bare bin names; it emits
`Promoter (<=1kb)`, `Exon (T1/G1, exon 2 of 4)`, `3' UTR`, `Downstream (<=300bp)` and
`Distal Intergenic`. EasyEpi now produces the same strings, including the sub-category
splits and the exon/intron numbering. An earlier pass of this module reported 99.88% on
*coarse* bins (a 6-bin fold of the same output); that number is superseded by this one.

Two things had to be right, and both are covered below:

1. both tools must read the **same gene model** (getting this wrong yields a misleading
   91%);
2. the port must follow ChIPseeker's assignment structures **literally**, not a
   re-derivation that merely looks equivalent (four such defects are listed under
   "What a literal port changed").

## Inputs

| What | Path | Identity |
|---|---|---|
| Peaks | `heart_peaks_clean/*.bed` | 12 narrowPeak sets (3 samples × 4 marks), chr1 only, 709–3,238 peaks each |
| Gene model, as shipped | `bench_annotation/annot/txdb_full.gtf` | md5 `1057e0510437725011527118a7f7b17f` |
| Gene model, repaired | `bench_annotation/annot/txdb_fixed.gtf` | md5 `c6f46990d238ac9695ef34c4e2cdf36d` |
| TxDb from repaired model | `bench_annotation/annot/TxDb.fixed.sqlite` | 412,044 transcripts |
| ChIPseeker package model | `TxDb.Hsapiens.UCSC.hg38.knownGene` 3.22.0 | 412,044 transcripts, but ChIPseeker uses 326,307 (see below) |

Peaks are chr1-subsampled narrowPeak, so this measures the annotation rules, not broad
marks; a 15-peak H3K27me3 narrowPeak set from the earlier ENCODE pass was excluded for
that reason and is not part of these numbers.

## Why the gene model had to be repaired

`txdb_full.gtf` cannot serve as a shared gene model, and the failure is silent:

* its `transcript` rows carry `transcript_id "ENST…"`, while its `exon`/`CDS`/`UTR` rows
  carry the UCSC numeric **gene** ids (`transcript_id "1"`); **0 of 412,044** exon ids
  resolve to a transcript row.
* a TxDb built from it reports **824,088** transcripts for 412,044 transcript rows --
  every numeric id becomes a transcript of its own.
* `txdb_hg38.gffread.gtf` *does* link the namespaces (412,044/412,044) but omits
  `gene_id` on 196,282 UTR/exon rows, so `makeTxDbFromGFF` aborts with
  `'genes$gene_id' must be a character vector with no NAs`.

`repair_gtf.py` therefore rebuilds the shared model from `txdb_hg38.gffread.gtf`,
backfilling a missing `gene_id` from the `transcript_id → gene_id` map learned from the
same file's transcript rows. It touches no coordinate and adds/drops no row
(3,701,168 rows, 199,111 `gene_id` values filled, 0 remaining).

Separately, ChIPseeker's `getGene(..., by="transcript")` calls `transcriptsBy()`, which
groups by gene and silently returns only the 326,307 transcripts that carry a `gene_id`
-- 85,737 are dropped. On the repaired model this number is 0 dropped, so the check in
`annotate_chipseeker_txdb.R` reports a truncation if it ever recurs.

## Method

1. `repair_gtf.py txdb_hg38.gffread.gtf txdb_fixed.gtf`
2. `annotate_chipseeker_txdb.R <bed> <out.tsv> TxDb.fixed.sqlite transcript` -- the
   ChIPseeker side, `tssRegion=c(-3000,3000)`, `level="transcript"`.
3. `easyepi annotate --peaks <bed> --gtf txdb_fixed.gtf --out <dir>` -- the EasyEpi side.
4. `matrix.py` -- compares every EasyEpi gene model × ChIPseeker gene model pair and
   writes `results/matrix.json`.
5. `verbatim.py` -- compares the `annotation` strings **without folding**: identical text,
   or it is a mismatch. Writes `results/verbatim_100.json`.

Comparison is positional row-by-row. The BED 0-based → GRanges 1-based conversion happens
once, in `run_annotate`, so both tools are fed the same 1-based start.

Result-changing parameters: `tssRegion`, `downstreamDistance`, `level=transcript`, the
gene model, and the peak sets. The BED offset is a coordinate contract, not a tuning knob:
if it is wrong the distance is off by one everywhere.

## Result

Matched gene model, verbatim string equality (`results/verbatim_100.json`):

| dataset | peaks | label match | % | dist match | % |
|---|---|---|---|---|---|
| CS13-12383_H3K27ac | 2,974 | 2,974 | 100.00 | 2,974 | 100.00 |
| CS13-12383_H3K27me3 | 709 | 709 | 100.00 | 709 | 100.00 |
| CS13-12383_H3K4me1 | 763 | 763 | 100.00 | 763 | 100.00 |
| CS13-12383_H3K4me3 | 1,888 | 1,888 | 100.00 | 1,888 | 100.00 |
| CS14-12408_H3K27ac | 3,003 | 3,003 | 100.00 | 3,003 | 100.00 |
| CS14-12408_H3K27me3 | 1,212 | 1,212 | 100.00 | 1,212 | 100.00 |
| CS14-12408_H3K4me1 | 2,201 | 2,201 | 100.00 | 2,201 | 100.00 |
| CS14-12408_H3K4me3 | 2,178 | 2,178 | 100.00 | 2,178 | 100.00 |
| CS16-12997_H3K27ac | 3,238 | 3,238 | 100.00 | 3,238 | 100.00 |
| CS16-12997_H3K27me3 | 909 | 909 | 100.00 | 909 | 100.00 |
| CS16-12997_H3K4me1 | 2,484 | 2,484 | 100.00 | 2,484 | 100.00 |
| CS16-12997_H3K4me3 | 2,087 | 2,087 | 100.00 | 2,087 | 100.00 |
| **total** | **23,646** | **23,646** | **100.0000** | **23,646** | **100.0000** |

Difference is exactly zero: no residual disagreement to explain away.

### Parameter sweep

The same comparison at seven `(tssRegion, downstreamDistance)` settings, 600 peaks each,
all exact on both label and distance:

| tssRegion | downstreamDistance | label | distance |
|---|---|---|---|
| ±3000 | 300 (default) | 600/600 | 600/600 |
| ±5000 | 300 | 600/600 | 600/600 |
| ±3000 | 5500 | 600/600 | 600/600 |
| ±5000 | 5500 | 600/600 | 600/600 |
| ±3000 | 7500 | 600/600 | 600/600 |
| ±1500 | 300 | 600/600 | 600/600 |
| ±3000 | 1000 | 600/600 | 600/600 |

Promoter sub-categories are `ceiling(max|tssRegion|/1000)` wide, so ±3000 gives the
familiar three buckets and ±5000 gives five. Below 2000 R skips sub-categorisation
entirely and emits the bare `Promoter`, which is also reproduced.

## What a literal port changed

Four defects were found by reading ChIPseeker's source and confirming each against
measured output. Each one was invisible to a coarse-bin comparison:

1. **Downstream runs after the priority loop.** In `rev(genomicAnnotationPriority)` the
   order is Intergenic, Downstream, Intron, Exon, 3UTR, 5UTR, Promoter -- and *both*
   Intergenic and Downstream fall into R's trailing `else` branch, so neither writes a
   label inside the loop. The Downstream rename happens afterwards. Treating Downstream
   as an in-loop stage was structurally wrong.
2. **`getFirstHitIndex` keeps the lowest array index.** `findOverlaps` returns hits
   ordered by subject index, and `GRegion <- unlist(...)` gives that index its meaning.
   Picking by coordinate, or by TxID, both selected the wrong transcript.
3. **`exonsBy()` returns minus-strand exons in descending order.** `exon_rank` follows
   that array (`exon_rank == seq_along(g)` is TRUE), so ranking a sorted-ascending array
   with `N - i` produced the right *number* but the wrong *row* when two exons overlapped
   one peak. The row list now keeps R's array order and a separate `(start, index)` view
   drives candidate lookup.
4. **Coordinate basis.** `annotatePeak` sees 1-based GRanges (`readPeakFile` converts
   BED). The +1 belongs at the reading boundary (`run_annotate`), not inside the
   classifier -- placing it inside double-shifted callers that already held GRanges
   coordinates and cost 1 bp on every distance.

## What the mismatched arms show

Same runs, deliberately wrong pairings -- these quantify how much of a naive comparison
is gene-model disagreement rather than rule disagreement:

| pair | datasets | label agree | dist ≤1 bp |
|---|---|---|---|
| `ee_fixed` × `cs_gtf_fixed` | 12 | **100.0000%** | **100.0000%** |
| `ee_txdb_full` × `cs_ucsc` | 12 | 91.24% | 64.94% |
| `ee_fixed` × `cs_gtf_full` | 1 | 57.63% | 100.00% |

The last row is the clearest single lesson: the distances already agree perfectly, but
the labels collapse, because the two tools take their exon/intron tables from different
files. A distance check alone would have declared that pairing a success.

## Scripts

| File | Role |
|---|---|
| `scripts/repair_gtf.py` | Backfill `gene_id`; produces the shared gene model |
| `scripts/annotate_chipseeker_txdb.R` | ChIPseeker against an explicit TxDb, with a truncation check |
| `scripts/annotate_chipseeker.R` | Original script pinned to the UCSC package TxDb (kept for the mismatched arm) |
| `scripts/matrix.py` | All-pairs gene-model comparison → `results/matrix.json` |
| `scripts/verbatim.py` | Unfolded label/distance equality → `results/verbatim_100.json` |
| `scripts/compare_heart.py` | Per-dataset agreement against the package TxDb |
| `scripts/diagnose_promoter.py` | GTF-based check of disputed Promoter calls |
| `scripts/flank_defect.R` | Quantifies ChIPseeker's `flank_gene_distances` zeroing defect and prints ten instances |
| `scripts/flank_reference.R` | Per-dataset flank reference in long format, both `dist_raw` and `dist_fixed` |
| `scripts/build_gencode_slice.sh` | chr1 slice of GENCODE, the only input that carries explicit `UTR` rows |

## `sameStrand` and `level`

`annotatePeak` has two switches that change which feature a peak is measured
against. Both are ported and verified on a 400-peak set built from the same
heart peaks but carrying real `+` (232) / `-` (168) strands:

| level | sameStrand | labels | distanceToTSS |
|---|---|---|---|
| transcript | FALSE | 400/400 | 400/400 |
| transcript | TRUE | 400/400 | 400/400 |
| gene | FALSE | 400/400 | 400/400 |
| gene | TRUE | 400/400 | 400/400 |

This set exists because it is the only way to reach those code paths: the
narrowPeak files used everywhere else in this benchmark carry `signalValue` in
column 6, and ChIPseeker's `readPeakFile` turns that into `*` for every peak --
so `sameStrand` never does anything and the default numbers above cannot
distinguish the branches.

Four behaviours had to be measured rather than inferred, each found by
comparing against R on this set:

1. `IRanges::follow`/`precede` are **strand aware**. For a minus-strand query
   they swap: `follow` returns the feature on the right and `precede` the one
   on the left. Getting this wrong cost 168 of 400 distances (204/400).
2. `distanceToTSS` reads `start(peaks)` for `psD` and `end(peaks)` for `peD`,
   and the **sign** comes from the *feature* strand -- so a minus peak and a
   minus feature both flip, independently.
3. A TSS falling **inside** the peak sets the distance to `0`
   (`ignoreOverlap`), regardless of which side would otherwise win.
4. `level = "gene"` replaces `features` wholesale (`getGene(TxDb, by = level)`,
   line 80), so the **TSS vector itself** changes, not just the follow() input.
   `annotation` is already final by line 74, so the categories stay the same;
   only the distance and `geneId` move.

`genes(TxDb)` reproduces exactly: 123,262 groups → 121,093 genes on the
GTF-derived TxDb, with 2,169 groups dropped for spanning seqlevels and 238 for
spanning strands -- the same counts R reports.

## Versions

Read from the runtime, not assumed: Python 3.11.16, easyepi 0.6.0, R 4.5.3, ChIPseeker
1.46.1, TxDb.Hsapiens.UCSC.hg38.knownGene 3.22.0, GenomicFeatures 1.62.0, txdbmaker
1.6.2, RSQLite 3.53.3. macs3: unavailable on the run host (peaks were called earlier
with it and are inputs here).

## Reporting: `annotation_stat.tsv`

EasyEpi writes one distribution, and it is ChIPseeker's own table -- the six-bin fold
that used to sit beside it has been removed, so there is no second vocabulary to
reconcile.

| File | Units | Categories | Source |
|---|---|---|---|
| `annotation_stat.tsv` | percent | ~9 (data dependent) | port of ChIPseeker's `getGenomicAnnoStat` |

Both table names are configurable: `--annotated-name` (default `annotated.tsv`) and
`--stat-name` (default `annotation_stat.tsv`).

`annotation_stat.tsv` is a literal port and reproduces ChIPseeker exactly:

* `Promoter (<=1kb)` / `(1-2kb)` / `(2-3kb)` stay distinct rows -- the count follows
  `ceiling(max|tssRegion|/1000)`, so ±5 kb gives five Promoter rows;
* `Exon (…, exon 1 of N)` becomes `1st Exon`, other exons `Other Exon`; likewise
  `1st Intron` / `Other Intron`;
* `Downstream (<=300bp)` becomes `Downstream (<=300)` -- **without** the `bp` suffix.
  That inconsistency is ChIPseeker's own (`paste0("Downstream (<=", dsd, ")")` in
  `getGenomicAnnoStat` versus `"Downstream (<=300bp)"` in the per-peak label) and is
  reproduced rather than corrected;
* frequencies are `table()/sum*100`, and zero-count categories are dropped;
* rows follow `lvs`: Promoter labels (sorted) then
  `5' UTR, 3' UTR, 1st Exon, Other Exon, 1st Intron, Other Intron, Downstream, <iglab>`.

Verified against `getAnnoStat()` on the same 600 peaks and the same TxDb: 9/9 rows
identical at ±3 kb and 11/11 at ±5 kb. The files are **byte-identical** (`diff -q`
clean, matching md5), not merely numerically equal.

One known locale dependency: `getGenomicAnnoStat` sorts the Promoter labels with
`sort()`. Under `en_US.UTF-8` that orders `(<=1kb), (1-2kb), (2-3kb)` by number, which
is what is reproduced; under `LC_ALL=C` ChIPseeker itself would emit a different order,
and no single implementation matches both.

## CLI parameters now wired

These previously parsed and were then ignored; they now reach the classifier:

| flag | default | effect |
|---|---|---|
| `--promoter-upstream` | 3000 | `tssRegion[1]`; also drives the Promoter category count |
| `--promoter-downstream` | 3000 | `tssRegion[2]` |
| `--downstream-distance` | 300 | `getOption("ChIPseeker.downstreamDistance")` |
| `--no-bins-summary` | off | skips `annotation_stat.tsv` |
| `--annotated-name` | `annotated.tsv` | per-peak table filename |
| `--stat-name` | `annotation_stat.tsv` | feature distribution filename |

## Tests

`tests/test_annotate_bins.py` asserts the exact strings ChIPseeker 1.46.1 returns for a
synthetic GTF -- `Promoter (<=1kb)`, `Exon (T1/G1, exon 2 of 4)`, `3' UTR`,
`Intron (T1/G1, intron 1 of 3)`, `Distal Intergenic` -- together with their distances, plus
the minus-strand TSS, the coordinate basis, and the `annotation_stat` category set,
percentage folding, missing-`bp` Downstream label and `ignore_downstream` option.
The expected values were read off ChIPseeker's own output for that fixture, not
paraphrased. Full suite: 670 passed; this file: 16 passed.

### Strand and level cases

`tests/test_annotate_bins.py` gained four cases that pin the behaviours above:
a minus peak swapping follow/precede, a TSS inside the peak forcing distance 0,
`sameStrand` gating sub-feature hits, and `level="gene"` keeping the category
while moving the distance. They are unit-level (synthetic `Features`), so they
run without R; the numbers in the table come from the CLI against R.

## Column-level agreement

`annotated.tsv` mirrors `as.data.frame(annotatePeak(...))`: 14 columns at
`--level transcript`, 13 at `--level gene` (no `transcriptId` there). Measured
against R on the 400-peak stranded set, per column:

| column | transcript | gene |
|---|---|---|
| `start` `end` `width` `strand` | 400/400 | 400/400 |
| `annotation` | 400/400 | 400/400 |
| `geneChr` | 400/400 | 400/400 |
| `geneStart` `geneEnd` `geneLength` `geneStrand` `geneId` | 399/400 | 400/400 |
| `transcriptId` | 399/400 | (absent) |
| `distanceToTSS` | 400/400 | 400/400 |

`start` is the **1-based GRanges coordinate**, which is what `start(peak)`
returns -- a BED row beginning at 713700 is written as 713701.

### The single remaining peak

One peak (`chr1:11846801-11853300`) differs in the six `gene*`/`transcriptId`
columns, and the cause is not a bug in the port. Three transcripts have a TSS
inside it, and `annotatePeak` keeps the one whose `subjectHits` index is
smallest:

```
subjectHits 211788  TxID:19547   <- kept
            211789  TxID:19548
            339265  TxID:1524
```

`subjectHits` is a position in `getGene(TxDb, by = "transcript")`, whose row
order comes from the TxDb's own `transcripts()` query: it is neither sorted by
TxID, nor by genomic position, nor grouped by seqlevel (23,613 out-of-order
pairs and 40,303 seqlevel changes on this TxDb). Rebuilding from the GTF gives
`TxID:1524` the lower index, so the tie resolves the other way. That ordering
is not recoverable from a GTF file, so this peak cannot be matched by any port
that reads a GTF instead of the TxDb.

`annotation` and `distanceToTSS` are unaffected: every candidate is a TSS
inside the peak, so the category is `Promoter` and the distance is 0 no matter
which one wins.

## `flank_*` columns: where EasyEpi deliberately differs

`--add-flank-info` writes `flank_txIds`, `flank_geneIds` and
`flank_gene_distances`, the three columns `annotatePeak(addFlankGeneInfo=TRUE,
flankDistance=...)` appends. The first two match ChIPseeker exactly. The third
does not, and the reason is a defect in the reference, not in the port.

### What the reference does

`ChIPseeker/R/getFlankingGene.R`, function `getAllFlankingGene`. In 1.46.1 this
is lines 21-23; in **1.49.5, the current devel, it is lines 111-113 and is
unchanged**:

```r
overlapHit <- findOverlaps(peak.gr2, unstrand(featureHit))
hitInfo$distance <- NA
hitInfo$distance[subjectHits(overlapHit)] <- 0      # queryHits is dropped
```

`featureHit` is `features[sh]`: the hit features **expanded per (peak, feature)
pair**, so a feature that sits in several peaks' windows occupies several
elements. `overlapHit` returns `(peak, featureHit element)` pairs, but line 113
keeps only the element half. "Peak A overlaps element k" is therefore written
to **row k**, and row k may belong to a different peak.

The `ii <- is.na(hitInfo$distance); hitInfo$distance[ii] <- dd[ii]` pass that
follows cannot undo it: a wrongly zeroed row is no longer `NA`.

### Measured on all 12 datasets

Every dataset in `heart_peaks_clean/` was run through the EasyEpi CLI with
`--add-flank-info --flank-distance 5000` and through `flank_reference.R`, which
reports each row twice: `dist_raw` as the reference writes it, and `dist_fixed`
as the same row scores when `queryHits` is kept. 23,646 peaks, 197,028
(peak, feature) rows. The run is repeated at `--level gene` (121,093 gene
features, 48,241 rows) because `level` swaps the feature set, and the defect
follows it. Numbers in this section and in the CLI tables below are reproduced
in `results/flank_12datasets.json`.

| dataset | peaks | rows | rows the reference zeros | wrongly zeroed | peaks affected |
|---|---|---|---|---|---|
| CS13-12383_H3K27ac | 2974 | 22,856 | 18,821 | 2,202 (9.63%) | 614 (20.6%) |
| CS13-12383_H3K27me3 | 709 | 6,689 | 5,633 | 841 (12.57%) | 131 (18.5%) |
| CS13-12383_H3K4me1 | 763 | 5,552 | 4,469 | 579 (10.43%) | 159 (20.8%) |
| CS13-12383_H3K4me3 | 1888 | 16,879 | 13,262 | 1,685 (9.98%) | 291 (15.4%) |
| CS14-12408_H3K27ac | 3003 | 24,253 | 19,678 | 2,674 (11.03%) | 741 (24.7%) |
| CS14-12408_H3K27me3 | 1212 | 10,316 | 8,900 | 1,256 (12.18%) | 288 (23.8%) |
| CS14-12408_H3K4me1 | 2201 | 17,696 | 14,786 | 2,668 (15.08%) | 608 (27.6%) |
| CS14-12408_H3K4me3 | 2178 | 19,683 | 15,851 | 2,129 (10.82%) | 385 (17.7%) |
| CS16-12997_H3K27ac | 3238 | 25,931 | 21,240 | 2,664 (10.27%) | 749 (23.1%) |
| CS16-12997_H3K27me3 | 909 | 8,142 | 6,900 | 880 (10.81%) | 175 (19.3%) |
| CS16-12997_H3K4me1 | 2484 | 20,412 | 16,979 | 3,151 (15.44%) | 655 (26.4%) |
| CS16-12997_H3K4me3 | 2087 | 18,619 | 14,794 | 1,939 (10.41%) | 359 (17.2%) |
| **total** | **23,646** | **197,028** | **161,313** | **22,668 (14.05% of zeros, 11.50% of rows)** | **5,155 (21.8%)** |

At `--level gene`, same peaks, gene-level features:

| level | features | rows | rows the reference zeros | wrongly zeroed | peaks affected |
|---|---|---|---|---|---|
| transcript | 412,044 | 197,028 | 161,313 | 22,668 (11.50%) | 5,155 (21.8%) |
| gene | 121,093 | 48,241 | 38,654 | 4,863 (10.08%) | 3,022 (12.78%) |

Gene level mislabels slightly less because there are 3.4x fewer features to
collide: the same gene lands in several windows less often than the same
transcript does. Both levels keep the one-directional property -- no genuine
overlap is ever missed.

The share is stable across marks (H3K27ac / H3K27me3 / H3K4me1 / H3K4me3) and
samples: 9.6-15.4% of rows, 15-28% of peaks. It is a property of how densely
the windows overlap, not of one dataset.

The error is one-directional: the reference never misses a real overlap, it
only marks extra rows. Downstream code that reads `flank_gene_distances == 0`
as "the peak lies inside this feature" therefore over-counts by ~14% and never
under-counts.

### CLI agreement across the 12 datasets

`flank_txIds` is compared as a set and each distance is matched per transcript,
because the element order inside a `;` list differs for the reason recorded
under "The single remaining peak" below.

| dataset | peaks | txId set equal | distances == corrected | distances == reference-as-written |
|---|---|---|---|---|
| CS13-12383_H3K27ac | 2974 | 2974 | 2974 | 2360 |
| CS13-12383_H3K27me3 | 709 | 709 | 709 | 578 |
| CS13-12383_H3K4me1 | 763 | 763 | 763 | 604 |
| CS13-12383_H3K4me3 | 1888 | 1888 | 1888 | 1597 |
| CS14-12408_H3K27ac | 3003 | 3003 | 3003 | 2262 |
| CS14-12408_H3K27me3 | 1212 | 1212 | 1212 | 924 |
| CS14-12408_H3K4me1 | 2201 | 2201 | 2201 | 1593 |
| CS14-12408_H3K4me3 | 2178 | 2178 | 2178 | 1793 |
| CS16-12997_H3K27ac | 3238 | 3238 | 3238 | 2489 |
| CS16-12997_H3K27me3 | 909 | 909 | 909 | 734 |
| CS16-12997_H3K4me1 | 2484 | 2484 | 2484 | 1829 |
| CS16-12997_H3K4me3 | 2087 | 2087 | 2087 | 1728 |
| **total** | **23,646** | **23,646 (100%)** | **23,646 (100%)** | 18,491 (78.2%) |

Per entry over all 197,028 pairs:

| outcome | entries | share |
|---|---|---|
| identical to both the reference and the correction | 174,360 | 88.50% |
| identical to the correction only | 22,668 | 11.50% |
| identical to the reference only | 0 | 0.00% |
| **identical to neither** | **0** | **0.00%** |

The two zeros are the point. EasyEpi never emits a distance that matches
neither scoring, and never emits one that matches only the defective scoring.
The 22,668 `correction only` entries equal, exactly, the number of rows the
reference mislabels -- the two scripts agree independently. The 23,646 - 18,491
= 5,155 peaks where the distances differ from the reference-as-written also
equal, exactly, the 5,155 affected peaks counted above.

At `--level gene` the CLI run is repeated over the same 12 datasets. R drops
`flank_txIds` at this level, so only `flank_geneIds` and
`flank_gene_distances` are written:

| dataset | peaks | geneId set equal | distances == corrected | distances == reference-as-written |
|---|---|---|---|---|
| CS13-12383_H3K27ac | 2974 | 2974 | 2974 | 2601 |
| CS13-12383_H3K27me3 | 709 | 709 | 709 | 632 |
| CS13-12383_H3K4me1 | 763 | 763 | 763 | 679 |
| CS13-12383_H3K4me3 | 1888 | 1888 | 1888 | 1731 |
| CS14-12408_H3K27ac | 3003 | 3003 | 3003 | 2558 |
| CS14-12408_H3K27me3 | 1212 | 1212 | 1212 | 1023 |
| CS14-12408_H3K4me1 | 2201 | 2201 | 2201 | 1840 |
| CS14-12408_H3K4me3 | 2178 | 2178 | 2178 | 1965 |
| CS16-12997_H3K27ac | 3238 | 3238 | 3238 | 2811 |
| CS16-12997_H3K27me3 | 909 | 909 | 909 | 798 |
| CS16-12997_H3K4me1 | 2484 | 2484 | 2484 | 2102 |
| CS16-12997_H3K4me3 | 2087 | 2087 | 2087 | 1884 |
| **total** | **23,646** | **23,646 (100%)** | **23,646 (100%)** | 20,624 (87.2%) |

Per entry over 48,241 gene pairs:

| outcome | entries | share |
|---|---|---|
| identical to both | 43,378 | 89.92% |
| identical to the correction only | 4,863 | 10.08% |
| identical to the reference only | 0 | 0.00% |
| **identical to neither** | **0** | **0.00%** |

4,863 equals the gene-level reference count exactly, and a further 2,990 peaks
carry `NA` in both -- no gene within 5 kb -- and agree there too. So the same
two zeros hold at gene level: no value matches only the defective scoring, and
none matches neither.

### Ten concrete instances

`R dist` is what ChIPseeker writes; `correct` is the same feature scored
against the peak it is actually listed under; `overlaps` names the peaks the
feature really intersects.

| peak | peak coords | transcript | transcript coords | R dist | correct | overlaps |
|---|---|---|---|---|---|---|
| 2 | `chr1:753301-753500` | `TxID:18228` | `chr1:756046-778647` | **0** | 25147 | 3 |
| 6 | `chr1:879301-879600` | `TxID:18310` | `chr1:868052-876876` | **0** | -2425 | 5 |
| 9 | `chr1:934301-937400` | `TxID:328` | `chr1:939275-944259` | **0** | -1875 | 10 |
| 10 | `chr1:940401-941100` | `TxID:320` | `chr1:925150-935793` | **0** | 15251 | 9 |
| 11 | `chr1:948301-949700` | `TxID:318` | `chr1:923923-944574` | **0** | 24378 | 9, 10 |
| 12 | `chr1:955001-957000` | `TxID:333` | `chr1:960584-965719` | **0** | -3584 | 13 |
| 13 | `chr1:960301-960600` | `TxID:18349` | `chr1:944203-959256` | **0** | -1045 | 11, 12 |
| 14 | `chr1:968201-968400` | `TxID:333` | `chr1:960584-965719` | **0** | 7617 | 13 |
| 16 | `chr1:1008901-1010100` | `TxID:18364` | `chr1:1004058-1006442` | **0** | -2459 | 15 |
| 17 | `chr1:1014801-1015700` | `TxID:348` | `chr1:1020120-1056116` | **0** | -4420 | 18 |

Peak 2 is the clearest: `TxID:18228` starts 2546 bp past the peak's end, so the
two do not intersect at all, yet the reference reports distance 0 because the
transcript also sits in peak 3's window and peak 3 genuinely overlaps it.

A minimal two-peak reproduction is reproducible without any TxDb:

```r
peaks  <- GRanges(c("c","c"), IRanges(c(100,1000), c(200,1100)))
genes  <- GRanges("c", IRanges(1050,1060), tx_name="F")
wide   <- peaks; start(ranges(wide)) <- start(ranges(wide)) - 5000
                 end(ranges(wide))   <- end(ranges(wide))   + 5000
hit    <- findOverlaps(wide, unstrand(genes))
featureHit <- genes[subjectHits(hit)]      # F appears twice, once per peak
overlapHit <- findOverlaps(peaks, unstrand(featureHit))
# queryHits(overlapHit) == 1 only: peak B is never hit
d <- rep(NA, 2); d[subjectHits(overlapHit)] <- 0
# d == c(0, 0) -- peak B is zeroed although F is 4050 bp away
```

### What EasyEpi does

`flank_info()` keeps `queryHits`, matching on the `(peak, feature)` pair:

```python
key_row = paste(peak_index, row_index)      # which row is which peak's element
key_hit = paste(queryHits, subjectHits)     # which pairs really overlap
rows in key_hit are scored 0
```

So `flank_txIds` and `flank_geneIds` are identical to ChIPseeker, while
`flank_gene_distances` carries the distance to the peak it is listed under.
Across all 12 datasets, at both levels, the id set matches on 23,646/23,646
peaks and every distance matches the corrected scoring; the peaks where the
output differs from the reference (5,155 transcript, 3,022 gene) are exactly
the ones the reference mislabels.

The per-element order inside each `;` list also differs: the reference emits
`getGene()` row order, which comes from the TxDb's SQLite query and is not
derivable from a GTF (see "The single remaining peak" above).

Both are recorded rather than reproduced. Reproducing the zeroing defect would
mean shipping "2546 bp away" as "distance 0" to downstream consumers.

## The other UTR path: explicit `UTR` rows

`gtf_features.py` can take UTRs two ways: from `five_prime_utr` / `three_prime_utr`
rows when the GTF has them, or, failing that, by subtracting CDS from exon
(`_derived_utr`). Every GTF in this benchmark until now carried **no UTR rows at
all**, so only the derived branch had ever run. That is the opposite of what an
earlier note here claimed.

GENCODE v50 (`/data1/.../reference/human/gencode.v50.primary_assembly.annotation.gtf.gz`)
does carry them: 1,518,503 `UTR` rows over 646,577 transcripts. A chr1 20 Mb
slice (`build_gencode_slice.sh`, 153,070 lines, 8,087 transcripts, 19,726 `UTR`
rows) is enough to exercise the branch:

| kind | python | R (TxDb) | |
|---|---|---|---|
| tss | 8,087 | 8,087 | 100% |
| exon | 70,854 | 70,854 | 100% |
| intron | 62,767 | 62,767 | 100% |
| utr3 | 11,554 | 11,554 | 100% |
| utr5 | 8,116 | 8,116 | 100% |

Two details are worth keeping:

* GENCODE tags these rows `UTR`, not `five_prime_utr` / `three_prime_utr`, so
  the parser does **not** read them: Python still reaches the derived branch.
  The table above therefore says the two derivations agree, not that the
  explicit branch is exercised by the parser.
* R's side reports utr3 + utr5 = 19,670 against 19,726 `UTR` rows. The 56-row
  difference is the stop codon, which `makeTxDbFromGFF` trims off the 3' UTR --
  the same trim `_trim_stop_codon` performs, which is why the two agree.

Only the slice was run. The full GENCODE file was not processed: R needs 17 s to
build a TxDb from 153k lines, so the 11.25M-line original is a scale test, not a
logic test, and the logic is what this table covers.

## Peaks on contigs the GTF never mentions

Real peak sets are not confined to chr1-22/X. MACS3 run genome-wide on four
ENCODE BAMs produced peaks on 41-55 seqlevels. That is a case the chr1-only
arms could not reach, so it was measured separately.

Input: the 282 peaks that fall on the 45 seqlevels the reference GTF does not
carry - 25 `chrUn_*`, 19 `*_random` and `chrM` - taken from all four
genome-wide MACS3 runs.

| | peaks in | rows out | lost |
|---|---|---|---|
| EasyEpi | 282 | **282** | 0 |
| ChIPseeker 1.46.1 | 282 | **208** | **74 (26.2%)** |

ChIPseeker does not annotate those peaks and does not label them
`Distal Intergenic` either; it **drops them**. Its own warning names the cause:

```
Each of the 2 combined objects has sequence levels not in the other:
  - in 'x': chrUn_KI270330v1, chrUn_KI270438v1, ...
```

`annotatePeak` intersects the peak set's seqlevels with the TxDb's and discards
whatever is not in the intersection. The loss is whole-contig, not
peak-by-peak: all 18 peaks on `chrUn_KI270438v1` vanished, as did all 16 on
`chr17_KI270729v1_random`. 16 of the 45 contigs in the input disappeared
entirely.

EasyEpi keeps them and reports `Distal Intergenic`, which is the honest answer:
there is no gene model on that contig.

### Why this matters beyond the missing rows

All 74 dropped peaks are `Distal Intergenic` in EasyEpi's output, so no
classification is *wrong* here - the data is simply gone. The danger is the
denominator. Any proportion computed from ChIPseeker's output (promoter share,
intron share) is over that 208-peak subset, not the 282 that were submitted.

The same mechanism fires on a much worse input: if the peak set names
chromosomes `1, 2, ...` while the TxDb names them `chr1, chr2, ...`, the
intersection is empty and **every peak is dropped** without an error. The
warning is printed, but the function still returns a valid, empty-ish object.

## `--keep-dups-all` and what duplicate filtering actually does

EasyEpi exposes `--keep-dup-all`, which maps to MACS3 `--keep-dup all`
(`peaks_cli.py`: `"auto" if not keep_dup_all else "all"`). Two measurements
were needed, because they disagree.

### Synthetic fixture: the flag does nothing, the coordinate stack does everything

A 260-read BAM with one real 200-read locus and one 60-read pile-up all at the
same coordinate, written twice - once with 59 reads flagged `0x400`, once with
none flagged:

| BAM | `--keep-dup` | peak produced |
|---|---|---|
| 59 reads flagged `0x400` | `auto` | `chr1:1000000-1000100` len=14 (the real locus) |
| 59 reads flagged `0x400` | `all` | `chr1:1002000-1002100` len=661 score=66.1 (the pile-up) |
| no flags | `auto` | `chr1:1000000-1000100` len=14 |
| no flags | `all` | `chr1:1002000-1002100` len=661 score=66.1 |

The two BAMs behave identically, so `0x400` is irrelevant. MACS3's `auto`
thresholds on the *number of reads sharing a coordinate*:

```
#1 calculate max duplicate tags in single position based on binomial distribution...
```

`all` skips that filter and turns the 60-deep pile-up into a strong false peak.

### Real ENCODE data: the flag changes almost nothing

`H3K4me3.bam`, 26,670,304 tags, `--nomodel --extsize 200 --qvalue 0.05`:

| | peaks | wall | tags after filtering |
|---|---|---|---|
| `auto` | 38,203 | 301 s | 26,670,279 (25 removed) |
| `all` | 38,203 | 311 s | no filtering line |

| comparison | result |
|---|---|
| peak count | identical |
| peak coordinates | **38,203 / 38,203 identical, zero in either-only** |
| peaks whose metrics differ | 3,137 (8.2%) |
| where they differ | the last 1-2 digits of `-log10(qvalue)` |

Examples:

```
chr11:129815817-129816405   auto q=8.63857   all q=8.63856
chr7:143620940-143621143    auto q=7.68728   all q=7.68727
chr16:13756654-13756890     auto q=4.16953   all q=4.16952
```

`auto` removed 25 tags out of 26.7M (0.0001%), because ENCODE alignments stack
shallowly: over the first 266,025 reads of this BAM, 6.15% of positions carry
more than one read and the deepest single position holds **3**. The binomial
threshold MACS3 derives admits depth 3, so essentially nothing is filtered.

### What the original assumption got right and wrong

The plan recorded "ENCODE filtered alignments are already deduplicated, so
`--keep-dup-all` makes no difference". Measured:

* *already deduplicated* - **false in form**. Every read in the two ENCODE BAMs
  checked has `0x400` clear: 0 of 2,000,000 flagged.
* *makes no difference* - **true in effect** on this data, but for a different
  reason: the coordinate stacks are too shallow to trip the threshold, not
  because the alignments arrive deduplicated.

The distinction matters downstream. On a library that is genuinely
duplicate-rich, or on unfiltered reads, `--keep-dup-all` will inflate peaks the
way the synthetic fixture shows - the flag is not inert, it is just inert on
shallow data.

## Limits

* chr1-only narrowPeak inputs; no broad-mark bin is exercised.
* UTR equivalence is verified on chr1 20 Mb of GENCODE v50, not on the
  full 11.25M-line file. The parser ignores GENCODE's `UTR` tag; the
  agreement shown is between two derivations, not a read of those rows.
* `flank_gene_distances` intentionally differs from ChIPseeker on the peaks
  affected by the reference's `subjectHits`/`queryHits` mix-up: the feature
  set and the ids match, the distances do not. See the `flank_*` section for
  the 12-dataset measurement and the minimal reproduction.
* One genome build (hg38) and one gene model.
* `5UTR`/`3UTR` are compared as the literal `5' UTR` / `3' UTR` strings ChIPseeker emits.
  No output file uses the old six-bin fold any more: `coarse_bin()` survives as a helper
  for callers that want it, but `annotation_stat.tsv` is the only distribution written.
* The repaired GTF is derived, not an upstream release -- `repair_gtf.py` and its md5 are
  the provenance chain.
* The `downstreamDistance > 1000` branch is covered by the parameter sweep (5500, 7500);
  `dsd` values that are not multiples of 1000 add a trailing-gap label, which the sweep
  exercises at 5500 and 7500.
* The off-contig arm uses 282 peaks from four genome-wide MACS3 runs; it tests
  contig handling, not genome-wide annotation throughput. The chr1-22 arms were
  deliberately not re-run at genome scale, since `annotate()` keys on `chrom`
  and the code path is identical.
* `--keep-dup-all` is compared on one BAM (`H3K4me3`) plus a synthetic fixture.
  The ENCODE BAMs are shallow enough that the flag is effectively inert; a
  duplicate-rich library is not represented.
* Peak calls used for the off-contig and keep-dup arms come from this machine's
  MACS3 3.0.4, not from the ENCODE pipeline that produced the 12-dataset peak
  sets, so the two sets are not directly comparable to each other.
