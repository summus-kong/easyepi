# Changelog

## 0.6.0 - unreleased

### Added

- **`easyepi concordance`**: replicate concordance between samples, as a
  subcommand. `--peaks` takes one peak file per sample and reports pairwise
  correlation *and* region-set Jaccard; `--matrix` takes the `diff --counts`
  layout and reports correlation only, because a count matrix has no genomic
  coordinates. `--design` is required.

  Two deliberate non-features: there is no pass/fail threshold and no
  `qc_verdict` block. `grouping_supported` compares mean within-group against
  mean between-group agreement and is a group-level statement -- measured on
  synthetic designs with real group structure, swapping one label in six still
  reports `True` -- so the command prints that caveat on every run. A value that
  could not be computed is `null` in the JSON and an empty field in the TSV,
  never `0`.

- **TSS enrichment** (`easyepi qc --tss`, `easyepi/tss_enrichment.py`): the
  ENCODE signal-to-background score at transcription start sites, piled up from
  read 5' ends in an apportioned window (`--tss-flank`, default 2000 bp).
  `--tss` accepts a GTF/GFF directly (`easyepi/annotation_tss.py` derives the
  positions at transcript or gene level) as well as a `chrom/tss` table, so a
  user does not have to build and maintain a second file that mirrors their
  annotation.

  Four states are reported rather than one nullable number: `measured`,
  `background_zero`, `no_reads` and `no_tss`. A sparse library whose flanking
  pileup is zero is *measured and failed*, not *never measured* -- those must not
  look alike downstream, so `qc --assay atac` scores that case as `0.0` (which
  can breach a floor threshold) instead of omitting the key.

- Replicate concordance module (`easyepi/concordance.py`) with `json_safe`, which
  converts NaN/inf to `null`. `json.dumps` otherwise emits a bare `NaN`, which is
  not JSON and makes a whole report unreadable to a strict parser.
- **`easyepi diff --method deseq2`**: negative-binomial count model via PyDESeq2
  0.5.4, replacing the rank/t-test path for designs where inference — not exploration
  — is wanted. It supports covariates (`--design-factors`), independent filtering, LFC
  shrinkage (`--lfc-shrink`) and records the model provenance (PyDESeq2 version,
  design formula, contrast, replicate counts, dispersion fit) beside the results table.

  The existing `mannwhitney` / `ttest` / `wilcoxon` methods are unchanged and remain
  the default; `--method deseq2` is opt-in.

### Fixed (found while building it)

- **PyDESeq2's `contrast` reports the reciprocal of what R DESeq2 does.** For
  `contrast=[factor, first, second]` it returns log2(second / first), so passing the
  operands in the natural reading order inverts every `log2FoldChange` while leaving
  p-values and significance untouched. Established by building data where the
  direction is known — group B doubled, group A halved — and reading off which operand
  order gives +1; the first version of the adapter reported ±the wrong sign and the
  test in `tests/test_deseq2_adapter.py` now pins the correct one.
- `ref_level` is deprecated in PyDESeq2 0.5 and has no effect. The adapter no longer
  passes it rather than exposing a parameter that is silently ignored.

### Input contract

`--method deseq2` refuses, with the reason rather than a rounded value: non-integer
counts (normalised signal), negative counts, NaN/inf, fewer than 4 samples, and fewer
than 2 replicates per group. A design factor with a single level is rejected too —
it cannot be estimated, and dropping it silently would hide a design mistake.

## 0.5.3

### Fixed

- **3 bp UTR slivers were emitted as real UTR intervals.** GENCODE defines some CDS
  blocks so a start or stop codon straddles an exon boundary, leaving a 1-3 bp piece
  outside the CDS. Building a TxDb normalises those away; deriving UTRs arithmetically
  does not, so the slivers became UTR features ChIPseeker has no counterpart for.

  Found by comparing transcript *sets* rather than peaks: 2,281 transcripts had a
  3'UTR here and not in ChIPseeker, and every one of them totalled at most 3 bp --
  2,281 in the 1-3 bp band, zero in 4-10, zero in 11-50, maximum 3 bp. Pieces shorter
  than 4 bp are now dropped, a threshold that clips no genuine UTR because nothing
  between 4 and 50 bp was ever affected.

  **Agreement with ChIPseeker goes 99.99% -> 99.9954%**; disagreeing peaks fall from
  4 to 2, and the 3'UTR interval count moves from 892,775 to 886,722 against
  ChIPseeker's 886,766.

- Corrects the 0.5.2 attribution, which called the residual 4 a data question. Three
  of them were caused by this, and the direction was determined by querying
  `threeUTRsByTranscript` at the exact peak coordinates: 0 hits, so ChIPseeker had
  UTR intervals this code lacked, not the reverse.

### Known gap

2 peaks (0.0046%) remain -- one `5'UTR -> Exon`, one `Exon -> 3'UTR`, one case in each
direction and no residual pattern. Matching them further would be fitting individual
peaks rather than correcting a rule.

## 0.5.2

### Fixed

- **3'UTR covered only the block abutting the stop codon.** `_extract_utr` kept the
  untranslated tail of the exon containing the stop codon and discarded every exon
  after it, on the reasoning that a downstream exon is "still an exon". That is true
  of the feature type but not of the annotation: a 36-exon gene routinely has several
  exons past translation end and every one of them is untranslated sequence.

  **Agreement with ChIPseeker goes 99.76% -> 99.99%**; disagreeing peaks fall from
  105 to 4. The UTR interval counts change accordingly: 3'UTR 617,876 -> 892,775 and
  5'UTR 487,261 -> 625,728.

- Corrects the 0.5.1 attribution. That release blamed the residual 105 on UTR
  provenance -- that ChIPseeker reads `threeUTRsByTranscript` while this derives UTRs
  arithmetically. Two checks disproved it: for a sampled transcript GENCODE's `UTR`
  rows and `exons - CDS` are the identical interval set, and each failing peak
  overlapped an Exon interval and zero UTR intervals. The sets disagreed about which
  regions exist, not about how to label them.

Three tests that had encoded the old behaviour were corrected, not deleted. They were
written from the same incorrect reasoning as the code and had been green against it.

### Known gap

4 peaks (0.01%) remain, all `Exon -> 3'UTR`, where a transcript's CDS status differs
between annotation builds. That is a data question rather than a rule question.

## 0.5.1

### Fixed

- **The overlap sweep searched an unsorted list.** Features were sorted by start, then
  the *end* list was passed to `bisect.bisect_left` to find where a peak's scan should
  begin. `ends` follows the start order and is therefore not sorted, so the returned
  position was arbitrary and overlapping features were skipped -- dropping peaks from
  every bin without any error. The scan is now bounded by `bisect_left(starts, qe)`,
  which is sound because start is the sort key, and pruned with a running maximum of
  ends.

  **Agreement with ChIPseeker goes 93.78% -> 99.76%**; the number of disagreeing peaks
  falls from 2,703 to 105. Three backends remain identical to each other at 100.0000%.

- Corrects the 0.5.0 attribution. That release blamed the residual 6.22% on intron
  modelling; every one of the 1,786 `Intron -> Intergenic` peaks actually does overlap
  an intron interval this code builds, and the cause was the scan, not the features.

### Known gap

105 peaks (0.24%) still differ, all UTR-versus-Exon boundaries: 3'UTR->Exon 91,
5'UTR->Exon 10, Exon->3'UTR 4. ChIPseeker reads `fiveUTRsByTranscript` and
`threeUTRsByTranscript` from the TxDb; this code derives UTRs by subtracting CDS from
exons, and the two disagree for transcripts with a gapped CDS or a retained intron.

## 0.5.0

### Fixed

- **Promoter is now decided by distance to the nearest TSS, not by whether the peak
  midpoint falls inside a +/-window interval.** ChIPseeker computes the distance from
  the peak edge to the nearest resized TSS and marks `tssRegion[1] <= distance <=
  tssRegion[2]` as Promoter, applied after the interval bins so it overrides them.
  Testing that rule against ChIPseeker's own `distanceToTSS` column matched on
  43,465/43,465 peaks with zero exceptions. Adopting it recovers 5,731 peaks that had
  been missed -- all within the window by distance but outside it by midpoint.
- The interval bins classify on **whole-peak overlap**, matching ChIPseeker's
  `findOverlaps`, where a peak whose body reaches into a region counts.

**Agreement with ChIPseeker rises from 80.39% to 93.78%**, crossing the 93% threshold
set for this work. The Promoter bin goes from 22,282 to 27,623 peaks against
ChIPseeker's 28,414. All three backends (`native`, `pybedtools`, `pyranges`) still
agree 100.0000% with each other.

### Known gap

6.22% of peaks still differ, dominated by `Intron -> Intergenic` (1,786 peaks, 4.11%).
Those peaks are a median 13.7 kb from the nearest TSS, so they lie inside a transcript
ChIPseeker modelled and this implementation did not: ChIPseeker builds introns from
its `intronsByTranscript` table rather than from gaps between the exon records read
here, and 4.5% of GENCODE transcripts have a single exon. Closing it means deriving
introns from transcript spans.

## 0.4.0

### Changed

- **One TSS per transcript instead of one per gene.** `_read_gtf_txs` keyed its table
  on `gene_id` and used the `gene` line's span, so every gene contributed a single
  TSS. ChIPseeker derives a TSS from every transcript, and on a 78,941-gene /
  646,577-transcript annotation that difference dominated the output: the Promoter
  bin was 29% instead of 65% and agreement with ChIPseeker was 58%. Records are now
  emitted per transcript, each carrying `transcript_id`, the gene id and the gene
  span, so downstream output can still report a gene.

  The parser also stopped relying on `(?:gene_id|transcript_id)`, which matched
  `gene_id` first and therefore attached every exon and CDS to the gene rather than
  the transcript.

  **Peak annotations change again.** On 43,465 H3K4me3 peaks the bins move from
  Promoter 29.2% / Exon 12.6% / Intron 33.6% / Intergenic 19.7% to
  Promoter 51.3% / Exon 5.3% / Intron 23.8% / Intergenic 16.3%.

### Verified

- ChIPseeker agreement rises from **57.66% to 80.39%** exact and from **59.19% to
  81.61%** coarse, with the largest per-bin delta falling from 36.2 to 14.1
  percentage points. Still below the 93% / 98% thresholds, so this version is
  **still not ChIPseeker-concordant**; the residual is dominated by
  `Promoter -> Intron/Exon` cases where the two disagree about which transcript
  starts nearest a peak.
- ChIPseeker's `Downstream` bin holds 19 of 43,465 peaks (0.04%), so folding it into
  `Intergenic` is confirmed harmless rather than assumed.
- **All three annotation backends agree exactly.** `native`, `pybedtools` (0.12.1,
  bedtools 2.31.1) and `pyranges` (0.1.7) classify all 43,465 peaks identically --
  100.0000% pairwise. Timing is comparable (44 / 53 / 48 s including GTF parsing),
  so the backend choice is about available tooling, not correctness.

## 0.3.0

### Fixed

- **The 5'UTR and 3'UTR bins could never be filled.** `GENOMIC_BINS` listed `Exon`
  before the UTRs, and the classifier returns the first bin whose interval contains
  the peak midpoint. An exon interval by definition contains its own UTRs, so every
  UTR peak was captured by the enclosing `Exon` and both UTR bins reported zero on
  every dataset. The priority is now `Promoter > 5'UTR > 3'UTR > Exon > Intron >
  Intergenic`, matching ChIPseeker. **Peak annotations change**: UTR counts rise from
  zero and `Exon` falls accordingly.
- **3'UTR swallowed every exon after the coding region.** `_extract_utr` compared
  each exon against the global `min(cds_start)` and `max(cds_end)`, so any exon
  between the first and last CDS block counted as 3'UTR. On a GENCODE model that is
  most internal exons. A UTR is now "exonic and non-coding": the piece adjacent to
  the coding region, with later exons of the same transcript reported as `Exon`, as
  ChIPseeker does.
- **Minus-strand UTRs were on the wrong side.** 5' and 3' were assigned from genomic
  coordinate alone, so every gene transcribed from the minus strand had its UTRs
  swapped. Roughly half a mammalian genome is minus-strand, so this was wrong in a
  way that partly averaged out instead of looking obviously broken.
- **TSS distance ignored strand.** `_nearest_tss_distance` returned `midpoint - tss`
  regardless of orientation, so a peak upstream of a minus-strand TSS got the wrong
  sign. It is now signed "negative upstream" on both strands.

### Verified — and found wanting

- **ChIPseeker benchmark run and failed.** 43,465 H3K4me3 peaks, one TxDb built from
  the same GENCODE v50 GTF for both sides, same 3000 bp window: per-peak exact
  agreement **57.66%**, coarse agreement **59.19%**, largest per-bin delta **36.2 pp**.
  Against thresholds of 93% / 98% / 2 pp this version is **not** ChIPseeker-concordant.
- **Cause found, fix deliberately deferred.** `_read_gtf_txs` keys on `gene_id` and
  takes the `gene` line's span, so each gene contributes one TSS; ChIPseeker derives a
  TSS from every transcript. The easyepi Promoter set is a strict subset of
  ChIPseeker's (0 disagreements in the other direction) and the two TSS-distance
  columns are uncorrelated (r = 0.007). Deriving one TSS per transcript roughly
  doubles the Promoter fraction, so it changes every annotation and needs its own
  before/after record — it is not part of this release.

### Added

- `annotate --promoter-upstream` / `--promoter-downstream`, default 3000 each (the
  ChIPseeker `tssRegion` default). The window was hardcoded.
- `annotate` writes `annotation_bins.tsv` beside the per-peak table;
  `--no-bins-summary` suppresses it.
- `tests/test_annotate_bins.py`, 19 tests over a hand-computed gene model. The
  classification logic previously had no coverage: only `_extract_utr` and
  `_nearest_tss_distance` were referenced anywhere in the suite.

## 0.2.0

- `easyepi[activity]` now works on Python 3.12 as well: `easyepi/activity.py` drives
  both decoupler generations — 1.9.x through `dc.run_ulm` and 2.x through
  `dc.mt.ulm`. The two APIs disagree on three things the adapter now handles: the
  entry point, the matrix orientation (2.x takes samples-by-features and 2.x's
  `mt/_run.py` already applies Benjamini-Hochberg, so `padj` is no longer corrected
  twice), and the return shape for AnnData input. Verified against decoupler 2.1.6:
  its t-values match 1.9.2 to 7.5e-07 on identical input.
- `easyepi[activity]` also declares `anndata`, which decoupler imports while loading;
  without it the extra installed and then raised ModuleNotFoundError: anndata.
- Figures are written as **PDF (vector) by default** and are produced without
  being asked: a run that should have plotted no longer looks identical to one that
  did. `--figure-format` still selects png/svg/tiff/jpg, `--no-figures` writes only
  tables and JSON, and every plotting entry point (`save_figure`, `pubfig`,
  `motifs.logo`, `viz`) reads one `DEFAULT_FIGURE_FORMAT`.
- Fixed figures silently disappearing on the motif-enrichment path: `np_log10`
  returned a bare ndarray and the following `.iloc` calls raised inside the
  best-effort render wrapper, which logged a warning and produced no file.
- The assay vocabulary is defined once (`easyepi/assays.py`) and read by `bam2bw`,
  `peaks`, `qc`, `assay_metrics` and `assay_commands`. `--assay` now accepts 31
  spellings everywhere instead of two different subsets, and `qc --assay` validates
  its value instead of silently emitting an empty metric block.
- Methylation QC reads beta values from the bisulfite BAM (MM/ML or Bismark XM), as
  the WGBS page always claimed; `qc --assay methylation --bam x.bam` previously
  reported read flags only.
- `differential-peaks --peak-dir` finds the per-sample peak layout that `easyepi
  peaks` writes, as well as the flat layout.
- Documentation: every one of the 81 `easyepi ...` examples in `docs/` and
  `README.md` is parsed by the live argparse tree in CI, and the four reference
  pages that existed but were unreachable from the nav are now published.
- Added `easyepi tf-annotate`: annotate peaks with transcription factors by scanning
  known PWMs (JASPAR/HOMER/MEME) in peak sequences fetched from an indexed genome
  FASTA. Thresholds come from the PWM score distribution under a binned uniform
  background (the FIMO construction), not from a fitted or invented cutoff.
- Fixed `Motif.to_log_odds` turning the whole PWM into NaN when a column had zero
  counts, which made the motif silently report "no site" everywhere.
- Fixed `footprint` hiding its significance test: the CLI now reports the
  target-vs-background Mann-Whitney and the centre-vs-shoulder Wilcoxon p-values,
  writes them to `footprint.json`, and errors when no target site falls in the window.
- Added `tests/conftest.py` selecting the Agg backend so plotting tests no longer
  depend on a display being present.
- Per-motif site p-values are now computed from one precomputed tail table instead of
  recomputing the convolution once per hit.
- Motifs too short to reach the requested per-sequence FPR are reported in
  `skipped_motifs` instead of appearing as "no binding site".
- Added `easyepi histone-states` and `easyepi differential-peaks` CLI entries for the
  pure-Python ChromHMM and DiffBind equivalents.

- Added `easyepi nucleosome-positioning`: nucleosome calls from paired-end fragment

  dyad density, with robust (median + MAD) local-maximum calling, per-call

  fuzziness and neighbour spacing. Single-end BAMs are refused explicitly, because

  a single-end read reports one cut site and no fragment span, so the dyad is

  unidentifiable.

- Added `easyepi viz.genomeviews`: MA plot, peak-set overlap matrix, TSS heatmap and

  stacked locus track. All four return Figures and export only when `save` is given.

- Peak overlap defaults to a base-pair Jaccard index computed on merged intervals.

  An interval-count metric reports a broad domain as fully covering many narrow

  peaks, which overstates similarity; the bp metric is symmetric and width-robust.

- Fixed `overlap_fraction` clipping only one interval end, which reported a

  candidate's whole length as shared for a partial overlap.

- Added `easyepi gene-panel`: one column per gene, one row per BigWig, with the
  gene model (exons and strand arrow) underneath. Reads gene bodies and exons from
  a GTF/GFF, or takes a BED of loci.
- Added `easyepi.viz.panels`: gene panels, UpSet intersection plots and a
  hierarchical-clustering sample-correlation heatmap.
- Added `easyepi.gc_background`: GC- and width-matched control windows from the
  genome FASTA, with interval-distance exclusion and per-candidate use-once
  sampling. Unmatched targets are reported rather than dropped.
- Peak-set overlap now defaults to a base-pair Jaccard index on merged intervals.

Defects fixed:

- UpSet membership tested whole intervals instead of elementary union segments, so
  a partially overlapping pair mislabelled the non-shared remainder and lost
  regions unique to the second set.
- GC matching sampled randomly instead of enumerating eligible windows, which
  reported a false "no match" when the only valid window was rare.
- GC matching compared absolute coordinates for the distance test, so any two
  distant windows counted as far apart and `exclude` did nothing at `min_distance=0`.
- `exclude` regions are now a hard exclusion independent of `min_distance`.
- Sample correlation left a NaN row/column when a sample had no variance.
- Gene panels had no genomic tick labels, and a single shared y limit per row
  flattened low-signal loci; limits are now robust per panel unless shared.

- `gene-panel`, `histone-states` and `integrate` now accept comma-separated
  values (e.g. `--bigwig H3K4me3,H3K27ac`, `--gene AGRN,TP73`), and `gene-panel`
  gained `--bigwig-dir` so bare labels resolve to `<NAME>.bw`/`<NAME>.bigWig`.
  Repeated flags keep working. Ambiguous or missing labels raise instead of
  producing a silently short track list.
- `histone-states` exposes `--n-init` and records the restart count and final
  log-likelihood in `states.json`; the option existed in the learner but was not
  reachable from the CLI.

- `differential-peaks --peaks` accepts three forms: `SAMPLE=BAM:PEAK`, `SAMPLE=BAM`
  (peak file found via `--peak-dir`) and bare `SAMPLE` (both found via `--bam-dir`
  and `--peak-dir`), with comma-separated labels. Unresolvable or ambiguous samples
  raise instead of silently changing the replicate count.
- `qc` gained `--input-dir`, which scans one folder and groups `*.bam`, `*.bw`,
  `*.bigWig`, peak and methylation files into samples by name, and both `--bam`/
  `--bigwig`/`--peak` now accept comma-separated lists.

Defects fixed:

- `qc` paired `--bigwig[i]`/`--peak[i]` with `--bam[i]` positionally, so a sample
  without a track, or a different flag order, silently attached the wrong file.
  Tracks are now matched to samples by name.
- `qc` wrote a report with an empty `samples` object and reported success when no
  input was recognised; it now fails with the paths it looked for.

- FRiP now reports `interpretable: false` with an explanatory warning when the peak
  set spans a small fraction of the mapped genome, because a regional numerator
  over a genome-wide denominator is not a FRiP. The peak span fraction and the
  read universe size are recorded alongside it.
- `maintenance/verify_linux.py` now checks the subcommands added since the last
  release (`tf-annotate`, `nucleosome-positioning`, `gene-panel`,
  `differential-peaks`, `histone-states`, `footprint`, `qc`, `integrate`).
- Added `scripts/verify_macs3_concordance.py`: byte-level comparison of the MACS3
  wrapper against a direct `macs3 callpeak` invocation.

## 0.1.0 release candidate - 2026-10-16

- Flattened the maintained Python module layout into universal and assay-specific entries.
- Standardized active source, comments, CLI messages and documentation in English.
- Replaced misleading peak adapters with official MACS3 commands and corrected ChromHMM stage inputs.
- Added decoupler ULM activity on measured matrices and Biopython-thresholded known motif enrichment.
- Corrected Tn5 negative-strand coordinates, consensus overlap support, differential pairing, BH missingness and methylation validation.
- Replaced invented normalization defaults and hand-written k-means with explicit factors and scikit-learn.
- Made BAM/BigWig and other specialist dependencies optional; no R backend is used.
- Disabled placeholder or unrelated implementations; documented unsupported MNase positioning and missing mature count/DMR models.
- Made figure exports opt-in and preserved public Figure-returning functions.
- Added regression tests, English usage/API/development documentation and a reliability audit.
- Removed the historical build, cache and output directories from the working tree; distribution manifests continue to exclude them.
- Added `[project.urls]` metadata and aligned the declared dependency extras with what the code actually imports.
- Resolved the `plot_heatmap` name collision: the table-based function is exported from `easyepi` / `easyepi.plot`, and the array-based one is now `plot_matrix_heatmap` (the old `easyepi.viz.plot_heatmap` remains as a deprecated alias).
- Fixed `_nearest_tss_distance` returning a global TSS index for a per-chromosome list, which selected the wrong gene on multi-chromosome input.
- Fixed `easyepi.specific.mnase.profile` failing to serialise the array-valued `occupancy` field to JSON.
