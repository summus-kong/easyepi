# Changelog

## 0.6.0 - unreleased

### Added (reproducibility)

- `easyepi/provenance.py`: input and output SHA-256, direct dependency
  versions and external tool versions, written into `status.json` under
  `provenance`. Tools resolve through `PATH` and then through explicit
  `search_paths`, so a conda env that the run does not activate still yields a
  version instead of a false `not_found`. Output hashing is opt-in
  (`RunRecorder.hash_outputs()`) because it costs real IO on a large BAM.
- A missing input, an absent tool and an unreadable path each record a reason
  (`not_found`, `skipped_too_large`, `not_a_file`, `unreadable`) instead of
  raising, so a manifest can be written after the fact.

### Added

- **Library metrics: duplication, complexity, mitochondrial content.**
  `assay_metrics.library_summary` reports `duplicate_fraction` with a
  `duplicate_status` that distinguishes a measured zero from a BAM that was never
  duplicate-marked, `library_complexity` (distinct templates over reads inspected),
  and `mitochondrial_fraction` gated on a new `--mito-pattern`. Registered as the
  `library` group for `atac-qc`, `chip-qc`, `histone-qc` and `mnase-qc`.
  The mitochondrial count is taken through the BAM index over the whole file: a
  scan bounded by `--max-reads` stops long before chrM, which is the last contig of
  a coordinate-sorted BAM, and the earlier scan-based version therefore reported
  `0.0` for a contaminated library. `mitochondrial_reads`,
  `mitochondrial_mapped_records` and `mitochondrial_denominator` are all reported so
  the fraction can be checked against the numbers it came from.
- **`--mito-pattern` and `--tss`/`--tss-flank` reach every QC entry point.**
  The five assay aliases build their own argument parser, so a flag registered only
  on `qc` was silently absent there and the caller's `getattr(..., None)` turned the
  absence into a default instead of an error. Both flag groups now come from shared
  helpers used by all six entry points, and
  `tests/test_qc_entry_point_parity.py` fails if they drift apart again.

- **`easyepi timecourse` ranks features by time-point specificity (TAU).**
  Answers two separate questions and keeps them apart: how *concentrated* a
  feature's profile is on one time point (TAU, Yanai et al. 2005) and whether it
  *trends* monotonically over time (Spearman). TAU is descriptive and carries no
  p-value anywhere, because treating it as a test statistic is the usual way the
  number gets misread. Counts are collapsed to per-time means before TAU so
  uneven replicate numbers do not weight a time point above its measurement, and
  the wide table writes **one `spec_<stage>` column per time point** beside the
  global `tau`; `max(spec) == tau` by construction, so the wide table, the
  `<out>.by_stage.tsv` long form and the JSON summary cannot tell contradictory
  stories. An all-zero feature reports a blank, not `0` -- "not measurable" is
  not "measured, no specificity" -- and it no longer aborts the run: the
  replicate diagnostic collected `tau_specificity`'s `None` and then called
  `min()` on it, so one silent feature killed the whole run with a `TypeError`
  before the blank row could be written. Three preconditions are enforced rather
  than assumed: non-negative input (an absolute value would manufacture
  specificity), at least three distinct time points (with two, tau is a
  deterministic function of the ratio and cannot rank), and a design whose
  `sample_id`/`time` match the matrix (time points are never inferred from
  column order). `tests/test_timecourse_stage_columns.py` pins the column count
  at 3, 4, 5 and 7 stages.

- **QC figures now include a TSS metagene, not just the enrichment score.**
  `tss_enrichment` had always returned the full per-base pileup in `profile`, but
  no figure read it: `_assay_metric_frame` keeps scalars plus two special-cased
  lists, so the nested `tss_enrichment_detail` was dropped and
  `plot_tss_enrichment` drew only the scalar ratio. The score says *whether* the
  signal is TSS-enriched; the shape says whether that enrichment is a sharp
  TSS-proximal peak over a flat background (a good ATAC library) or a broad hump
  that happens to peak (fragmented or over-tagmented DNA), and two samples with
  the same score can differ entirely in shape. `plot_tss_profile` draws one curve
  per sample with its x axis derived from `profile_center_index`, so a peak at the
  TSS lands at 0 rather than at whatever the list midpoint happens to be.
  `profile` keeps its documented length of `2*flank+1`: an earlier attempt binned
  it to shrink the JSON and was reverted, because `tests/test_tss_enrichment.py`
  asserts that length and a resampled series would have silently moved the
  TSS. `tests/test_tss_profile_figure.py` pins the plumbing, the length and the
  axis. Verified end to end on a synthetic ATAC fixture (score 50.57, profile
  4001 values, `tss_profile.png` written): 960 tests pass before and after.

- **`easyepi peaks` is MACS3, and that is now pinned rather than assumed.**
  The command imports `MACS3.Commands.callpeak_cmd` and calls `run()` on a
  Namespace built from MACS3's own argument parser, so "output equivalence" is
  the definition rather than a property to measure. What is measured instead:
  the run emits MACS3's own peak file (coordinates identical to a direct
  `macs3 callpeak`), `--macs-arg` reaches the parsed Namespace, a rejected option
  raises naming its tokens, and a missing parser refuses rather than accepting
  tokens it cannot apply. An earlier changelog entry here described a
  hand-written Poisson caller compared against MACS3 as a peer; that caller had
  already been retired and the comparison was withdrawn. See
  `docs/dev/TOOL_COMPARISONS.md` and `tests/test_peaks_is_macs3.py`.


- **`--macs-arg` can no longer be silently dropped.** `build_namespace` builds
  MACS3's argv and parses it with MACS3's own argparse, which supports every
  option; when that parser cannot be read it falls back to a hand-written
  Namespace with a fixed attribute set. Extra tokens used to be appended to the
  argv, recorded in the provenance, and then ignored -- the run reported a
  passthrough that never reached the calculation. It now raises, naming the
  tokens and the reason, because refusing is the only outcome that does not
  misreport what ran. `tests/test_peaks_is_macs3.py` also pins, end to end, that
  a real `easyepi peaks` run writes the same peak coordinates as a direct
  `macs3 callpeak` invocation with the same options -- the command is MACS3, so
  agreement is a definition rather than a coincidence.
- **`run_tobias_footprint` named a subcommand TOBIAS does not have.** It built
  `TOBIAS Footprint --cores 4`; 0.17.5 rejects `Footprint` at argument parsing
  and `FootprintScores --help` states that `ScoreBigwig` replaces it, with usage
  `ScoreBigwig --signal <bw> --regions <bed> --output <bw>` and no `--cores`.
  Two defects in one call, and nothing invoked it, so it could never have worked.
  Now `ScoreBigwig` with the three documented options;
  `tests/test_tobias_command.py` pins the subcommand against TOBIAS's accepted
  list, the absence of `--cores`, and invokes the real binary.
- **`easyepi.motifs.scan_sequences` agrees with FIMO on where the site is.** On a
  12-sequence fixture (planted consensus, its reverse complement, and noise) the
  best window matches in 11/12 cases and the ranking of best-window scores is
  Spearman rho = 1.000; the one disagreement is a noise-only sequence within 1 bp
  of a tie. The absolute scores are deliberately not compared -- on four fixed
  7-mers the difference changes sign with base composition (GATATCA +0.40,
  GGGGGGG -1.14), and a constant offset, a pseudocount and background smoothing
  were each ruled out, so FIMO's `score` is its own statistic. No downstream
  decision uses it: a motif call comes from a threshold or a ranking.
  `tests/test_motifs_fimo_agreement.py` and `tests/oracle/compare_to_fimo.py`.
  Details, including the two harness bugs that would have made the comparison
  vacuous: `docs/dev/TOOL_COMPARISONS.md`.

### Changed

- **Coverage metadata names the counting unit.** `bam2bw` now writes a
  `counting` block into the `.bw.json` with `unit`, `fragments`,
  `unpaired_reads` and per-chromosome counts, so a fragment-counted track can
  be told apart from a read-counted one and the two are not compared. The unit
  reports what was counted, not what the BAM header claims: a file whose
  records carry no mate reports `read`. The cut-site path (`--assay atac`)
  counts read ends, so it writes no block at all rather than naming a unit it
  did not use.
- **`easyepi dmr` no longer needs R.** The per-site test for designs with two
  or more replicates per group is now a pure numpy/scipy reimplementation of
  `DSS::DMLtest(smoothing = TRUE)`, including DSS's shrunken dispersion
  (`est.prior.BSseq.logN` + `dispersion.shrinkage.BSseq`). Verified against
  DSS 2.58.0 intermediate by intermediate on a 300-site, two-chromosome,
  irregularly spaced fixture: smoothed levels 5e-16, weighted levels 7.5e-12,
  `phi` 1.8e-11, variances 1.3e-13, `diff` 7.5e-12, `se` 1.5e-12, `pval`
  1.1e-10, and the identical set of sites called at p < 0.05. The oracle
  (`tests/oracle/dss_oracle.R`) asserts internally that its own step-by-step
  chain equals `DMLtest` to 1e-12 before any comparison is made.
  `easyepi/dmr_r.py` is kept and still correct, but nothing on the product
  path imports it. Details and the six transcription traps:
  `docs/dev/DSS_REPLACEMENT_PLAN.md`.

- **The variance kernel is 2.7x faster, with identical numbers.** Profiling put
  `rho ** lag` at 60% of the per-site kernel body (96.4 of 116.6 us/site at ~101
  CpG per 500 bp window) for an exponent that only takes the integers
  `0 .. span//50 + 1`. A lookup table replaces it exactly -- verified with
  `array_equal`, not a tolerance. `nitem_bin` was still a Python two-pointer
  loop; sorted positions make it two binary searches, 20.6 ms -> 1.0 ms per 20k
  sites. Total: 379 us/site -> 139 us/site, i.e. 5.9 h -> 2.2 h extrapolated for
  a 28M-CpG genome against two groups. All five DSS oracle comparisons were
  re-run afterwards and every intermediate is unchanged. Numbers and the stage
  breakdown: `docs/dev/DSS_REPLACEMENT_PLAN.md`.

- **`easyepi dmr --no-smoothing` exposes DSS's unsmoothed site test.**
  `DMLtest(smoothing = FALSE)` is a separate code path in DSS, not a span of
  zero: `compute.mean.noSmooth` drops the window but keeps the pseudocount,
  `compute.waldStat.noSmooth` divides per-replicate sums by `rowSums(n)^2`,
  and the variance floor is `1e-05` one-sided instead of `1e-04` two-sided.
  Both floors are kept as separate functions rather than one shared helper,
  because reusing the smoothed one would change every low-coverage p-value.
  The flag only reaches the beta-binomial branch: with fewer than two
  replicates per group the pooled test still runs and still says so.

- **All five `DMLtest` modes are now compared against DSS numerically**, not
  just the default one: `smoothing` x `equal.disp` (four combinations) plus
  the one-replicate case that `DSS:::DMLtest.noSmooth` resolves by sharing a
  dispersion. Every mode agrees on every intermediate -- smoothed and
  unsmoothed levels 5e-16, `phi` 1.9e-11, variances 1.5e-12, `se` 1.5e-11,
  `pval` 1.2e-10 -- and on the set of sites called at p < 0.05.
  `tests/oracle/compare_modes.py` prints the whole table.

- **`DMLtest` and `DSS:::DMLtest.noSmooth` disagree about legal inputs.**
  The public function rejects `smoothing = FALSE` with one replicate per
  group ("There is no biological replicates in at least one condition"),
  while the internal one handles it. The oracle takes a `direct` route to
  reach that branch; Python reports `equal_disp` as what was used, not as
  what the caller asked for.
- **The three DSS smoothing kernels are reimplemented in numpy**, so the
  package compiles nothing at install or run time and no longer carries a
  `ctypes` build probe or a C source file. `tests/test_dss_kernels.py` checks
  each one against a literal transcription of the C loops across window widths
  1-500 and uniform, irregular and duplicate site spacing. One behavioural
  difference is real and preserved: `compute_var_smooth.c` treats the window's
  right edge as exclusive while `filter.c` treats the same half-widths as
  inclusive.

- **The replicate-aware DMR path now passes real genomic coordinates.** The
  old R bridge handed DSS a fabricated `chr1` with positions `1..n_sites`,
  which put every smoothing window on a one-base-per-site axis; `hasCoverage`
  and the window kernel are both genomic, so the fabricated axis changed which
  sites were tested. It also applies DSS's own `hasCoverage` filter (both
  groups, `ws = 200`, threshold `2 * n_samples`) before computing anything,
  because the dropped sites move every later window. Dropped sites come back
  as `NaN` with `tested == False` and are excluded before BH correction, and
  `stats["dropped_by_dss_coverage_filter"]` counts them.

- **`annotated.tsv` gained `flank_txIds` / `flank_geneIds` / `flank_gene_distances`**
  behind `--add-flank-info` and `--flank-distance` (default 5000), mirroring
  `annotatePeak(addFlankGeneInfo=TRUE)`. The id lists match ChIPseeker; the
  distances deliberately do not. The reference writes
  `hitInfo$distance[subjectHits(overlapHit)] <- 0`, dropping `queryHits`, so a
  feature that sits in several peaks' windows has its overlap attributed to the
  wrong peak's row. On 3003 real peaks that mislabels 11% of the entries and
  touches 24.7% of the peaks (measured by `benchmark/a10/scripts/flank_defect.R`;
  present in 1.46.1 and unchanged in 1.49.5). EasyEpi keeps `queryHits` and
  reports the distance to the peak the feature is actually listed under, rather
  than reproducing a wrong answer.

- **`annotated.tsv` now carries annotatePeak's full column set.** 14 columns at
  `--level transcript` (`seqnames,start,end,width,strand,annotation,geneChr,
  geneStart,geneEnd,geneLength,geneStrand,geneId,transcriptId,distanceToTSS`)
  and 13 at `--level gene`, where `transcriptId` does not exist. Column names
  follow the pandas-safe spellings (`chrom` for `seqnames`, `tss_distance` for
  `distanceToTSS`) and the values are the reference's own.

- **`start` is now the 1-based GRanges coordinate**, not BED's 0-based offset:
  `annotatePeak` writes `start(peak)`, so a BED row beginning at 713700 appears
  as 713701. `width` stays `end - start`.
- **`easyepi annotate` now reproduces ChIPseeker's `sameStrand` and `level`.**
  `--same-strand` and `--level {transcript,gene}` mirror `annotatePeak`'s
  arguments. Verified against ChIPseeker 1.46.1 on a 400-peak set carrying
  real `+`/`-` strands: all four (level, sameStrand) combinations agree on
  400/400 labels and 400/400 `distanceToTSS`, end to end through the CLI.

- **Strand is read from BED column 6 only when it is `+` or `-`.** In a
  narrowPeak that column holds `signalValue`, so reading it positionally
  turns every peak into `*` -- which is exactly what ChIPseeker's own
  `readPeakFile` does, and why the strand paths went untested until now.

- **`easyepi annotate` now emits ChIPseeker's exact labels and feature table.**
  `annotated.tsv` carries the reference strings verbatim --
  `Promoter (<=1kb)`, `Exon (T1/G1, exon 2 of 4)`, `3' UTR`,
  `Downstream (<=300bp)` -- and the new `annotation_stat.tsv` is a port of
  `getAnnoStat(annotatePeak(...))`: percentages over the reference's own
  category set (Promoter sub-categories kept, exons/introns split into
  1st/Other, `Downstream` its own row). Verified byte-identical to ChIPseeker
  1.46.1 on 600 peaks at +/-3 kb and +/-5 kb.

  `annotation_bins.tsv` (six coarse bins, fractions) is **removed**: keeping a
  second, non-reference vocabulary alongside the exact one invited
  disagreement. Labels and distances match ChIPseeker on 100% of 23,646 peaks
  across 12 data sets.

- `annotate`'s `--promoter-upstream`, `--promoter-downstream` and
  `--no-bins-summary` were parsed and then ignored; they now reach the
  classifier. `--downstream-distance` was added.

- `annotate` gains `--annotated-name` and `--stat-name` so the two output
  tables can be named for downstream joins.

### Added (repository completeness and tooling)

- **Six modules the package imported were never tracked.** `annotate/__init__.py`
  imports `chipseeker_rules` and `tx_features`, `dmr.py` imports
  `dmr_replicate_aware`, and `cli.py` imports `qc_html`; none of them were in the
  repository, so a clone failed at import while the README described code that
  could not be checked out. `gtf_reader` and `dmr_r` are tracked as well even
  though nothing imports them today -- `gtf_reader` records why the transcript
  reader was rewritten once already: the GTF-only version ran 740 s on GENCODE
  v50 and silently turned every GFF3 exon into its own transcript, which made
  Intron and 5'UTR come out as exactly zero.

- **`scripts/txdb_to_gtf.R`** exports a TxDb with one consistent transcript-id
  namespace. `rtracklayer::export()` writes the TxDb's own tables and does not
  expand the per-transcript GRangesLists, so the result has exon and CDS rows but
  zero UTR rows, and a peak inside a UTR can only fall through to Promoter or
  Intron. It also leaves two namespaces in one file -- `transcripts()` returns an
  unnamed GRanges while every `*ByTranscript()` accessor names its elements by
  internal numeric tx_id -- so a naive export writes transcript rows as ENST...
  and exon rows as "1", where nothing links up. That is the failure the
  annotation benchmark inherited: zero linkable exon ids across 412,044
  transcripts. The exporter keys everything by tx_name.

- **ENCODE acquisition scripts** (`scripts/encode_*.py`, `fetch_small_bams.py`,
  `verify_pilot_bams.py`) and `benchmark/annotation/` for the R side of the
  annotation comparison. `fetch_small_bams.py` retries per attempt so one dropped
  connection does not lose the whole set; `verify_pilot_bams.py` runs before any
  analysis touches the files.

- **Invariant tests that were only described in comments**:
  `test_coordinate_contract.py` (which convention each layer uses, 1-based closed
  vs 0-based half-open), `test_boundary_oracles.py` (where the promoter window is
  clamped), `test_gff3_multi_parent.py` (a GFF3 child with several parents reaches
  all of them), `test_histone_state_window.py`, `test_dmr_replicate_aware.py`.

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

- **A single-file HTML QC report** (`easyepi/qc_html.py`, written automatically
  beside the JSON; `--no-html` turns it off). Figures are inlined as base64 so the
  report has no external references and survives being emailed or reopened months
  later. Every `<img>` carries its intrinsic pixel `width`/`height`, read from the
  image header: a percentage width alone resolves to zero in a viewer whose
  containing block is not laid out yet, and the figure then renders at 0x0 while
  the file is fine.

  Three rendering rules, each pinned by a test: the verdict badge shows only what
  the threshold layer produced (no re-judging), an unmeasured metric renders as a
  dash and never as `0`, and `tss_enrichment = 0.0` is shown as absent when its
  status says `background_zero` -- that value is a sentinel so a floor threshold
  can fire, not a score.

- **Nucleosome repeat length** (`assay_metrics.nucleosome_repeat_length`, reported
  by `qc --assay mnase`, plus fragment-length and autocorrelation figures). The
  NRL is the *spacing* of the mono/di/tri ladder, not the mononucleosome fragment
  length, so it is read from the lag that maximises the autocorrelation of the
  length histogram rather than from the position of its first peak; confusing the
  two is the usual way this number is reported 30-50 bp too small. An estimate is
  withheld (`error`, no `nrl_bp`) when the peak autocorrelation is below 0.05 or
  the uncertainty spans the whole search window: fragments all at one length
  previously produced "138.1 bp" from pure noise, autocorrelation -0.009.

- **`qc` rendered its figures before writing the JSON they read.** Every plot
  helper in `easyepi._plots` loads the written payload, so the first pass found no
  data and produced nothing; a second pass from `_cmd_qc` then wrote PDF-only
  files with the bare default. A QC run therefore carried no images at all. The
  figure block now runs after the JSON is on disk, and the duplicate render in
  `_cmd_qc` is gone, so one pass writes the figures the report embeds.

- **`_assay_metric_frame` matched only files named `*_qc.json`.** `qc --out`
  writes whatever filename the user chose (`atac.json`), so every assay figure
  built from that frame was silently skipped.

- **MNase figures were gated on the module name.** `qc --assay mnase` renders
  through the same entry point with module `"qc"`, so a `"mnase" in module_key`
  gate never opened and the command that produces MNase metrics never drew them.
  The gate now asks the payload whether MNase periodicity output is present.

- **A shared `figures/` directory leaked across runs.** The report scanned the
  whole directory, so a MNase report embedded an ATAC run's charts when both wrote
  to the same output folder. The report is now given the exact list this run
  produced.

- **PDF was the only format written, so nothing could be embedded.** A PDF cannot
  be shown in `<img>`, and the guard meant to add a raster format tested
  `figure_formats` for `None` instead of testing the resolved list, so it did
  nothing whenever that list was non-empty.

### Fixed (signal tracks, DMR model statements)

- **`bam2bigwig` normalisation and duplicate handling were implicit.** Three
  choices changed the track scale while the argument for them was unstated: the
  denominator a track is normalised by, whether duplicate reads are kept, and how
  a proper pair contributes. `tests/test_bam2bigwig_denominator.py` now pins the
  first -- a track built from one chromosome must not silently use that
  chromosome's mapped count as the whole-genome depth -- and
  `test_paired_end_parity.py` pins the third. `--keep-dup-all` is documented as
  what it is: MACS3 `--keep-dup all`, which thresholds on the number of reads
  sharing a coordinate, not on the BAM duplicate flag.

- **The DMR run summary described the wrong model.** Its `method` and
  `limitation` fields were hard-coded to the pooled binomial text, so a
  replicate-aware run wrote "binomial-based, not beta-binomial" into
  `<out>.json` about its own beta-binomial output. Both fields are now
  conditional on the model that actually ran, and
  `tests/test_dmr_summary_is_honest.py` pins the branches apart.
- **The DMR page stated `--pooled` as if it modelled replicates.** It does not:
  the pooled-binomial site test sums each group's per-site counts before testing,
  so it sees no between-replicate variance and is not a beta-binomial model.
  Measured p-values are therefore anti-conservative on heterogeneous replicates.
  The flag stays, because reproducing an older run is a legitimate use, and the
  command page now says which one it is.

- **`diffpeaks` described its test only by family.** It is a negative-binomial GLM
  Wald test with a single common dispersion, estimated by moments from
  within-group residuals and shrunk toward zero. That is the same family as
  DESeq2's default and not the same procedure; the module now states the
  difference where the statistic is produced, so a reader does not attribute one
  result to the other.

- **`easyepi diff --method deseq2` needed a hand-assembled environment.**
  `pyproject.toml` gains a `deseq2` extra (`pydeseq2>=0.5`) and lists it in
  `all`, so a clean install can run the documented model instead of being told at
  runtime to pip-install it. The sign convention is part of the contract and is
  stated where it is produced: log2(A / B), so a feature that doubles in group B
  reports about -1.

### Fixed (TSS enrichment and duplicate handling)

- **TSS enrichment counted both mates of a pair.** A proper pair has two mates
  that are the two ends of the *same* fragment, so piling both 5' ends up doubles
  the score relative to single-end data and makes the reported number depend on
  library layout rather than on the biology. Read1 only is now kept, which is the
  convention `coverage_parallel` already used, so the TSS track and the coverage
  track agree on what one fragment contributes.

- **Two TSS enrichment definitions shared one threshold.** `tss_enrichment` is
  measured from BAM read 5' ends in +-25 bp over the outer 100 bp;
  `tss_enrichment_bigwig` from +-100 bp over the 1900-2000 bp ring. They are the
  same quantity measured two ways, not two pieces of evidence, and sharing the
  shipped ENCODE floor of 2.0 was wrong in both directions: applied to the BAM
  definition the floor is no floor at all (one library scored 1220.08, so nothing
  would ever breach it), and letting a sample report both lets the better of two
  scales decide the verdict. The BigWig rule is now skipped whenever the BAM rule
  was evaluated, and `qc_payload` records a `tss_enrichment_disagreement` block
  naming both window definitions when the two differ by more than 3x.

- **`qc --annotation` was documented and ignored.** TSS enrichment accepted only
  `--tss`, so the obvious invocation -- pointing at the annotation the QC page
  itself describes -- silently produced no score. `--annotation` is now read as a
  TSS source too; `--tss` still wins when both are given.

### Fixed (found while building it)

- **Two TSS enrichment readers were writing the same key, 28.8x apart.** The package
  had a BigWig reader (`qc_extended.tss_enrichment`, +-100 bp over the 1900-2000 bp
  ring) and a BAM reader (`tss_enrichment.tss_enrichment`, read 5' ends in +-25 bp
  over the outer 100 bp). Both published `tss_enrichment`, and
  `build_sample_payload` updates the payload from each in turn, so whichever ran
  last decided the number silently. Measured on one library with both inputs
  supplied: 1220.08 versus 42.38, ratio 28.8x. Both were then judged against a
  single floor of 2.0, which only means anything on the BigWig scale -- applied to
  the BAM definition nothing could ever breach it.

  The readers now publish distinct keys (`tss_enrichment` for the BAM definition,
  `tss_enrichment_bigwig` for the other), `tss_enrichment_source` records which
  produced the value, and a sample that supplies both inputs reports
  `tss_enrichment_disagreement` with the ratio and the reason rather than leaving
  the conflict to be discovered. The BigWig rule is skipped as `not_applicable`
  when the BAM rule was evaluated, so one finding is never scored twice.

  Found by building a library, computing both definitions on it, and comparing;
  an earlier check used a transcript-only annotation, which made the BigWig reader
  find no TSS at all and hid the collision instead of ruling it out.

- **`qc --annotation` produced no TSS enrichment.** `--annotation` is the
  documented flag and the one in the generated command page, but TSS enrichment
  read only `--tss`, so the obvious invocation silently produced no score. It now
  accepts either, with `--tss` winning when both are given because a hand-curated
  table is a deliberate statement and a GTF is a convenience. A bad annotation is
  reported and skipped rather than aborting a QC run whose other metrics are fine.

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
- **`easyepi footprint`'s score is TOBIAS's footprint score divided by the
  flank level.** TOBIAS 0.17.5's `tobias_footprint_array` reports
  `flank_mean - fp_mean` in raw-signal units while `footprint_profile` reports
  `(shoulder - centre) / shoulder`, so the two are related by
  `TOBIAS = flank_level * easyepi`. Measured ratio on seven synthetic sites:
  **1.000000** (min = max = mean), asserted at `rel=1e-6` in
  `tests/test_tobias_footprint_score.py`. The TOBIAS side runs the real binary.
- Note the definitional constraint this uncovered: TOBIAS accumulates only
  **negative** values into its footprint sum, so a footprint that is merely a
  dip in an all-positive track leaves its score at the flank level exactly. A
  comparison has to use an ATACorrect-style (negative) footprint;
  `test_purely_positive_dip_leaves_tobias_at_the_flank_level` pins that so it is
  documentation rather than a rediscovered surprise.

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
