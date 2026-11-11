# Changelog

## 0.2.0 - unreleased

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
