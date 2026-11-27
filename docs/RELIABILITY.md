# Reliability audit

Audit date: 2026-10-16. Scope: the active `easyepiv1` Python package, public module paths, input contracts and local installation. This is a software release-candidate review, not biological validation or a claim of publication readiness.

## Changes made

| Area | Defect | Resolution |
| --- | --- | --- |
| Layout | Redundant `easyepi.python` nesting and stale imports | Independent modules now live in `easyepi.universal` and `easyepi.specific`; shared helpers moved to package root |
| Language | Chinese source/help/documentation mixed with English | Active source, CLI messages and maintained documentation use English |
| Runtime | Python 3.10 declared despite unconditional `tomllib` | Minimum Python 3.11 declared |
| Packaging | Stale builds, R artifacts and obsolete nested paths could enter distributions | Explicit package discovery, manifest exclusions and isolated clean-source build verification |
| Peaks | ChIP narrow/broad entries invoked ATAC native caller | Official MACS3 commands distinguish assay, broad/narrow peaks and BAM/BAMPE |
| MNase | Positioning invoked ATAC peaks | Misleading entry rejects execution pending an established positioning backend |
| Counts | BigWig integral called read count | Explicit BAM alignment-count vs BigWig integrated-signal modes and units |
| Consensus | Only exactly matching peak coordinates counted | Interval support across distinct samples; adjacent retained segments merged |
| Tn5 | Negative-strand shift used alignment start | Negative-strand cuts use reference end minus five; flags/MAPQ checked |
| Normalization | Default library sizes invented; invalid values replaced silently | Explicit measured factors, validation, corrected RPGC coverage denominator |
| Differential | Pair IDs validated without sample reordering | Reorder paired matrices by subject; reject duplicate/missing design IDs |
| Missingness | Paired values dropped independently; failed tests could appear significant | Drop complete pairs, preserve untested NaNs and exclude them from BH test count |
| Region coverage | Default random 20,000-peak cap | Retain all regions unless an explicit limit is requested |
| Methylation | Values clipped, fractional coordinates truncated, count coverage set to one | Reject invalid values; derive measured coverage; preserve unknown coverage |
| Motif enrichment | Every finite scan score considered a motif hit | Biopython PWM score thresholds, Fisher exact test and BH |
| Motif scan | Ignored forward-only option and reverse-hit coordinate convention | Honor strand option and return forward-sequence coordinates |
| External tools | Nonzero exit ignored; FIMO stdout discarded; MEME objective invalid | Raise on external failure, retain FIMO hits, use classic objective with ZOOPS model |
| Motif workflow | Silent HOMER-to-native fallback and false equivalence claims | The native de novo engine is gone (0.7.1): `homer_like` and `de_novo_homer` were removed after measuring 3.90 bp mean shared consensus with HOMER, so there is no fallback to be silent about. De novo requires HOMER and raises without it; a nonzero HOMER exit raises. Known-motif enrichment and PWM scanning stay and need no external program |
| TF activity | Sequence-only motif averages described as chromVAR-like activity | decoupler ULM on measured feature matrices and a declared network |
| ChromHMM | Headered CellMark, wrong BinarizeBed positional args, unused threads/seed | Official headerless table, peak mode, bin size, explicit paths, seed and processor options |
| ChromHMM stages | Imported stage main did not set its stage; standalone training could not read prior binary output | Callable stage wrappers and `--binary-dir` |
| State labels | Renaming anonymous columns called biological annotation | Coverage summaries plus optional manually reviewed labels |
| Integration | Hand-written k-means and misleading cluster count on identical profiles | scikit-learn KMeans; reject fewer distinct profiles than requested clusters |
| Plotting | Implicit PNGs, selected PDF also emitted PNG, NaN p-values filled with zero | No export by default, explicit formats only, missing p-values preserved |
| QC | Missing reader returned a successful-looking file-only summary | Missing BAM/BigWig dependencies fail explicitly when needed |

## Validation evidence and limits

Regression tests cover strand semantics, MACS3 commands, official ChromHMM arguments, paired sample order, missing observations, BH correction, methylation validation, JASPAR parsing, calibrated motif thresholds, FIMO output, external failure propagation, plotting defaults and public imports. Existing API tests are retained. See `VERIFICATION.md` for the final run results.

The decoupler adapter was checked against the installed 1.9.2 implementation using synthetic measured matrices on the available development host. That host is outside the supported Linux runtime and required a numba environment adjustment for a local smoke check. No Windows-specific workaround was added to the package. A Linux CI workflow and environment specification are supplied; results remain pending until executed on Linux.

ChromHMM command construction was checked against the [official manual](https://compbio.mit.edu/ChromHMM/ChromHMM_manual.pdf).

External-tool concordance is now measured for five of the six tools, with the numbers in `docs/dev/TOOL_COMPARISONS.md`. MACS3 output is byte-identical because `easyepi peaks` invokes MACS3 rather than reimplementing it. FIMO agrees on 11 of 12 best-motif sites (Spearman 1.000) and never on absolute score, because the two tools define the score differently. The TOBIAS ports reproduce the reference at machine precision at every layer that can be isolated -- scoring (1.8e-15), cut sites and k-mers (18/18), the sampled count matrices (0.000e+00) and `fast_rolling_math` (1e-12) -- but whole-track ATACorrect agreement is 0.999992 forward and 0.989637 reverse, and the reverse residual is unexplained; it is recorded as a known limitation rather than claimed as equivalent. ChromHMM was trained end to end on real ENCODE marks. HOMER remains PARTIAL: its installed genomes are absent, so only the FASTA path is exercised.

A real cohort benchmark now exists for annotation: twelve peak sets (three donors x four histone marks) compared against ChIPseeker on a shared gene model agree on every annotation string and every distance, 23,646 peaks in total. That establishes implementation equivalence *given the same gene model*; it does not establish that either tool picks the right TSS set, and it does not extend to any other module. The remaining comparisons still rest on one dataset, one chromosome or a synthetic fixture, so command tests and reader mocks likewise still do not demonstrate biological validity for them. See `docs/dev/BENCHMARKS.md` for the numbers, the tool versions and the commands.

## Remaining development priorities

1. **Count-based differential analysis:** integrate an explicitly versioned PyDESeq2 adapter with raw fragment-count input, covariates, contrasts, independent filtering, dispersion QC and real-data concordance. The present rank/t-tests are exploratory signal comparisons.
2. **ATAC/ChIP QC:** add FRiP, TSS enrichment, fragment-length distributions, duplicate/mitochondrial fractions, library complexity and replicate concordance with documented filtering and control rules.
3. **Fragment quantification:** provide an established fragment-count backend; current BAM mode counts alignments, not paired fragments.
4. **Motif background and networks:** implement validated GC/length/accessibility-matched background selection and reliable motif-hit-to-peak/TF network conversion. ULM itself does not supply chromVAR's bias correction.
5. **Methylation:** add a maintained Python coverage-aware site/region model or an explicitly supported external mature tool. Current site tests do not detect DMRs or model sequencing coverage.
6. **MNase:** `easyepi nucleosome` reports occupancy and a nucleosome-free fraction from a signal track, and `nucleosome-positioning` calls dyads from paired-end fragments. Neither reproduces NucleoATAC or MINCE, so positions are not interchangeable with those tools. There is no reserved positioning module; the stub that used to reject execution was deleted once `nucleosome-positioning` existed.
7. **Experimental designs:** support time trends, batches and repeated subjects through appropriate established models; do not treat a time course as arbitrary two-group tests.
8. **Annotation:** improve GFF3/isoform support, blacklist handling, coordinate/build compatibility checks and scalable interval indexing. Methylation annotation currently expects gene-level GTF.
9. **Scale:** replace dense per-base BigWig extraction/large all-pairs operations where justified. Binned depth over 20 Mb now has measured scaling (0.902 s at one worker, 4.24x at eight, output identical throughout, work memory about 38 MB independent of worker count). 100k+ peaks and many-sample numbers are still missing.
10. **Reproducibility:** add hashed inputs, complete tool/database versions, environment locks, Linux CI and programmatic argument provenance. Shared provenance currently records process argv and selected dependency versions.
11. **Visualization:** extend assay-specific diagnostic functions using established definitions. Generic table plots alone do not constitute complete module QC.
12. **Publication validation:** prespecify datasets, comparisons and metrics; reproduce original tool outputs, quantify runtime/memory and evaluate biological recovery on held-out data. The software's organizing interface alone does not establish methodological novelty.

## Legacy compatibility

Low-level native peak and motif utilities remain importable because the user requested preservation of underlying functions. They are legacy exploratory helpers and are not recommended as confirmatory backends. Disabled metadata-only peak/nucleosome CLI commands raise errors instead of producing apparent successful analyses. `nucleosome-positioning` is the only positioning entry point.

## Filesystem cleanup

Historical build directories, caches, output folders and the old helper-only `easyepi/python`, `easyepi/backends` and `easyepi/scripts` directories have been removed from the working tree with `maintenance/clean_legacy.ps1 -Apply`. Those paths were never part of the maintained source layout and remain explicitly excluded from installable artifacts by `pyproject.toml` and `MANIFEST.in`.

The only retained local-only tree is `.verification/`, which holds the project-local dependency runtimes, wheel/sdist builds and installed-package checks recorded above. It is git-ignored, pruned from artifacts and is not part of the distributable release.
