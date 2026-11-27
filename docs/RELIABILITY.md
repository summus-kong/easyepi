# Reliability audit

Audit date: 2026-11-27, against commit `05ace7f`. Scope: the active package
tree, public module paths, input contracts and local installation. This is a
software release-candidate review, not biological validation or a claim of
publication readiness.

This page is a **risk register**: what is known to be unfinished, wrong or
unverified. Priorities live in [ROADMAP.md](ROADMAP.md), datasets and thresholds
in [VALIDATION_PLAN.md](VALIDATION_PLAN.md), and the claim matrix in
[INNOVATION.md](INNOVATION.md). An earlier revision of this file carried its own
prioritised to-do list; that list is superseded by the roadmap.

**Corrections to the 2026-10-16 revision** (kept as a dated note rather than
silently rewritten): the earlier text said native motif helpers remain available
for compatibility and that `motifs --backend auto` may fall back to a built-in
engine. Both were removed on 2026-11-27 — `easyepi/motifs/homer_like.py`,
`de_novo_homer` and `discover_motifs` are gone, and no fallback exists, because
the substitute shared a mean 3.90 bp of consensus with HOMER's own top motifs.

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

## Open gaps, by consequence

The prioritised work items moved to [ROADMAP.md](ROADMAP.md). What stays here is
the register of what is currently *not* established, ordered by what it would
invalidate:

| Gap | Consequence if a result is reported without it |
| --- | --- |
| No cohort-level comparison against DSS for methylation | A DMR call is unvalidated; the five oracle tests skip because R is absent on the analysis host |
| No accessibility-matched motif background | An ATAC motif enrichment p-value is computed against a null that does not match the data; the ranking may reflect composition rather than binding |
| No background-sensitivity analysis | A motif ranking cannot be attributed to biology rather than to the choice of null |
| No cohort-level DESeq2 comparison | Differential results are directionally checked on synthetic data only |
| Linux CI has never executed | "CI passes" cannot be asserted |
| The TOBIAS reverse-strand residual is unexplained | Whole-track ATACorrect agreement is 0.989637, not 1.0; equivalence is not claimed |
| 100k+ peak and many-sample resource numbers are missing | Scale claims cannot be made |
| Motif enrichment is not TF activity, and ULM is not chromVAR | Reporting one as the other is a statistical error, not a documentation slip |

## Legacy compatibility

Low-level native **peak** utilities remain importable because the user requested
preservation of underlying functions. They are legacy exploratory helpers and are
not recommended as confirmatory backends.

Native **motif** discovery utilities do not remain: `homer_like`, `de_novo_homer`
and `discover_motifs` were removed on 2026-11-27. They were measured against HOMER
on real peak sequences and did not reproduce it — the k-mer counter shared 3.90 bp
of mean consensus with HOMER's top ten (the overlap of unrelated short motifs) and
the EM variant returned `AAAAAAAA` and `CCAGCCTG` as its best two motifs while
claiming all 600 sequences as sites. De novo discovery now delegates to HOMER and
raises when HOMER is absent; PWM scanning and known-motif enrichment remain and
need no external program. Disabled metadata-only peak/nucleosome CLI commands raise errors instead of producing apparent successful analyses. `nucleosome-positioning` is the only positioning entry point.

## Filesystem cleanup

Historical build directories, caches, output folders and the old helper-only `easyepi/python`, `easyepi/backends` and `easyepi/scripts` directories have been removed from the working tree with `maintenance/clean_legacy.ps1 -Apply`. Those paths were never part of the maintained source layout and remain explicitly excluded from installable artifacts by `pyproject.toml` and `MANIFEST.in`.

The only retained local-only tree is `.verification/`, which holds the project-local dependency runtimes, wheel/sdist builds and installed-package checks recorded above. It is git-ignored, pruned from artifacts and is not part of the distributable release.
