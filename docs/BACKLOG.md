# Development backlog

Status of each item is one of: `TODO` / `IN PROGRESS` / `DONE` / `BLOCKED`.
Every item names the module it will live in, so an entry is actionable rather than
aspirational. Measured claims go in `VERIFICATION.md`, not here.

Last reviewed: 2026-11-12 (v0.6.0).

## Cross-cutting constraints

* No item may break `tests/` — 542 tests pass at v0.5.3 and that number only grows.
* Any behavioural change gets a version bump and a `CHANGELOG.md` entry.
* Anything measurable gets a number in `VERIFICATION.md`, including failures.
* Optional dependencies are installed with micromamba from `conda-forge`/`bioconda`;
  `pip` cannot reach PyPI from the analysis host.

---

## 1. Count-based differential analysis (PyDESeq2)

**Status:** DONE (v0.6.0)
**Module:** `easyepi/differential_count.py` (new), wired into `easyepi/diff.py`
**Why it matters:** the current tests are rank-sum/t-tests on normalised signal. They
cannot model counts, covariates, batches or independent filtering, so a "significant"
call is exploratory rather than inferential.

| Sub-task | Status | Notes |
| --- | --- | --- |
| Install `pydeseq2` via micromamba | DONE | 0.5.4 from conda-forge |
| Adapter: raw count matrix to `DeseqDataSet` | DONE | `easyepi/differential/deseq2.py` |
| Covariates and design formula passthrough | DONE | `--design-factors`; single-level factors rejected |
| Contrast specification | DONE | explicit operands; sign measured, see below |
| Independent filtering | DONE | count of untested features reported |
| LFC shrinkage (`--lfc-shrink`) | DONE | recorded in the info block |
| Version pin recorded in provenance | DONE | `deseq2.tsv.deseq2.json` |

**Measured note:** PyDESeq2's `contrast=[factor, first, second]` reports
log2(second / first), the reciprocal of R DESeq2's convention. The first version of the
adapter silently inverted every log2FoldChange — p-values unchanged, so any
significance-only check would have passed it. Settled by building data with a known
direction (group B doubled, group A halved): a correct adapter reports ≈+1 and ≈−1 and
zeros for the null features. Pinned by `tests/test_deseq2_adapter.py`.
| Dispersion QC as a plot | TODO | the fit type is reported, not yet plotted |
| Multi-factor interactions | TODO | additive terms only; an interaction needs an explicit formula |

## 2. ATAC/ChIP QC completeness

**Status:** TODO
**Module:** `easyepi/assay_metrics.py`, `easyepi/tss_enrichment.py` (new)

| Sub-task | Status | Notes |
| --- | --- | --- |
| TSS enrichment score | TODO | needs `genome.fa` + TSS table; both exist on the host |
| Fragment-length distribution | TODO | derive from BAM CIGAR/TLEN; report nucleosome-free/mono/di peaks |
| Duplicate fraction | TODO | pysam `get_duplicate_flag`; distinct from library complexity |
| Mitochondrial fraction | TODO | needs a chrom-name pattern; make it a flag, not a guess |
| Library complexity (preseq-style estimate) | TODO | at minimum: distinct fragments / total fragments |
| Replicate concordance | TODO | pairwise correlation + Jaccard on peaks; needs a group column |
| Documented filtering rules | TODO | thresholds must be in `QC.md`, not only in code |

## 3. Fragment counting backend

**Status:** TODO
**Module:** `easyepi/fragment_counts.py` (new)

Current BAM mode counts alignments. For paired-end ATAC that over-counts: two mates of
one fragment are two alignments.

| Sub-task | Status | Notes |
| --- | --- | --- |
| Paired-fragment counting | TODO | proper pair + same chrom + TLEN-based dedup |
| Single-end fallback with explicit warning | TODO | never silently treat SE as fragments |
| Report which mode was used | TODO | in metrics JSON and the manifest |

## 4. TF activity bias correction

**Status:** TODO
**Module:** `easyepi/activity.py`

ULM gives activity relative to a supplied background, not to GC/accessibility bias.
chromVAR's approach (background peaks matched on GC + counts) is the established fix.

| Sub-task | Status | Notes |
| --- | --- | --- |
| GC + accessibility matched background peaks | TODO | reuse `easyepi/gc_background.py` |
| Bias-corrected deviation z-scores | TODO | document what is and is not chromVAR-equivalent |
| Replicate-aware variance | TODO | requires the group column from item 2 |

## 5. Methylation model

**Status:** TODO
**Module:** `easyepi/dmr.py`

Current DMR test is a Wald test per site with a minimum-coverage guard. It is
anti-conservative at high coverage and unstable at low coverage.

| Sub-task | Status | Notes |
| --- | --- | --- |
| Beta-binomial site model | TODO | over-dispersion estimated, not assumed |
| Smoothing across neighbouring CpGs | TODO | the main power gain in real DMR callers |
| Permutation-based region significance | TODO | avoids relying on the asymptotic approximation |
| Keep the current test as an explicit fast mode | TODO | do not delete; mark it exploratory |

## 6. MNase

**Status:** TODO
**Module:** `easyepi/nucleosome_positioning.py`

| Sub-task | Status | Notes |
| --- | --- | --- |
| Establish what `nucleosome-positioning` currently claims | TODO | it must not imply a validated caller |
| Report the heuristic explicitly in output | TODO | or connect a real caller if one is installable |

## 7. Experimental design

**Status:** TODO
**Module:** `easyepi/differential/`, `easyepi/project.py`

| Sub-task | Status | Notes |
| --- | --- | --- |
| Covariate columns in the project spec | TODO | validated, with clear errors |
| Batch as a design term | TODO | not a post-hoc correction |
| Time course handled explicitly | TODO | refuse to treat it as arbitrary two-group |
| Repeated subjects | TODO | needs a blocking factor |

## 8. Annotation robustness

**Status:** TODO
**Modules:** `easyepi/annotate/`, `easyepi/formats.py`

| Sub-task | Status | Notes |
| --- | --- | --- |
| GFF3 isoforms without transcript lines | TODO | currently degrades to one unit per gene |
| Genome build compatibility check | TODO | refuse hg19 GTF with hg38 peaks |
| Chromosome naming reconciliation | TODO | `chr1` vs `1` must not silently produce zero overlaps |
| Interval index for large peak sets | TODO | `pyranges` backend exists; benchmark at 100k+ |

## 9. Scale and performance

**Status:** TODO

| Sub-task | Status | Notes |
| --- | --- | --- |
| Benchmark 100k peaks x N samples | TODO | record memory and wall time in `VERIFICATION.md` |
| Avoid dense per-base extraction where avoidable | TODO | chunked accumulation already exists in `coverage_parallel.py` |
| Streaming for large tables | TODO | only if the benchmark shows it is needed |

## 10. Reproducibility

**Status:** TODO
**Module:** `easyepi/provenance.py` (new)

| Sub-task | Status | Notes |
| --- | --- | --- |
| SHA-256 of every input file | TODO | recorded in the run manifest |
| Output file hashes | TODO | lets a downstream step verify inputs |
| Direct dependency versions | TODO | package versions, not just Python |
| External tool versions | TODO | `macs3 --version`, `bedtools --version`, etc. |
| Command-line provenance | DONE | argv already recorded |
| Environment lock file | TODO | export from micromamba |

## 11. Assay-specific diagnostics

**Status:** IN PROGRESS
**Module:** `easyepi/_plots.py`, `easyepi/viz/`

Generic table plots alone are not module QC.

| Sub-task | Status | Notes |
| --- | --- | --- |
| ATAC: fragment-size and TSS-enrichment profiles | TODO | depends on items 2 and 3 |
| ChIP: profile and heatmap around TSS | DONE | `viz` metagene |
| MNase: nucleosome repeat length | TODO | |
| Methylation: beta distribution per sample | DONE | read from BAM or table |

## 12. External tool benchmarks on a real cohort

**Status:** TODO

Command construction is verified; output equivalence is not.

| Tool | Sub-task | Status |
| --- | --- | --- |
| MACS3 | compare peak sets across parameter settings on real data | TODO |
| ChromHMM | run the official binary end to end, compare states | TODO |
| MEME/FIMO | compare hits against the web service output | TODO |
| HOMER | compare motif enrichment against HOMER's own output | TODO |
| TOBIAS | compare footprint scores | TODO |

## 13. Release

**Status:** TODO
**Why it matters:** the only finished-but-unusable part of the project.

| Sub-task | Status | Notes |
| --- | --- | --- |
| Choose the first public version number | TODO | 0.5.3 has never been published |
| `twine check` on the built artifacts | TODO | wheel job validates metadata; twine is the publisher's check |
| Test the install from a clean env, no extras | TODO | verifies the minimal dependency set |
| Publish to TestPyPI first | TODO | |
| Then PyPI | TODO | requires the user's credentials |
