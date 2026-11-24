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
* Commits and `git push` happen on the local workstation only; the analysis host is
  not a push portal. Remote work is copied back and checksum-verified before it is
  committed. The local workstation is **only a relay**: tests and measurements run
  on the analysis host. See `DEVELOPMENT.md`, "Where commits are made and pushed".
* No item may be pushed before `tests/` passes on the analysis host and
  `scripts/gen_docs.py` exits 0.

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

**Status:** DONE (v0.6.0)
**Module:** `easyepi/assay_metrics.py`, `easyepi/tss_enrichment.py` (new)

| Sub-task | Status | Notes |
| --- | --- | --- |
| TSS enrichment score | DONE (v0.6.0) | `easyepi/tss_enrichment.py`; `qc --tss` accepts a GTF directly. Four states, so a zero-background library reads as measured-and-failed rather than never measured |
| Fragment-length distribution | DONE (v0.6.0) | `assay_metrics.fragment_length_summary`; proper pairs only, each template counted once, nucleosome-free / mono / di nucleosome fractions reported. Single-end input returns `not_applicable` rather than a fabricated length |
| Duplicate fraction | DONE (v0.6.0) | `assay_metrics.library_summary`; `get_duplicate_flag` when present, and `duplicate_status` says `not_marked` rather than printing a 0.0 that would claim a clean unmarked library |
| Mitochondrial fraction | DONE (v0.6.0) | `--mito-pattern`, never guessed; counted through the index over the whole file because chrM sorts last and a bounded scan missed it. Reports reads, mapped records and the denominator population; a non-matching pattern yields `null`, not 0.0 |
| Library complexity (preseq-style estimate) | DONE (v0.6.0) | `library_complexity` = distinct templates / reads inspected, matched pairs counted once. Not a preseq extrapolation: it reports observed complexity and makes no claim about the unseen portion of the library |
| Replicate concordance | DONE (v0.6.0) | `easyepi concordance`; `--peaks` gives correlation + Jaccard, `--matrix` correlation only. No threshold by design |
| Documented filtering rules | DONE (v0.6.0) | `QC.md` has "Which reads are discarded": the baseline filter, the `--min-mapq` defaults including why `dmr` uses 10, the fragment-length window, one-record-per-fragment, blacklist opt-in, and the fact that `is_qcfail` is deliberately not tested. Derived from the six reader loops, not from intent |

## 3. Fragment counting backend

**Status:** DONE (v0.6.0)
**Module:** `easyepi/coverage_parallel.py`, `easyepi/bam2bigwig.py`

The body path used to add one span per alignment, so a paired-end ATAC fragment
contributed twice. It now counts one span per template.

| Sub-task | Status | Notes |
| --- | --- | --- |
| Paired-fragment counting | DONE (v0.6.0) | mate two skipped via `is_read1`, span taken from TLEN to the mate's outer edge. `test_paired_end_parity` pins that PE coverage equals the SE body |
| Single-end fallback with explicit warning | DONE (v0.6.0) | unpaired records are counted at read extent and reported as `read`, so SE is never silently called a fragment |
| Report which mode was used | DONE (v0.6.0) | `counting` block in the BigWig JSON: `unit`, `fragments`, `unpaired_reads` and per-chromosome counts. Reported per chromosome, never summed across chromosomes, because chunks are padded and clipped |

Not covered: the cut-site path (`--assay atac`, `shift=tn5`) counts read ends, so
`counting` is absent there by design rather than reporting a wrong unit.

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

**Status:** PARTIAL (v0.6.0) — every field is recorded; the lock export is
environment-specific and stays out of the package
**Module:** `easyepi/provenance.py`, written into `status.json` under
`provenance` by `RunRecorder.finish()`

| Sub-task | Status | Notes |
| --- | --- | --- |
| SHA-256 of every input file | DONE | `provenance.digest_files()`; streams in 1 MiB
blocks so a 2.6 GB BAM does not have to fit in memory |
| Output file hashes | DONE | opt-in via `RunRecorder.hash_outputs()`; hashing every
output by default would double the wall time of a run that nobody asked to verify |
| Direct dependency versions | DONE | `provenance.package_versions()`; the list lives in
one place instead of being spelled out at the call site |
| External tool versions | DONE | `provenance.external_tools()`; resolved through `PATH`
then through explicit `search_paths`, because a conda env bin directory is not on `PATH` |
| Command-line provenance | DONE | argv already recorded |
| Environment lock file | TODO | out of scope: `micromamba env export` is
environment-specific and cannot be produced by this package |

A missing file, an absent tool and an unreadable path each produce a recorded
reason (`not_found`, `skipped_too_large`, `not_a_file`, `unreadable`) rather than
an exception, because a manifest written after the fact is more useful than a
run that fails while writing it. `describe_gaps()` turns those into sentences the
recorder logs, so a manifest that established nothing does not look like one that
had nothing to establish.

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
