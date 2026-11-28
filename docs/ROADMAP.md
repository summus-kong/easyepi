# Product roadmap

Last reviewed: 2026-11-27, against `05ace7f`.

This page answers one question: what is being built next, and how will we know it
is done. It deliberately contains no measured results — those belong in
[VERIFICATION.md](VERIFICATION.md) and [dev/BENCHMARKS.md](dev/BENCHMARKS.md) — and
no risk register, which stays in [RELIABILITY.md](RELIABILITY.md).

Three documents, three jobs:

| Document | Answers | Changes when |
| --- | --- | --- |
| **this page** | What is next, for whom, and what "done" means | a task is added, reprioritised or finished |
| [VALIDATION_PLAN.md](VALIDATION_PLAN.md) | Which dataset and which reference tool prove it | a comparison is run or its threshold moves |
| [INNOVATION.md](INNOVATION.md) | Which manuscript claim each result supports | a claim gains or loses evidence |

## How a task is defined here

Every task carries the same eight fields, because a task missing any of them is
not actionable:

`Priority` · `Problem` · `Approach` · `Depends on` · `Done when` · `Test data` ·
`Risk` · `Status`

"Improve motif enrichment" is not a task. "Report the GC and length distributions
of target and background, and fail when they differ" is.

Priority is about *when*, not *how much it would be nice to have*:

- **P0** — a published result would be wrong without it, or the package cannot be
  installed and run as documented.
- **P1** — results are usable but a reviewer would ask for it before accepting a
  claim.
- **P2** — genuinely useful, no current claim depends on it.

Status is one of `TODO`, `IN PROGRESS`, `DONE`, `BLOCKED`, `DEFERRED`. `BLOCKED`
names what blocks it. `DEFERRED` names why it is not worth doing now.

---

## Phase 1 — Release correctness (P0, blocking publication)

No analysis claim matters if the documented install does not work. This phase is
short and is a prerequisite for everything else.

### 1.1 Dependencies, extras and external tools agree with the docs

- **Priority** P0 · **Status** DONE (`05ace7f`)
- **Problem** `all` did not include `hmmlearn` or `matplotlib` while
  README and PUBLISHING described it as the union of every extra; README annotated
  `hmmlearn` on `[bam,viz]`; `install.md` claimed Python 3.9 against `>=3.11`;
  MACS3 and HOMER were not declared anywhere machine-readable.
- **Approach** Union the extras; state the external programs in
  `[tool.easyepi.external]` rather than faking them as pip extras, which cannot
  install a conda package.
- **Done when** A test compares `all` against every extra by requirement name, and
  a test compares CI's install list against `pyproject.toml`.
- **Test data** `pyproject.toml`, `.github/workflows/linux.yml`.
- **Risk** An extra that only aliases another (as `bigwig` does) makes a naive
  label comparison fail; the test compares resolved package names instead.

### 1.2 Linux CI has actually run

- **Priority** P0 · **Status** BLOCKED (workflow must execute on GitHub Actions)
- **Problem** The workflow is supplied as configuration and has never executed in
  this project's history, so "CI passes" is not a claim anyone can make yet. Its
  `EXTRAS` also missed `tobias`, which made TOBIAS port tests skip.
- **Approach** Push and read the run; treat a green local suite as no evidence
  about the workflow.
- **Done when** A run is linked from [VERIFICATION.md](VERIFICATION.md) with its
  commit, Python versions and skip count.
- **Risk** The TOBIAS and DSS-related tests may legitimately skip on a runner
  without those tools; the workflow's own skip check must then be relaxed *with a
  stated reason*, not by deleting the check.

### 1.3 Reliability audit is stale and mis-scoped

- **Priority** P0 · **Status** TODO
- **Problem** `RELIABILITY.md` is dated 2026-10-16 and scopes itself to the
  `easyepiv1` tree, while the active tree is `easyepiv1_relay` and the latest
  commit is 2026-11-27. Two paragraphs are now factually wrong: they say native
  motif helpers remain available for compatibility, and that `motifs --backend
  auto` may fall back to a built-in engine. Both were removed on 2026-11-27.
- **Approach** Re-date and re-scope the audit; correct the two paragraphs; keep the
  document as a *risk register* only, moving priorities to this page.
- **Done when** No statement in `RELIABILITY.md` contradicts the code, and a test
  asserts that the removed modules are absent.
- **Test data** `tests/test_release_metadata.py` (already asserts `homer_like` is
  gone through `easyepi.__init__`).
- **Risk** Re-writing the audit invites rewriting history. Existing entries that
  were wrong are corrected with a dated note, not deleted.

---

## Phase 2 — Core reliability of the numbers (P0/P1)

These decide whether a reported number is defensible. Ordered by how much a
reviewer would object if they were missing.

### 2.1 Differential peaks and design-aware models

- **Priority** P0 · **Status** PARTIAL
- **Problem** The self-contained path compares normalised signal with rank/t
  tests. That cannot model counts, batch, repeated subjects or a time trend, and
  two-group tests are not a design. A PyDESeq2 adapter exists but the two paths
  are not described as answering different questions.
- **Approach** Keep the count model as the inferential path; label the signal path
  exploratory in its own output. Require an explicit contrast and report the
  design matrix, the number of features filtered, and outliers.
- **Depends on** nothing.
- **Done when** Every differential output states which model produced it, its
  design, its contrast and the filtering; a cohort-level comparison against
  DESeq2 exists per [VALIDATION_PLAN.md](VALIDATION_PLAN.md) §2.1.
- **Test data** Synthetic counts with a known direction (already used to pin the
  log2FC sign), then a real multi-donor ENCODE cohort.
- **Risk** A cohort with confounded donors and groups cannot separate batch from
  biology; the comparison must state what it cannot distinguish.

### 2.2 Methylation DMR

- **Priority** P0 · **Status** PARTIAL
- **Problem** The pooled binomial model ignores between-replicate variance. A
  beta-binomial path exists but has no multi-cohort comparison against DSS.
- **Approach** Treat beta-binomial as the default once it is validated; keep the
  binomial model and say in its output that it is not replicate-aware.
- **Depends on** R + DSS on the analysis host (currently absent, so five oracle
  tests skip).
- **Done when** Effect size, significance and region overlap are compared against
  DSS on at least two independent cohorts; the skips are either run or recorded as
  an explicit gap.
- **Test data** `tests/test_beta_binomial_matches_dss.py` (skips without R), plus a
  real WGBS cohort.
- **Risk** Installing DSS on the analysis host may not be possible; if so the gap
  is stated, never silently dropped.

### 2.3 ATAC/ChIP QC completeness

- **Priority** P1 · **Status** PARTIAL
- **Problem** FRiP, TSS enrichment, duplicate fraction, library complexity and
  mitochondrial fraction exist; there is no single documented filtering rule per
  assay, and no replicate-concordance metric.
- **Approach** Publish thresholds in `QC.md` with the source for each, and report
  the threshold used in every verdict so a non-default run is visible.
- **Done when** Each QC verdict names its threshold and the rule it implements;
  a sample that fails a rule cannot report "pass".
- **Test data** ENCODE ATAC and ChIP libraries with published QC values.
- **Risk** Thresholds are assay- and depth-dependent; a single default presented
  as universal would be worse than no default.

---

## Phase 3 — Motif chain (P0/P1/P2)

Motif analysis is a chain, not a feature: input and background, then site-level
scanning, then enrichment, then discovery, then TF activity, then figures. Each
step has different evidence needs and none substitutes for another.

### 3.1 Input and background preparation

- **Priority** P0 · **Status** PARTIAL (accessibility-matched mode DONE; the
  sensitivity comparison is not)
- **Problem** `easyepi background` produced GC- and length-matched controls only.
  For ATAC data a composition-matched genomic background is not the same null as
  an accessibility-matched one: it compares peaks against mostly-inaccessible
  sequence, which no plausible motif enrichment can be read against.
- **Approach** Accessible-universe sampling was chosen over shift-and-match.
  Shift-and-match needs an accessibility definition independent of the targets —
  the hard part — and a shifted control shares the target's local chromatin
  neighbourhood, so it is not independent. `--accessible BED` now restricts
  controls to a supplied universe; a control must lie **wholly** inside one merged
  accessible interval, and GC and length matching still apply.
- **Done when** (a) the report names the construction — DONE; (b) a target set
  whose universe already covers it is refused as circular — DONE; (c) the
  sensitivity comparison between unmatched, matched and accessibility-matched
  backgrounds is reported — **not done**.
- **Test data** 300 real H3K27ac peaks on chr1 against an independent 2000-interval
  accessible set: 297/300 matched, GC 0.520 vs 0.508 (KS p = 0.248), and the
  circular configuration (`--accessible` = the peaks) is refused with exit 1.
- **Risk** An accessibility background built from the same peaks it controls for is
  circular. The implementation refuses a universe that wholly contains every
  target, but a universe that merely *overlaps* the targets heavily is still
  partly circular and is not detected. The report prints the target fraction
  inside the universe so a reader can judge.
- **Remaining note** This was the highest-value item in the motif chain. What is
  left is not the construction but the evidence that the construction changes a
  ranking — see [VALIDATION_PLAN.md](VALIDATION_PLAN.md) §3.3.

### 3.2 Site-level known-motif scanning

- **Priority** P0 · **Status** DONE (`05ace7f`)
- **Problem** `scan_sequences` returned one best score per sequence, so no site
  table with coordinates and strand could be produced.
- **Approach** `scan_sites` returns every above-threshold site with sequence index,
  start, end, strand, score and optional p-value. Thresholds come from the exact
  convolution of per-column scores up to 10 bp.
- **Done when** Coordinates, both strands, ambiguous bases, short and empty inputs
  are covered by tests; a palindromic motif is not counted twice.
- **Test data** `tests/test_motif_parsing.py`; a real 10 bp HOMER motif on 300 real
  200 bp peak sequences.
- **Risk** Reported coordinates are offsets into the supplied sequence. Calling
  them genomic coordinates without the BED offset would be a silent error; the
  docstring and payload say offsets.

### 3.3 Known-motif enrichment

- **Priority** P1 · **Status** PARTIAL
- **Problem** The statistical unit is a sequence, the threshold comes from a
  uniform background, and the multiple-testing scope is across motifs only. The
  per-motif site threshold is not a cross-motif FDR, and the output does not
  currently separate the two.
- **Approach** State the unit, the threshold provenance and the correction scope in
  the payload; report rankings under every available background and the change
  between them.
- **Depends on** 3.1 for the accessibility arm — satisfied.
- **Done when** Two backgrounds that differ only in matching produce a reported
  change in ranking, and the output names which background was used — **DONE**, and
  the change is large (rho 0.343 between unmatched and GC-matched; the unmatched
  null's top hit falls to rank 113).
- **Risk** A ranking that is stable across backgrounds may reflect strong biology
  or a weak test; stability alone is not correctness. And the measurement shows the
  ranking is *not* stable, so a published ranking without its background is
  uninterpretable — that is now a documented requirement, not a caveat.

### 3.4 De novo discovery via external tools

- **Priority** P1 · **Status** PARTIAL
- **Problem** HOMER is wired and reports the tool that produced a result. MEME and
  FIMO wrappers exist under `universal.motif_analysis` but not behind the same
  reported interface, so the choice of engine is not visible in one place.
- **Approach** One backend selector that records tool name, version, arguments and
  raw output location. Do not reimplement HOMER, and do not claim equivalence.
- **Done when** Running each backend writes a provenance block naming the executable
  and version; a missing tool raises instead of substituting.
- **Test data** HOMER 4.11 and MEME 5.5.9 on the analysis host (both installed).
- **Risk** HOMER exits 0 when its helper programs are missing; that specific
  failure mode is already handled and must stay handled for any new backend.

### 3.5 TF annotation and activity

- **Priority** P2 · **Status** PARTIAL
- **Problem** Motif similarity, motif enrichment, TF activity and footprinting are
  four different pieces of evidence. `tf-annotate` reports motif presence;
  `activity` runs decoupler ULM on a user-supplied network. Nothing currently
  prevents a reader from treating the first as the third.
- **Approach** Keep them separate in output and in documentation; record the motif
  database version; state that ULM is not chromVAR and does not supply
  background-deviation correction.
- **Done when** No output file mixes a motif score with an activity score, and the
  ambiguity between homologous motifs is reported where a TF name is produced.
- **Risk** A TF name derived from a similar motif is a guess; describing it as an
  identification would be a publishable error.

### 3.6 Figures and reports

- **Priority** P2 · **Status** TODO
- **Problem** No logo or enrichment figure exists for the motif chain, and reports
  do not separate scanning, enrichment, activity and footprint results.
- **Approach** Plot from the tables the chain already writes, so a figure cannot
  show a number the table lacks.
- **Done when** Every panel traces to a table cell, and the report labels which
  step produced each section.
- **Risk** A logo drawn from a PWM that failed its threshold would look like a
  positive result.

---

## Phase 3.5 — Architecture (P1, done as part of the motif work)

### 3.5.1 Group and signal class are separate axes

- **Priority** P1 · **Status** DONE
- **Problem** `canonical_assay` answered "which group", and the Tn5 shift was then
  decided by re-matching the raw spelling. CUT&Tag and CUT&RUN both map to `chip`,
  so `resolve_defaults(canonical_assay("cuttag"))` returned `shift: none` while
  `resolve_defaults("cuttag")` returned `tn5` — a 4-5 bp displacement of every cut
  site, silently, for any caller that normalised first.
- **Approach** Two functions with two contracts. `canonical_assay(x)` keeps
  returning the group, so MACS3 arguments are unchanged. `describe_assay(x)`
  returns spelling, group and signal class, and anything deciding *how to process*
  reads both fields from that one object.
- **Done when** The answer does not depend on call order, and CUT&RUN never
  receives a Tn5 offset in any spelling.
- **Test data** `tests/test_assay_vocabulary.py` — 8 new tests covering both
  assays, every accepted spelling, and the ATAC/WGBS classes.
- **Risk** A caller that passes only a group label (`"chip"`) still cannot recover
  the signal class, and that is correct: the group genuinely does not say. The
  docs state it rather than adding a fallback that would guess.

## Phase 4 — Explanation and scale (P2)

### 4.1 Cross-omics interpretation

- **Priority** P2 · **Status** PARTIAL
- **Problem** `integrate` links peaks to genes by nearest TSS and clusters tracks.
  That is descriptive; it is not regulatory inference.
- **Done when** Outputs say "nearest TSS", not "regulatory link", everywhere.

### 4.2 Scale and resource measurement

- **Priority** P2 · **Status** PARTIAL
- **Problem** Binned depth over 20 Mb has measured scaling; 100k+ peaks and
  many-sample runs have no numbers.
- **Done when** Runtime and peak memory are reported for the sizes a real study
  reaches, with the machine described.
- **Risk** A benchmark on one host understates the variance a cluster user sees.

---

## Not planned

Stated so they stop being implied:

- Reimplementing HOMER, MEME, MACS3, DSS or chromVAR. The measured attempts at a
  HOMER substitute are why `homer_like` was removed rather than deprecated.
- A GUI. The interface is files in, files out.
- Native Windows support. The local workstation is a relay for committing, not an
  execution environment.
- Calling motif enrichment "TF activity", or ULM "chromVAR".
