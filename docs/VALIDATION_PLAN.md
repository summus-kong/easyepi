# Validation plan

Last reviewed: 2026-11-27, against `05ace7f`.

This page says, for each analysis: which data proves it, which reference tool it
is compared against, which metric decides, and what threshold counts as agreement.
It contains no priorities — those are in [ROADMAP.md](ROADMAP.md) — and no results,
which go in [VERIFICATION.md](VERIFICATION.md).

Three levels of evidence, and they are not interchangeable:

| Level | Question it answers | What it cannot show |
| --- | --- | --- |
| **Format and contract** | Does the parser read the real file; are coordinates, strands and units what they claim? | Anything about biology |
| **Agreement with a reference** | Given the same input, does it produce the same answer as the tool it stands beside? | That either tool is right |
| **Biological recovery** | On held-out data with known truth, is the expected signal recovered? | That the method generalises beyond that dataset |

A test suite establishes the first level. This plan exists because the project
repeatedly mistook the first level for the second: the removed `homer_like` module
passed its own round-trip tests while sharing 3.90 bp of consensus with HOMER.

## Standing rules

1. **Compare against actual output of the reference, never against its
   documentation.** Every oracle in `tests/oracle/` asserts the reference agrees
   with itself before the comparison is made.
2. **Record the version of everything.** Tool name, version, arguments, database
   release, and the commit under test. A result without these is an anecdote.
3. **State the residual.** Where agreement is not exact, the number and its
   unexplained part are both written down. The TOBIAS reverse-strand residual is
   the model for this: 0.989637 is recorded as a limitation, not smoothed away.
4. **A skip is a result.** Five DSS tests skip on the analysis host because R is
   absent. That is a gap in coverage, recorded as a gap.
5. **One dataset is a fixture, not a cohort.** A single chromosome or a single
   donor establishes the code path, not the method.

---

## §1 Release correctness

### 1.1 Install surface

| Item | Value |
| --- | --- |
| Data | `pyproject.toml`, `.github/workflows/linux.yml`, a clean venv |
| Reference | the documented install line in README |
| Metric | every extra resolvable; every documented command importable |
| Threshold | zero missing requirements, zero import errors |
| Covered | `tests/test_release_metadata.py`, `scripts/verify_installed.py` |
| Not covered | the GitHub Actions run, which has not executed |

### 1.2 Package imports from an installed wheel

| Item | Value |
| --- | --- |
| Data | built wheel, installed into a clean venv, run from a directory that is not the checkout |
| Metric | import succeeds and resolves to the installed path, not the source tree |
| Threshold | exact |
| Covered | `scripts/verify_installed.py`, the `wheel` job in CI |
| Not covered | that the CI job runs |

An existing global editable install once masked missing submodules; the check
exists because of that, and `python -S` was needed to see it.

---

## §2 Core analyses

### 2.1 Differential peaks

| Item | Value |
| --- | --- |
| Data | (a) synthetic counts with a known direction per feature; (b) a real multi-donor ENCODE cohort |
| Reference | DESeq2 via PyDESeq2 for counts; the signal path has **no** reference and is labelled exploratory |
| Metric | sign and magnitude of log2FC; rank correlation of p-values; count of features called at FDR 0.05 with and without independent filtering |
| Threshold | sign exact on the synthetic set; Spearman ρ ≥ 0.90 and ≥ 80% overlap at FDR 0.05 on the cohort |
| Covered | `tests/test_deseq2_adapter.py` — direction, sign flip on group swap, batch term, provenance |
| Not covered | a cohort run; batch/repeated-measure models; anything about causal interpretation |
| Known trap | PyDESeq2's `contrast` reports log2(second/first), the reciprocal of R DESeq2. A significance-only check passes an inverted adapter, which is why the synthetic set carries a direction. |

### 2.2 Methylation DMR

| Item | Value |
| --- | --- |
| Data | two independent real WGBS cohorts with replicates |
| Reference | DSS `DMLtest` |
| Metric | per-site effect-size correlation; overlap of significant sites; overlap of merged regions at the same FDR |
| Threshold | ρ ≥ 0.95 on effect size; ≥ 80% site overlap; region overlap reported but not thresholded, because merging amplifies small differences |
| Covered | `tests/test_beta_binomial_matches_dss.py` — **skips**: R and DSS are absent on the analysis host |
| Not covered | everything. This is the largest gap in the project's evidence. |
| Unblocking | R + DSS on the analysis host, or a cohort run reported as an explicit gap |

### 2.3 ATAC/ChIP QC

| Item | Value |
| --- | --- |
| Data | ENCODE libraries with published QC values, plus one deliberately poor library |
| Reference | the published metric definition, computed independently |
| Metric | FRiP, TSS enrichment score, fragment-length modes, duplicate fraction, mitochondrial fraction |
| Threshold | each metric within 5% of the reference; a poor library must fail the documented rule |
| Covered | metric units and denominators; the "interpretable: false" path for regional peak sets |
| Not covered | cross-lab stability of thresholds; replicate concordance (no metric yet) |
| Known trap | a scan bounded by `--max-reads` stops before chrM, which made an early version report 0.0 mitochondrial fraction for a contaminated library. The count now goes through the BAM index. |

---

## §3 Motif chain

Ordered as the chain runs, because a later step inherits an earlier step's null.

### 3.1 Background construction

| Item | Value |
| --- | --- |
| Data | real ATAC and H3K27ac peak sets on a real genome; synthetic genomes with a known GC layout as a truth check |
| Reference | none — this is a construction, so the check is internal consistency plus sensitivity |
| Metric | (a) do target and background GC and length distributions agree? (b) is every matched control inside the tolerance the report claims? (c) how much does a downstream ranking move between background designs? |
| Threshold | (a) KS p > 0.01 on both GC and length; (b) zero controls outside tolerance, verified by independent recomputation; (c) reported, not thresholded |
| Covered | `tests/test_gc_background.py` — 23 tests including a control recomputed against its own interval GC; `tests/test_background_report.py` — 11 tests |
| Not covered | the sensitivity comparison **across** background designs — the accessibility arm exists but its effect on a downstream ranking has not been measured |
| Accessibility | `--accessible BED` restricts controls to a supplied universe; a control must lie wholly inside one merged interval. A universe that wholly contains the targets is refused as circular. Measured on 300 real peaks: 297/300 matched, and the circular case exits 1 |
| Measured | 300 real H3K27ac peaks, chr1, 1000 bp grid: target GC 0.520 vs background 0.513, KS p = 0.563; 100% matched after the interval-GC fix; 2.5 s for 1.245M windows |

Two failures this section already caught, recorded because they are the class of
error the plan exists to find:

- 144 of 197 controls sat outside their own stated 0.05 tolerance, because the
  target was matched on one window's GC while the control was measured over its
  whole interval.
- 73% match rate on a real peak set, with every 7.7 kb peak failing, because the
  eligibility test rejected any target wider than a grid window.

Neither was visible from easyepi's own tests; both appeared the first time a real
peak set was used together with a report of the distributions.

### 3.2 Site-level scanning

| Item | Value |
| --- | --- |
| Data | synthetic sequences with planted sites on both strands, ambiguous bases, short and empty inputs; real 200 bp peak sequences |
| Reference | FIMO 5.5.9 on the same PWM library |
| Metric | number and position of sites above a shared threshold; strand assignment; score correlation |
| Threshold | identical site sets where the threshold is comparable; ρ ≥ 0.98 on scores |
| Covered | `tests/test_motif_parsing.py` — coordinates, both strands, N-containing windows, empty input, palindromic dedup, agreement with `scan_sequences` |
| Not covered | a full-library FIMO comparison at cohort scale |
| Measured | FIMO 5.5.9 accepted easyepi-written MEME after the format fix: 15 sites for one motif, 11324 across a 469-motif library. easyepi scan: 300 sequences × 469 motifs in 24 s, FIMO in 12 s |
| Known trap | `pwm_score_distribution` bins the score axis and reports bin edges, so a perfect match (7.8242266) fell below its own threshold (7.8248567) and every site was dropped. The exact convolution is used for motifs up to 10 bp. |

### 3.3 Known-motif enrichment

| Item | Value |
| --- | --- |
| Implementation | **AME 5.5.9** (MEME Suite), called from `easyepi.motifs.ame_enrichment`. The test (`--ame-method`) and scoring (`--ame-scoring`) are recorded in every result |
| Data | the same real peak sets as 3.1, with three backgrounds: unmatched genomic, GC/length-matched, and accessibility-matched |
| Reference | HOMER `knownResults.txt`, plus FIMO, under the same sequence sets |
| Metric | rank correlation of enrichment across motifs; overlap of the top 20; stability of a motif's rank between backgrounds |
| Threshold | **not met, and withdrawn**: agreement with HOMER was the original gate (ρ ≥ 0.80) and four implementations failed it for reasons that turned out not to be implementation defects — see below |
| Covered | delegation to the published tool; argument validation; control-set requirement; the exact invocation is written into the output |
| Not covered | a second cohort; a calibrated p-value that does not depend on the partition choice |

**Why the HOMER-agreement gate was withdrawn.** Four implementations were compared on
identical input (300 peaks, 300 background sequences, 449 motifs):

| Implementation | Wall time | rho(hit rate) vs HOMER | rho(log P) vs HOMER |
| --- | --- | --- | --- |
| local scan + binomial | 77 s | 0.109 | 0.237 |
| port of HOMER's own hypergeometric + threshold search | 88 s | 0.120 | 0.206 |
| FIMO scan + the same hypergeometric | 89 s | 0.169 | 0.233 |
| **AME 5.5.9 (fisher / avg)** | **5 s** | — | 0.396 (n=144) |

FIMO against the local scan directly, with the same test applied to both, reaches
only rho 0.292 on hit rates and 0.435 on log P, with 6/20 top-20 overlap. Two
independent standard scanners therefore disagree with each other about as much as
either disagrees with HOMER. The ranking is set by scanner and partition choice, not
by which implementation is correct, so "match HOMER" was not a reachable target and
no amount of local reimplementation was going to reach it.

The concrete reason HOMER cannot be reproduced by a standard scan was also measured.
For Sp1 on this input HOMER reports 30/300 target hits and **0/300** background hits;
a direct PWM scan of the same background sequences cannot produce 0 hits at any
threshold, and the best pair it can reach (169, 30) scores *more* significant
(log P 34.89) than HOMER's own answer (17.21). HOMER scores background sequence
through its own k-mer subsystem, so the disagreement starts in the scanner, before
any statistic is applied. The source was read to confirm this
(`Motif2.cpp:1300` threshold search, `Motif2.cpp:4707` hypergeometric, `Motif2.cpp:900`
background rescaling), and a direct port of `loghypergeoD` was checked against scipy
on 3000 random tables (max |log difference| 1.7e-12) — the formula was never the
problem.
| **Sensitivity (measured)** | `scripts/bench_background_sensitivity.py` compares unmatched, GC/length-matched and accessibility-matched backgrounds on the same targets and motifs |
Measured on 300 real H3K27ac peaks (chr1), 120 HOMER motifs, three backgrounds:

| Background | GC target / bg | KS p (GC) | n matched |
| --- | --- | --- | --- |
| genomic (unmatched) | 0.520 / 0.396 | 7e-30 | 300/300 |
| GC + length matched | 0.520 / 0.507 | 0.342 | 300/300 |
| accessible + GC + length | 0.520 / 0.509 | 0.174 | 297/300 |

Spearman rho of the enrichment ranking between configurations, and the top-10
overlap:

| Comparison | rho | top-10 overlap |
| --- | --- | --- |
| genomic vs GC+length | **0.343** | 2/10 |
| genomic vs accessible | **0.422** | 3/10 |
| GC+length vs accessible | **0.196** | — |

Largest single-motif moves, competing-null to matched:

```
AP-2gamma(AP2)    rank   1 -> 113   (+112)
ETS:E-box         rank  12 -> 117   (+105)
DUX(Homeobox)     rank 113 ->  16   ( -97)
```

The unmatched null's top hit, AP-2gamma at log P 12.97, falls to rank 113 once the
background is GC-matched. That is what "enrichment" measured when the null was
mostly-AT-rich sequence: composition, not binding.
| Known limit | AME partitions sequences by affinity score; HOMER searches a threshold. They answer overlapping but different questions, so a p-value from one is not interchangeable with the other even at the same nominal threshold. The 6 tests x 4 scorings AME offers do not agree with each other either, which is why both are recorded rather than fixed. |
| Known limit | the per-motif site threshold is not a cross-motif FDR, and the two are corrected by different mechanisms. Reporting one as the other would be a statistical error. |

### 3.4 De novo discovery

| Item | Value |
| --- | --- |
| Data | real peak sets with a known composition (e.g. a H3K27ac set where AP-1 is expected) |
| Reference | HOMER 4.11 and MEME 5.5.9, run natively |
| Metric | motif recovery; whether expected factors appear; whether the tool name and version are recorded |
| Threshold | the expected factor appears in the top 5 of both; provenance block complete |
| Covered | HOMER end to end on real peaks (16 de novo motifs, populated `knownResults.txt`) |
| Not covered | MEME/STREME behind the same reported interface; a multi-dataset recovery rate |
| Known trap | HOMER exits 0 when its helper programs are missing, writing a 0-byte motif file. That was reported as success until fixed. Any new backend needs the same check. |

### 3.5 TF activity

| Item | Value |
| --- | --- |
| Data | a measured, normalised matrix and a user-supplied TF–target network |
| Reference | decoupler's published ULM implementation |
| Metric | agreement with the same method run outside easyepi |
| Threshold | exact within float tolerance |
| Covered | `tests/test_activity_backends.py`; rejection of sequence-only input and constant networks |
| Not covered | any comparison to chromVAR; the project does not implement background-deviation correction and must not imply it does |

---

## §4 What this plan cannot establish

- That any motif is a functional binding site. Motif identity is shared between
  paralogues, chromatin context is ignored, and no footprint is computed.
- That a peak-to-gene link is regulatory. Linkage is nearest TSS.
- That a p-value is calibrated for the assay. The nulls are sequence-composition
  nulls unless an accessibility background is supplied.
- That results generalise across labs, antibodies and sequencing depths. Every
  cohort used here is one cohort.

Agreement with HOMER, FIMO or DSS supports "produces a consistent answer under the
stated inputs", not "is biologically correct". The manuscript claims in
[INNOVATION.md](INNOVATION.md) are written to that limit.
