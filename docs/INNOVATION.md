# Publication evidence

Last reviewed: 2026-11-27, against `05ace7f`.

This page holds the claim–evidence–limit matrix a manuscript would be built from.
Each row is one thing the software could assert, the evidence that exists for it,
and the boundary beyond which the evidence does not reach. A claim with no
evidence row does not go in the paper.

The roadmap that produces this evidence is [ROADMAP.md](ROADMAP.md); the datasets
and thresholds are in [VALIDATION_PLAN.md](VALIDATION_PLAN.md); measured numbers
are in [VERIFICATION.md](VERIFICATION.md) and [dev/BENCHMARKS.md](dev/BENCHMARKS.md).

## The contribution, stated honestly

The software is an independent-module Python interface with shared sample metadata,
provenance, explicit optional dependencies and reusable plotting. It **combines
established methods**; it does not propose a new statistical method. Any manuscript
that claims a new method would be overstating what exists.

What is defensible is narrower and worth stating precisely: every analysis step is
a separate command reading files and writing files, each step records the tool and
version that produced it, and steps that are *not* validated say so in their own
output rather than in a footnote.

## Claim matrix

Status: `SUPPORTED` — evidence exists and is reproducible here.
`PARTIAL` — evidence exists but covers one dataset, one chromosome or one tool.
`NONE` — the claim cannot be made yet.

### Interface and reproducibility

| # | Claim | Evidence | Limit | Status |
| --- | --- | --- | --- | --- |
| I1 | The MACS3 wrapper does not alter parameters or output | Byte-identical peak files vs a direct `macs3 callpeak`, three configurations (BAM narrow, bedGraph narrow, bedGraph broad), `scripts/verify_macs3_concordance.py` | Says nothing about biological correctness; says nothing about any pure-Python helper | SUPPORTED |
| I2 | One command per analysis step, each reading and writing files | CLI surface check, 31 subcommands, `maintenance/verify_linux.py`; `tests/test_cli_help_coverage.py` requires help on every option | Structure only; no claim that the steps compose into a valid study design | SUPPORTED |
| I3 | Every run records the provenance of its inputs and tools | `easyepi/provenance.py`: input and output SHA-256, dependency versions, external tool versions and arguments in `status.json` | A tool absent from the run is recorded as absent; a tool whose banner cannot be parsed records `status: ok` with the raw banner | SUPPORTED |
| I4 | The installed wheel works outside the source tree | `scripts/verify_installed.py`; a clean-venv install run from another directory | Only on the hosts where it was run | PARTIAL |
| I5 | Linux CI executes and passes | none — the workflow has never run | **This is why "CI passes" cannot appear in a paper** | NONE |

### Concordance with established tools

| # | Claim | Evidence | Limit | Status |
| --- | --- | --- | --- | --- |
| C1 | Annotation agrees with ChIPseeker | 12 peak sets (3 donors × 4 histone marks), 23,646 peaks, every annotation string and distance identical on a shared gene model | Establishes equivalence *given the same gene model*; does not establish that either tool picks the right TSS set; does not extend to another module | SUPPORTED |
| C2 | The PWM scanner agrees with FIMO | 11 of 12 best-motif sites, Spearman 1.000 on the artificial set; FIMO 5.5.9 accepts easyepi-written MEME libraries (15 sites on one motif, 11324 across 469) | Never agrees on absolute score, because the two tools define the score differently; one dataset | PARTIAL |
| C3 | The TOBIAS ports reproduce the reference | Scoring 1.8e-15; cut sites and k-mers 18/18; sampled count matrices 0.000e+00; `fast_rolling_math` 1e-12 | Whole-track ATACorrect is 0.999992 forward but 0.989637 reverse, and the reverse residual is **unexplained** — recorded as a known limitation, not claimed as equivalence | PARTIAL |
| C4 | ChromHMM's model is reproduced structurally | End-to-end training on real ENCODE marks recovers the planted two-mark structure | No state call is compared to a published annotation; emissions depend on the seed | PARTIAL |
| C5 | HOMER's enrichment is reproduced | none — HOMER is used directly; the project no longer implements a substitute | `homer_like` was removed after measuring 3.90 bp mean shared consensus with HOMER's top ten | NONE |
| C6 | DMR statistics agree with DSS | none — the comparison tests skip because R and DSS are absent on the analysis host | The single largest evidence gap in the project | NONE |
| C7 | Differential expression agrees with DESeq2 | Sign, sign-flip on group swap, batch term and provenance pass on synthetic counts | No cohort-level comparison; no independent-filtering comparison | PARTIAL |

### Motif chain

| # | Claim | Evidence | Limit | Status |
| --- | --- | --- | --- | --- |
| M1 | Backgrounds are matched and the match is checked | `easyepi background` reports both GC and length distributions, a KS test on each, the matched fraction and the seed. `--accessible BED` additionally restricts controls to a supplied accessible universe, refusing a universe that wholly contains the targets. 300 real H3K27ac peaks, chr1: GC 0.520 vs 0.507 (KS p = 0.342), 100% matched; with an independent 2000-interval accessible set: 297/300 matched, GC 0.520 vs 0.508 (KS p = 0.248) | A universe that heavily *overlaps* the targets without wholly containing them is still partly circular and is not detected. Nothing yet shows that the choice of background changes an enrichment ranking | PARTIAL |
| M2 | Site-level scanning reports coordinates, strand and significance | `scan_sites` with 9 tests: coordinates, both strands, N-containing windows, short and empty input, palindromic dedup, agreement with `scan_sequences`; 115 sites on 300 real sequences in 0.33 s | Coordinates are offsets into the supplied sequence, not genomic positions; no full-library FIMO comparison at cohort scale | PARTIAL |
| M3 | Enrichment ranking is **not** stable across backgrounds, which is why the background must be stated | `scripts/bench_background_sensitivity.py`, 300 real H3K27ac peaks, 120 HOMER motifs: genomic vs GC/length-matched rho = **0.343** (top-10 overlap 2/10); the unmatched null's top hit AP-2gamma (log P 12.97) falls to rank 113 when matched | One peak set, one chromosome, one motif library. It shows the null dominates the ranking, not what the "right" ranking is | PARTIAL |
| M4 | De novo discovery reproduces a reference tool | HOMER 4.11 run end to end on real peaks: 16 de novo motifs and a populated `knownResults.txt` | One peak set; MEME/STREME not behind the same reported interface; no recovery rate | PARTIAL |
| M5 | TF activity is inferred by a published method | decoupler ULM on a user-supplied network, checked against the installed implementation | ULM is not chromVAR: no background-deviation correction. Motif presence is not binding, and motif enrichment is not activity | PARTIAL |

### Scope limits that apply to every row

| Limit | Consequence |
| --- | --- |
| Motif identity is shared between paralogues | A TF name from a similar motif is a hypothesis, not an identification |
| Nearest-TSS linkage is not regulatory inference | `integrate` output must not be described as a regulatory network |
| Nulls are sequence-composition nulls unless an accessibility background is supplied | p-values are conditional on that choice and change with it |
| Every cohort used is one cohort | No claim about generalisation across labs, antibodies or depths |
| MNase positions are not NucleoATAC or MINCE positions | Not interchangeable |
| Methylation site tests are not a DMR analysis | The pooled binomial model ignores between-replicate variance |

## What would change a row from PARTIAL to SUPPORTED

Each of these is a roadmap item; listed here so the evidence gap and the work are
the same list:

| Row | Needs |
| --- | --- |
| C2 | A cohort-scale site comparison against FIMO on a shared threshold |
| C3 | An explanation of the reverse-strand residual, or a statement that it is a reference-implementation artifact |
| C4 | A state call compared against a published ChromHMM annotation |
| C7 | A multi-donor DESeq2 comparison including independent filtering |
| M1 | An accessibility-matched background, plus the sensitivity comparison showing what it changes |
| M3 | The background sensitivity analysis in [VALIDATION_PLAN.md](VALIDATION_PLAN.md) §3.3 |
| M4 | A recovery rate across several peak sets, and MEME/STREME behind the same interface |
| I5 | A green GitHub Actions run, linked with its commit and skip count |

C5 and C6 have no route to SUPPORTED by reimplementation. For DSS the route is to
run the comparison; for HOMER the answer is to delegate, which is what the code
now does.

## Citation basis

Methods this software delegates to or reimplements, with the primary reference for
each. Full entries are in [REFERENCES.md](REFERENCES.md).

| Method | Reference |
| --- | --- |
| MACS peak calling | Zhang et al., 2008 |
| ChromHMM | Ernst & Kellis, 2012 |
| DSS methylation statistics | Feng et al., 2014 |
| DESeq2 | Love et al., 2014 |
| chromVAR (not implemented; cited to define what ULM is not) | Schep et al., 2017 |
| HOMER | Heinz et al., 2010 |
| MEME Suite | Bailey et al., 2009 |
| FIMO | Grant et al., 2011 |
