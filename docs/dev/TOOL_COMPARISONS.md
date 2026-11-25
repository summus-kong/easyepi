# External tool comparison: MACS3, FIMO, TOBIAS

Backlog item 12 asks whether `easyepi` agrees with the reference tools it calls.
This page records what was measured, and — just as important — what each
comparison cannot show.

Measured with **MACS3 3.0.4**, **MEME Suite 5.5.9 (FIMO)** and **TOBIAS 0.17.5**
in the `easyepi` conda environment, on the analysis host.

Reproduce:

```bash
python -m pytest tests/test_peaks_is_macs3.py       -q
python -m pytest tests/test_motifs_fimo_agreement.py -q
python -m pytest tests/test_tobias_command.py        -q
python tests/oracle/compare_to_fimo.py               # prints the full table
```

## MACS3 — `easyepi peaks` *is* MACS3

The first thing to establish is architectural, because it decides what a
comparison can even mean: `easyepi peaks` is not a reimplementation. It imports
`MACS3.Commands.callpeak_cmd` and calls `run()` on a Namespace built from MACS3's
own argument parser (`easyepi/macs_options.py` reads `prepare_argparser` out of
the `macs3` console script). `peaks_cli.py` states the policy in its docstring:
*"MACS3 is the only backend."* An earlier hand-written Poisson caller was
retired; the message that retired it records that it found 37% of MACS3's peaks
while emitting 1000 bp peaks against MACS3's 191 bp.

So "output equivalence" is not a property to measure — it is the definition. What
*is* worth pinning:

| property | test | result |
|---|---|---|
| the run produces MACS3's own peak file | `test_the_run_produces_macs3s_own_peak_file` | coordinates identical to a direct `macs3 callpeak` |
| `--macs-arg` reaches MACS3 | `test_passthrough_reaches_the_namespace_when_the_parser_is_available` | extra options appear in the parsed Namespace |
| a rejected option is reported, not ignored | `test_a_rejected_passthrough_names_the_tokens` | `ValueError` naming the tokens |
| no silent drop when the parser is unavailable | `test_passthrough_is_refused_when_the_parser_cannot_be_read` | refuses with the reason |

The last one fixed a real defect: the fallback Namespace is hand-written with a
fixed attribute set, so `--macs-arg` tokens were accepted, recorded in the
provenance, and then **silently ignored**. A run reported a `passthrough` that
never reached the calculation. It now refuses instead of pretending.

## FIMO — best window and ranking

Here the two are genuinely different implementations: `easyepi.motifs.read_pwm`
builds `log2(p_i / 0.25)` and `scan_sequences` reports the argmax window's sum;
FIMO scores every window and converts to a p-value.

**Best window agrees 11/12** across a 12-sequence fixture with a planted
consensus, its reverse complement, and noise at varying offsets. The single
disagreement is a noise-only sequence where both callers land within 1 bp of a
weak tie.

**Ranking agrees perfectly**: Spearman ρ = 1.000 on the best-window score.

### The absolute scores are *not* comparable, and this was measured

On four fixed 7-mers the difference changes sign with base composition:

| k-mer | easyepi | FIMO | delta |
|---|---|---|---|
| GATATCA (ML consensus) | 12.5785 | 12.1807 | +0.40 |
| GATATCG | 9.5798 | 9.4699 | +0.11 |
| AAAACGT | −10.6206 | −10.8976 | +0.26 |
| GGGGGGG | −13.1243 | −12.0060 | −1.14 |

A constant offset, a pseudocount (0.01 shifts by only 0.02), and a
background-smoothing blend were each ruled out. FIMO's `score` column is its own
statistic, not `read_pwm`'s log-odds sum. Establishing its exact definition would
mean reading MEME Suite's source; it is not pursued because **no downstream
decision uses the absolute score** — a motif call comes from a threshold or a
ranking.

Two harness bugs were found and fixed along the way, both of which would have
produced a meaningless comparison: a MEME file was handed to `read_pwm` (which
documents JASPAR input), and the JASPAR writer transposed the matrix the wrong way
so a 7-position motif was silently compared as 4 positions.
`test_the_matrix_is_not_truncated` now guards the second.

## TOBIAS — the wrapper named a subcommand that does not exist

`run_tobias_footprint` built `TOBIAS Footprint --cores 4 ...`. TOBIAS 0.17.5
rejects it at argument parsing:

```
TOBIAS: error: argument : invalid choice: 'Footprint'
(choose from 'ATACorrect', 'ScoreBigwig', 'FootprintScores', 'BINDetect', ...)
```

and `FootprintScores --help` states in its own banner that **ScoreBigwig replaces
the previous FootprintScores command**, with usage
`ScoreBigwig --signal <cutsites.bw> --regions <regions.bed> --output <output.bw>`
— no `--cores`. So two defects in one call: a non-existent subcommand and an
unsupported flag. Nothing tested it, so it could never have worked.

Fixed to `ScoreBigwig` with the three documented options.
`tests/test_tobias_command.py` pins the subcommand against TOBIAS's own accepted
list, pins the absence of `--cores`, and invokes the real binary so a future
rename fails loudly.

### The two scores are the same quantity in different units

`easyepi footprint` also reports a per-locus footprint score, so a second
implementation does exist and the numbers can be compared. Read from TOBIAS
0.17.5's own Cython source (`tobias/utils/signals.pyx`,
`tobias_footprint_array`):

```
fp_sum     = sum of the NEGATIVE values inside the footprint window
fp_mean    = fp_sum / footprint_w                 (<= 0)
flank_mean = (left_sum + right_sum) / (2 * flank_w)
footprint  = flank_mean - fp_mean                 additive, raw-signal units
```

and `easyepi.footprint.footprint_profile` computes

```
score = (shoulder_baseline - centre) / shoulder_baseline      fractional
```

With the same windows these are related by exactly

```
TOBIAS == flank_level * easyepi
```

**Measured: ratio = 1.000000 on all seven sites** (min = max = mean),
`tests/test_tobias_footprint_score.py::test_tobias_score_equals_flank_level_times_easyepi`.
The relation is an identity, not an approximation, so it is asserted at
`rel=1e-6`.

### Why an obvious fixture silently produced a constant

TOBIAS accumulates **only negative values** into `fp_sum`. A first probe used a
purely positive dip (flanks 10.0, centre 3.0) — the shape an ATAC track has
before correction — and every output value came back as exactly 10.0. That is
TOBIAS's definition working as written, not a failure: with no negative values
`fp_mean = 0` and `footprint = flank_mean`. A comparison therefore has to use an
ATACorrect-style track where the footprint is negative.
`test_purely_positive_dip_leaves_tobias_at_the_flank_level` pins this so the
constraint is documented rather than rediscovered.

## ChromHMM — input contract and a trained model

ChromHMM 1.28 is installed (`share/chromhmm-1.28-0/ChromHMM.jar`, Java 8) and both
stages run: `BinarizeBed` on real hg38 chromosome lengths from
`reference/human/genome.fa.fai`, then `LearnModel`.

### The trained model recovers the planted structure

Measured with two marks over **disjoint** regions of chr1 (H3K4me3 at 1–20 Mb,
H3K27ac at 30–60 Mb) and `-peaks`, asking for three states:

* 11 artifacts written (`model_3.txt`, `emissions_3.txt`, `transitions_3.txt`,
  `ct_3_segments.bed`, dense/expanded BEDs, PNG/SVG emissions and transitions)
* 9801 segments; **all three states used** (`E1`, `E2`, `E3`)
* the emission matrix separates the marks exactly as planted:

| state | H3K27ac | H3K4me3 |
|---|---|---|
| 1 | 2.9e-25 | 3.0e-27 |
| 2 | **0.999999999** | 4.0e-222 |
| 3 | 6.4e-211 | **0.999999999** |

So training is verified as *recovering the signal*, not against a golden number:
the assertions require two different states to dominate the two marks and a third
to emit neither.

### Two ways to get there were wrong at first

`LearnModel` aborts with

    On this data the INFORMATION initialization strategy can only support 2 states

if the fixture binarizes one mark and asks for three. Two corrections were
needed, both matching what `easyepi/specific/histone/chromhmm.py` emits: pass
`-peaks` (ChIP marks are peak calls, not raw signal) and give the model two marks
with disjoint occupancy. A single mark has only present/absent structure, which
is why the reported ceiling is 2.

### The input contract, because four rules bite silently

| rule | what happens when broken |
|---|---|
| cellmarkfiletable is the 3rd positional arg | with no args: usage banner, **exit 0**; with some: `-b` is read as the lengths file and throws `FileNotFoundException: -b` |
| file column is a **bare filename** | a path is concatenated with the bed dir, giving `<bed>/<bed>/x.bed` |
| column 2's mark name must match the file | `Warning did not find data for ... treating as missing`, zero bins |
| BED needs a strand column (6 cols) | `strand column expected, but not found` |

The bare-filename rule costs the most time: a caller who writes the natural thing
— the full path they already hold — gets a `FileNotFoundException` naming a
doubled path, which reads like a ChromHMM defect.

## HOMER — FASTA mode works, genome mode does not

`findMotifsGenome.pl`, `findMotifs.pl` and `homer2` are all present, and the
motif library under `share/homer/motifs` is complete. But
`share/homer/data/genomes` **does not exist**, so `findMotifsGenome.pl` has no
genome to scan and the genome mode is unusable in this environment.

`easyepi.motifs.external_tools.run_homer_motifs` therefore has exactly one
usable path here: FASTA mode, which requires an explicit background FASTA and
refuses without one (`ValueError`). With target and background it runs and
writes `homerResults` / `knownResults`.

Install for genome mode: `$ENV/share/homer/configureHomer.pl -install hg38`.

## What none of this establishes

- **Not a real cohort.** All three fixtures are synthetic with known ground truth.
  The backlog item is titled "on a real cohort"; that remains TODO.
- **Not MACS3's default model.** The comparison used `--nomodel --extsize 200` to
  match easyepi's fixed bins. MACS3's empirical fragment model is not covered.
- **Not motif-database competition.** One PWM per run, so FIMO never chooses
  between motifs.
- **Not enrichment statistics.** FIMO's q-value and easyepi's HOMER-like
  enrichment are different statistics on different inputs.
- **TOBIAS is verified as a score relation, not on real data.** The identity
  `TOBIAS == flank_level * easyepi` holds on a synthetic track where both windows
  are configured to match (`--flank-min/--flank-max/--fp-min/--fp-max`). TOBIAS's
  own window search, `ATACorrect` normalisation, and `BINDetect` differential
  testing are not covered.

- **ChromHMM's model is verified structurally, not against a reference
  segmentation.** Training recovers the planted two-mark structure, but no state
  call is compared to a published annotation, and the emission/transition values
  depend on the seed.
- **HOMER's enrichment is not compared.** FASTA mode runs, but HOMER's
  hypergeometric enrichment is not reproduced or cross-checked by easyepi's own
  `homer_like` implementation.

## Status against backlog item 12

| tool | status |
|---|---|
| MACS3 | **DONE** — equivalence is architectural; passthrough and refusal are pinned |
| MEME/FIMO | **DONE** — best window 11/12, ranking ρ = 1.000; score scale documented as incomparable |
| TOBIAS | **DONE** - broken command fixed and pinned; score identity `TOBIAS = flank_level x easyepi` measured at ratio 1.000000 |
| ChromHMM | **DONE** — BinarizeBed and LearnModel run on real hg38 lengths; trained model recovers the planted two-mark structure; input contract pinned |
| HOMER | **PARTIAL** — FASTA mode runs and is pinned; genome mode blocked (no installed genome) |
| real cohort | TODO |
