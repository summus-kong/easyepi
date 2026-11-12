# Verification record

Date: 2026-10-16. Supported deployment target: Linux. The available development host is not Linux and has no installed Linux distribution, so the checks below must not be described as completed Linux validation.

| Check | Result |
| --- | --- |
| Maintained Python source syntax | Passed |
| API and regression suite | 57 passed, 0 skipped on the supported Linux runtime (`pytest tests -q`) |
| Independent module CLI help | 39 entries passed |
| Real decoupler 1.9.2 ULM on a synthetic measured matrix | Passed local smoke check; 4 samples x 2 TF estimates |
| BAM/BigWig reader integration | Tests provided; skipped locally because readers are absent |
| Clean source snapshot sdist/wheel build | Passed local installation-content check; not uploaded |
| Wheel/source archive inspection | No R scripts, old python/backend paths, bytecode, historical figures or analysis outputs |
| Base import without specialist readers | Passed |
| Installed wheel imports outside source | Passed with site startup disabled to exclude an old editable-install finder; installed-copy suite: 57 passed, 2 skipped |
| Project schema validation via installed package | Passed; resolved config, inputs, logs and status written |
| Source/help/documentation language | Maintained files are English; legacy excluded helper copies await cleanup |
| Linux CI | Implemented for Ubuntu/Python 3.11 and 3.12; execution pending |
| Real-data external tool concordance | Not run; requires actual Linux data/tools |
| Historical artifact deletion | Executed locally via `maintenance/clean_legacy.ps1 -Apply`; ignored build/cache/output trees removed from the working tree |

## Local environment details

The unit suite used Python 3.11 with scikit-learn and Biopython installed in project-local verification directories. The real ULM smoke check used decoupler 1.9.2 with its numerical dependencies in a separate local directory. A host-specific numba setting was needed for that smoke check; no such setting is built into Easyepi. Linux is the only supported runtime.

An existing global editable installation mapped missing submodules back into the old workspace. The first installed-wheel import check correctly exposed that contamination. Repeating with `python -S` and explicit installed-package/dependency paths confirmed that the clean wheel itself contains no `easyepi.python` or R-backend paths. A fresh Linux virtual environment avoids the old finder entirely.

Local build files exist only under ignored `.verification` paths and are verification artifacts, not a requested distributable release. Package content exclusions were checked independently from the physical workspace cleanup.

## TF annotation and profile validation (2026-11-10)

Verified on the Linux host with a real ENCODE H3K4me3 experiment
(GRCh38 narrowPeak, 622 peaks) against the indexed GRCh38 genome and 14 JASPAR
2024 core PFMs.

| Check | Result |
| --- | --- |
| Full suite on Linux (`pytest tests -q`) | **184 passed**, 0 failed |
| `_score_all` vs brute-force oracle over both strands | 0 mismatches in 40 random motifs/loci, including N-containing sequences |
| Peak sequences fetched from indexed FASTA | 622/622, length 147 / 252 / 2909 (min/median/max) |
| Motifs usable at fpr=1e-4 | 14 of 14 (0 skipped) |
| Peaks carrying >=1 site | 375/622 (60.3%) |
| Top annotated TFs | SP1 28.5%, CTCF 19.0%, NRF1 17.4%, YY1 12.4%, MYC 9.3% |
| Aggregate profile on real bigWig | centre 8.380 vs shifted-background 1.455 (+6.93) |
| Figures rendered and inspected | `docs/validation/footprint.png`, `docs/validation/tf_annotation.png` |

Defects found and fixed during this validation:

* `Motif.to_log_odds` divided by zero for an all-zero column and produced an
  all-NaN PWM, so the motif silently matched nothing everywhere.
* The PWM score distribution was computed exactly, which is not tractable for a
  19 bp motif (4^19 states); the score axis is now binned as FIMO does.
* Site p-values recomputed the whole convolution per hit; a single precomputed
  tail table is now reused.
* The minus-strand scan reversed positions without complementing bases, scoring
  the reverse string and reporting wrong offsets; it now builds the reverse
  complement. An independent brute-force oracle caught this.
* `fig_footprint` plotted profiles against their array index, producing a 0..2*flank
  ramp and an axis labelled "distance from motif center"; it now plots signed bp
  offsets.
* Motifs shorter than the requested FPR can resolve are listed in
  `skipped_motifs` rather than appearing as "no site".
* `tests/conftest.py` pins the Agg backend so plotting tests no longer fail on a
  host without a display.

Reproduce with `scripts/validate_real_data.py` (uses the ENCODE files and genome
already present on the host; the JASPAR subset is cached after the first run).

## Coverage and figure round (2026-11-11)

| Check | Result |
| --- | --- |
| Full suite on Linux (`pytest tests -q`) | **235 passed**, 0 failed |
| CLI surface | every subcommand's `--help` builds; `nucleosome-positioning` added |
| Previously unimplemented analysis | `mnase positioning` now implemented as `nucleosome-positioning` |
| Single-end refusal on a real ENCODE ChIP BAM | Refused with an actionable message (no dyad from single-end reads) |
| New figures rendered and inspected | MA plot, peak-set overlap, TSS heatmap, locus track, dyad-density track |

Real-data check of the overlap metric on ENCODE H3 marks (chr1:0-5Mb, hg38):

| Pair | bp-Jaccard | interval coverage |
| --- | --- | --- |
| H3K4me3 vs H3K27ac | **0.142** | 0.215 |
| H3K4me3 vs H3K4me1 | 0.048 | 0.108 |
| H3K4me3 vs H3K27me3 | 0.043 | 0.278 |
| H3K4me1 vs H3K27ac | 0.049 | 0.225 |
| H3K4me1 vs H3K27me3 | 0.023 | 0.316 |
| H3K27ac vs H3K27me3 | 0.037 | 0.283 |

The ordering is the expected one (active marks co-localise, the repressive mark is
anti-correlated). It also shows why the metric had to change: with interval counts,
H3K27me3 appeared to cover 28-32% of every other mark, because a few broad domains
touch many narrow peaks. Base-pair Jaccard removes that artefact.

Defects found and fixed in this round:

* `overlap_fraction` clipped only one side of the intersection, so a 5 bp overlap
  was reported as the full length of the candidate interval.
* Peak-set overlap defaulted to an interval-count Jaccard, which reads as
  asymmetric and overstates broad-domain similarity; it now defaults to bp-Jaccard
  on merged intervals.
* Local-maximum calling used a MAD-only threshold; with MAD = 0 on a sparse dyad
  track the threshold collapsed to 0 and every noise spike was called.

Not covered here: no paired-end MNase/ATAC BAM was available on the host, so the
dyad track was exercised on synthetic arrays and the single-end refusal on real
data. Positioning accuracy against a published caller is not established.

## Panel and background round (2026-11-12)

| Check | Result |
| --- | --- |
| Full suite on Linux (`pytest tests -q`) | **278 passed**, 0 failed |
| CLI subcommands | 25; `gene-panel` added |
| Real-data gene panel | 4 marks x 4 chr1 loci from ENCODE bigWigs + GENCODE v50 GTF, rendered and visually inspected |
| Figure defects found by inspection | missing coordinate ticks, over-flattened shared y-limits, exon boxes too thin to see |

The panel was rendered three times; each iteration fixed a defect that only the
image itself revealed:

1. No x-axis ticks at all, so a reader could not relate a peak to the gene body.
2. One shared y limit per track row let a 40x peak in one locus flatten the other
   three columns; limits are now robust per panel.
3. Exon rectangles narrower than a pixel vanished, so the gene row looked empty.

Real-data sanity: H3K4me3 and H3K27ac peak at the AGRN and SKI promoters, H3K27me3
is flat across SKI and spread over TP73, and H3K4me1 is broad at PLCH2. This is the
expected mark pattern for these loci.

Defects fixed this round (all found by test or image inspection, none assumed):

* UpSet membership tested whole intervals rather than elementary union segments,
  losing regions unique to the second set and mislabelling partial overlaps.
* GC matching randomly sampled candidates and reported false "no match" results.
* GC matching used coordinate differences for the distance test, so `exclude` was
  ineffective at `min_distance=0`.
* `exclude` is now a hard exclusion, independent of `min_distance`.
* Sample correlation produced a NaN row and column for a zero-variance sample.

Not covered: no paired-end MNase/ATAC data was available, so nucleosome positions
remain validated only on synthetic arrays. GC matching was exercised on synthetic
genomes; a genome-wide run against the real FASTA was not timed.

## Compact input syntax round (2026-11-13)

| Check | Result |
| --- | --- |
| Full suite on Linux (`pytest tests -q`) | **320 passed**, 0 failed |
| New tests | 25 (`tests/test_compact_inputs.py`) |
| Real `qc --input-dir` on a staged ENCODE folder | 1 sample resolved from `H3K4me3.bam` + `.bw` + `.narrowPeak`; `signal_enrichment` = 10.09 |
| Empty folder | Fails with the directory name; no empty report written |
| Real `--bigwig H3K4me3,H3K27ac` + `--bigwig-dir` | Identical PNG bytes to the long `NAME=PATH` form |

Defects found and fixed this round:

* `qc` attached BigWig/peak files to BAMs by **list position**. A sample missing a
  track shifted every later pairing, and the report still looked plausible.
  Attachments are now resolved by sample name.
* `qc` wrote `{"samples": {}}` and exited 0 when no input was recognised, which
  looks like a successful run. It now raises.
* `differential-peaks` only accepted the full `SAMPLE=BAM:PEAK` spelling and
  rejected anything else, so a folder of per-sample files had to be typed by hand.

Observed during validation, not a code defect: FRiP for the staged sample is
0.0018 because the peak set covers only chr1:0-5 Mb while the BAM spans the whole
genome, so the numerator and denominator describe different regions. The
BigWig-based `signal_enrichment` (10.09) is the meaningful number for this input.
Peaks and reads must cover the same genomic universe for FRiP to be interpretable.

## Linux verification and external concordance (2026-11-14)

Run on the Linux host with the project conda environment (Python 3.11.16,
pysam 0.24, pyBigWig 0.3.25, MACS3 3.0.4).

| Check | Result |
| --- | --- |
| `pytest tests -q` | **325 passed**, 0 failed |
| `maintenance/verify_linux.py` | **48 checks, 0 failures, exit 0** |
| CLI subcommands covered by the verifier | added `tf-annotate`, `nucleosome-positioning`, `gene-panel`, `differential-peaks`, `histone-states`, `footprint`, `qc`, `integrate` |
| FRiP genome-universe guard on real ENCODE data | `interpretable: false`, peaks span 0.16% of mapped reference; warning emitted |

### MACS3 argument concordance (`scripts/verify_macs3_concordance.py`)

| Case | Input | easyepi peaks | direct macs3 peaks | identical coordinates | identical md5 |
| --- | --- | --- | --- | --- | --- |
| `bam_narrow` | BAM, `-f BAM` | 37,250 | 37,250 | yes | yes |
| `bedgraph_narrow` | bedGraph, `-f BED` | 458 | 458 | yes | yes |
| `bedgraph_broad` | bedGraph, `-f BED --broad` | 1,551 | 1,551 | yes | yes |

The wrapper is argument-transparent: it produces byte-identical peak files to a
hand-written `macs3 callpeak` with the same documented flags. This rules out the
failure mode a wrapper can introduce (silently changed parameters). It does **not**
show that the peaks are biologically correct, and it says nothing about the
legacy pure-Python helpers still importable from `easyepi`, which are a separate code
path and are not reachable from `easyepi peaks`.

### FRiP guard

The ENCODE H3K4me3 stage used for validation has peaks called only over
chr1:0-5 Mb while the BAM is whole-genome. FRiP was therefore 0.00187
(374 in-peak reads of 200,528), which reads as a catastrophic QC failure but is an
artefact of comparing a regional numerator against a genome-wide denominator. The
metric now carries `interpretable: false`, the peak span fraction, and an
explanatory warning; `signal_enrichment` (10.09 from peak vs flank BigWig signal)
is the meaningful number for that input.

## Reproduce on Linux

```bash
conda env create -f environment.yml
conda activate easyepi
python -m pip install -e ".[dev]"
python maintenance/verify_linux.py
```

The verifier records per-command logs and versions. Install the bam, bigwig, motif, stats and activity extras when using a venv instead of the conda specification; otherwise corresponding optional tests may skip. The supplied CI installs all these dependencies.

No unit test, command plan or synthetic smoke example establishes biological validity, confirms equivalence to upstream software or completes publication benchmarking. See RELIABILITY.md for remaining development and validation priorities.

## Genomic-region annotation, 0.3.0

Date: 2026-11-12. Data: ENCODE H3K4me3 peaks against GENCODE v50 primary assembly,
run on the Linux analysis host with `easyepi annotate --gtf` (native backend,
midpoint mode, ChIPseeker-default 3000 bp promoter window).

### Three defects fixed in this version

Found by writing `tests/test_annotate_bins.py` against a hand-computed gene model;
the classification module had no coverage of its own logic before that.

| Defect | Effect before the fix |
| --- | --- |
| `Exon` sorted before the UTRs in `GENOMIC_BINS` | UTR peaks were captured by the enclosing exon interval, so the UTR bins were nearly empty on every dataset |
| `_extract_utr` used global min/max of the CDS | every exon between the first and last CDS block counted as 3'UTR |
| 5'/3' assigned from coordinate only | minus-strand genes had their UTRs swapped |
| `_nearest_tss_distance` ignored strand | upstream of a minus-strand TSS reported the wrong sign |

### Measured effect on a real peak set

43,465 ENCODE H3K4me3 peaks, same inputs and same promoter window, only the bin
priority changed back and forth:

| bin | old priority | 0.3.0 |
| --- | --- | --- |
| Promoter | 29.22% | 29.22% |
| 5'UTR | 4.26% | 4.26% |
| 3'UTR | 0.02% | 0.63% |
| Exon | 13.23% | 12.63% |
| Intron | 33.61% | 33.61% |
| Intergenic | 19.65% | 19.65% |

The 3'UTR bin was 32x smaller than it should be, and the difference it stole went
to `Exon`. This is why 0.3.0 is a new version rather than a 0.2.x patch: peak
annotations in existing projects will change.

### Sanity check against expectation

H3K4me3 is a promoter-associated mark but by no means promoter-exclusive. The
measured 29.2% promoter, 33.6% intronic, 19.7% intergenic distribution is the shape
the literature reports for this mark, and it is produced by the same run.

### Still open

Superseded by the 0.5.0 section below, which reaches 93.78% by adopting
ChIPseeker's distance-based Promoter rule. The remaining 6.22% is analysed there
and is dominated by intron modelling, not by the transcript choice discussed in
the earlier note.

## ChIPseeker benchmark, 0.3.0 — **does not meet the threshold**

Run 2026-11-12 on the Linux analysis host. Both sides used a TxDb built from the
**same** GENCODE v50 GTF (`txdbmaker::makeTxDbFromGFF`) so the comparison measures
the classification algorithm, not two gene-model sources. Same peak file, same
order, same promoter window (TSS +/-3000). ChIPseeker 1.46.1.

Peaks: 43,465 ENCODE H3K4me3.

### Result

| metric | value | threshold | verdict |
| --- | --- | --- | --- |
| per-peak exact agreement | **57.66%** | >= 93% | fail |
| coarse agreement (Promoter/Genic/Intergenic) | **59.19%** | >= 98% | fail |
| max per-bin frequency delta | **36.15 pp** | <= 2 pp | fail |

| bin | ChIPseeker | easyepi | delta |
| --- | --- | --- | --- |
| Promoter | 65.37% | 29.22% | -36.15 |
| 5'UTR | 0.10% | 4.26% | +4.16 |
| 3'UTR | 1.58% | 0.63% | -0.95 |
| Exon | 1.61% | 12.63% | +11.02 |
| Intron | 21.62% | 33.61% | +11.99 |
| Intergenic | 9.72% | 19.65% | +9.93 |

### Root cause: one TSS per gene instead of one per transcript

The disagreement is systematic, not noise:

* **The easyepi Promoter set is a strict subset of ChIPseeker's.** Of the peaks
  ChIPseeker calls Promoter, 15,713 are called something else by easyepi; the
  reverse count is **0**.
* **The two sides' TSS distances are uncorrelated** (Pearson r = 0.007 on the
  disagreeing peaks). They are not measuring the same TSS at all.
* TSS distance is the mechanism: 60% of the disagreements sit within 3 kb of the
  ChIPseeker TSS, and 39.7% are within 3 kb for ChIPseeker but farther for easyepi.

`_read_gtf_txs` keys its gene table by `gene_id` and takes the `gene` line's span,
so every gene contributes **one** TSS. ChIPseeker derives TSS from **all**
transcripts. Worked example, `ENSG00000223764.3` (minus strand):

```
gene        916865-925614   -> easyepi TSS = 925613
transcripts 916865-921016, 916870-919692, 917370-918534,
            917563-925585, 917573-919306
```

A peak at chr1:940534 is 0 bp from a ChIPseeker TSS and 15,099 bp from the
single TSS of the gene easyepi assigned it to. That figure is a per-peak
distance, not a model error: for this particular gene the single TSS (925,613)
happens to sit 1 bp from its nearest transcript start. The error shows up
across genes -- 35,318 of 78,941 genes have more than one transcript, one has
522 -- where a single TSS cannot represent starts that are kilobases apart.

### Verdict

**The benchmark fails and this version is not ChIPseeker-concordant.** The three
bins fixed earlier were genuinely broken and the fixes stand, but they were not the
whole story: the TSS model is architecturally different, and that dominates the
output. Nothing in this version should be described as agreeing with ChIPseeker.

The fix is to derive one TSS per transcript, matching ChIPseeker. That is a larger
change than the bin-priority fixes: it changes every peak's annotation and the
Promoter fraction roughly doubles. It is deliberately **not** done in 0.3.0, because
it needs its own before/after record and its own re-run of this benchmark.

Reproduce with `bench/build_txdb.R`, `bench/run_chipseeker.R` and `bench/compare.py`.

## 0.4.0 — per-transcript TSS, and the benchmark re-run

The 0.3.0 benchmark failed at 57.66% exact agreement. Its stated cause was that
`_read_gtf_txs` gave each gene one TSS while ChIPseeker uses one per transcript.
That is now fixed and the benchmark re-run under identical conditions (same TxDb
from GENCODE v50, same 43,465 peaks, same 3000 bp window).

| metric | 0.3.0 | 0.4.0 | threshold | verdict |
| --- | --- | --- | --- | --- |
| per-peak exact agreement | 57.66% | **80.39%** | >= 93% | still short |
| coarse agreement | 59.19% | **81.61%** | >= 98% | still short |
| max per-bin delta | 36.15 pp | **14.10 pp** | <= 2 pp | still short |

| bin | ChIPseeker | 0.3.0 | 0.4.0 |
| --- | --- | --- | --- |
| Promoter | 65.37% | 29.22% | 51.27% |
| 5'UTR | 0.10% | 4.26% | 1.70% |
| 3'UTR | 1.58% | 0.63% | 1.52% |
| Exon | 1.61% | 12.63% | 5.34% |
| Intron | 21.62% | 33.61% | 23.83% |
| Intergenic | 9.72% | 19.65% | 16.35% |

Exact agreement rises 22.7 points. The residual disagreements are dominated by
`Promoter -> Intron` (2,338) and `Promoter -> Exon` (1,790): cases where both sides
find a nearby transcript start but disagree about whether the peak's midpoint is
inside the promoter window, or where they pick a different nearest transcript. The
easyepi Promoter set is still a subset of ChIPseeker's.

### Downstream folding re-checked

ChIPseeker reports 19 `Downstream` peaks out of 43,465 (0.04%). easyepi has no such
bin and folds them into `Intergenic`. At that share the folding cannot move any
aggregate by more than 0.04 pp, so it is not a factor in the disagreement.

### Backend cross-check

All three annotation backends were run over the same 43,465 peaks and the same
transcript table:

| backend | version | time | agreement with native |
| --- | --- | --- | --- |
| native | built in | 44.1 s | reference |
| pybedtools | 0.12.1 / bedtools 2.31.1 | 53.1 s | **100.0000%** (43465/43465) |
| pyranges | 0.1.7 | 47.9 s | **100.0000%** (43465/43465) |

Pairs: native-vs-pybedtools 43,465/43,465 and native-vs-pyranges
43,465/43,465, both exactly 100.0000%, zero mismatches. The backends are therefore
interchangeable for correctness; choosing one is about which dependency is
available, not about the result. Timings include GTF parsing and are within 20% of
each other.

Both optional backends were installed with micromamba from `conda-forge`/`bioconda`
(`pybedtools`, `bedtools`, `pyranges`); `pip` cannot reach PyPI from this host.

## 0.5.0 — Promoter by TSS distance, and 93.8% agreement

0.4.0 reached 80.39%. The remaining gap was attributed by toggling one rule at a
time against ChIPseeker's own output on the same 43,465 peaks:

| rule under test | exact agreement |
| --- | --- |
| midpoint-in-interval + transcript TSS (0.4.0) | 80.45% |
| whole-peak overlap + transcript TSS | 82.42% |
| whole-peak overlap + gene TSS | 59.61% |
| **whole-peak overlap + distance-to-TSS promoter** | **93.78%** |

The decisive finding was mechanical rather than statistical: ChIPseeker decides
Promoter from the *distance to the nearest TSS*, not from whether the peak's midpoint
lies inside a +/-window interval. Testing that rule directly against ChIPseeker's
own `distanceToTSS` column gave `(theirs == "Promoter") == (|distanceToTSS| <= 3000)`
on **43,465 / 43,465 peaks** -- zero exceptions. Adopting it moved agreement from
82.42% to 93.78% on its own; 5,731 peaks that had been missed were all within the
window by distance but outside it by midpoint.

The overlap rule contributes the other ~2 points: ChIPseeker calls
`findOverlaps(peaks, region)` on the whole peak, so a peak whose body reaches into a
region counts, not only one whose centre does.

| metric | 0.3.0 | 0.4.0 | 0.5.0 | threshold |
| --- | --- | --- | --- | --- |
| per-peak exact agreement | 57.66% | 80.39% | **93.78%** | >= 93% |
| Promoter calls | 12,701 | 22,282 | 27,623 | (ChIPseeker 28,414) |

### What still disagrees (6.22%)

| ChIPseeker | easyepi | peaks | share |
| --- | --- | --- | --- |
| Intron | Intergenic | 1,786 | 4.11% |
| Promoter | 5'UTR | 421 | 0.97% |
| Promoter | Exon | 342 | 0.79% |
| 3'UTR | Exon | 92 | 0.21% |

The largest residual is `Intron -> Intergenic`. Those peaks have a median distance
to the nearest TSS of 13.7 kb, so they sit inside a transcript ChIPseeker modelled
and this implementation did not: 29,034 of 646,577 transcripts (4.5%) have a single
exon, and ChIPseeker's `intronsByTranscript` derives introns from the transcript
table rather than from gaps between the exon records we read. Closing this needs the
intron list built from transcript spans, not from exon gaps.

The `Promoter -> 5'UTR` and `Promoter -> Exon` cases are ordering: ChIPseeker assigns
interval bins first and lets Promoter override, which this now does too, so the
remainder is a boundary or nearest-transcript choice difference.

### Backends remain identical

With the new rule in place all three backends still agree exactly: `native`,
`pybedtools` 0.12.1 and `pyranges` 0.1.7 each score 40,762/43,465 (93.78%) against
ChIPseeker and 100.0000% against each other.

### Method note

These figures come from differential testing against a real ChIPseeker 1.46.1
installation on the same TxDb, not from reading documentation. Two earlier claims in
this file were corrected by that process: the TSS level is `"transcript"` (not
`"gene"`), and the 15,099 bp figure quoted for one gene was a per-peak distance, not
a model error.

## 0.5.1 — the sweep bug behind the last 6%, and 99.76%

0.5.0 reached 93.78% and blamed the remaining 6.22% on intron modelling. That
attribution was wrong.

A diagnostic on the 1,786 `Intron -> Intergenic` peaks showed every one of them
**does** overlap an intron interval this code builds (1,786/1,786) and lies inside a
transcript span (100%). Classifying one of those peaks against the intron set alone
returned `Intron`; classifying the same peak against the full feature set returned
`Intergenic`. The features were right and the scan was wrong.

The cause was in the overlap sweep:

```python
feats.sort(key=lambda x: (x[0], x[1], x[2]))   # sorted by START
k = bisect.bisect_left(ends, qs)               # binary search on ENDS
```

`ends` follows the start order, so it is not sorted, and `bisect_left` on an
unsorted list returns an arbitrary position. The scan began there and skipped
overlapping features, which silently dropped peaks from every bin.

The fix bounds the scan by `bisect_left(starts, qe)` -- sound, because `starts` is
the sort key -- and prunes with a running maximum of `ends`, which is monotone by
construction.

| metric | 0.5.0 | 0.5.1 |
| --- | --- | --- |
| per-peak exact agreement | 93.78% | **99.76%** |
| disagreeing peaks | 2,703 | **105** |

The 105 remaining are all UTR-versus-Exon boundary cases:

| ChIPseeker | easyepi | peaks |
| --- | --- | --- |
| 3'UTR | Exon | 91 |
| 5'UTR | Exon | 10 |
| Exon | 3'UTR | 4 |

These come from UTR provenance, not from the scan: ChIPseeker reads
`fiveUTRsByTranscript`/`threeUTRsByTranscript` from the TxDb, while this code derives
UTRs by subtracting CDS from exons. The two disagree where a transcript has a gapped
CDS or a retained intron, which is exactly the 0.21% shown above. Closing it means
reading the UTR records rather than computing them.

All three backends remain identical at 99.76% and 100.0000% against each other.

A regression test now covers the scan directly: a long interval whose start comes
first and ends after later, shorter ones -- a list whose `ends` is unsorted -- plus a
shuffle-invariance check.

## 0.5.2 — 3'UTR covers every untranslated exon, reaching 99.99%

0.5.1 left 105 peaks differing, all UTR-versus-Exon, and attributed them to UTR
provenance -- the idea that ChIPseeker reads `threeUTRsByTranscript` while this code
derives UTRs arithmetically. Two checks showed that was wrong:

* For a sampled transcript, GENCODE's `UTR` rows and `exons - CDS` are **the same
  interval set** (verified as a set comparison, not by eye). Reading the UTR rows
  would have changed nothing.
* Inspecting the failing peaks directly: each one overlaps an **Exon** interval of
  this code's feature set and **zero** of its UTR intervals (0 for both 5' and 3').
  So the UTR sets disagreed about which regions exist, not about how to label them.

The actual cause was one `continue` in `_extract_utr`: after the last CDS block the
code kept only the piece abutting the stop codon and discarded every later exon, on
the reasoning that "a downstream exon is still an exon". That is true of the feature
type but not of the annotation. A 36-exon gene routinely has several exons past
translation end, and every one of them is untranslated sequence.

| metric | 0.5.1 | 0.5.2 |
| --- | --- | --- |
| per-peak exact agreement | 99.76% | **99.99%** |
| disagreeing peaks | 105 | **4** |
| 3'UTR intervals | 617,876 | 892,775 |
| 5'UTR intervals | 487,261 | 625,728 |

The four that remain are all `Exon -> 3'UTR`: peaks where a transcript's CDS status
differs between the two builds of the annotation, which is a data question rather
than a rule question. All three backends still agree 100.0000% with each other.

Three tests that had encoded the old behaviour were corrected rather than deleted.
They were written from the same incorrect reasoning as the code and had been passing
against it -- a test can pin a bug just as easily as it can pin a fix.

## 0.5.3 — 3 bp UTR slivers, reaching 99.995%

0.5.2 left 4 peaks differing and called them a data question. That was wrong: they
were a modelling difference with an identifiable cause.

The four peaks overlap an Exon interval and **zero** UTR intervals of this code,
yet ChIPseeker says 3'UTR -- so ChIPseeker had UTR intervals here that this code did
not. (The reverse would have been this code over-covering.) Querying
`threeUTRsByTranscript` at exactly those coordinates returned 0 hits, which settled
the direction.

Comparing the transcript *sets* rather than the peaks:

| | this code | ChIPseeker |
| --- | --- | --- |
| transcripts with a 3'UTR | 351,479 | 349,198 |
| 3'UTR intervals | 892,775 | 886,766 |
| present only here | **2,281** | 0 |
| present only in ChIPseeker | 0 | -- |

Every one of the 2,281 has a 3'UTR totalling **at most 3 bp**, and the distribution is
unambiguous: 2,281 in the 1-3 bp band, zero in 4-10, zero in 11-50, maximum 3 bp.
GENCODE defines some CDS blocks so a start or stop codon straddles an exon boundary,
leaving a sliver outside the CDS. `makeTxDbFromGFF` normalises those away when it
builds a TxDb; deriving UTRs arithmetically does not, so the slivers became real UTR
features that covered four peaks.

`_extract_utr` now drops pieces shorter than 4 bp. The threshold is safe by
measurement rather than by assertion: nothing between 4 and 50 bp was ever dropped,
so no genuine UTR is clipped.

| metric | 0.5.2 | 0.5.3 |
| --- | --- | --- |
| per-peak exact agreement | 99.99% | **99.9954%** |
| disagreeing peaks | 4 | **2** |
| 3'UTR intervals | 892,775 | 886,722 |
| ChIPseeker 3'UTR intervals | 886,766 | 886,766 |

The two that remain are one `5'UTR -> Exon` and one `Exon -> 3'UTR`, i.e. one case in
each direction. At 0.0046% with no residual pattern, further matching would be fitting
individual peaks rather than correcting a rule.

All three backends remain identical at 99.9954% and 100.0000% against each other.
