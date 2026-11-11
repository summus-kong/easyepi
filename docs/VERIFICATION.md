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
