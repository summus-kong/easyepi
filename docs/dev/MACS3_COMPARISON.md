# MACS3 comparison

Backlog item 12 asks whether `easyepi peaks` agrees with the reference callers it
claims to replace. This page records the MACS3 half: what was compared, what the
numbers were, and what the numbers do **not** establish.

Reproduce with:

```bash
python tests/oracle/compare_to_macs3.py     # full table
python -m pytest tests/test_peaks_macs3_agreement.py -q   # the assertions
```

Measured with **MACS3 3.0.4** and `easyepi` at commit `a961902`, on the analysis
host, in the `easyepi` conda environment.

## The two callers are different algorithms

A low overlap would not by itself prove a defect, and an identity assertion would
be wrong. The models differ in a way that matters:

| | easyepi | MACS3 |
|---|---|---|
| background | median of the **whole** binned signal | **local** lambda from genome-wide + 5 kb windows |
| test | Poisson upper tail per bin, merged | Poisson tail from local lambda |
| resolution | fixed bin (200 bp default) | base pair, with `--extsize` on the fragment model |
| multiple testing | BH over merged peaks | q-value from the same tail |

So the comparison asks the questions that survive the model difference:
recovery of regions known to be there, invention of regions that are not, and
how far apart the edges land.

## Fixture

A 400 kb `chr1` with **6 planted peaks** of known centre and half-width, over
uniform background, written as a paired-end BAM (141,928 read starts). Planted
enrichments span 8x to 30x so the comparison covers both a marginal and an
obvious peak.

## Results

| measure | easyepi | MACS3 |
|---|---|---|
| peaks called | 6 | 7 |
| mean fraction of each planted region recovered | **1.000** | 0.976 |
| peaks lying mostly outside a planted region | **0 of 6** | 1 of 7 |
| interval Jaccard between the two sets | \multicolumn — | **0.831** |
| easyepi peaks with a MACS3 partner | \multicolumn — | **6/6** |
| MACS3 peaks with an easyepi partner | \multicolumn — | 6/7 |
| median \|left edge difference\| on shared peaks | \multicolumn — | **95 bp** |
| median \|right edge difference\| | \multicolumn — | **90 bp** |

Bin size is 200 bp, so both edge differences are sub-bin: the callers agree about
where a region starts and ends to within one bin.

## What this establishes

1. **Recovery is complete.** Every planted region is found by both callers. easyepi
   recovers 100% of each region's span; MACS3 recovers ~97.6% consistently, which
   is its `--extsize 200` fragment model trimming the edges slightly, not a
   disagreement about location.
2. **easyepi did not invent anything.** Zero of its six peaks fell mostly outside a
   planted region. MACS3 reported one extra peak in background; with six
   background-free regions over 400 kb that is a plausible local-lambda call, and
   it is reported here rather than hidden because it is the one place the two
   callers differ in kind.
3. **Edges agree sub-bin.** Median 95/90 bp on a 200 bp bin.

## What this does **not** establish

- **Not equivalence.** `easyepi` uses a global median background; MACS3 uses a
  local lambda. On a genome with copy-number structure or uneven mappability the
  two will diverge, and the global-median model is the weaker one there. This
  fixture is uniform by construction, which is what makes the comparison clean and
  also what makes it narrow.
- **Not a real-cohort result.** The backlog item is titled "on a real cohort".
  This is a synthetic fixture with known ground truth. A real ENCODE-scale BAM
  would exercise the local-background difference, and that has not been run.
- **MACS3 invoked with `--nomodel --extsize 200`.** Those flags make the fragment
  model comparable to easyepi's fixed-bin signal. MACS3's default model would
  build a strand-shifted empirical fragment length, and the numbers above do not
  describe that configuration.
- **MACS3's own output was not validated.** It is treated as the reference; this
  page compares the two, it does not audit MACS3.

## Status

MACS3 row of backlog item 12: **partially done**. Synthetic-fixture agreement is
measured and pinned by tests; the real-cohort comparison remains TODO.
