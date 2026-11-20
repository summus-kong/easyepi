# MNase-seq

MNase digests linker DNA and leaves nucleosome-protected fragments. Two different
questions need two different commands, and one of them cannot answer as much as it
used to claim.

## Occupancy and nucleosome-free fraction (signal track)

```bash
easyepi bam2bw mnase.bam -o tracks/ --assay mnase --jobs 8
easyepi nucleosome --bigwig tracks/mnase.bw \
    --chrom chr1 --start 0 --end 5000000 --bin-size 10 --out nuc/
```

Writes `nucleosome_summary.json` with the z-scored occupancy mean, the
nucleosome-free fraction, and a nucleosome-repeat estimate that carries its own
verdict (see below).

## Dyad positions (paired-end fragments)

```bash
easyepi nucleosome-positioning --bam mnase_paired.bam \
    --chrom chr1 --start 0 --end 2000000 \
    --bin-size 10 --min-length 120 --max-length 200 \
    --out positioning/ --figure-format pdf
```

Calls come from the **fragment midpoint** (the dyad), with a robust median + MAD
threshold and a per-call fuzziness measure. The BED interval is the measured
dyad cluster, so widths differ between a well-positioned call and a fuzzy one.
Column 6 is `.`: a dyad has no strand, because the midpoint of a fragment does
not carry a direction.

## What the repeat estimate does and does not mean

`nucleosome` and `nucleosome-positioning` both report a nucleosome repeat, and
both attach `reliable` to it. Read that flag before quoting the number.

Measured on `GSE247785` (SRR26830867, GRCh38, 48.5 M read pairs, 98.3 % aligned,
median fragment 161 bp):

| quantity | value |
|---|---|
| calls in chr1:0–5 Mb | 15,231 |
| median spacing between neighbouring dyads | 210 bp |
| median fuzziness | 370 bp |
| repeat from a whole-region PSD | **1280 bp — `reliable: false`** |

The 1280 bp figure is not a nucleosome repeat. Two mechanisms produce it:

1. `welch()` runs with `nperseg=256`, so the PSD has only 128 frequency points
   and cannot resolve any period longer than 256 bins (2560 bp at 10 bp bins).
   The peak therefore lands on the lowest bins and is **identical for every
   window size from 20 kb to 5 Mb**.
2. raw dyad density carries a 1/f trend, and `argmax(psd)` follows that trend
   instead of a repeat peak.

A detrended autocorrelation of the same 5 Mb track found no ~185 bp peak either
(peak 120–130 bp, ACF amplitude 0.035–0.14). That is expected rather than
surprising: whole-region dyad spacing is heterogeneous because phasing restarts
across genes and chromatin domains. **A nucleosome repeat has to be measured on
aligned regions or a metagene profile.** When the estimate falls outside
120–400 bp the module says so and refuses to print a number as if it were a
repeat.

## Nucleosome-free fraction

The NFR fraction is the lowest `nfr_quantile` (default 0.2) of the raw signal,
not a fixed z-score cut. MNase coverage is strongly right-skewed — on the same
dataset 38.6 % of 10 bp bins in chr1:0–5 Mb are exactly zero while the tail
reaches 17× — so after z-scoring the minimum is only about −0.9. A `z < −1` rule
therefore reported **0 % NFR on every real input**. Ties are included: if the
20th percentile is 0, every zero bin counts, and the fraction reports the tied
population rather than exactly the requested quantile.

## Why single-end is refused

A single-end read reports one cut site and no fragment span, so the dyad cannot be
located. The command raises rather than returning plausible-looking coordinates, and
points to `easyepi nucleosome` for the occupancy-based alternative.

This is also why ENCODE's 2011 MNase-seq (K562/GM12878) is not usable as a dyad
input: those runs are 35 bp single-end, so the pipeline has nothing to take a
midpoint of.

## Assumptions

- the length window 120–200 bp selects mononucleosomes; widen it for soft digestion
- `--min-spacing` (default 100 bp) prevents one broad nucleosome from being called
  several times
- `--fuzziness-radius` (default 500 bp) caps the dyad-cluster search. On the
  dataset above about 28 % of calls still hit that cap, so raise it when
  neighbouring nucleosomes sit closer than 500 bp
