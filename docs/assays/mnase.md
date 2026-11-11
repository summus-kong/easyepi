# MNase-seq

MNase digests linker DNA and leaves nucleosome-protected fragments. Two different
questions need two different commands.

## Occupancy and periodicity (signal track)

```bash
easyepi bam2bw mnase.bam -o tracks/ --assay mnase --jobs 8
easyepi nucleosome --bam mnase.bam --out nuc/ --figure-format pdf
```

Reports the fragment-length distribution and the periodic spacing of the dyad signal.

## Dyad positions (paired-end fragments)

```bash
easyepi nucleosome-positioning --bam mnase_paired.bam \
    --chrom chr1 --start 0 --end 2000000 \
    --bin-size 10 --min-length 120 --max-length 200 \
    --out positioning/ --figure-format pdf
```

Calls come from the **fragment midpoint** (the dyad), with a robust
median + MAD threshold and a per-call fuzziness measure.

## Why single-end is refused

A single-end read reports one cut site and no fragment span, so the dyad cannot be
located. The command raises rather than returning plausible-looking coordinates, and
points to `easyepi nucleosome` for the occupancy-based alternative.

## Assumptions

- the length window 120–200 bp selects mononucleosomes; widen it for soft digestion
- `--min-spacing` (default 100 bp) prevents one broad nucleosome from being called
  several times
