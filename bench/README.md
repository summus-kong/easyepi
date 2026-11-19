# ChIPseeker concordance benchmark

Compares `easyepi annotate` against `ChIPseeker::annotatePeak` on the same peaks,
the same order and the same promoter window.

## Why a shared TxDb matters

ChIPseeker's usual input is `TxDb.Hsapiens.UCSC.hg38.knownGene` (RefSeq/UCSC).
GENCODE v50 has different TSS positions and a different transcript set, so
easyepi-against-GENCODE compared with ChIPseeker-against-RefSeq would measure the
**annotation sources**, not the classification algorithm. `build_txdb.R` builds the
TxDb from the same GTF easyepi reads, removing that confound.

## Requirements

R (>= 4.4) with `ChIPseeker`, `txdbmaker`, `GenomicFeatures`, `GenomeInfoDb`,
`GenomeInfoDbData`, installed via micromamba from `conda-forge`/`bioconda`:

```bash
micromamba install -n <env> -c conda-forge -c bioconda \
    bioconductor-chipseeker bioconductor-txdbmaker bioconductor-genomeinfodbdata
```

`makeTxDbFromGFF` moved out of `GenomicFeatures` (defunct since 1.61.1) into
`txdbmaker`; calling it on the old package raises a `.Defunct` error.

## Running

```bash
Rscript build_txdb.R   <gencode.gtf.gz> <gencode.sqlite>     # ~6 min
Rscript run_chipseeker.R <peaks> <gencode.sqlite> cs.tsv     # ~3 min
easyepi annotate --peaks <peaks> --gtf <gencode.gtf.gz> --out ee/
python compare.py --chipseeker cs.tsv --easyepi ee/annotated.tsv
```

## Conventions handled in the comparison

* ChIPseeker returns 1-based closed intervals; easyepi writes 0-based half-open BED.
  The comparison normalises before checking that rows correspond, so a systematic
  offset cannot pass as agreement.
* ChIPseeker's `Downstream` bin has no easyepi equivalent. Both are folded into
  `Intergenic`; otherwise `Downstream` would show up as a fake disagreement.
* Only `Promoter`, `5'UTR`, `3'UTR`, `Exon` and `Intron` are compared directly;
  everything else is `Intergenic` on both sides.
