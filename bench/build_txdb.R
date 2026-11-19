#!/usr/bin/env Rscript
# Build a TxDb from the same GENCODE GTF that easyepi annotates against.
#
# Both sides of the benchmark must read one gene model. ChIPseeker's usual input is
# TxDb.Hsapiens.UCSC.hg38.knownGene (RefSeq/UCSC); GENCODE v50 has different TSS
# positions and transcript sets, so comparing easyepi-against-GENCODE with
# ChIPseeker-against-RefSeq would measure the annotation sources, not the
# classification algorithm. Building the TxDb from the same file removes that.
#
# GenomicFeatures >= 1.61.1 defunct'd makeTxDbFromGFF; it lives in txdbmaker now.
suppressMessages({
  library(txdbmaker)
  library(GenomicFeatures)
})

args <- commandArgs(trailingOnly = TRUE)
gtf <- args[[1]]
out <- args[[2]]

cat("GTF:", gtf, "\n")
cat("txdbmaker:", as.character(packageVersion("txdbmaker")), "\n")

txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf",
                                   organism = "Homo sapiens", taxonomyId = 9606)
saveDb(txdb, out)

cat("saved:", out, "\n")
cat("transcripts:", length(transcripts(txdb)), "\n")
cat("exons      :", length(exons(txdb)), "\n")
cat("CDS        :", length(cds(txdb)), "\n")
cat("genes      :", length(genes(txdb)), "\n")
