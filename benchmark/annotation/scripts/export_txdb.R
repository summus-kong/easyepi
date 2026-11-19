#!/usr/bin/env Rscript
# Export the TxDb used by ChIPseeker as a GTF, so EasyEpi annotates against the SAME
# transcript set.
#
# Why this matters for the benchmark
# ----------------------------------
# The claim under test is "EasyEpi reproduces ChIPseeker's region labels". If EasyEpi
# reads GENCODE v50 while ChIPseeker reads TxDb.Hsapiens.UCSC.hg38.knownGene, any
# disagreement could come from the annotation rather than the implementation, and the
# comparison would be uninterpretable.
#
# Exporting the TxDb removes that confound: both tools then see identical gene models,
# so a disagreement is an implementation difference and nothing else.
#
# Only `gene` and `transcript` features are written. ChIPseeker's region classes
# (promoter / 5'UTR / exon / intron / downstream / distal intergenic) are derived from
# the transcript model plus TSS, so exon and UTR features are included as well to keep
# the two tools able to classify the same way.

suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicFeatures)
})

args <- commandArgs(trailingOnly = TRUE)
out <- args[1]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# `asGFF` writes a GFF3 that GTF readers generally accept; writeGFF3 is explicit.
gff <- asGFF(txdb)
rtracklayer::export(gff, out, format = "gff3")
cat(sprintf("exported %d features -> %s\n", length(gff), out))

# Report the transcript count so the benchmark README can state the annotation basis.
cat(sprintf("transcripts in TxDb: %d\n", length(transcripts(txdb))))
cat(sprintf("genes in TxDb      : %d\n", length(genes(txdb))))
