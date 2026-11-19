#!/usr/bin/env Rscript
# Dump ChIPseeker's own nearest-feature index and distance, per peak.
#
# Why: the bin assignment disagreed even after porting the rules, and the only
# remaining unknown is `distance` itself. Calling
# getNearestFeatureIndicesAndDistances directly isolates that function from the
# priority loop, so a per-peak diff says exactly where the two implementations
# diverge instead of only reporting an aggregate percentage.

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicFeatures)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
bed <- args[1]
out <- args[2]
level <- args[3]          # "gene" or "transcript"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile(bed, as = "GRanges")

features <- if (level == "transcript") {
  ChIPseeker:::getGene(txdb, by = "transcript")
} else {
  ChIPseeker:::getGene(txdb, by = "gene")
}

# ChIPseeker's default arguments: sameStrand = FALSE, overlap = "TSS".
res <- ChIPseeker:::getNearestFeatureIndicesAndDistances(
  peaks, features,
  sameStrand = FALSE, ignoreOverlap = FALSE,
  ignoreUpstream = FALSE, ignoreDownstream = FALSE,
  overlap = "TSS"
)

# `res` omits peaks where every bound was NA, so re-key by peak coordinates.
key <- paste(as.character(seqnames(res$peak)), start(res$peak), end(res$peak), sep = ":")
df <- data.frame(
  key = key,
  index = res$index,
  distance = res$distance,
  tx_start = start(features[res$index]),
  tx_end = end(features[res$index]),
  tx_strand = as.character(strand(features[res$index])),
  stringsAsFactors = FALSE
)
write.table(df, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d peaks with a resolved nearest feature (of %d)\n", level, nrow(df), length(peaks)))
