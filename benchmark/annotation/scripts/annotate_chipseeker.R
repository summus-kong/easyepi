#!/usr/bin/env Rscript
# Annotate peaks with ChIPseeker against the same TxDb EasyEpi would use.
#
# Why ChIPseeker is the reference and not a second implementation: it is the tool the
# repository itself claims agreement with, so it is the only reference that can settle
# whether that claim generalises to independent data.
#
# Both region label and distance-to-TSS are written. A high label agreement with poor
# distance agreement would mean the labels match for an uninteresting reason (e.g.
# everything lands in promoters in this particular data), so the driver compares both.

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicRanges)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
bed <- args[1]
out <- args[2]
label <- args[3]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile(bed, as = "GRanges")

anno <- annotatePeak(
  peaks,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  level = "gene",
  verbose = FALSE
)

df <- as.data.frame(anno)
cat(sprintf("%s: annotated %d peaks\n", label, nrow(df)))

# `annotation` is ChIPseeker's region label; `distanceToTSS` is signed bp to nearest TSS.
res <- data.frame(
  seqnames = as.character(df$seqnames),
  start = as.integer(df$start),
  end = as.integer(df$end),
  annotation = as.character(df$annotation),
  distanceToTSS = as.numeric(df$distanceToTSS),
  stringsAsFactors = FALSE
)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)

tab <- sort(table(res$annotation), decreasing = TRUE)
cat("top region labels:\n")
print(head(tab, 8))

summary <- list(
  label = label,
  input = bed,
  peaks = nrow(res),
  txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
  tss_region = c(-3000, 3000),
  label_counts = as.list(as.integer(tab)),
  names = names(tab)
)
write_json(summary, paste0(out, ".json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("written -> %s\n", out))
