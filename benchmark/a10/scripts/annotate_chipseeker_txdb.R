#!/usr/bin/env Rscript
# Annotate peaks with ChIPseeker against a TxDb built FROM the same GTF EasyEpi reads.
#
# Why a second script rather than a flag on the first one: the original
# annotate_chipseeker.R pins TxDb.Hsapiens.UCSC.hg38.knownGene, and the whole point
# of this variant is that the two tools must resolve features from ONE gene model.
# Keeping them as separate files makes the gene model an explicit, auditable input
# instead of a hidden default.
#
# What was wrong with the first comparison (measured, not assumed):
#   * TxDb.Hsapiens.UCSC.hg38.knownGene has 412,044 transcripts, but ChIPseeker's
#     transcriptsBy()-based feature loader drops the 85,737 that carry no gene_id,
#     leaving 326,307 features / 263,047 distinct TSS.
#   * The GTF EasyEpi reads yields 341,049 distinct TSS, of which 78,002 are absent
#     from ChIPseeker's restricted set and none are absent from the GTF-derived set.
#   * So ChIPseeker was answering with a strictly smaller TSS catalogue than EasyEpi.
#     Every "EasyEpi says Promoter, ChIPseeker says Intron/Exon/Intergenic" dispute
#     is explained by that gap, not by an EasyEpi rule error.
#
# Usage: annotate_chipseeker_txdb.R <bed> <out.tsv> <txdb.sqlite> [level] [label]

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: annotate_chipseeker_txdb.R <bed> <out.tsv> <txdb.sqlite> [level] [label]")
}
bed <- args[1]
out <- args[2]
txdb_path <- args[3]
level <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "transcript"
label <- if (length(args) >= 5 && nzchar(args[5])) args[5] else level

txdb <- loadDb(txdb_path)
peaks <- readPeakFile(bed, as = "GRanges")

# Warn rather than silently annotate against a truncated catalogue: if the DB
# loader ever drops transcripts again, the disagreement would look like a rule bug.
n_tx <- length(transcripts(txdb))
n_feat <- length(ChIPseeker:::getGene(txdb, by = level))
if (n_feat != n_tx) {
  cat(sprintf("WARNING: %s features %d < %d transcripts; catalogue is truncated\n",
              level, n_feat, n_tx))
}

anno <- annotatePeak(
  peaks,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  level = level,
  verbose = FALSE
)

df <- as.data.frame(anno)
cat(sprintf("%s: annotated %d peaks\n", label, nrow(df)))

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
  txdb = txdb_path,
  transcripts_in_txdb = n_tx,
  features_used = n_feat,
  peaks = nrow(res),
  tss_region = c(-3000, 3000),
  level = level,
  label_counts = as.list(as.integer(tab)),
  names = names(tab)
)
write_json(summary, paste0(out, ".json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("written -> %s\n", out))
