#!/usr/bin/env Rscript
# Export TxDb features as plain tables, so the Python port consumes exactly the
# objects ChIPseeker sees.
#
# Why not a GTF round-trip: a GTF export and re-parse introduced TSS positions
# that do not exist in the TxDb (e.g. chr1:199107), which pushed peaks into the
# promoter window and made every downstream comparison wrong. Writing the
# features directly removes format conventions, strand handling and ID remapping
# from the equation -- the only thing left to differ is the algorithm.
#
# Output (tab separated, no header quoting, 1-based inclusive like GRanges):
#   tx.tsv       tx_id chrom start end strand
#   tss.tsv      tx_id chrom pos strand
#   feature.tsv  tx_id kind chrom start end      (kind: exon|intron|utr5|utr3)

suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicFeatures)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
level <- if (length(args) >= 2) args[2] else "transcript"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# `getGene` is what annotatePeak calls; replicate it exactly.
features <- if (level == "transcript") {
  ChIPseeker:::getGene(txdb, by = "transcript")
} else {
  ChIPseeker:::getGene(txdb, by = "gene")
}
cat("level:", level, " features:", length(features), "\n")

# Identity of each feature. The child GRangesLists are named by `tx_name`
# (ENST...), so tx.tsv uses the same namespace -- re-deriving a composite key
# made the join fail silently and every child feature dropped out.
tx_names <- as.character(mcols(features)$tx_name)
if (all(is.na(tx_names))) {
  # Gene-level features have no tx_name; fall back to a coordinate composite.
  tx_names <- sprintf("%s:%d-%d:%s", as.character(seqnames(features)),
                      start(features), end(features),
                      as.character(strand(features)))
}
tx_key <- tx_names

tx_df <- data.frame(
  key = tx_key,
  chrom = as.character(seqnames(features)),
  start = start(features),
  end = end(features),
  strand = as.character(strand(features)),
  stringsAsFactors = FALSE
)
write.table(tx_df, file.path(outdir, "tx.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

r <- resize(features, width = 1)
tss_df <- data.frame(
  key = tx_key,
  chrom = as.character(seqnames(r)),
  pos = start(r),
  strand = as.character(strand(r)),
  stringsAsFactors = FALSE
)
write.table(tss_df, file.path(outdir, "tss.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("tss rows:", nrow(tss_df), "\n")

# Child features. These come from their own TxDb accessors, not from
# `getGene`: ChIPseeker's get_exonList() calls exonsBy(TxDb) -- which defaults
# to by = "gene" -- while get_intronList() calls intronsByTranscript(TxDb).
# Their names are therefore gene ids / transcript ids, unrelated to the feature
# keys above, so they are written with a synthetic group name and matched on
# coordinates only.
emit <- function(grl, kind) {
  if (is.null(grl) || length(grl) == 0) return(NULL)
  flat <- unlist(grl, use.names = FALSE)
  data.frame(
    key = rep(names(grl), times = lengths(grl)),
    kind = kind,
    chrom = as.character(seqnames(flat)),
    start = start(flat),
    end = end(flat),
    stringsAsFactors = FALSE
  )
}

parts <- list(
  emit(intronsByTranscript(txdb, use.names = TRUE), "intron"),
  emit(exonsBy(txdb, use.names = TRUE), "exon"),
  emit(fiveUTRsByTranscript(txdb, use.names = TRUE), "utr5"),
  emit(threeUTRsByTranscript(txdb, use.names = TRUE), "utr3")
)
feat <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
write.table(feat, file.path(outdir, "feature.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("feature rows:", nrow(feat), " by kind:",
    paste(names(table(feat$kind)), table(feat$kind), sep = "=", collapse = " "), "\n")
