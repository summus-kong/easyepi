#!/usr/bin/env Rscript
# Export the same feature tables as export_features.R, but starting from a GTF.
#
# Why: `export_features.R` reads a TxDb, while `gtf_features.py` reads a GTF.
# To claim the two agree, both must see the same gene model, so this script
# builds a TxDb *from the GTF* with makeTxDbFromGFF and then applies the same
# accessors ChIPseeker uses. Any remaining difference is then the parser, not
# the annotation source.
#
# Output matches export_features.R exactly, so the two directories can be
# compared row by row.

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  if (!requireNamespace("txdbmaker", quietly = TRUE)) {
    stop("txdbmaker is required: makeTxDbFromGFF moved there in GenomicFeatures >= 1.61")
  }
})

args <- commandArgs(trailingOnly = TRUE)
gtf_path <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("building TxDb from", gtf_path, "...\n")
# GenomicFeatures >= 1.61 defuncts its own makeTxDbFromGFF and points at
# txdbmaker; the accessors below are unchanged.
txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf",
                                   circ_seqs = character())

# ChIPseeker's getGene(by = "transcript"): one row per transcript, with the
# transcription start site strand-resolved.
tx <- transcripts(txdb, columns = c("tx_name", "gene_id"))
tx_name <- as.character(mcols(tx)$tx_name)
gene_id <- as.character(mcols(tx)$gene_id)
strand_chr <- as.character(strand(tx))

tx_df <- data.frame(
  key = tx_name, chrom = as.character(seqnames(tx)),
  start = start(tx), end = end(tx), strand = strand_chr,
  stringsAsFactors = FALSE
)
write.table(tx_df, file.path(outdir, "tx.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# resize(width = 1) picks the left edge on '+' and the right edge on '-'.
tss_pos <- ifelse(strand_chr == "-", end(tx), start(tx))
tss_df <- data.frame(
  key = tx_name, chrom = as.character(seqnames(tx)),
  pos = as.integer(tss_pos), strand = strand_chr,
  stringsAsFactors = FALSE
)
write.table(tss_df, file.path(outdir, "tss.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Child features, keyed by the TxDb group name. Two accessors differ from the
# transcript table on purpose: exonsBy() defaults to by = "gene" and
# intronsByTranscript() is per transcript, exactly as ChIPseeker calls them.
names(tx) <- tx_name
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
parts <- parts[!vapply(parts, is.null, logical(1))]
feat <- do.call(rbind, parts)
write.table(feat, file.path(outdir, "feature.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat(sprintf("features: %d  feature rows: %d\n", nrow(tx_df), nrow(feat)))
cat("by kind:", paste(names(table(feat$kind)), table(feat$kind),
                      sep = "=", collapse = " "), "\n")
