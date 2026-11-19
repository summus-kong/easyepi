#!/usr/bin/env Rscript
# Export a TxDb to GTF with one consistent transcript-id namespace.
#
# Why not `rtracklayer::export(txdb, "gtf")`
# ------------------------------------------
# That call writes the rows stored in the TxDb's own tables and does not expand
# the per-transcript GRangesLists, so the result has exon/CDS rows but zero UTR
# rows. ChIPseeker queries fiveUTRsByTranscript/threeUTRsByTranscript, which
# return real intervals, so a peak inside a UTR could only fall through to
# Promoter or Intron.
#
# The namespace trap
# ------------------
# `transcripts()` returns a plain GRanges with **no names**, while every
# `*ByTranscript()` accessor names its GRangesList elements by the TxDb's
# internal numeric tx_id ("1", "2", ...). Writing both through unchanged
# produces a file where transcript rows say ENST... and exon rows say "1":
# nothing links up, and every downstream join silently fails or doubles.
# This script therefore keys everything by tx_name, mapping the internal id once.
#
# Coordinates: TxDb/GRanges are 1-based closed and GTF is 1-based closed, so the
# values pass through unchanged. No +-1 anywhere.

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: txdb_to_gtf.R <txdb.sqlite|package> <out.gtf>")
}
txdb_path <- args[1]
out_path <- args[2]

if (grepl("\\.sqlite$|\\.db$", txdb_path)) {
  txdb <- loadDb(txdb_path)
} else {
  suppressPackageStartupMessages(library(txdb_path, character.only = TRUE))
  txdb <- get(txdb_path)
}

message("reading transcripts ...")
tx <- transcripts(txdb, columns = c("tx_name", "gene_id"))
tx_name <- as.character(GenomicRanges::mcols(tx)$tx_name)
gene_id <- as.character(GenomicRanges::mcols(tx)$gene_id)

# `transcripts()` does not expose `tx_id`, yet every `*ByTranscript()` accessor
# names its elements by that internal id ("1", "2", ...). The two are aligned by
# construction — row i of `transcripts()` *is* tx_id i — so the mapping is taken
# positionally rather than looked up in a column that does not exist.
tx_id <- as.character(seq_along(tx))

blank <- is.na(tx_name) | tx_name == ""
tx_name[blank] <- tx_id[blank]

# internal tx_id -> tx_name, and tx_name -> gene_id. These are the only two
# lookups the writer needs.
name_of_id <- stats::setNames(tx_name, tx_id)
gene_of_name <- stats::setNames(
  ifelse(is.na(gene_id) | gene_id == "", tx_name, gene_id),
  tx_name
)

emit <- function(grl, feature, out_fh) {
  if (length(grl) == 0) return(0L)
  widths <- lengths(grl)
  keep <- widths > 0
  if (!any(keep)) return(0L)
  grl <- grl[keep]
  widths <- widths[keep]

  raw <- rep(names(grl), widths)
  ids <- unname(name_of_id[raw])
  ids[is.na(ids)] <- raw[is.na(ids)]

  flat <- unlist(grl, use.names = FALSE)
  gids <- unname(gene_of_name[ids])
  gids[is.na(gids)] <- ids[is.na(gids)]

  df <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(flat)),
    source = "TxDb",
    feature = feature,
    start = as.integer(GenomicRanges::start(flat)),
    end = as.integer(GenomicRanges::end(flat)),
    score = ".",
    strand = as.character(GenomicRanges::strand(flat)),
    frame = ".",
    attributes = sprintf('gene_id "%s"; transcript_id "%s";', gids, ids),
    stringsAsFactors = FALSE
  )
  write.table(df, file = out_fh, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE, append = TRUE)
  nrow(df)
}

cat('##description: TxDb export with exon/CDS/UTR rows\n', file = out_path)
cat('##format: gtf\n', file = out_path, append = TRUE)

total <- 0L

# The transcript rows themselves: one per transcript, named by tx_name so the
# child rows can reference them.
tx_named <- tx
names(tx_named) <- tx_name
total <- total + emit(split(tx_named, names(tx_named)), "transcript", out_path)

for (spec in list(
  list(fun = exonsBy, arg = "tx", feat = "exon"),
  list(fun = cdsBy, arg = "tx", feat = "CDS"),
  list(fun = fiveUTRsByTranscript, arg = NULL, feat = "five_prime_utr"),
  list(fun = threeUTRsByTranscript, arg = NULL, feat = "three_prime_utr")
)) {
  grl <- tryCatch(
    if (is.null(spec$arg)) spec$fun(txdb) else spec$fun(txdb, by = spec$arg),
    error = function(e) {
      message("  skipped ", spec$feat, ": ", conditionMessage(e))
      NULL
    })
  if (is.null(grl)) next
  n <- emit(grl, spec$feat, out_path)
  message(sprintf("  %-15s %d rows", spec$feat, n))
  total <- total + n
}

message(sprintf("written -> %s (%d rows, %.1f MB)",
                out_path, total, file.size(out_path) / 1024^2))
