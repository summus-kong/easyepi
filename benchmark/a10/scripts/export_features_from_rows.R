#!/usr/bin/env Rscript
# Export feature tables straight from GTF rows, without makeTxDbFromGFF.
#
# Why not a TxDb: makeTxDbFromGFF synthesises a full-length exon for every
# transcript row, so on a 2,548-transcript slice it yields 17,575 exons where
# the file only has 15,027 exon rows. Comparing that against a parser that
# reads explicit exon rows can never reach parity, and the mismatch says
# nothing about the parser. Building the tables from the same rows removes the
# difference in *input semantics* and leaves only the rules under test:
#
#   tss     <- transcript row, resized to 1 bp, strand aware
#   exon    <- explicit exon rows, grouped by gene, NOT merged
#   intron  <- gaps between a transcript's own exon spans
#   3' UTR  <- explicit UTR rows, stop-codon bases trimmed
#   5' UTR  <- explicit UTR rows
#
# Usage: export_features_from_gtf.R <gtf> <outdir> <transcript|gene>

suppressPackageStartupMessages({
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: export_features_from_gtf.R <gtf> <outdir> [transcript|gene]")
}
gtf_path <- args[1]
outdir <- args[2]
level <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "transcript"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("reading ", gtf_path, " ...")
df <- read.table(gtf_path, header = FALSE, sep = "\t", quote = "",
                 comment.char = "#", stringsAsFactors = FALSE,
                 colClasses = c("character", "character", "character",
                                "integer", "integer", "character", "character",
                                "character", "character"))
names(df) <- c("chrom", "src", "feature", "start", "end", "score", "strand",
               "frame", "attr")
df <- as.data.frame(df, stringsAsFactors = FALSE)
message("rows read: ", nrow(df), "  features: ",
        paste(unique(df$feature), collapse = ","))

# transcript_id is the key ChIPseeker joins on.
attr_of <- function(attr, field) {
  pat <- paste0(field, ' "([^"]+)"')
  out <- rep(NA_character_, length(attr))
  hit <- regmatches(attr, regexec(pat, attr))
  for (i in seq_along(hit)) {
    if (length(hit[[i]]) >= 2) out[i] <- hit[[i]][2]
  }
  out
}

df$tx_id <- attr_of(df$attr, "transcript_id")
df$gene_id <- attr_of(df$attr, "gene_id")

tx <- as.data.frame(df[df$feature == "transcript" & !is.na(df$tx_id), ],
                    stringsAsFactors = FALSE)
if (nrow(tx) == 0) {
  # Some GTF variants omit explicit transcript rows; derive them from exons.
  span <- do.call(rbind, lapply(split(df[!is.na(df$tx_id), ], df$tx_id[!is.na(df$tx_id)]),
                                function(d) data.frame(tx_id = d$tx_id[1],
                                                       chrom = d$chrom[1],
                                                       start = min(d$start),
                                                       end = max(d$end),
                                                       strand = d$strand[1],
                                                       gene_id = d$gene_id[1])))
  tx <- as.data.frame(span, stringsAsFactors = FALSE)
}
tx <- as.data.frame(tx[!duplicated(tx$tx_id), , drop = FALSE],
                    stringsAsFactors = FALSE)
tx$idx <- seq_len(nrow(tx)) - 1L
tx_idx <- seq_len(nrow(tx)) - 1L
names(tx_idx) <- as.character(tx$tx_id)

# resize(width = 1) is strand aware: left edge on +, right edge on -.
tss_pos <- ifelse(tx$strand == "-", tx$end, tx$start)

tss_tbl <- data.frame(key = paste0("T", tx$idx),
                      chrom = tx$chrom, pos = tss_pos,
                      strand = tx$strand, stringsAsFactors = FALSE)
write.table(tss_tbl, file.path(outdir, "tx.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
write.table(tss_tbl, file.path(outdir, "tss.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# Map a transcript id to its 0-based feature owner index; NA when unknown.
owner_of <- function(ids) {
  v <- tx_idx[as.character(ids)]
  as.integer(ifelse(is.na(v), NA_integer_, v))
}

rows <- list()

# ---- exon: explicit rows, grouped by gene, intervals kept separate ----------
ex <- df[df$feature == "exon" & !is.na(df$tx_id), ]
ex$gene <- ifelse(is.na(ex$gene_id), ex$tx_id, ex$gene_id)
ex <- ex[order(ex$chrom, ex$gene, ex$start, ex$end), ]
ex_owner <- ave(seq_len(nrow(ex)), paste(ex$chrom, ex$gene), FUN = seq_along) - 1L
rows[["exon"]] <- data.frame(key = paste0("T", ex_owner), kind = "exon",
                             chrom = ex$chrom, start = ex$start, end = ex$end)

# ---- intron: gaps between a transcript's own exon spans --------------------
intron_rows <- list()
for (tid in unique(ex$tx_id)) {
  d <- ex[ex$tx_id == tid, ]
  d <- d[order(d$start, d$end), ]
  d <- d[!duplicated(paste(d$start, d$end)), ]
  if (nrow(d) < 2) next
  for (i in seq_len(nrow(d) - 1)) {
    s <- d$end[i] + 1
    e <- d$start[i + 1] - 1
    if (e >= s) {
      intron_rows[[length(intron_rows) + 1]] <-
        data.frame(key = paste0("T", owner_of(tid)), kind = "intron",
                   chrom = d$chrom[i], start = s, end = e)
    }
  }
}
rows[["intron"]] <- do.call(rbind, intron_rows)

# ---- UTRs: explicit rows, 3' side trimmed of stop-codon bases --------------
stops <- df[df$feature == "stop_codon" & !is.na(df$tx_id), ]

for (kind in c("five_prime_utr", "three_prime_utr")) {
  u <- df[df$feature == kind & !is.na(df$tx_id), ]
  if (nrow(u) == 0) next
  out <- list()
  for (i in seq_len(nrow(u))) {
    s <- u$start[i]; e <- u$end[i]; tid <- u$tx_id[i]
    if (kind == "three_prime_utr") {
      st <- stops[stops$tx_id == tid, ]
      for (j in seq_len(nrow(st))) {
        ss <- st$start[j]; se <- st$end[j]
        if (se < s || ss > e) next
        if (u$strand[i] != "-") {
          if (se < e) { s <- se + 1 } else { s <- e + 1 }
        } else {
          if (ss > s) { e <- ss - 1 } else { e <- s - 1 }
        }
      }
    }
    if (e >= s) {
      out[[length(out) + 1]] <-
        data.frame(key = paste0("T", owner_of(tid)),
                   kind = if (kind == "five_prime_utr") "utr5" else "utr3",
                   chrom = u$chrom[i], start = s, end = e)
    }
  }
  nm <- if (kind == "five_prime_utr") "utr5" else "utr3"
  rows[[nm]] <- do.call(rbind, out)
}

feat <- do.call(rbind, rows)
feat <- feat[order(feat$kind, feat$chrom, feat$start, feat$end), ]
write.table(feat, file.path(outdir, "feature.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

cat(sprintf("transcripts: %d  feature rows: %d\n", nrow(tx), nrow(feat)))
cat("by kind:", paste(names(table(feat$kind)), table(feat$kind),
                      sep = "=", collapse = " "), "\n")
