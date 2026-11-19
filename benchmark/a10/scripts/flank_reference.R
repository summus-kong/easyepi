#!/usr/bin/env Rscript
# Build the flank_* reference for every dataset, in two flavours:
#
#   dist_raw    what ChIPseeker actually writes
#               (hitInfo$distance[subjectHits(overlapHit)] <- 0)
#   dist_fixed  the same rows scored on the (peak, feature) pair, i.e. what the
#               line above would produce if it kept queryHits(overlapHit)
#
# Output is long format (dataset, peakIdx, key, dist_raw, dist_fixed) so a
# caller can compare per feature instead of per semicolon-joined string --
# ChIPseeker's and EasyEpi's list orders differ for reasons documented in
# benchmark/a10/README.md.
#
# `level` selects the feature set the same way annotatePeak(level=) does:
#   transcript -> getGene(TxDb, by = "transcript"); key is tx_name
#   gene       -> getGene(TxDb, by = "gene");       key is the gene name
# At gene level there is no flank_txIds column, so `key` is the gene id.
#
# Usage: flank_reference.R <txdb.sqlite> <peaksDir> <outDir> [distance] [level]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: flank_reference.R <txdb> <peaksDir> <outDir> [distance] [level]")
txdb_path <- args[1]; peaks_dir <- args[2]; out_dir <- args[3]
distance <- if (length(args) >= 4) as.integer(args[4]) else 5000L
level    <- if (length(args) >= 5) args[5] else "transcript"
stopifnot(level %in% c("transcript", "gene"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
    library(ChIPseeker); library(GenomicRanges); library(txdbmaker)
})
txdb <- loadDb(txdb_path)
feat <- ChIPseeker:::getGene(txdb, by = level)
cat("level:", level, " features:", length(feat), "\n")

files <- sort(list.files(peaks_dir, pattern = "\\.bed$", full.names = TRUE))
summary_rows <- list()

for (fp in files) {
    name <- sub("\\.bed$", "", basename(fp))
    df <- read.delim(fp, header = FALSE, stringsAsFactors = FALSE)
    usable_strand <- ncol(df) >= 6 && all(df$V6 %in% c("+", "-"))
    gr <- if (usable_strand) {
        GRanges(df$V1, IRanges(df$V2 + 1L, df$V3), strand = df$V6)
    } else {
        GRanges(df$V1, IRanges(df$V2 + 1L, df$V3))
    }
    N <- length(gr)

    ## ---- the function's internals, once ----
    peak.gr <- gr; peak.gr2 <- gr
    start(ranges(peak.gr)) <- start(ranges(peak.gr)) - distance
    end(ranges(peak.gr))   <- end(ranges(peak.gr))   + distance
    hit <- findOverlaps(peak.gr, unstrand(feat))
    qh <- queryHits(hit); sh <- subjectHits(hit)
    featureHit <- feat[sh]; names(featureHit) <- NULL
    n_rows <- length(sh)

    overlapHit <- findOverlaps(peak.gr2, unstrand(featureHit))

    ## reference's scoring
    d_raw <- rep(NA_real_, n_rows)
    d_raw[subjectHits(overlapHit)] <- 0
    ## corrected scoring
    key_row <- paste(qh, seq_len(n_rows))
    key_hit <- paste(queryHits(overlapHit), subjectHits(overlapHit))
    d_fix <- rep(NA_real_, n_rows)
    d_fix[key_row %in% key_hit] <- 0

    ## non-overlapping rows get the signed min(|psD|,|peD|), feature strand
    psD <- ifelse(as.character(strand(featureHit)) == "+",
                  start(peak.gr2[qh]) - start(featureHit),
                  end(featureHit) - end(peak.gr2[qh]))
    peD <- ifelse(as.character(strand(featureHit)) == "+",
                  end(peak.gr2[qh]) - start(featureHit),
                  end(featureHit) - start(peak.gr2[qh]))
    dd <- ifelse(abs(psD) > abs(peD), peD, psD)

    d_raw[is.na(d_raw)] <- dd[is.na(d_raw)]
    d_fix[is.na(d_fix)] <- dd[is.na(d_fix)]

    if (level == "transcript") {
        key <- as.character(featureHit$tx_name)
    } else {
        cn <- names(as.data.frame(featureHit))
        key <- as.character(as.data.frame(featureHit)$gene_id)
    }

    out <- data.frame(dataset = name,
                      peakIdx = qh,
                      key     = key,
                      dist_raw   = d_raw,
                      dist_fixed = d_fix,
                      stringsAsFactors = FALSE)
    write.table(unique(out), file.path(out_dir, paste0(name, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    n_zero_raw <- sum(d_raw == 0, na.rm = TRUE)
    wrong <- sum(d_raw == 0 & d_fix != 0, na.rm = TRUE)
    badp  <- length(unique(qh[d_raw == 0 & d_fix != 0]))
    summary_rows[[name]] <- data.frame(
        dataset = name, peaks = N, rows = n_rows,
        zero_raw = n_zero_raw, zero_fixed = sum(d_fix == 0, na.rm = TRUE),
        wrong = wrong, peaks_affected = badp,
        stringsAsFactors = FALSE)
    cat(sprintf("%-24s peaks=%5d rows=%6d zero_ref=%6d wrong=%5d (%.2f%%) peaks_affected=%4d (%.2f%%)\n",
                name, N, n_rows, n_zero_raw, wrong, 100 * wrong / n_rows,
                badp, 100 * badp / N))
}

summ <- do.call(rbind, summary_rows)
write.table(summ, file.path(out_dir, "_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== totals ===\n")
cat(sprintf("datasets %d  peaks %d  rows %d\n",
            nrow(summ), sum(summ$peaks), sum(summ$rows)))
cat(sprintf("reference zeros %d, of which wrong %d (%.2f%% of zeros, %.2f%% of rows)\n",
            sum(summ$zero_raw), sum(summ$wrong),
            100 * sum(summ$wrong) / sum(summ$zero_raw),
            100 * sum(summ$wrong) / sum(summ$rows)))
cat(sprintf("peaks affected %d (%.2f%%)\n",
            sum(summ$peaks_affected), 100 * sum(summ$peaks_affected) / sum(summ$peaks)))
