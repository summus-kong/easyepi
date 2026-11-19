#!/usr/bin/env Rscript
# Quantify ChIPseeker's flank_gene_distances defect.
#
# getAllFlankingGene() writes
#     hitInfo$distance[subjectHits(overlapHit)] <- 0
# which throws away queryHits(overlapHit). `featureHit` is features[sh], i.e.
# the hit features expanded per (peak, feature) pair, so one feature occupies
# several elements. "peak A overlaps element k" then lands on row k, which may
# belong to a different peak.
#
# This script rebuilds the function's internals once and scores each row twice:
# once the way the reference does (subjectHits only) and once matching on the
# (peak, element) pair. No per-peak loop is needed, so it finishes in seconds
# where a lapply over thousands of peaks would take minutes.
#
# Usage:
#   Rscript flank_defect.R <txdb.sqlite> <peaks.bed|narrowPeak> [flankDistance]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("usage: flank_defect.R <txdb.sqlite> <peaks> [distance]")
}
txdb_path <- args[1]
peak_path <- args[2]
distance  <- if (length(args) >= 3) as.integer(args[3]) else 5000L

suppressPackageStartupMessages({
    library(ChIPseeker); library(GenomicRanges); library(txdbmaker)
})

txdb <- loadDb(txdb_path)
feat <- ChIPseeker:::getGene(txdb, by = "transcript")

df <- read.delim(peak_path, header = FALSE, stringsAsFactors = FALSE)
# Column 6 of a narrowPeak file is signalValue, not strand; only accept + - .
# Only + and - are usable; anything else (signalValue, ".", "*") means unstranded.
usable_strand <- ncol(df) >= 6 && all(df$V6 %in% c("+", "-"))
gr <- if (usable_strand) {
    GRanges(df$V1, IRanges(df$V2 + 1L, df$V3), strand = df$V6)
} else {
    GRanges(df$V1, IRanges(df$V2 + 1L, df$V3))
}
N <- length(gr)

## ---- reproduce the function's internals once -------------------------------
peak.gr  <- gr
peak.gr2 <- gr                                  # R keeps the un-widened copy
start(ranges(peak.gr)) <- start(ranges(peak.gr)) - distance
end(ranges(peak.gr))   <- end(ranges(peak.gr))   + distance

hit <- findOverlaps(peak.gr, unstrand(feat))
qh  <- queryHits(hit)
sh  <- subjectHits(hit)
featureHit <- feat[sh]
names(featureHit) <- NULL
n_rows <- length(sh)

overlapHit <- findOverlaps(peak.gr2, unstrand(featureHit))

## the reference's way: only the subject half survives
d_ref <- rep(NA_real_, n_rows)
d_ref[subjectHits(overlapHit)] <- 0

## correct way: match on (peak, feature-element)
key_row <- paste(qh, seq_len(n_rows))
key_hit <- paste(queryHits(overlapHit), subjectHits(overlapHit))
d_ok <- rep(NA_real_, n_rows)
d_ok[key_row %in% key_hit] <- 0

## ---- report ---------------------------------------------------------------
wrong <- which(!is.na(d_ref) & is.na(d_ok))
missed <- which(is.na(d_ref) & !is.na(d_ok))
zeroed <- which(!is.na(d_ref))

cat("peaks                        :", N, "\n")
cat("(peak, feature) rows         :", n_rows, "\n")
cat("rows the reference zeros     :", length(zeroed), "\n")
cat("  of those, genuine          :", length(zeroed) - length(wrong), "\n")
cat("  of those, WRONG            :", length(wrong),
    sprintf("(%.2f%% of the zeros)\n", 100 * length(wrong) / length(zeroed)),
    sep = " ")
cat("wrong entries, all rows      :", length(wrong),
    sprintf("(%.2f%%)\n", 100 * length(wrong) / n_rows), sep = " ")
cat("genuine overlaps missed      :", length(missed), "\n")
bad_peaks <- unique(qh[wrong])
cat("peaks with >=1 wrong entry   :", length(bad_peaks),
    sprintf("(%.1f%% of %d)\n", 100 * length(bad_peaks) / N, N), sep = " ")

## ---- ten concrete instances ----------------------------------------------
psD <- ifelse(as.character(strand(featureHit)) == "+",
              start(peak.gr2[qh]) - start(featureHit),
              end(featureHit) - end(peak.gr2[qh]))
peD <- ifelse(as.character(strand(featureHit)) == "+",
              end(peak.gr2[qh]) - start(featureHit),
              end(featureHit) - start(peak.gr2[qh]))
dd <- ifelse(abs(psD) > abs(peD), peD, psD)

ov_all <- findOverlaps(gr, unstrand(featureHit))
ov_by_elem <- split(queryHits(ov_all), subjectHits(ov_all))

cat("\n%-5s %-24s %-13s %-24s %8s %9s  %s\n",
    "peak", "peak coords", "transcript", "transcript coords",
    "ref", "correct", "really overlaps")
cat(strrep("-", 122), "\n", sep = "")
shown <- 0L
seen <- integer(0)
for (k in wrong) {
    p <- qh[k]
    if (p %in% seen) next
    seen <- c(seen, p)
    tx <- as.character(featureHit$tx_name[k])
    j <- match(tx, as.character(feat$tx_name))
    cat(sprintf("%-5d %-24s %-13s %-24s %8d %9d  %s\n",
                p,
                sprintf("%s:%d-%d", as.character(seqnames(gr))[p],
                        start(gr)[p], end(gr)[p]),
                tx,
                sprintf("%s:%d-%d", as.character(seqnames(feat))[j],
                        start(feat)[j], end(feat)[j]),
                0L, as.integer(dd[k]),
                paste(head(ov_by_elem[[as.character(k)]], 4), collapse = ",")))
    shown <- shown + 1L
    if (shown >= 10L) break
}
