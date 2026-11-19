#!/usr/bin/env Rscript
# Run ChIPseeker on the same peaks with the same promoter window easyepi uses, so the
# comparison is like-for-like.
#
# Defaults deliberately left at ChIPseeker's:
#   tssRegion = c(-3000, 3000)   -> easyepi's --promoter-downstream/upstream 3000
#   genomicAnnotationPriority    -> Promoter, 5UTR, 3UTR, Exon, Intron,
#                                   Downstream, Intergenic
# easyepi has no Downstream bin (it folds those peaks into Intergenic), so the
# comparison maps Downstream -> Intergenic rather than letting it show up as a
# disagreement.
suppressMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
peaks_bed <- args[[1]]
txdb_path <- args[[2]]
out_tsv <- args[[3]]

cat("peaks:", peaks_bed, "\n")
cat("txdb :", txdb_path, "\n")

txdb <- loadDb(txdb_path)

peaks <- read.table(peaks_bed, header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "")
# easyepi accepts narrowPeak/broadPeak/BED; take the first three columns as intervals.
gr <- GRanges(seqnames = peaks$V1,
              ranges = IRanges(start = peaks$V2 + 1L, end = peaks$V3))

cat("peaks read:", length(gr), "\n")

anno <- ChIPseeker::annotatePeak(gr, tssRegion = c(-3000, 3000), TxDb = txdb,
                                 verbose = FALSE)
df <- as.data.frame(anno)

# The annotation column is "Promoter (<=1kb)", "5' UTR", "Exon", ... Keep the label
# only, because the distance suffix varies and is not what is being compared.
df$bin <- sub(" \\(.*\\)$", "", df$annotation)
df$bin <- trimws(df$bin)

# Fold categories easyepi does not model, so a missing bin is not counted as a
# disagreement. easyepi's five non-intergenic bins are Promoter, 5'UTR, 3'UTR, Exon,
# Intron; everything else it reports is Intergenic.
df$bin_mapped <- ifelse(df$bin %in% c("Promoter", "5' UTR", "3' UTR", "Exon",
                                      "Intron"),
                        df$bin, "Intergenic")
df$bin_mapped[df$bin_mapped == "5' UTR"] <- "5'UTR"
df$bin_mapped[df$bin_mapped == "3' UTR"] <- "3'UTR"

write.table(df[, c("seqnames", "start", "end", "bin", "bin_mapped", "distanceToTSS")],
            out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("written:", out_tsv, "\n")
cat("raw bins:\n"); print(table(df$bin))
cat("mapped bins:\n"); print(table(df$bin_mapped))
