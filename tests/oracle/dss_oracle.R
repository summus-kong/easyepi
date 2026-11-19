#!/usr/bin/env Rscript
# Oracle for the pure-Python beta-binomial test: DSS reference values,
# intermediates included, so the port can be checked step by step and not just
# on its final p-value.
#
# The product no longer needs R (easyepi.dmr_beta_binomial replaces the DSS
# subprocess), but DSS remains the reference this implementation is checked
# against. Everything below calls the same DSS functions DMLtest calls, in the
# same order, and then also runs DMLtest itself and asserts the two agree before
# writing anything -- so a green comparison downstream is a comparison with DSS,
# not with this script's idea of DSS.
#
# Usage: Rscript dss_oracle.R <counts.txt> <out.tsv> [smoothing] [equal_disp]
#   counts.txt: n_sites, n_samples, chr (space separated), pos, meth row-major,
#               total row-major, design
#   smoothing:  TRUE (default) or FALSE -- DMLtest(smoothing = ...)
#   equal_disp: FALSE (default) or TRUE  -- DMLtest(equal.disp = ...)
#
# The output has the same columns in every mode, so the harness compares all of
# them regardless of which switch produced the reference.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: dss_oracle.R <counts.txt> <out.tsv> [smoothing] [equal_disp]")
}
data_path <- args[1]
out_path <- args[2]
SMOOTHING <- if (length(args) >= 3) as.logical(args[3]) else TRUE
EQUAL_DISP <- if (length(args) >= 4) as.logical(args[4]) else FALSE
# `DIRECT`: call DSS:::DMLtest.noSmooth instead of the public DMLtest. Needed
# for one replicate + smoothing=FALSE, which DMLtest refuses up front even
# though DMLtest.noSmooth handles it by sharing one dispersion.
DIRECT <- if (length(args) >= 5) as.logical(args[5]) else FALSE
if (is.na(SMOOTHING) || is.na(EQUAL_DISP) || is.na(DIRECT))
  stop("smoothing/equal_disp/direct unparsable")

con <- file(data_path, "r")
on.exit(close(con), add = TRUE)
n_sites <- as.integer(readLines(con, n = 1))
n_samples <- as.integer(readLines(con, n = 1))
if (length(n_sites) != 1 || length(n_samples) != 1) stop("bad header")
chr <- scan(con, what = character(), n = n_sites, quiet = TRUE)
pos <- scan(con, what = integer(), n = n_sites, quiet = TRUE)
meth <- matrix(scan(con, what = integer(), n = n_sites * n_samples, quiet = TRUE),
               nrow = n_samples, byrow = TRUE)
total <- matrix(scan(con, what = integer(), n = n_sites * n_samples, quiet = TRUE),
                nrow = n_samples, byrow = TRUE)
design <- scan(con, what = integer(), n = n_samples, quiet = TRUE)

suppressPackageStartupMessages({
  library(DSS)
  library(bsseq)
  library(GenomicRanges)
  library(IRanges)
})

SPAN <- 500

loci <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
sample_names <- paste0("s", seq_len(n_samples))
bs <- BSseq(M = t(meth), Cov = t(total), gr = loci,
            sampleNames = sample_names, rmZeroCov = FALSE)

group1 <- sample_names[design == 0]
group2 <- sample_names[design == 1]

# ---- DSS's own site filter, applied exactly where DMLtest applies it: on the
# full locus list, before anything else is computed.
cov1 <- as.array(getCoverage(bs[, group1], type = "Cov"))
cov2 <- as.array(getCoverage(bs[, group2], type = "Cov"))
keep <- DSS:::hasCoverage(cov1, pos) & DSS:::hasCoverage(cov2, pos)

bs1 <- bs[, group1][keep]
bs2 <- bs[, group2][keep]
x1 <- as.array(getCoverage(bs1, type = "M"))
n1 <- as.array(getCoverage(bs1, type = "Cov"))
x2 <- as.array(getCoverage(bs2, type = "M"))
n2 <- as.array(getCoverage(bs2, type = "Cov"))
allchr <- as.character(seqnames(bs1))
allpos <- start(bs1)
nreps1 <- ncol(x1)
nreps2 <- ncol(x2)

mean_at <- function(x, n) {
  if (SMOOTHING) DSS:::compute.mean.Smooth(x, n, allchr, allpos, SPAN)
  else           DSS:::compute.mean.noSmooth(x, n)
}

# DMLtest.noSmooth shares one dispersion when either group has a single
# replicate; DMLtest.Smooth does not (it errors instead).
SHARE <- EQUAL_DISP || (!SMOOTHING && (nreps1 == 1 || nreps2 == 1))

raw1 <- mean_at(x1, n1)
raw2 <- mean_at(x2, n2)

if (SHARE) {
  prior <- DSS:::est.prior.BSseq.logN(cbind(x1, x2), cbind(n1, n2))
  phi1 <- phi2 <- DSS:::est.dispersion.BSseq(
    cbind(x1, x2), cbind(n1, n2), cbind(raw1, raw2), 1)
} else {
  prior1 <- DSS:::est.prior.BSseq.logN(x1, n1)
  prior2 <- DSS:::est.prior.BSseq.logN(x2, n2)
  phi1 <- DSS:::est.dispersion.BSseq(x1, n1, raw1, 1)
  phi2 <- DSS:::est.dispersion.BSseq(x2, n2, raw2, 1)
}

group_var <- function(x, n, phi) {
  wt <- 1 / (1 + (n - 1) * phi)
  wt <- wt / mean(wt, na.rm = TRUE)
  lvl <- mean_at(x * wt, n * wt)
  e <- lvl[, 1]
  if (SMOOTHING) {
    v <- DSS:::compute.var.Smooth(e, n * wt, phi, SPAN, allchr, allpos)
  } else {
    # compute.waldStat.noSmooth, verbatim:
    #   rowSums(n1 * p * (1-p) * (1 + (n1-1)*phi1)) / rowSums(n1)^2
    v <- rowSums((n * wt) * e * (1 - e) * (1 + ((n * wt) - 1) * phi)) /
         (rowSums(n * wt))^2
  }
  list(wt = wt, level = e, var = v)
}

g1 <- group_var(x1, n1, phi1)
g2 <- group_var(x2, n2, phi2)

diff <- g1$level - g2$level
vv <- g1$var + g2$var
if (SMOOTHING) {
  vv[vv < 1e-05] <- 1e-04
  vv[vv > 1 - 1e-05] <- 1 - 1e-04
} else {
  # compute.waldStat.noSmooth floors to 1e-05, one-sided, and does not cap.
  vv[vv < 1e-05] <- 1e-05
}
se <- sqrt(vv)
stat <- diff / se
pval <- 2 * pnorm(-abs(stat))

# ---- DMLtest itself: prove this chain is the same computation.
if (DIRECT) {
  fit <- DSS:::DMLtest.noSmooth(bs1, bs2, EQUAL_DISP, 1)
} else {
  fit <- DMLtest(bs, group1 = group1, group2 = group2,
                 smoothing = SMOOTHING, smoothing.span = SPAN,
                 equal.disp = EQUAL_DISP, ncores = 1)
}
stopifnot(nrow(fit) == length(diff))
stopifnot(all(fit$pos == allpos))
tol <- 1e-12
stopifnot(max(abs(fit$diff - diff)) <= tol)
stopifnot(max(abs(fit$diff.se - se)) <= tol)
stopifnot(max(abs(fit$pval - pval)) <= tol)
stopifnot(max(abs(fit$phi1 - phi1)) <= tol)
stopifnot(max(abs(fit$phi2 - phi2)) <= tol)

out <- data.frame(
  rn = which(keep), chr = allchr, pos = allpos,
  e1 = raw1[, 1], phi1 = phi1, e1w = g1$level, var1 = g1$var,
  e2 = raw2[, 1], phi2 = phi2, e2w = g2$level, var2 = g2$var,
  diff = diff, se = se, stat = stat, pval = pval
)
prior_line <- if (SHARE) c(prior, prior) else c(prior1, prior2)
writeLines(c(
  paste("#prior", sprintf("%.17g", prior_line), sep = "\t"),
  paste("#n_sites", n_sites, "#n_tested", nrow(out), sep = "\t"),
  paste("#smoothing", SMOOTHING, "#equal_disp", EQUAL_DISP,
        "#direct", DIRECT, sep = "\t")
))
write.table(out, out_path, sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE, na = "NA")
cat("DSS", as.character(packageVersion("DSS")),
    "smoothing =", SMOOTHING, " equal.disp =", EQUAL_DISP, "->", out_path, "\n")
