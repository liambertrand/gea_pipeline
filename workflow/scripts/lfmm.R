#!/usr/bin/env Rscript
# LFMM2 (ridge, via LEA::lfmm2) genome-wide, all selected predictors jointly, K latent factors.
# Run on ALL SNPs (not only RDA outliers) so the genomic-inflation calibration
# and FDR are valid; the two methods are intersected afterwards.
#
# Rscript lfmm.R --config config.yaml --geno imputed.rds --env env_selected.tsv \
#   --k K.txt --outdir results/sp/lfmm
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
suppressPackageStartupMessages(library(LEA))

a <- parse_args(); need(a, c("config", "geno", "env", "k", "outdir"))
cfg <- read_config(a$config)$lfmm
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

G <- readRDS(a$geno)
env <- read_samples_tab(a$env)
stopifnot(identical(env$sample_id, rownames(G)))
X <- as.matrix(env[, -1])
kinfo <- fread(a$k, header = FALSE)
K <- as.integer(as_auto(cfg$K, kinfo[V1 == "K", as.integer(V2)]))
msg("LFMM2: ", ncol(G), " SNPs, predictors ", paste(colnames(X), collapse = ", "), ", K = ", K)

Y <- G; storage.mode(Y) <- "double"
mod <- LEA::lfmm2(input = Y, env = X, K = K, lambda = as.numeric(cfg$lambda))
tst <- LEA::lfmm2.test(object = mod, input = Y, env = X, linear = TRUE,
                       genomic.control = identical(cfg$calibrate, "gif"))
# LEA returns one row per predictor (or a vector for one predictor): make SNPs x predictors
as_snp_by_env <- function(m) {
  m <- as.matrix(m)
  if (nrow(m) != ncol(Y)) m <- t(m)
  stopifnot(nrow(m) == ncol(Y), ncol(m) == ncol(X))
  m
}
P <- as_snp_by_env(tst$pvalues)
Z <- as_snp_by_env(tst$zscores)
colnames(P) <- colnames(Z) <- colnames(X)
gif <- data.table(variable = colnames(X), gif = rep_len(as.numeric(tst$gif), ncol(X)))
msg("GIF: ", paste(sprintf("%s=%.2f", gif$variable, gif$gif), collapse = " "))
fwrite(gif, file.path(a$outdir, "gif.tsv"), sep = "\t")

Q <- apply(P, 2, safe_q)
out <- split_id(colnames(G))
for (v in colnames(X)) {
  out[, paste0("p_", v) := P[, v]]
  out[, paste0("q_", v) := Q[, v]]
  out[, paste0("z_", v) := Z[, v]]
}
Qn <- Q; Qn[is.na(Qn)] <- 1
best <- max.col(-Qn, ties.method = "first")        # predictor with the smallest q
out[, lfmm_min_p := apply(P, 1, min, na.rm = TRUE)]
out[, lfmm_min_q := Qn[cbind(seq_len(nrow(Qn)), best)]]
out[, lfmm_best_env := colnames(X)[best]]
out[, lfmm_sig_env := apply(Q < 0.05, 1, function(r) paste(colnames(X)[which(r)], collapse = ";"))]
fwrite(out, file.path(a$outdir, "lfmm_snps.tsv.gz"), sep = "\t")
msg("SNPs with q < 0.05 for >= 1 predictor: ", sum(out$lfmm_min_q < 0.05, na.rm = TRUE))

pdf(file.path(a$outdir, "lfmm_plots.pdf"), width = 10, height = 3 + 1.5 * ncol(X))
par(mfrow = c(ncol(X), 2), mar = c(3, 4, 2, 1))
for (v in colnames(X)) {
  hist(P[, v], breaks = 50, main = sprintf("%s calibrated p (GIF %.2f)", v, gif[variable == v, gif]), xlab = "")
  plot(-log10(P[, v]), pch = 20, cex = 0.3, col = ifelse(Q[, v] < 0.05, "firebrick", "grey50"),
       ylab = "-log10 p", xlab = "", main = v)
}
dev.off()
msg("done")
