#!/usr/bin/env Rscript
# Population structure on LD-pruned SNPs (PCA + sNMF), choice of K,
# and ancestry-aware imputation of the full genotype matrix.
#
# Rscript structure.R --config config.yaml --pruned pruned.raw --full all.raw \
#   --outdir results/sp/structure --geno-out results/sp/geno/imputed.rds [--threads 4]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
suppressPackageStartupMessages(library(LEA))

a <- parse_args(); need(a, c("config", "pruned", "full", "outdir", "geno_out"))
cfg <- read_config(a$config)$structure
threads <- as.integer(if (is.null(a$threads)) 1 else a$threads)
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

# ---------------- PCA on pruned SNPs
pr <- read_raw(a$pruned)
msg("pruned matrix: ", nrow(pr), " individuals x ", ncol(pr), " SNPs")
prm <- apply(pr, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
prm <- prm[, apply(prm, 2, var) > 0, drop = FALSE]
pca <- prcomp(prm, center = TRUE, scale. = TRUE)
n_keep <- min(20, ncol(pca$x))
pve <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(n_keep)]

# ---------------- sNMF for K = 1..K_max
geno_file <- file.path(a$outdir, "pruned.geno")
g9 <- pr; g9[is.na(g9)] <- 9L
LEA::write.geno(g9, geno_file)
K_max <- min(as.integer(cfg$K_max), nrow(pr) - 1)
proj <- LEA::snmf(geno_file, K = 1:K_max, entropy = TRUE, repetitions = as.integer(cfg$snmf_reps),
                  alpha = as.numeric(cfg$snmf_alpha), CPU = threads, project = "new", seed = 1)
ce <- sapply(1:K_max, function(k) mean(LEA::cross.entropy(proj, K = k)))
K <- as_auto(cfg$K, which.min(ce))
K <- max(1L, as.integer(K))
msg("cross-entropy by K: ", paste(round(ce, 4), collapse = " "), " -> K = ", K)
best <- which.min(LEA::cross.entropy(proj, K = K))
Q <- LEA::Q(proj, K = K, run = best)
colnames(Q) <- paste0("Q", seq_len(K))

n_pcs <- as.integer(as_auto(cfg$n_pcs, K))
n_pcs <- max(1L, min(n_pcs, n_keep))
writeLines(c(paste0("K\t", K), paste0("n_pcs\t", n_pcs)), file.path(a$outdir, "K.txt"))

pcs <- data.table(sample_id = rownames(pr), pca$x[, seq_len(n_keep), drop = FALSE])
fwrite(pcs, file.path(a$outdir, "pcs.tsv"), sep = "\t")
fwrite(data.table(sample_id = rownames(pr), Q), file.path(a$outdir, "ancestry_Q.tsv"), sep = "\t")
fwrite(data.table(K = 1:K_max, cross_entropy = ce), file.path(a$outdir, "cross_entropy.tsv"), sep = "\t")

# ---------------- plots
pdf(file.path(a$outdir, "structure.pdf"), width = 9, height = 6)
par(mfrow = c(2, 2))
plot(pve, type = "b", xlab = "PC", ylab = "Proportion of variance", main = "PCA scree (pruned SNPs)")
cl <- max.col(Q)
plot(pca$x[, 1], pca$x[, 2], col = cl, pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA coloured by sNMF cluster")
plot(1:K_max, ce, type = "b", xlab = "K", ylab = "Cross-entropy", main = "sNMF"); abline(v = K, lty = 2)
o <- order(cl, -apply(Q, 1, max))
barplot(t(Q[o, , drop = FALSE]), col = seq_len(K), border = NA, space = 0, main = paste("Ancestry, K =", K))
dev.off()

# ---------------- ancestry-aware imputation of full matrix
G <- read_raw(a$full)
stopifnot(identical(rownames(G), rownames(pr)))
msg("full matrix: ", nrow(G), " individuals x ", ncol(G), " SNPs; missing = ",
    round(mean(is.na(G)) * 100, 2), "%")

mode_fill <- function(M) {
  # most common genotype per column, NA-aware, vectorised
  c0 <- colSums(M == 0L, na.rm = TRUE); c1 <- colSums(M == 1L, na.rm = TRUE); c2 <- colSums(M == 2L, na.rm = TRUE)
  mode <- max.col(cbind(c0, c1, c2), ties.method = "first") - 1L
  idx <- which(is.na(M), arr.ind = TRUE)
  if (nrow(idx)) M[idx] <- mode[idx[, 2]]
  M
}
if (identical(cfg$impute, "median")) {
  med <- apply(G, 2, median, na.rm = TRUE)
  idx <- which(is.na(G), arr.ind = TRUE); G[idx] <- as.integer(round(med[idx[, 2]]))
} else {
  for (k in unique(cl)) {
    rows <- which(cl == k)
    G[rows, ] <- mode_fill(G[rows, , drop = FALSE])
  }
  G <- mode_fill(G)  # anything left (column missing in a whole cluster)
}
# drop SNPs that became monomorphic
v <- apply(G, 2, var)
G <- G[, v > 0, drop = FALSE]
msg("after imputation: ", ncol(G), " polymorphic SNPs")
saveRDS(G, a$geno_out, compress = FALSE)
msg("done")
