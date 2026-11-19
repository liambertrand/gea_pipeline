#!/usr/bin/env Rscript
# Redundancy analysis following Capblancq & Forester (2021):
#  1. variance partitioning  (environment | geography | structure)
#  2. main (partial) RDA  G ~ env + Condition(structure [+ geography])
#  3. rdadapt outlier test on the first K constrained axes
#
# Rscript rda.R --config config.yaml --geno imputed.rds --env env_selected.tsv \
#   --geo geo.tsv --pcs pcs.tsv --k K.txt --outdir results/sp/rda [--threads 4]
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
suppressPackageStartupMessages(library(vegan))

a <- parse_args(); need(a, c("config", "geno", "env", "geo", "pcs", "k", "outdir"))
cfg <- read_config(a$config)$rda
threads <- as.integer(if (is.null(a$threads)) 1 else a$threads)
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

G <- readRDS(a$geno)
env <- read_samples_tab(a$env); geo <- read_samples_tab(a$geo); pcs <- read_samples_tab(a$pcs)
kinfo <- fread(a$k, header = FALSE); n_pcs <- kinfo[V1 == "n_pcs", as.integer(V2)]
stopifnot(identical(env$sample_id, rownames(G)), identical(geo$sample_id, rownames(G)),
          identical(pcs$sample_id, rownames(G)))
X <- as.matrix(env[, -1]); GEO <- as.matrix(geo[, -1])
STR <- as.matrix(pcs[, 2:(1 + n_pcs)])
msg("G: ", nrow(G), " x ", ncol(G), "; env: ", paste(colnames(X), collapse = ", "),
    "; structure PCs: ", n_pcs)

# ---------------- 1. variance partitioning
vp_rows <- list()
if (isTRUE(cfg$variance_partitioning)) {
  fit_part <- function(name, Y, Z) {
    m <- if (is.null(Z)) rda(G, Y) else rda(G, Y, Z)
    r2 <- RsquareAdj(m)
    p <- NA_real_
    if (as.integer(cfg$permutations) > 0) p <- anova.cca(m, permutations = as.integer(cfg$permutations), parallel = threads)$`Pr(>F)`[1]
    data.table(model = name, inertia = m$CCA$tot.chi, prop_total = m$CCA$tot.chi / m$tot.chi,
               R2 = r2$r.squared, R2_adj = r2$adj.r.squared, p = p)
  }
  vp_rows <- rbindlist(list(
    fit_part("full: env + geo + structure", cbind(X, GEO, STR), NULL),
    fit_part("pure env | geo + structure", X, cbind(GEO, STR)),
    fit_part("pure structure | env + geo", STR, cbind(X, GEO)),
    fit_part("pure geo | env + structure", GEO, cbind(X, STR))
  ))
  full <- vp_rows$prop_total[1]
  vp_rows[, prop_of_explained := prop_total / full]
  vp_rows[, confounded := c(full - sum(prop_total[2:4]), NA, NA, NA)]
  msg("variance partitioning:\n", paste(capture.output(print(vp_rows)), collapse = "\n"))
} else {
  vp_rows <- data.table(model = "skipped")
}
fwrite(vp_rows, file.path(a$outdir, "variance_partitioning.tsv"), sep = "\t")

# ---------------- 2. main model
cond <- NULL
co <- unlist(cfg$condition_on)
if ("structure" %in% co) cond <- cbind(cond, STR)
if ("geo" %in% co) cond <- cbind(cond, GEO)
mod <- if (is.null(cond)) rda(G, X) else rda(G, X, cond)
ev <- mod$CCA$eig
msg("constrained axes: ", length(ev), "; proportion of constrained inertia: ",
    paste(round(ev / sum(ev), 3), collapse = " "))

K <- as_auto(cfg$K, NA)
if (is.na(K)) K <- sum(ev > mean(ev))
K <- as.integer(max(2, min(K, length(ev))))
if (length(ev) < 2) K <- 1L
msg("rdadapt K = ", K)

perm_txt <- "model permutation test skipped"
if (as.integer(cfg$permutations) > 0) {
  an <- anova.cca(mod, permutations = as.integer(cfg$permutations), parallel = threads)
  perm_txt <- capture.output(print(an))
}

# ---------------- 3. rdadapt
res <- rdadapt(mod, K)
ids <- colnames(G)
loads <- mod$CCA$v[, seq_len(K), drop = FALSE]
colnames(loads) <- paste0("RDA", seq_len(K), "_loading")
# which predictor each SNP tracks most closely (raw correlation)
R <- cor(G, X)
top <- max.col(abs(R), ties.method = "first")
out <- cbind(split_id(ids), res, as.data.table(loads),
             data.table(rda_top_env = colnames(X)[top], rda_top_env_r = R[cbind(seq_len(nrow(R)), top)]))
fwrite(out, file.path(a$outdir, "rda_snps.tsv.gz"), sep = "\t")
msg("RDA outliers at q < 0.05: ", sum(out$q_rda < 0.05, na.rm = TRUE))

sink(file.path(a$outdir, "rda_summary.txt"))
cat("Predictors:", paste(colnames(X), collapse = ", "), "\n")
cat("Conditioned on:", if (length(co)) paste(co, collapse = " + ") else "nothing", "\n")
cat("Structure PCs:", n_pcs, "\n")
cat("Constrained proportion of total inertia:", round(mod$CCA$tot.chi / mod$tot.chi, 4), "\n")
cat("Adjusted R2:", round(RsquareAdj(mod)$adj.r.squared, 4), "\n")
cat("Axis proportions (constrained):", paste(round(ev / sum(ev), 3), collapse = " "), "\n")
cat("rdadapt K:", K, "\n")
cat("Outliers q<0.05:", sum(out$q_rda < 0.05, na.rm = TRUE), "\n\n")
cat(perm_txt, sep = "\n")
sink()

# ---------------- plots
pdf(file.path(a$outdir, "rda_plots.pdf"), width = 9, height = 8)
barplot(ev / sum(ev), main = "Constrained axes", ylab = "Proportion", col = ifelse(seq_along(ev) <= K, "grey30", "grey80"))
if (length(ev) >= 2) {
  sc_site <- scores(mod, display = "sites", choices = 1:2, scaling = 1)
  sc_bp <- scores(mod, display = "bp", choices = 1:2, scaling = 1)
  sc_snp <- scores(mod, display = "species", choices = 1:2, scaling = 1)
  out_idx <- which(out$q_rda < 0.05)
  plot(sc_snp, pch = 20, col = "grey75", cex = 0.4, main = "RDA triplot (outlier SNPs in red)",
       xlim = range(c(sc_snp[, 1], sc_site[, 1])), ylim = range(c(sc_snp[, 2], sc_site[, 2])))
  points(sc_snp[out_idx, , drop = FALSE], pch = 20, col = "firebrick", cex = 0.6)
  points(sc_site, pch = 21, bg = "steelblue", cex = 0.8)
  mul <- ordiArrowMul(sc_bp)
  arrows(0, 0, sc_bp[, 1] * mul, sc_bp[, 2] * mul, length = 0.08, col = "black")
  text(sc_bp * mul * 1.1, rownames(sc_bp), cex = 0.8)
}
hist(out$p_rda, breaks = 50, main = "rdadapt p-values (should be ~uniform + spike near 0)", xlab = "p")
dev.off()
msg("done")
