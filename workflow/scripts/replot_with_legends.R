#!/usr/bin/env Rscript
# =====================================================================
# replot_with_legends.R
# Redraws the GEA pipeline's figures with full legends, labelled axes and
# readable variable names. Large plots (millions of SNPs) are written as PNG
# (a few MB) instead of PDF (hundreds of MB); small plots as PDF + PNG.
#
# Uses only files the pipeline already wrote; re-runs no analysis.
#
# Usage (from the gea-pipeline folder, with the pipeline's own R):
#   ENV=$(ls -d .snakemake/conda/*_ | head -1)
#   $ENV/bin/Rscript workflow/scripts/replot_with_legends.R results_viridis/Cviridis
#
# Output: <results dir>/figures/  (existing pipeline figures are not touched)
# =====================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
D <- if (length(args)) args[1] else "results_viridis/Cviridis"
SP <- basename(normalizePath(D))
OUT <- file.path(D, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
f <- function(...) file.path(D, ...)
have <- function(...) file.exists(f(...))

# ---------------------------------------------------------------- style
CAT  <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
INK  <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; GREY1 <- "#b5b3ad"; GREY2 <- "#d6d4ce"
TIER <- c(tier1 = "#eb6834", rda_only = "#2a78d6", lfmm_only = "#1baf7a")
TIER_LAB <- c(tier1 = "Both RDA and LFMM2 (tier 1)", rda_only = "RDA only", lfmm_only = "LFMM2 only")
VAR <- c(bio1 = "bio1 Annual mean temperature", bio2 = "bio2 Mean diurnal range",
         bio3 = "bio3 Isothermality", bio4 = "bio4 Temperature seasonality",
         bio5 = "bio5 Max temp. warmest month", bio6 = "bio6 Min temp. coldest month",
         bio7 = "bio7 Temperature annual range", bio8 = "bio8 Mean temp. wettest quarter",
         bio9 = "bio9 Mean temp. driest quarter", bio10 = "bio10 Mean temp. warmest quarter",
         bio11 = "bio11 Mean temp. coldest quarter", bio12 = "bio12 Annual precipitation",
         bio13 = "bio13 Precip. wettest month", bio14 = "bio14 Precip. driest month",
         bio15 = "bio15 Precipitation seasonality", bio16 = "bio16 Precip. wettest quarter",
         bio17 = "bio17 Precip. driest quarter", bio18 = "bio18 Precip. warmest quarter",
         bio19 = "bio19 Precip. coldest quarter", lc_developed = "Developed land (NLCD)",
         lc_grass = "Grassland (NLCD)", lc_shrub = "Shrubland (NLCD)", lc_forest = "Forest (NLCD)",
         lc_agri = "Agriculture (NLCD)", lc_wetland = "Wetland (NLCD)", lc_barren = "Barren (NLCD)",
         lc_water = "Water (NLCD)")
vname <- function(v) ifelse(v %in% names(VAR), VAR[v], v)

theme <- function() par(family = "sans", col.axis = INK2, col.lab = INK, fg = INK2, col.main = INK,
                        font.main = 2, las = 1, tcl = -0.25, mgp = c(2.4, 0.5, 0), bty = "l")
open_png <- function(name, w, h, res = 150) {
  path <- file.path(OUT, paste0(name, ".png"))
  ok <- tryCatch({ png(path, width = w, height = h, units = "in", res = res, type = "cairo"); TRUE },
                 error = function(e) FALSE)
  if (!ok) png(path, width = w, height = h, units = "in", res = res)
  theme()
}
# small figures: draw once to PDF and once to PNG
both <- function(name, w, h, draw) {
  pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h); theme(); draw(); dev.off()
  open_png(name, w, h); draw(); dev.off()
  msg("wrote figures/", name, ".pdf/.png")
}
legend_right <- function(...) legend("topright", inset = c(-0.02, 0), xpd = NA, bty = "n", text.col = INK, ...)

# chromosome layout: order contigs by length, cumulative offsets
chrom_layout <- function(dt) {
  L <- dt[, .(len = max(pos)), by = chrom][order(-len)]
  L[, offset := c(0, cumsum(as.numeric(len))[-.N])]
  L[, mid := offset + len / 2]
  L[, idx := seq_len(.N)]
  L
}
chrom_label <- function(x) { x <- sub("^scaffold[-_]?", "", x, ignore.case = TRUE); ifelse(nchar(x) > 10, substr(x, 1, 10), x) }

# =====================================================================
# 1. Genome-wide Manhattan plots (RDA and LFMM2), coloured by tier
# =====================================================================
if (have("rda/rda_snps.tsv.gz") && have("lfmm/lfmm_snps.tsv.gz")) {
  msg("reading per-SNP results (this takes a minute)")
  R <- fread(f("rda/rda_snps.tsv.gz"), select = c("snp_id", "chrom", "pos", "p_rda", "q_rda"))
  Lf <- fread(f("lfmm/lfmm_snps.tsv.gz"))
  S <- merge(R, Lf[, .(snp_id, lfmm_min_p, lfmm_min_q, lfmm_best_env)], by = "snp_id")
  S[, chrom := as.character(chrom)]
  qthr <- 0.05
  S[, tier := fifelse(q_rda < qthr & lfmm_min_q < qthr, "tier1",
               fifelse(q_rda < qthr, "rda_only", fifelse(lfmm_min_q < qthr, "lfmm_only", "none")))]
  S[is.na(tier), tier := "none"]
  CL <- chrom_layout(S)
  S <- merge(S, CL[, .(chrom, offset, idx)], by = "chrom")
  S[, x := pos + offset]
  ntier <- S[, .N, by = tier]; nt <- setNames(ntier$N, ntier$tier); nt[setdiff(names(TIER), names(nt))] <- 0
  genes <- if (have("genes.tsv")) fread(f("genes.tsv")) else NULL
  if (!is.null(genes)) genes[, chrom := as.character(chrom)]

  # background SNPs: thin the non-significant ones with large p for speed and file size
  set.seed(1)
  bg <- S[tier == "none" & (pmin(p_rda, lfmm_min_p) < 0.05 | runif(.N) < 0.15)]
  fg <- S[tier != "none"]
  open_png("manhattan_with_legend", 13, 9)
  layout(matrix(1:2, 2)); par(mar = c(3.5, 4.5, 3.5, 15), xpd = FALSE)
  for (what in c("rda", "lfmm")) {
    yb <- if (what == "rda") -log10(bg$p_rda) else -log10(bg$lfmm_min_p)
    yf <- if (what == "rda") -log10(fg$p_rda) else -log10(fg$lfmm_min_p)
    qcol <- if (what == "rda") S$q_rda else S$lfmm_min_q
    pcol <- if (what == "rda") S$p_rda else S$lfmm_min_p
    thr <- suppressWarnings(max(pcol[qcol < qthr], na.rm = TRUE))
    ymax <- max(c(yb, yf), na.rm = TRUE) * 1.05
    plot(NA, xlim = range(S$x), ylim = c(0, ymax), xaxt = "n", xlab = "", ylab = expression(-log[10](italic(p))),
         main = if (what == "rda") "A. RDA (rdadapt): SNPs associated with the multivariate environment"
                else "B. LFMM2: strongest association across the 6 predictors")
    points(bg$x, yb, pch = 16, cex = 0.25, col = ifelse(bg$idx %% 2 == 0, GREY1, GREY2))
    o <- order(match(fg$tier, c("lfmm_only", "rda_only", "tier1")))
    points(fg$x[o], yf[o], pch = 16, cex = 0.45, col = TIER[fg$tier[o]])
    if (is.finite(thr)) abline(h = -log10(thr), lty = 2, col = INK2)
    axis(1, at = CL$mid, labels = chrom_label(CL$chrom), las = 2, cex.axis = 0.6, tick = FALSE)
    if (!is.null(genes) && what == "lfmm") {
      top <- head(genes[called == TRUE][order(p_enrich)], 10)
      if (nrow(top)) {
        top <- merge(top, CL[, .(chrom, offset)], by = "chrom")
        lab <- ifelse(!is.na(top$name) & nzchar(top$name), top$name, top$gene_id)
        text(top$offset + (top$start + top$end) / 2, ymax * 0.97, lab, srt = 90, adj = 1, cex = 0.55, col = INK)
      }
    }
    par(xpd = NA)
    legend(par("usr")[2] + diff(par("usr")[1:2]) * 0.01, par("usr")[4],
           legend = c(sprintf("%s (n = %s)", TIER_LAB, format(nt[names(TIER)], big.mark = ",")),
                      "Not significant (alternating chromosomes)", "FDR q = 0.05 threshold"),
           col = c(TIER, GREY1, INK2), pch = c(16, 16, 16, 16, NA), lty = c(NA, NA, NA, NA, 2),
           pt.cex = 1.2, bty = "n", cex = 0.8, title = "SNP category", title.adj = 0)
    par(xpd = FALSE)
  }
  mtext(sprintf("%s: %s SNPs tested; top 10 called genes labelled in B", SP, format(nrow(S), big.mark = ",")),
        side = 1, line = 2.5, cex = 0.7, col = INK2)
  dev.off(); msg("wrote figures/manhattan_with_legend.png")

  # ---------------- per-predictor LFMM2 Manhattans + p-value histograms
  pv <- grep("^p_", names(Lf), value = TRUE); vars <- sub("^p_", "", pv)
  gif <- if (have("lfmm/gif.tsv")) fread(f("lfmm/gif.tsv")) else data.table(variable = vars, gif = NA)
  Lf[, chrom := as.character(chrom)]
  Lf <- merge(Lf, CL[, .(chrom, offset, idx)], by = "chrom"); Lf[, x := pos + offset]
  open_png("lfmm_per_predictor_with_legend", 13, 2.3 * length(vars) + 1)
  layout(matrix(seq_len(2 * length(vars)), ncol = 2, byrow = TRUE), widths = c(4, 1.2))
  for (i in seq_along(vars)) {
    v <- vars[i]; p <- Lf[[pv[i]]]; q <- Lf[[paste0("q_", v)]]
    p[is.na(p)] <- 1; q[is.na(q)] <- 1
    keep <- q < qthr | p < 0.05 | runif(length(p)) < 0.15
    par(mar = c(if (i == length(vars)) 4 else 1.5, 4.5, 2, 1))
    sig <- keep & q < qthr
    plot(Lf$x[keep & !sig], -log10(p[keep & !sig]), pch = 16, cex = 0.25,
         col = ifelse(Lf$idx[keep & !sig] %% 2 == 0, GREY1, GREY2), xaxt = "n", xlab = "",
         ylab = expression(-log[10](italic(p))), ylim = c(0, max(-log10(p[keep]), na.rm = TRUE) * 1.05),
         main = sprintf("%s: %s SNPs with q < 0.05 (GIF %.2f)", vname(v), format(sum(q < qthr, na.rm = TRUE), big.mark = ","),
                        gif[variable == v, gif][1]), cex.main = 0.95)
    points(Lf$x[sig], -log10(p[sig]), pch = 16, cex = 0.45, col = CAT[(i - 1) %% 8 + 1])
    if (i == length(vars)) axis(1, at = CL$mid, labels = chrom_label(CL$chrom), las = 2, cex.axis = 0.55, tick = FALSE)
    if (i == 1) legend("topright", c("q < 0.05 for this predictor", "Not significant"), pch = 16,
                       col = c(CAT[1], GREY1), bty = "n", cex = 0.75)
    par(mar = c(if (i == length(vars)) 4 else 1.5, 3, 2, 1))
    hist(p, breaks = 40, col = GREY1, border = "white", main = "Calibrated p-values", xlab = if (i == length(vars)) "p" else "",
         ylab = "", cex.main = 0.8, yaxt = "n")
  }
  dev.off(); msg("wrote figures/lfmm_per_predictor_with_legend.png")
  rm(Lf, R); gc(verbose = FALSE)

  # ---------------- RDA SNP loadings (RDA1 x RDA2)
  RL <- fread(f("rda/rda_snps.tsv.gz"), select = c("snp_id", "q_rda", "RDA1_loading", "RDA2_loading", "rda_top_env"))
  if (all(c("RDA1_loading", "RDA2_loading") %in% names(RL))) {
    out <- RL[q_rda < qthr]; set.seed(2); bgL <- RL[q_rda >= qthr][sample(.N, min(.N, 200000))]
    envs <- sort(unique(out$rda_top_env)); ecol <- setNames(CAT[(seq_along(envs) - 1) %% 8 + 1], envs)
    open_png("rda_loadings_with_legend", 10, 8)
    par(mar = c(4.5, 4.5, 4, 15))
    plot(bgL$RDA1_loading, bgL$RDA2_loading, pch = 16, cex = 0.3, col = GREY2, xlab = "SNP loading on RDA axis 1",
         ylab = "SNP loading on RDA axis 2",
         xlim = range(c(bgL$RDA1_loading, out$RDA1_loading)), ylim = range(c(bgL$RDA2_loading, out$RDA2_loading)),
         main = "RDA SNP loadings: outlier SNPs coloured by the predictor they correlate with most")
    abline(h = 0, v = 0, col = GRID)
    points(out$RDA1_loading, out$RDA2_loading, pch = 16, cex = 0.5, col = ecol[out$rda_top_env])
    par(xpd = NA)
    cnt <- out[, .N, by = rda_top_env]; cn <- setNames(cnt$N, cnt$rda_top_env)
    legend(par("usr")[2] + diff(par("usr")[1:2]) * 0.02, par("usr")[4],
           c(sprintf("%s (n = %s)", vname(envs), format(cn[envs], big.mark = ",")), "Not an RDA outlier (random 200,000 shown)"),
           pch = 16, col = c(ecol, GREY2), pt.cex = 1.3, bty = "n", cex = 0.75, title = "RDA outlier SNPs (q < 0.05)", title.adj = 0)
    dev.off(); msg("wrote figures/rda_loadings_with_legend.png")
  }
  rm(RL); gc(verbose = FALSE)
} else msg("SKIP Manhattan/RDA/LFMM figures: rda/rda_snps.tsv.gz or lfmm/lfmm_snps.tsv.gz not found")

# =====================================================================
# 2. Variance partitioning
# =====================================================================
if (have("rda/variance_partitioning.tsv")) {
  vp <- fread(f("rda/variance_partitioning.tsv"))
  if ("prop_total" %in% names(vp) && nrow(vp) >= 4) {
    comp <- c("Environment (pure)" = vp$prop_total[2], "Geography (pure)" = vp$prop_total[4],
              "Population structure (pure)" = vp$prop_total[3],
              "Shared / confounded" = max(0, vp$prop_total[1] - sum(vp$prop_total[2:4])),
              "Unexplained" = 1 - vp$prop_total[1])
    cols <- c(CAT[3], CAT[1], CAT[2], GREY1, GRID)
    both("variance_partitioning_with_legend", 9, 4.5, function() {
      par(mar = c(4.5, 2, 4, 2))
      bp <- barplot(matrix(comp * 100), horiz = TRUE, col = cols, border = "white", xlab = "% of total genetic variance",
                    main = "Variance partitioning (partial RDAs)", xlim = c(0, 100), space = 0.6)
      x0 <- c(0, cumsum(comp * 100))[-length(comp) - 1] + comp * 100 / 2
      lab <- sprintf("%.1f%%", comp * 100); lab[comp < 0.03] <- ""
      text(x0, bp, lab, cex = 0.8, col = INK)
      legend("bottom", inset = c(0, -0.62), xpd = NA, names(comp), fill = cols, border = NA, bty = "n", ncol = 3, cex = 0.8)
    })
  }
}

# =====================================================================
# 3. Population structure
# =====================================================================
if (have("structure/ancestry_Q.tsv") && have("structure/pcs.tsv")) {
  Q <- fread(f("structure/ancestry_Q.tsv")); P <- fread(f("structure/pcs.tsv"))
  CE <- if (have("structure/cross_entropy.tsv")) fread(f("structure/cross_entropy.tsv")) else NULL
  qc <- grep("^Q", names(Q), value = TRUE); K <- length(qc); kc <- CAT[(seq_len(K) - 1) %% 8 + 1]
  Qm <- as.matrix(Q[, ..qc]); cl <- max.col(Qm); Q[, cl := cl]
  M <- merge(P, Q[, .(sample_id, cl)], by = "sample_id")
  pcv <- sapply(P[, -1], var); pcv <- pcv / sum(pcv) * 100
  both("structure_with_legends", 12, 9, function() {
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 9))
    plot(seq_along(pcv), pcv, type = "b", pch = 16, col = CAT[1], xlab = "Principal component",
         ylab = "% variance (of first 20 PCs)", main = "A. PCA scree plot (LD-pruned SNPs)")
    plot(M$PC1, M$PC2, pch = 21, bg = kc[M$cl], col = "white", cex = 1.4, xlab = "PC1", ylab = "PC2",
         main = "B. PCA, coloured by main ancestry cluster")
    tb <- table(factor(cl, levels = seq_len(K)))
    legend("topright", inset = c(-0.32, 0), xpd = NA, sprintf("Cluster %d (n = %d)", seq_len(K), tb), pt.bg = kc,
           pch = 21, col = "white", pt.cex = 1.4, bty = "n", cex = 0.8, title = "sNMF cluster")
    if (!is.null(CE)) {
      kb <- CE$K[which.min(CE$cross_entropy)]
      plot(CE$K, CE$cross_entropy, type = "b", pch = 16, col = INK2, xlab = "Number of ancestral populations (K)",
           ylab = "Cross-entropy (lower = better)", main = "C. Choice of K (sNMF)")
      points(kb, min(CE$cross_entropy), pch = 16, cex = 2, col = CAT[2])
      legend("topright", sprintf("Chosen K = %d", kb), pch = 16, col = CAT[2], pt.cex = 1.6, bty = "n", cex = 0.8)
    } else plot.new()
    o <- order(cl, -apply(Qm, 1, max))
    barplot(t(Qm[o, , drop = FALSE]), col = kc, border = NA, space = 0, ylab = "Ancestry proportion",
            xlab = sprintf("Individuals (n = %d), grouped by main cluster", nrow(Qm)), main = sprintf("D. Ancestry proportions, K = %d", K))
    legend("topright", inset = c(-0.3, 0), xpd = NA, sprintf("Cluster %d", seq_len(K)), fill = kc, border = NA,
           bty = "n", cex = 0.8, title = "Ancestral cluster")
  })
}

# =====================================================================
# 4. Environmental predictors
# =====================================================================
if (have("env/env_all_candidates_raw.tsv")) {
  E <- fread(f("env/env_all_candidates_raw.tsv")); E[, sample_id := NULL]
  E <- E[, which(sapply(E, function(x) sd(x, na.rm = TRUE) > 0)), with = FALSE]
  sel <- if (have("env/pruning_log.tsv")) fread(f("env/pruning_log.tsv"))[step == "selected", variable] else character(0)
  Rm <- cor(E); o <- hclust(as.dist(1 - abs(Rm)), "average")$order; Rm <- Rm[o, o]
  pal <- colorRampPalette(c("#184f95", "#6da7ec", "#f0efec", "#f19c9b", "#e34948"))(41)
  labs <- paste0(ifelse(colnames(Rm) %in% sel, "* ", ""), vname(colnames(Rm)))
  both("env_correlations_with_legend", 11, 10, function() {
    layout(matrix(1:2, 1), widths = c(6, 0.8)); par(mar = c(13, 13, 4, 1))
    image(seq_len(ncol(Rm)), seq_len(ncol(Rm)), Rm, zlim = c(-1, 1), col = pal, axes = FALSE, xlab = "", ylab = "",
          main = "Correlations among candidate predictors\n(* = kept for RDA and LFMM2: |r| < 0.7, VIF < 5)")
    abline(v = seq_len(ncol(Rm) - 1) + 0.5, h = seq_len(ncol(Rm) - 1) + 0.5, col = "white")
    axis(1, seq_len(ncol(Rm)), labs, las = 2, cex.axis = 0.7, tick = FALSE, font = 1)
    axis(2, seq_len(ncol(Rm)), labs, las = 2, cex.axis = 0.7, tick = FALSE)
    par(mar = c(13, 0.5, 4, 3.5))
    z <- seq(-1, 1, length.out = 41)
    image(1, z, matrix(z, 1), col = pal, axes = FALSE, xlab = "", ylab = "")
    axis(4, at = seq(-1, 1, 0.5), las = 1, cex.axis = 0.8)
    mtext("Pearson r", side = 4, line = 2.3, las = 0, cex = 0.8)
  })
}

# =====================================================================
# 5. Gene-level summaries
# =====================================================================
if (have("genes.tsv")) {
  G <- fread(f("genes.tsv"))
  G[, env_tier1 := as.character(env_tier1)]
  called <- G[called == TRUE][order(p_enrich)]
  msg(nrow(called), " called genes")
  parse_env <- function(s) {
    s <- s[!is.na(s) & nzchar(s)]; if (!length(s)) return(integer(0))
    tok <- unlist(strsplit(s, "[;|]")); tok <- tok[nzchar(tok)]
    nm <- sub("\\(.*$", "", tok); n <- suppressWarnings(as.integer(sub("^.*\\((\\d+)\\)$", "\\1", tok))); n[is.na(n)] <- 1L
    tapply(n, nm, sum)
  }
  if (nrow(called)) {
    top <- head(called, 40)
    lab <- ifelse(!is.na(top$name) & nzchar(top$name) & top$name != top$gene_id, paste0(top$name, "  (", top$gene_id, ")"), top$gene_id)
    drv <- vapply(top$env_tier1, function(s) { x <- parse_env(s); if (length(x)) names(x)[which.max(x)] else "unassigned" }, "")
    envs <- unique(drv); ecol <- setNames(CAT[(seq_along(envs) - 1) %% 8 + 1], envs); ecol["unassigned"] <- GREY1
    both("called_genes_with_legend", 10, 2 + 0.22 * nrow(top), function() {
      par(mar = c(4.5, 16, 4, 12))
      bp <- barplot(rev(top$n_tier1), horiz = TRUE, names.arg = rev(lab), las = 1, cex.names = 0.65, border = NA,
                    col = rev(ecol[drv]), xlab = "Tier-1 SNPs in gene \u00b1 5 kb (significant in RDA and LFMM2)",
                    main = sprintf("Environment-associated genes (top %d of %d called), coloured by main predictor", nrow(top), nrow(called)))
      text(rev(top$n_tier1), bp, sprintf("  q = %.1e", rev(top$q_enrich)), adj = 0, cex = 0.55, col = INK2, xpd = NA)
      legend("topright", inset = c(-0.38, 0), xpd = NA, vname(names(ecol)), fill = ecol, border = NA, bty = "n",
             cex = 0.75, title = "Main predictor")
    })

    # gene x predictor heatmap with colour key
    vars <- if (have("env/env_selected.tsv")) setdiff(names(fread(f("env/env_selected.tsv"), nrows = 1)), "sample_id") else character(0)
    hg <- head(called, 60)
    Mx <- t(vapply(hg$env_tier1, function(s) { x <- parse_env(s); v <- setNames(integer(length(vars)), vars)
                                               k <- intersect(names(x), vars); v[k] <- as.integer(x[k]); v }, integer(length(vars))))
    if (length(vars) == 1) Mx <- matrix(Mx, ncol = 1, dimnames = list(NULL, vars))
    rownames(Mx) <- ifelse(!is.na(hg$name) & nzchar(hg$name), hg$name, hg$gene_id)
    Mx <- Mx[, colSums(Mx) > 0, drop = FALSE]
    if (ncol(Mx) && nrow(Mx)) {
      if (nrow(Mx) > 2) Mx <- Mx[hclust(dist(Mx > 0), "average")$order, , drop = FALSE]
      mx <- max(Mx); nb <- min(mx, 8)
      brks <- c(-0.5, 0.5, seq(0.5, mx + 0.5, length.out = nb + 1)[-1])
      cols <- c("#f7f7f5", colorRampPalette(c("#cde2fb", "#3987e5", "#0d366b"))(length(brks) - 2))
      both("gene_env_heatmap_with_legend", 3.5 + 0.9 * ncol(Mx), 2.5 + 0.2 * nrow(Mx), function() {
        layout(matrix(1:2, 1), widths = c(5, 1.3)); par(mar = c(10, 9, 4, 1))
        image(seq_len(ncol(Mx)), seq_len(nrow(Mx)), t(Mx), breaks = brks, col = cols, axes = FALSE, xlab = "", ylab = "",
              main = "Tier-1 SNPs per gene and predictor")
        abline(v = seq_len(ncol(Mx) - 1) + 0.5, h = seq_len(nrow(Mx) - 1) + 0.5, col = "white", lwd = 1.5)
        axis(1, seq_len(ncol(Mx)), vname(colnames(Mx)), las = 2, tick = FALSE, cex.axis = 0.75)
        axis(2, seq_len(nrow(Mx)), rownames(Mx), las = 1, tick = FALSE, cex.axis = 0.6)
        nz <- which(Mx > 0, arr.ind = TRUE)
        text(nz[, 2], nz[, 1], Mx[nz], cex = 0.55, col = ifelse(Mx[nz] > mx / 2, "white", INK))
        par(mar = c(10, 0, 4, 0)); plot.new()
        pb <- brks[-1]; lo <- ceiling(pb[-length(pb)]); hi <- floor(pb[-1])
        lab <- c("0", ifelse(lo == hi, as.character(lo), paste0(lo, "-", hi)))
        legend("left", lab, fill = cols, border = GREY1, bty = "n", cex = 0.75, title = "Tier-1 SNPs")
      })
    }
  }
}

# =====================================================================
# 6. Allele frequency vs. environment for the lead SNP of top genes
# =====================================================================
if (have("plots/lead_snps.tsv") && have("geno/imputed.rds") && have("env/env_all_candidates_raw.tsv")) {
  lead <- fread(f("plots/lead_snps.tsv"))
  if (nrow(lead)) {
    msg("reading genotype matrix for allele-frequency plots")
    Gm <- readRDS(f("geno/imputed.rds"))
    raw <- fread(f("env/env_all_candidates_raw.tsv"), colClasses = list(character = "sample_id"))
    raw <- raw[match(rownames(Gm), sample_id)]
    # group individuals sharing identical coordinates/environment into sites
    site <- if (have("env/samples_clean.csv")) {
      sc <- fread(f("env/samples_clean.csv"), colClasses = list(character = "sample_id"))
      sc <- sc[match(rownames(Gm), sample_id)]
      if ("site" %in% names(sc)) as.character(sc$site)
      else if (all(c("lat", "lon") %in% names(sc))) paste(round(sc$lat, 3), round(sc$lon, 3))
      else do.call(paste, raw[, -1])
    } else do.call(paste, raw[, -1])
    gn <- if (have("genes.tsv")) fread(f("genes.tsv"))[, .(gene_id, name)] else data.table(gene_id = character(), name = character())
    lead <- merge(lead, gn, by = "gene_id", all.x = TRUE, sort = FALSE)
    lead <- lead[snp_id %in% colnames(Gm) & driver %in% names(raw)]
    np <- min(12, nrow(lead))
    if (np) {
      both("allele_freq_vs_env_with_legend", 12, 3 + 2.6 * ceiling(np / 4), function() {
        par(mfrow = c(ceiling(np / 4), 4), mar = c(4.2, 4.2, 3.2, 0.8), oma = c(3, 0, 2.5, 0))
        for (i in seq_len(np)) {
          L <- lead[i]; g <- Gm[, L$snp_id]; e <- raw[[L$driver]]
          d <- data.table(site, g, e)[, .(freq = sum(g) / (2 * .N), env = mean(e), n = .N), by = site]
          plot(d$env, d$freq, type = "n", ylim = c(0, 1), xlab = vname(L$driver), ylab = "Allele frequency",
               main = ifelse(!is.na(L$name) & nzchar(L$name), L$name, L$gene_id), cex.main = 0.95, cex.lab = 0.8)
          abline(h = c(0.25, 0.5, 0.75), col = GRID)
          fit <- tryCatch(glm(cbind(g, 2 - g) ~ e, family = binomial), error = function(err) NULL)
          if (!is.null(fit)) { xs <- seq(min(e), max(e), length.out = 100)
            lines(xs, predict(fit, data.frame(e = xs), type = "response"), col = INK, lwd = 2) }
          points(d$env, d$freq, pch = 21, bg = CAT[1], col = "white", cex = 0.6 + 1.8 * sqrt(d$n / max(d$n)))
          mtext(sprintf("%s  LFMM2 q = %.1e", L$snp_id, L$q_lfmm), side = 3, line = 0.2, cex = 0.5, col = INK2)
        }
        mtext("Lead SNP of each top environment-associated gene: allele frequency along its main predictor",
              outer = TRUE, line = 0.8, font = 2, cex = 0.9)
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
        legend("bottom", horiz = TRUE, bty = "n", cex = 0.85, xpd = NA,
               c("Sampling site (size = number of snakes)", "Logistic fit to individual genotypes"),
               pch = c(21, NA), pt.bg = c(CAT[1], NA), col = c("white", INK), lty = c(NA, 1), lwd = c(NA, 2), pt.cex = 1.6)
      })
    }
  }
} else msg("SKIP allele-frequency figure: needs plots/lead_snps.tsv and geno/imputed.rds")

msg("done; figures are in ", OUT)
