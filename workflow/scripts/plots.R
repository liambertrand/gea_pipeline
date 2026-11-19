#!/usr/bin/env Rscript
# Per-species figures:
#   allele_freq_vs_env.pdf  lead SNP of each top gene: allele frequency vs its driver
#   regional_zoom.pdf       RDA and LFMM2 signal across each top gene +/- flank, with gene models
#   gene_env_heatmap.pdf    called genes x predictors, filled by number of tier1 SNPs
# plus lead_snps.tsv and gene_env_matrix.tsv (the numbers behind the figures).
#
# Rscript plots.R --config cfg.yaml --species sp --genes genes.tsv --snps candidate_snps.tsv \
#   --rda rda_snps.tsv.gz --lfmm lfmm_snps.tsv.gz --geno imputed.rds --env env_selected.tsv \
#   --raw env_all_candidates_raw.tsv --samples samples.csv --gff annotation.gff3 --outdir plots
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))

a <- parse_args()
need(a, c("config", "species", "genes", "snps", "rda", "lfmm", "geno", "env", "raw", "samples", "gff", "outdir"))
C <- read_config(a$config)
P <- C$plots
if (is.null(P)) P <- list()
top_n <- as.integer(if (is.null(P$top_genes)) 12 else P$top_genes)
flank <- as.integer(if (is.null(P$zoom_flank_bp)) 100000 else P$zoom_flank_bp)
pop_col <- if (is.null(P$pop_col)) "site" else P$pop_col
max_heat <- as.integer(if (is.null(P$heatmap_max_genes)) 60 else P$heatmap_max_genes)
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

genes <- fread(a$genes)
genes[, env_tier1 := as.character(env_tier1)]
cand <- fread(a$snps)
sel_vars <- setdiff(names(read_samples_tab(a$env)), "sample_id")
ecol <- env_colours(sel_vars)

# Genes to show: called genes first, then genes with the most tier1 SNPs
show <- genes[n_tier1 > 0][order(-as.integer(called), p_enrich, -n_tier1)]
show <- head(show, top_n)
msg(a$species, ": ", sum(genes$called), " called genes; plotting ", nrow(show))

empty_pdf <- function(path, text) {
  pdf(path, width = 7, height = 4); plot.new(); text(0.5, 0.5, text, col = PAL$ink2); dev.off()
}

# ---------------------------------------------------------------------
# Lead SNP per gene: tier1 SNP inside the gene window with the smallest LFMM q
# (ties broken by RDA q); driver = the gene's most frequent tier1 predictor
# ---------------------------------------------------------------------
cand[, genes := as.character(genes)]
cand[is.na(genes), genes := ""]
lead <- rbindlist(lapply(seq_len(nrow(show)), function(i) {
  g <- show$gene_id[i]
  s <- cand[tier == "tier1" & vapply(strsplit(genes, ";", fixed = TRUE), function(x) g %in% x, TRUE)]
  if (!nrow(s)) return(NULL)
  s <- s[order(lfmm_min_q, q_rda)][1]
  drv <- top_env_name(show$env_tier1[i])
  if (is.na(drv)) drv <- s$assoc_env
  drv <- strsplit(drv, ";", fixed = TRUE)[[1]][1]
  data.table(gene_id = g, snp_id = s$snp_id, chrom = s$chrom, pos = s$pos,
             q_rda = s$q_rda, q_lfmm = s$lfmm_min_q, driver = drv)
}))
fwrite(if (is.null(lead)) data.table() else lead, file.path(a$outdir, "lead_snps.tsv"), sep = "\t")

# ---------------------------------------------------------------------
# 1. Allele frequency vs environment
# ---------------------------------------------------------------------
G <- readRDS(a$geno)
raw <- read_samples_tab(a$raw)
stopifnot(identical(raw$sample_id, rownames(G)))
S <- fread(a$samples, colClasses = setNames(list(C$env$sample_id_col), "character"))
setnames(S, C$env$sample_id_col, "sample_id")
S <- S[match(rownames(G), sample_id)]
use_pop <- nzchar(pop_col) && pop_col %in% names(S)
if (nzchar(pop_col) && !use_pop) msg("column '", pop_col, "' not in samples table; plotting individuals")
pop <- if (use_pop) as.character(S[[pop_col]]) else NULL

afile <- file.path(a$outdir, "allele_freq_vs_env.pdf")
if (is.null(lead) || !nrow(lead)) {
  empty_pdf(afile, "No tier1 SNPs in genes: nothing to plot")
} else {
  pdf(afile, width = 11, height = 8.5)
  per_page <- 12
  for (start in seq(1, nrow(lead), by = per_page)) {
    rows <- start:min(nrow(lead), start + per_page - 1)
    par(mfrow = c(3, 4), mar = c(3.6, 3.6, 2.6, 0.8), oma = c(0, 0, 2, 0)); base_theme()
    for (i in rows) {
      L <- lead[i]
      if (!L$snp_id %in% colnames(G) || !L$driver %in% names(raw)) { plot.new(); next }
      g <- G[, L$snp_id]; e <- raw[[L$driver]]
      col <- ecol[L$driver]; if (is.na(col)) col <- PAL$cat[1]
      fit <- tryCatch(glm(cbind(g, 2 - g) ~ e, family = binomial), error = function(err) NULL)
      if (use_pop) {
        d <- data.table(pop, g, e)[, .(freq = sum(g) / (2 * .N), env = mean(e), n = .N), by = pop]
        plot(d$env, d$freq, type = "n", ylim = c(0, 1), xlab = L$driver, ylab = "Allele frequency",
             main = L$gene_id, bty = "l")
        abline(h = c(0.25, 0.5, 0.75), col = PAL$grid, lwd = 0.8)
        points(d$env, d$freq, pch = 21, bg = col, col = "white", lwd = 1.2, cex = 0.7 + 1.6 * sqrt(d$n / max(d$n)))
      } else {
        plot(e, g / 2 + runif(length(g), -0.04, 0.04), type = "n", ylim = c(-0.05, 1.05), xlab = L$driver,
             ylab = "Genotype / 2", main = L$gene_id, bty = "l")
        abline(h = c(0, 0.5, 1), col = PAL$grid, lwd = 0.8)
        points(e, g / 2 + runif(length(g), -0.04, 0.04), pch = 16, col = adjustcolor(col, 0.45), cex = 0.8)
      }
      if (!is.null(fit)) {
        xs <- seq(min(e), max(e), length.out = 100)
        lines(xs, predict(fit, data.frame(e = xs), type = "response"), col = PAL$ink, lwd = 2)
      }
      mtext(sprintf("%s  q(LFMM)=%.1e", L$snp_id, L$q_lfmm), side = 3, line = 0.2, cex = 0.6, col = PAL$ink2)
    }
    mtext(sprintf("%s: lead SNP of each top gene vs. its main predictor (%s)", a$species,
                  if (use_pop) "points = populations, size = n" else "points = individuals"),
          outer = TRUE, cex = 0.95, col = PAL$ink)
  }
  dev.off()
}

# ---------------------------------------------------------------------
# 2. Regional zoom plots
# ---------------------------------------------------------------------
rda <- fread(a$rda, select = c("snp_id", "chrom", "pos", "p_rda"))
lf <- fread(a$lfmm, select = c("snp_id", "lfmm_min_p"))
R <- merge(rda, lf, by = "snp_id")
R[, chrom := as.character(chrom)]
R <- merge(R, cand[, .(snp_id, tier)], by = "snp_id", all.x = TRUE)
R[is.na(tier), tier := "none"]
gff <- read_gff_genes(a$gff, C$gene_feature, C$gene_id_attr)
qthr <- as.numeric(C$candidates$q_threshold)

zfile <- file.path(a$outdir, "regional_zoom.pdf")
if (!nrow(show)) {
  empty_pdf(zfile, "No genes with tier1 SNPs")
} else {
  pdf(zfile, width = 10, height = 7)
  for (i in seq_len(nrow(show))) {
    gi <- show[i]
    lo <- max(1, gi$start - flank); hi <- gi$end + flank
    r <- R[chrom == gi$chrom & pos >= lo & pos <= hi][order(tier != "none")]  # highlighted on top
    gg <- gff[chrom == gi$chrom & end >= lo & start <= hi]
    layout(matrix(1:3, ncol = 1), heights = c(3, 3, 1.6))
    for (what in c("rda", "lfmm")) {
      par(mar = c(0.5, 4.5, if (what == "rda") 3 else 0.5, 1)); base_theme()
      y <- if (what == "rda") -log10(r$p_rda) else -log10(r$lfmm_min_p)
      plot(r$pos, y, type = "n", xlim = c(lo, hi), ylim = c(0, max(c(y, 2), na.rm = TRUE) * 1.05),
           xaxt = "n", xlab = "", ylab = if (what == "rda") "RDA  -log10 p" else "LFMM2  -log10 p", bty = "l")
      rect(gi$start, par("usr")[3], gi$end, par("usr")[4], col = adjustcolor(PAL$grid, 0.6), border = NA)
      points(r$pos, y, pch = 21, bg = PAL$tier[r$tier], col = "white", lwd = 0.6,
             cex = ifelse(r$tier == "none", 0.7, 1.1))
      if (what == "rda") {
        title(main = sprintf("%s  %s:%s-%s  (%d tier1 SNPs, q = %.1e)  %s", gi$gene_id, gi$chrom,
                             format(gi$start, big.mark = ","), format(gi$end, big.mark = ","),
                             gi$n_tier1, gi$q_enrich, ifelse(is.na(gi$env_tier1), "", gi$env_tier1)),
              cex.main = 0.95)
        legend("topright", c("tier1 (both)", "RDA only", "LFMM2 only", "not significant"), pt.bg = PAL$tier[c(4, 2, 3, 1)],
               pch = 21, col = "white", bty = "n", cex = 0.8, text.col = PAL$ink2)
      }
    }
    # gene track
    par(mar = c(3.5, 4.5, 0.3, 1)); base_theme()
    plot(NA, xlim = c(lo, hi), ylim = c(0, 1), yaxt = "n", ylab = "", xlab = paste(gi$chrom, "position (bp)"), bty = "n")
    if (nrow(gg)) {
      lane <- rep(c(0.65, 0.3), length.out = nrow(gg))
      focal <- gg$gene_id == gi$gene_id
      rect(pmax(gg$start, lo), lane - 0.1, pmin(gg$end, hi), lane + 0.1,
           col = ifelse(focal, PAL$tier["tier1"], PAL$bg), border = NA)
      arrows_x <- ifelse(gg$strand == "-", pmax(gg$start, lo), pmin(gg$end, hi))
      text((pmax(gg$start, lo) + pmin(gg$end, hi)) / 2, lane + 0.22,
           ifelse(is.na(gg$name), gg$gene_id, gg$name), cex = 0.6, col = ifelse(focal, PAL$ink, PAL$ink2))
      text(arrows_x, lane, ifelse(gg$strand == "-", "<", ">"), cex = 0.6, col = PAL$ink2)
    }
  }
  dev.off()
}

# ---------------------------------------------------------------------
# 3. Gene x predictor heatmap
# ---------------------------------------------------------------------
hg <- genes[called == TRUE][order(p_enrich)]
if (nrow(hg) > max_heat) { msg("heatmap limited to the top ", max_heat, " genes"); hg <- head(hg, max_heat) }
hfile <- file.path(a$outdir, "gene_env_heatmap.pdf")
if (!nrow(hg)) {
  fwrite(data.table(), file.path(a$outdir, "gene_env_matrix.tsv"), sep = "\t")
  empty_pdf(hfile, "No called genes")
} else {
  M <- t(vapply(hg$env_tier1, function(s) {
    x <- parse_env_counts(s); v <- setNames(integer(length(sel_vars)), sel_vars)
    k <- intersect(names(x), sel_vars); v[k] <- as.integer(x[k]); v
  }, integer(length(sel_vars))))
  if (length(sel_vars) == 1) M <- matrix(M, ncol = 1, dimnames = list(NULL, sel_vars))
  rownames(M) <- hg$gene_id
  fwrite(data.table(gene_id = rownames(M), M), file.path(a$outdir, "gene_env_matrix.tsv"), sep = "\t")
  M <- M[, colSums(M) > 0, drop = FALSE]
  if (!ncol(M)) {
    empty_pdf(hfile, "Called genes have no predictor attribution")
  } else {
    # cluster rows so genes sharing a driver sit together
    if (nrow(M) > 2) M <- M[hclust(dist(M > 0), "average")$order, , drop = FALSE]
    brks <- c(-0.5, 0.5, seq(1.5, max(M) + 0.5, length.out = min(max(M), 7)))
    cols <- c(PAL$seq[1], seq_ramp(length(PAL$seq))[-1][seq_len(length(brks) - 2)])
    pdf(hfile, width = 3 + 0.7 * ncol(M), height = 2 + 0.22 * nrow(M))
    par(mar = c(5, 8, 3, 1)); base_theme()
    image(seq_len(ncol(M)), seq_len(nrow(M)), t(M), breaks = brks, col = cols, axes = FALSE,
          xlab = "", ylab = "", main = paste(a$species, ": tier1 SNPs per gene and predictor"))
    abline(v = seq_len(ncol(M) - 1) + 0.5, h = seq_len(nrow(M) - 1) + 0.5, col = "white", lwd = 2)
    axis(1, seq_len(ncol(M)), colnames(M), las = 2, tick = FALSE, cex.axis = 0.8)
    axis(2, seq_len(nrow(M)), rownames(M), las = 1, tick = FALSE, cex.axis = 0.65)
    nz <- which(M > 0, arr.ind = TRUE)
    text(nz[, 2], nz[, 1], M[nz], cex = 0.6, col = ifelse(M[nz] >= max(M) / 2, "white", PAL$ink))
    dev.off()
  }
}
msg("done")
