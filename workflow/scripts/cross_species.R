#!/usr/bin/env Rscript
# Compare environment-associated genes across species (each analysed on its own).
#  * gene x species table of calls and associated predictors
#  * pairwise overlap tests (hypergeometric, universe = genes testable in both)
#  * genes called in >= 2 species
# If species were mapped to different references, provide config `orthology:`
# (TSV with columns species, gene_id, orthogroup) and comparison is by orthogroup.
#
# Figures: overlap.pdf, upset.pdf, pairwise_heatmap.pdf, shared_gene_drivers.pdf,
#          variance_partitioning.pdf (if --varpart given)
#
# Rscript cross_species.R --config config.yaml --species sppA,sppB --genes a.tsv b.tsv \
#   [--varpart a_vp.tsv b_vp.tsv] --outdir out
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
a <- parse_args(); need(a, c("config", "species", "genes", "outdir"))
C <- read_config(a$config)
sp <- strsplit(a$species, ",")[[1]]
files <- a$genes
stopifnot(length(sp) == length(files))
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

orth <- NULL
if (!is.null(C$orthology) && nzchar(C$orthology)) orth <- fread(C$orthology)

L <- rbindlist(lapply(seq_along(sp), function(i) {
  g <- fread(files[i])[, .(gene_id, n_snps, n_tier1, called, env_tier1)]
  g[, env_tier1 := as.character(env_tier1)]
  g[is.na(env_tier1), env_tier1 := ""]
  g[, species := sp[i]]
  if (!is.null(orth)) {
    g <- merge(g, orth[species == sp[i], .(gene_id, orthogroup)], by = "gene_id")
    g[, unit := orthogroup]
  } else g[, unit := gene_id]
  # collapse to one row per unit (an orthogroup may hold several paralogs)
  g[, .(tested = any(n_snps > 0), called = any(called), n_tier1 = sum(n_tier1),
        env_tier1 = paste(unique(env_tier1[nzchar(env_tier1)]), collapse = "|")), by = .(unit, species)]
}))

wide <- dcast(L, unit ~ species, value.var = "called", fill = FALSE)
env_w <- dcast(L[called == TRUE], unit ~ species, value.var = "env_tier1", fill = "")
setnames(env_w, sp[sp %in% names(env_w)], paste0("env_", sp[sp %in% names(env_w)]))
wide[, n_species_called := rowSums(.SD), .SDcols = sp]
out <- merge(wide[n_species_called > 0], env_w, by = "unit", all.x = TRUE)
setorder(out, -n_species_called, unit)
setnames(out, "unit", if (is.null(orth)) "gene_id" else "orthogroup")
fwrite(out, file.path(a$outdir, "gene_overlap.tsv"), sep = "\t")

# pairwise tests
pairs <- combn(sp, 2, simplify = FALSE)
pt <- rbindlist(lapply(pairs, function(p) {
  A <- L[species == p[1]]; B <- L[species == p[2]]
  U <- intersect(A[tested == TRUE, unit], B[tested == TRUE, unit])
  a1 <- intersect(A[called == TRUE, unit], U); b1 <- intersect(B[called == TRUE, unit], U)
  k <- length(intersect(a1, b1))
  pval <- phyper(k - 1, length(a1), length(U) - length(a1), length(b1), lower.tail = FALSE)
  data.table(sp1 = p[1], sp2 = p[2], universe = length(U), called_sp1 = length(a1), called_sp2 = length(b1),
             shared = k, expected = length(a1) * length(b1) / max(length(U), 1),
             fold = k / max(length(a1) * length(b1) / max(length(U), 1), 1e-9),
             jaccard = k / max(length(union(a1, b1)), 1), p_hyper = pval)
}))
pt[, q_hyper := p.adjust(p_hyper, "BH")]
fwrite(pt, file.path(a$outdir, "pairwise_tests.tsv"), sep = "\t")
print(pt)
msg("genes called in >= 2 species: ", sum(out$n_species_called >= 2))

# combination counts (UpSet-style bar chart in base R)
combo <- out[, .(combo = paste(sp[unlist(.SD)], collapse = " & ")), by = seq_len(nrow(out)), .SDcols = sp][, .N, by = combo]
setorder(combo, -N)
pdf(file.path(a$outdir, "overlap.pdf"), width = 9, height = 6)
par(mar = c(4, 12, 2, 1))
if (nrow(combo)) barplot(rev(combo$N), names.arg = rev(combo$combo), horiz = TRUE, las = 1, cex.names = 0.7,
        main = "Called genes by species combination", xlab = "genes")
dev.off()

# =====================================================================
# Figures
# =====================================================================
nsp <- length(sp)
sp_col <- setNames(PAL$cat[(seq_len(nsp) - 1) %% length(PAL$cat) + 1], sp)
empty_pdf <- function(path, text) { pdf(path, 7, 4); plot.new(); text(0.5, 0.5, text, col = PAL$ink2); dev.off() }

# ---------------- UpSet plot (base R): intersection sizes over a dot matrix
up <- file.path(a$outdir, "upset.pdf")
if (!nrow(out)) {
  empty_pdf(up, "No genes called in any species")
} else {
  memb <- as.matrix(out[, ..sp])
  key <- apply(memb, 1, function(r) paste(as.integer(r), collapse = ""))
  tab <- sort(table(key), decreasing = TRUE)
  combos <- do.call(rbind, lapply(strsplit(names(tab), ""), as.integer))
  colnames(combos) <- sp
  set_size <- colSums(memb)
  nb <- length(tab)
  pdf(up, width = max(6, 2.5 + 0.45 * nb), height = 3.5 + 0.35 * nsp)
  layout(matrix(c(0, 1, 3, 2), 2, byrow = TRUE), widths = c(2.2, max(3.5, 0.45 * nb)), heights = c(3, 0.35 * nsp + 0.6))
  # top: intersection bars
  par(mar = c(0.5, 3.5, 3, 1)); base_theme()
  xl <- c(0.2, 1.4 * nb + 0.2)
  bp <- barplot(as.integer(tab), col = PAL$ink2, border = NA, space = 0.4, ylab = "Genes", axes = TRUE,
                xlim = xl, xaxs = "i", ylim = c(0, max(tab) * 1.15), main = "Environment-associated genes shared across species")
  text(bp, as.integer(tab), as.integer(tab), pos = 3, cex = 0.7, col = PAL$ink2)
  # dot matrix
  par(mar = c(1, 3.5, 0.3, 1))
  plot(NA, xlim = xl, ylim = c(0.5, nsp + 0.5), axes = FALSE, xlab = "", ylab = "", xaxs = "i")
  for (j in seq_len(nb)) {
    on <- which(combos[j, ] == 1)
    points(rep(bp[j], nsp), seq_len(nsp), pch = 16, cex = 1.3, col = PAL$grid)
    if (length(on) > 1) segments(bp[j], min(on), bp[j], max(on), col = PAL$ink, lwd = 2)
    points(rep(bp[j], length(on)), on, pch = 21, bg = sp_col[on], col = "white", cex = 1.6)
  }
  # left: set sizes
  par(mar = c(1, 1, 0.3, 0.5))
  plot(NA, xlim = c(max(set_size) * 1.1, 0), ylim = c(0.5, nsp + 0.5), axes = FALSE, xlab = "", ylab = "")
  rect(set_size, seq_len(nsp) - 0.3, 0, seq_len(nsp) + 0.3, col = sp_col, border = NA)
  text(set_size, seq_len(nsp), paste0(sp, " (", set_size, ")  "), pos = 2, cex = 0.75, col = PAL$ink, xpd = NA)
  dev.off()
}

# ---------------- pairwise heatmap: fold enrichment (colour) + shared count
hm <- file.path(a$outdir, "pairwise_heatmap.pdf")
if (!nrow(pt)) {
  empty_pdf(hm, "Fewer than two species")
} else {
  F <- matrix(NA_real_, nsp, nsp, dimnames = list(sp, sp)); N <- F; Qm <- F
  for (i in seq_len(nrow(pt))) {
    F[pt$sp1[i], pt$sp2[i]] <- F[pt$sp2[i], pt$sp1[i]] <- log2(pmax(pt$fold[i], 1e-3))
    N[pt$sp1[i], pt$sp2[i]] <- N[pt$sp2[i], pt$sp1[i]] <- pt$shared[i]
    Qm[pt$sp1[i], pt$sp2[i]] <- Qm[pt$sp2[i], pt$sp1[i]] <- pt$q_hyper[i]
  }
  diag(N) <- vapply(sp, function(s) L[species == s & called == TRUE, .N], 1L)
  Fraw <- F
  lim <- min(3, max(1, abs(F), na.rm = TRUE))            # colour scale capped at 8x either way
  F <- pmax(pmin(F, lim), -lim)
  ncol_div <- 11
  div <- grDevices::colorRampPalette(c("#184f95", "#6da7ec", "#f0efec", "#f19c9b", "#e34948"))(ncol_div)
  pdf(hm, width = max(7, 2.5 + 0.9 * nsp), height = max(6, 2 + 0.9 * nsp))
  layout(matrix(1:2, 1), widths = c(5, 1))
  par(mar = c(5, 6, 3, 0.5)); base_theme()
  image(seq_len(nsp), seq_len(nsp), t(F[nsp:1, ]), zlim = c(-lim, lim), col = div, axes = FALSE, xlab = "", ylab = "",
        main = "Pairwise overlap of called genes")
  rect(seq_len(nsp) - 0.5, rev(seq_len(nsp)) - 0.5, seq_len(nsp) + 0.5, rev(seq_len(nsp)) + 0.5, col = PAL$grid, border = NA)
  abline(v = seq_len(nsp - 1) + 0.5, h = seq_len(nsp - 1) + 0.5, col = "white", lwd = 2)
  axis(1, seq_len(nsp), sp, las = 2, tick = FALSE); axis(2, seq_len(nsp), rev(sp), las = 1, tick = FALSE)
  for (i in seq_len(nsp)) for (j in seq_len(nsp)) {
    y <- nsp - i + 1
    lab <- if (i == j) paste0(N[i, j], "\ncalled") else
      paste0(N[i, j], " shared\n", sprintf("%.1fx", 2^Fraw[i, j]), if (!is.na(Qm[i, j]) && Qm[i, j] < 0.05) " *" else "")
    text(j, y, lab, cex = 0.65, col = PAL$ink)
  }
  par(mar = c(5, 0.5, 3, 3.5))
  z <- seq(-lim, lim, length.out = ncol_div)
  image(1, z, matrix(z, 1), col = div, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = pretty(z), labels = sprintf("%.2gx", 2^pretty(z)), las = 1, cex.axis = 0.7)
  mtext("fold vs. chance", side = 4, line = 2.4, cex = 0.7, col = PAL$ink2, las = 0)
  dev.off()
}

# ---------------- driver dot plot for genes called in >= 2 species
dr <- file.path(a$outdir, "shared_gene_drivers.pdf")
shared <- out[n_species_called >= 2]
if (!nrow(shared)) {
  empty_pdf(dr, "No genes called in two or more species")
} else {
  shared <- head(shared, 60)
  unit_col <- names(out)[1]
  D <- rbindlist(lapply(sp, function(s) {
    e <- if (paste0("env_", s) %in% names(shared)) shared[[paste0("env_", s)]] else rep("", nrow(shared))
    e[is.na(e)] <- ""
    data.table(unit = shared[[unit_col]], species = s, called = shared[[s]],
               driver = vapply(e, top_env_name, ""),
               n = vapply(e, function(x) { v <- parse_env_counts(x); if (length(v)) as.integer(sum(v)) else 0L }, 1L))
  }))
  D <- D[called == TRUE]
  D[is.na(driver), driver := "unassigned"]
  freq <- sort(table(D$driver), decreasing = TRUE)
  dcol <- env_colours(names(freq))
  pdf(dr, width = 3.5 + 0.8 * nsp, height = 2 + 0.25 * nrow(shared))
  layout(matrix(1:2, 1), widths = c(3 + 0.6 * nsp, 2))
  par(mar = c(5, 8, 3, 0.5)); base_theme()
  yy <- setNames(rev(seq_len(nrow(shared))), shared[[unit_col]])
  plot(NA, xlim = c(0.5, nsp + 0.5), ylim = c(0.5, nrow(shared) + 0.5), axes = FALSE, xlab = "", ylab = "",
       main = "Main predictor of genes shared across species")
  abline(h = yy, col = PAL$grid, lwd = 0.8)
  axis(1, seq_len(nsp), sp, las = 2, tick = FALSE)
  axis(2, yy, names(yy), las = 1, tick = FALSE, cex.axis = 0.65)
  sz <- 0.8 + 1.4 * sqrt(D$n / max(1, max(D$n)))
  points(match(D$species, sp), yy[D$unit], pch = 21, bg = dcol[D$driver], col = "white", lwd = 1, cex = sz)
  par(mar = c(5, 0, 3, 0.5)); plot.new()
  legend("topleft", names(dcol), pt.bg = dcol, pch = 21, col = "white", pt.cex = 1.6, bty = "n",
         title = "Predictor", cex = 0.8, text.col = PAL$ink2)
  legend("bottomleft", c("few", "many"), pch = 21, pt.cex = c(0.9, 2.2), pt.bg = PAL$bg, col = "white",
         bty = "n", title = "tier1 SNPs", cex = 0.8, text.col = PAL$ink2)
  dev.off()
}

# ---------------- variance partitioning across species
vpf <- file.path(a$outdir, "variance_partitioning.pdf")
if (is.null(a$varpart)) {
  empty_pdf(vpf, "No variance partitioning files given")
} else {
  vp_files <- a$varpart
  stopifnot(length(vp_files) == nsp)
  parts <- c("Environment", "Geography", "Structure", "Confounded")
  V <- sapply(seq_len(nsp), function(i) {
    v <- fread(vp_files[i])
    if (!"prop_total" %in% names(v) || nrow(v) < 4) return(rep(NA_real_, 4))
    c(v$prop_total[2], v$prop_total[4], v$prop_total[3], max(0, v$confounded[1], na.rm = TRUE))
  })
  V <- matrix(V, nrow = 4, dimnames = list(parts, sp))
  fwrite(data.table(component = parts, V), file.path(a$outdir, "variance_partitioning.tsv"), sep = "\t")
  if (all(is.na(V))) {
    empty_pdf(vpf, "Variance partitioning was switched off")
  } else {
    pcol <- PAL$cat[c(3, 1, 2)]; pcol <- c(pcol, PAL$bg)
    pdf(vpf, width = 3.5 + 0.9 * nsp, height = 5)
    par(mar = c(5, 4.5, 3, 9), xpd = NA); base_theme()
    Vp <- V * 100; Vp[is.na(Vp)] <- 0
    bp <- barplot(Vp, col = pcol, border = "white", lwd = 2, space = 0.5, las = 2,
                  ylab = "% of genetic variance explained", main = "Variance partitioning (RDA)",
                  ylim = c(0, max(colSums(Vp)) * 1.12))
    text(bp, colSums(Vp), sprintf("%.1f%%", colSums(Vp)), pos = 3, cex = 0.7, col = PAL$ink2)
    legend(par("usr")[2] + 0.2, max(colSums(Vp)), rev(parts), fill = rev(pcol), border = NA, bty = "n",
           cex = 0.8, text.col = PAL$ink2)
    dev.off()
  }
}
msg("figures written")
