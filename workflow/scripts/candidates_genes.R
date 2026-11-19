#!/usr/bin/env Rscript
# Combine RDA and LFMM2, then move from SNPs to windows and genes.
#   tier1      = q < threshold in RDA AND in LFMM2 (>= 1 predictor)
#   rda_only / lfmm_only = one method only (tier2)
# Windows and genes are tested for an excess of tier1 SNPs relative to the
# genome-wide tier1 rate (one-sided binomial), which accounts for the number
# of SNPs they contain (so long / SNP-dense genes are not favoured).
#
# Rscript candidates_genes.R --config config.yaml --rda rda_snps.tsv.gz \
#   --lfmm lfmm_snps.tsv.gz --gff annotation.gff3 --outdir results/sp
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))

a <- parse_args(); need(a, c("config", "rda", "lfmm", "gff", "outdir"))
C <- read_config(a$config)
qthr <- as.numeric(C$candidates$q_threshold)

rda <- fread(a$rda); lf <- fread(a$lfmm)
keep_lf <- c("snp_id", "lfmm_min_p", "lfmm_min_q", "lfmm_best_env", "lfmm_sig_env")
S <- merge(rda, lf[, ..keep_lf], by = "snp_id", sort = FALSE)
S[is.na(lfmm_sig_env), lfmm_sig_env := ""]
S[, chrom := as.character(chrom)]
S[, rda_sig := !is.na(q_rda) & q_rda < qthr]
S[, lfmm_sig := !is.na(lfmm_min_q) & lfmm_min_q < qthr]
S[, tier := fifelse(rda_sig & lfmm_sig, "tier1", fifelse(rda_sig, "rda_only", fifelse(lfmm_sig, "lfmm_only", "none")))]
# associated predictor: LFMM's significant predictors if any, else RDA's closest predictor
S[, assoc_env := fifelse(nzchar(lfmm_sig_env), lfmm_sig_env, rda_top_env)]
msg("SNPs: ", nrow(S), "; ", paste(names(table(S$tier)), table(S$tier), sep = "=", collapse = " "))
rate1 <- mean(S$tier == "tier1")
rate_any <- mean(S$tier != "none")

binom_p <- function(k, n, r) ifelse(k > 0 & n > 0, pbinom(k - 1, n, r, lower.tail = FALSE), 1)

# ---------------- windows
wsize <- as.integer(C$windows$size_bp); wstep <- as.integer(C$windows$step_bp)
W <- S[, .(start = seq(1L, max(pos), by = wstep)), by = chrom]
W[, end := start + wsize - 1L]
Sx <- S[, .(chrom, start = pos, end = pos, tier)]
setkey(W, chrom, start, end)
ov <- foverlaps(Sx, W, nomatch = 0L)
Wc <- ov[, .(n_snps = .N, n_tier1 = sum(tier == "tier1"), n_any = sum(tier != "none")), by = .(chrom, start, end)]
Wc[, p_tier1 := binom_p(n_tier1, n_snps, rate1)]
Wc[, q_tier1 := p.adjust(p_tier1, "BH")]
setorder(Wc, p_tier1)
fwrite(Wc, file.path(a$outdir, "windows.tsv.gz"), sep = "\t")
msg("windows with q < ", qthr, ": ", sum(Wc$q_tier1 < qthr))

# ---------------- genes
genes <- read_gff_genes(a$gff, C$gene_feature, C$gene_id_attr)
flank <- as.integer(C$genes$flank_bp)
Gx <- genes[, .(chrom, start = pmax(1L, start - flank), end = end + flank, gene_id)]
setkey(Gx, chrom, start, end)
Sy <- S[, .(snp_id, chrom, start = pos, end = pos, tier, q_rda, lfmm_min_q, assoc_env)]
og <- foverlaps(Sy, Gx, nomatch = 0L)
top_env <- function(e) {
  e <- unlist(strsplit(e[nzchar(e)], ";", fixed = TRUE))
  if (!length(e)) return("")
  t <- sort(table(e), decreasing = TRUE)
  paste(sprintf("%s(%d)", names(t), as.integer(t)), collapse = ";")
}
GT <- og[, .(n_snps = .N,
             n_tier1 = sum(tier == "tier1"),
             n_rda = sum(tier %in% c("tier1", "rda_only")),
             n_lfmm = sum(tier %in% c("tier1", "lfmm_only")),
             min_q_rda = min(q_rda, na.rm = TRUE),
             min_q_lfmm = min(lfmm_min_q, na.rm = TRUE),
             env_tier1 = top_env(assoc_env[tier == "tier1"]),
             env_any = top_env(assoc_env[tier != "none"])), by = gene_id]
GT <- merge(genes[, .(gene_id, name, chrom, start, end, strand)], GT, by = "gene_id", all.x = TRUE)
GT[is.na(n_snps), `:=`(n_snps = 0L, n_tier1 = 0L, n_rda = 0L, n_lfmm = 0L)]
GT[, p_enrich := binom_p(n_tier1, n_snps, rate1)]
GT[, q_enrich := p.adjust(p_enrich, "BH")]
GT[, called := q_enrich < as.numeric(C$genes$q_threshold) & n_tier1 >= as.integer(C$genes$min_tier1_snps)]
GT[, support := fifelse(called, "enriched_tier1",
                 fifelse(n_tier1 > 0, "has_tier1",
                 fifelse(n_rda + n_lfmm > 0, "tier2_only", "none")))]
GT <- GT[order(-as.integer(called), p_enrich, -n_tier1)]
fwrite(GT, file.path(a$outdir, "genes.tsv"), sep = "\t")
msg("genes called: ", sum(GT$called), " of ", nrow(GT), " (", sum(GT$n_snps > 0), " with SNPs)")

# ---------------- candidate SNP table with gene context
cand <- S[tier != "none"]
cg <- og[tier != "none", .(genes = paste(unique(gene_id), collapse = ";")), by = snp_id]
cand <- merge(cand, cg, by = "snp_id", all.x = TRUE)
cand[is.na(genes), genes := ""]
# nearest gene (TSS-agnostic, edge distance) for intergenic candidates
g_by_chr <- split(genes, genes$chrom)
cand[, c("nearest_gene", "dist_to_gene") := {
  gg <- g_by_chr[[chrom[1]]]
  if (is.null(gg)) list(NA_character_, NA_integer_) else {
    d <- outer(pos, gg$start, "-"); d2 <- outer(pos, gg$end, "-")
    inside <- d >= 0 & d2 <= 0
    dist <- ifelse(inside, 0L, pmin(abs(d), abs(d2)))
    j <- max.col(-dist, ties.method = "first")
    list(gg$gene_id[j], as.integer(dist[cbind(seq_along(pos), j)]))
  }
}, by = chrom]
setorder(cand, tier, q_rda)
cols <- c("snp_id", "chrom", "pos", "tier", "q_rda", "lfmm_min_q", "assoc_env", "lfmm_sig_env",
          "rda_top_env", "rda_top_env_r", "genes", "nearest_gene", "dist_to_gene")
fwrite(cand[, ..cols], file.path(a$outdir, "candidate_snps.tsv"), sep = "\t")

# ---------------- Manhattan plots
chr_len <- S[, .(len = max(pos)), by = chrom]
chr_len[, offset := c(0, cumsum(as.numeric(len))[-.N])]
S2 <- merge(S, chr_len[, .(chrom, offset)], by = "chrom")
S2[, x := pos + offset]
cols_chr <- ifelse(match(S2$chrom, chr_len$chrom) %% 2 == 0, "grey40", "grey70")
tcol <- PAL$tier
mids <- chr_len[, offset + len / 2]
pdf(file.path(a$outdir, "manhattan.pdf"), width = 12, height = 8)
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))
for (what in c("rda", "lfmm")) {
  y <- if (what == "rda") -log10(S2$p_rda) else -log10(S2$lfmm_min_p)
  col <- ifelse(S2$tier == "none", cols_chr, tcol[S2$tier])
  plot(S2$x, y, pch = 20, cex = 0.35, col = col, xaxt = "n", xlab = "", ylab = "-log10 p",
       main = if (what == "rda") "RDA (rdadapt)" else "LFMM2 (min over predictors)")
  axis(1, mids, chr_len$chrom, las = 2, cex.axis = 0.6)
  legend("topright", legend = names(tcol)[-1], col = tcol[-1], pch = 20, bty = "n", cex = 0.8)
}
dev.off()
msg("done")
