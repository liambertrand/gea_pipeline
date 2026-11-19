#!/usr/bin/env Rscript
# Simulation only: score the pipeline against the known truth.
# Rscript evaluate.R --snps candidate_snps.tsv --genes genes.tsv --all-snps rda_snps.tsv.gz \
#   --truth truth_snps.tsv --truth-genes truth_genes.tsv --out evaluation.tsv
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
a <- parse_args(); need(a, c("snps", "genes", "all_snps", "truth", "truth_genes", "out"))

cand <- fread(a$snps); genes <- fread(a$genes); all <- fread(a$all_snps, select = "snp_id")
truth <- fread(a$truth); tg <- fread(a$truth_genes)
causal <- intersect(truth[role == "causal", snp_id], all$snp_id)
true_any <- truth$snp_id
msg(length(causal), " of ", truth[role == "causal", .N], " causal SNPs passed filtering")

snp_row <- function(label, ids) {
  data.table(level = "snp", set = label, n_called = length(ids),
             true_positives = sum(ids %in% causal), n_true = length(causal),
             recall = if (length(causal)) mean(causal %in% ids) else NA_real_,
             fdr = if (length(ids)) mean(!ids %in% true_any) else NA_real_)
}
rows <- list(
  snp_row("tier1 (RDA & LFMM)", cand[tier == "tier1", snp_id]),
  snp_row("RDA (q<thr)", cand[tier %in% c("tier1", "rda_only"), snp_id]),
  snp_row("LFMM (q<thr)", cand[tier %in% c("tier1", "lfmm_only"), snp_id]),
  snp_row("RDA or LFMM", cand$snp_id)
)
gene_row <- function(label, ids) {
  data.table(level = "gene", set = label, n_called = length(ids),
             true_positives = sum(ids %in% tg$gene_id), n_true = nrow(tg),
             recall = mean(tg$gene_id %in% ids),
             fdr = if (length(ids)) mean(!ids %in% tg$gene_id) else NA_real_)
}
rows <- c(rows, list(
  gene_row("called (enriched tier1)", genes[called == TRUE, gene_id]),
  gene_row("any tier1 SNP", genes[n_tier1 > 0, gene_id]),
  gene_row("any candidate SNP", genes[support != "none", gene_id])
))
res <- rbindlist(rows)

# did we attribute the right driver? (driver name must appear in env_tier1)
drv <- merge(tg[, .(gene_id, driver, class)], genes[, .(gene_id, called, env_tier1)], by = "gene_id")
drv[, driver_pattern := fifelse(driver == "forest", "lc_forest", driver)]
drv[, driver_found := mapply(function(p, e) grepl(paste0("(^|;)", p, "\\("), e), driver_pattern, env_tier1)]
fwrite(drv, sub("\\.tsv$", "_genes.tsv", a$out), sep = "\t")
res <- rbind(res, data.table(level = "gene", set = "called genes with correct driver",
                             n_called = drv[called == TRUE, .N],
                             true_positives = drv[called == TRUE & driver_found == TRUE, .N],
                             n_true = nrow(tg), recall = NA_real_, fdr = NA_real_), fill = TRUE)
res[, recall := round(recall, 3)][, fdr := round(fdr, 3)]
print(res)
fwrite(res, a$out, sep = "\t")
