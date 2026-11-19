#!/usr/bin/env Rscript
# Per-species markdown summary.
# Rscript report.R --species sppA --dir results/sppA --out results/sppA/report.md
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))
a <- parse_args(); need(a, c("species", "dir", "out"))
d <- a$dir
md_table <- function(dt, digits = 3) {
  dt <- copy(dt)
  for (j in names(dt)) if (is.numeric(dt[[j]])) set(dt, j = j, value = signif(dt[[j]], digits))
  c(paste0("| ", paste(names(dt), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
    apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
k <- fread(file.path(d, "structure/K.txt"), header = FALSE)
envlog <- fread(file.path(d, "env/pruning_log.tsv"))
vp <- fread(file.path(d, "rda/variance_partitioning.tsv"))
gif <- fread(file.path(d, "lfmm/gif.tsv"))
snps <- fread(file.path(d, "candidate_snps.tsv"))
genes <- fread(file.path(d, "genes.tsv"))
ev <- file.path(d, "evaluation.tsv")

L <- c(paste("# GEA report:", a$species), "",
       paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
       "## Population structure", "",
       sprintf("sNMF K = %s; structure PCs used as covariates = %s", k$V2[k$V1 == "K"], k$V2[k$V1 == "n_pcs"]), "",
       "## Environmental predictors", "",
       paste("Selected:", paste(envlog[step == "selected", variable], collapse = ", ")), "",
       md_table(envlog), "",
       "## Variance partitioning (RDA)", "", md_table(vp), "",
       "```", readLines(file.path(d, "rda/rda_summary.txt")), "```", "",
       "## LFMM2 genomic inflation factors", "", md_table(gif), "",
       "## Candidate SNPs", "", md_table(snps[, .N, by = tier]), "",
       "## Top genes", "",
       md_table(head(genes[, .(gene_id, chrom, start, end, n_snps, n_tier1, p_enrich, q_enrich, called, env_tier1)], 30)), "")
if (file.exists(ev)) L <- c(L, "## Evaluation against simulated truth", "", md_table(fread(ev)), "")
L <- c(L, "## Files", "",
       "- `candidate_snps.tsv`: every SNP significant in RDA and/or LFMM2 with gene context",
       "- `genes.tsv`: all genes with SNP counts, enrichment test and associated predictors",
       "- `windows.tsv.gz`: sliding-window enrichment",
       "- `plots/allele_freq_vs_env.pdf`: lead SNP of each top gene vs. its main predictor",
       "- `plots/regional_zoom.pdf`: RDA and LFMM2 signal around each top gene, with gene models",
       "- `plots/gene_env_heatmap.pdf`: called genes x predictors (tier1 SNP counts)",
       "- `manhattan.pdf`, `rda/rda_plots.pdf`, `lfmm/lfmm_plots.pdf`, `structure/structure.pdf`, `env/env_correlations.pdf`")
writeLines(L, a$out)
