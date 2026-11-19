#!/usr/bin/env Rscript
# =====================================================================
# snp_annotation_classes.R
# Classifies every tested SNP as CDS / UTR-or-other-exon / intron /
# within 5 kb of a gene / intergenic, and compares candidate sets
# (tier 1, RDA only, LFMM2 only) against all tested SNPs.
# Also exports a gene -> homolog symbol table for GO / STRING / MGI.
#
# Usage (from gea-pipeline, with the pipeline's own R):
#   ENV=$(ls -d .snakemake/conda/*_ | head -1)
#   $ENV/bin/Rscript snp_annotation_classes.R results_viridis/Cviridis \
#       ../Ref_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf
# Output (in <results dir>/annotation/):
#   snp_class_summary.tsv   counts and % per class for each SNP set, with Fisher tests
#   tier1_snps_classified.tsv
#   gene_homologs.tsv       gene_id, Anolis/Python/Thamnophis homolog symbols, best symbol
# =====================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
D <- a[1]; GTF <- a[2]
OUT <- file.path(D, "annotation"); dir.create(OUT, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

msg("reading GTF")
g <- fread(cmd = paste("grep -v '^#'", shQuote(GTF)), sep = "\t", header = FALSE, quote = "",
           col.names = c("chrom", "src", "type", "start", "end", "score", "strand", "phase", "attr"))
g[, chrom := as.character(chrom)]
getattr <- function(x, key) { m <- regmatches(x, regexec(paste0(key, ' "([^"]*)"'), x)); vapply(m, function(v) if (length(v)) v[2] else NA_character_, "") }
msg("feature types: ", paste(names(table(g$type)), table(g$type), sep = "=", collapse = " "))

# ---- homolog table (gene lines)
gn <- g[type == "gene"]
gn[, gene_id := getattr(attr, "gene_id")]
for (k in c("Anolis_Homolog", "Python_Homolog", "Thamnophis_Homolog", "Anolis_Blast_Type", "Python_Blast_Type", "Thamnophis_Blast_Type"))
  set(gn, j = k, value = getattr(gn$attr, k))
sym <- function(h) { s <- sub("^[^_]+_[^_]+_", "", h); s[grepl("^LOC", s) | is.na(s)] <- NA; toupper(s) }
gn[, `:=`(sym_anolis = sym(Anolis_Homolog), sym_python = sym(Python_Homolog), sym_thamnophis = sym(Thamnophis_Homolog))]
# best symbol: prefer an RBB hit with a named symbol, then any named symbol
gn[, best_symbol := fifelse(Anolis_Blast_Type == "RBB" & !is.na(sym_anolis), sym_anolis,
                    fifelse(Python_Blast_Type == "RBB" & !is.na(sym_python), sym_python,
                    fifelse(Thamnophis_Blast_Type == "RBB" & !is.na(sym_thamnophis), sym_thamnophis,
                    fcoalesce(sym_anolis, sym_python, sym_thamnophis))))]
fwrite(gn[, .(gene_id, chrom, start, end, best_symbol, sym_anolis, Anolis_Blast_Type, sym_python, Python_Blast_Type,
              sym_thamnophis, Thamnophis_Blast_Type)], file.path(OUT, "gene_homologs.tsv"), sep = "\t")
msg(sum(!is.na(gn$best_symbol)), " of ", nrow(gn), " genes have a named homolog symbol")

# ---- intervals
cds  <- g[type == "CDS", .(chrom, start, end)]
exon <- g[type %in% c("exon", "five_prime_UTR", "three_prime_UTR", "UTR"), .(chrom, start, end)]
gene <- g[type == "gene", .(chrom, start, end)]
flank <- copy(gene)[, `:=`(start = pmax(1L, start - 5000L), end = end + 5000L)]
for (x in list(cds, exon, gene, flank)) setkey(x, chrom, start, end)

msg("reading tested SNPs")
S <- fread(file.path(D, "rda/rda_snps.tsv.gz"), select = c("snp_id", "chrom", "pos", "q_rda"))
L <- fread(file.path(D, "lfmm/lfmm_snps.tsv.gz"), select = c("snp_id", "lfmm_min_q"))
S <- merge(S, L, by = "snp_id"); rm(L)
S[, chrom := as.character(chrom)]
S[, tier := fifelse(q_rda < 0.05 & lfmm_min_q < 0.05, "tier1", fifelse(q_rda < 0.05, "rda_only", fifelse(lfmm_min_q < 0.05, "lfmm_only", "none")))]
S[, `:=`(start = pos, end = pos)]
hit <- function(iv) { o <- foverlaps(S[, .(snp_id, chrom, start, end)], iv, nomatch = 0L, which = TRUE); u <- logical(nrow(S)); u[unique(o$xid)] <- TRUE; u }
msg("overlapping SNPs with CDS / exons / genes / flanks")
in_cds <- hit(cds); in_exon <- hit(exon); in_gene <- hit(gene); in_flank <- hit(flank)
S[, class := fifelse(in_cds, "CDS", fifelse(in_exon, "UTR / non-coding exon", fifelse(in_gene, "intron",
             fifelse(in_flank, "within 5 kb of gene", "intergenic (>5 kb)"))))]
lev <- c("CDS", "UTR / non-coding exon", "intron", "within 5 kb of gene", "intergenic (>5 kb)")

# ---- summary + Fisher test of each class vs all other tested SNPs
sets <- list(`all tested` = rep(TRUE, nrow(S)), tier1 = S$tier == "tier1", `RDA only` = S$tier == "rda_only",
             `LFMM2 only` = S$tier == "lfmm_only", `not significant` = S$tier == "none")
res <- rbindlist(lapply(names(sets), function(nm) {
  sel <- sets[[nm]]
  rbindlist(lapply(lev, function(cl) {
    a1 <- sum(sel & S$class == cl); b1 <- sum(sel & S$class != cl)
    a0 <- sum(!sel & S$class == cl); b0 <- sum(!sel & S$class != cl)
    ft <- if (nm != "all tested") fisher.test(matrix(c(a1, b1, a0, b0), 2)) else NULL
    data.table(set = nm, class = cl, n = a1, pct = round(100 * a1 / sum(sel), 2),
               odds_ratio_vs_rest = if (is.null(ft)) NA_real_ else round(unname(ft$estimate), 3),
               p = if (is.null(ft)) NA_real_ else signif(ft$p.value, 3))
  }))
}))
fwrite(res, file.path(OUT, "snp_class_summary.tsv"), sep = "\t")
fwrite(S[tier == "tier1", .(snp_id, chrom, pos, class)], file.path(OUT, "tier1_snps_classified.tsv"), sep = "\t")
print(dcast(res, class ~ set, value.var = "pct"))
msg("done; see ", OUT)
