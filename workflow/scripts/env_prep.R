#!/usr/bin/env Rscript
# Build the environmental predictor matrix for one species:
#  * keep only genotyped individuals (order matches the genotype matrix)
#  * NLCD point classes -> collapsed habitat dummies (or use buffer proportions)
#  * drop near-constant predictors
#  * correlation pruning (|r| >= threshold) respecting a priority list
#  * iterative VIF pruning
#  * scale + write selected predictors and geographic coordinates
#
# Rscript env_prep.R --config config.yaml --samples samples.csv \
#   --genotyped results/sp/structure/pcs.tsv --outdir results/sp/env
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "common.R"))

a <- parse_args(); need(a, c("config", "samples", "genotyped", "outdir"))
cfg <- read_config(a$config)$env
dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

S <- fread(a$samples, colClasses = setNames(list(cfg$sample_id_col), "character"))
id_col <- cfg$sample_id_col
setnames(S, id_col, "sample_id")
S[, sample_id := as.character(sample_id)]
geno_ids <- read_samples_tab(a$genotyped)$sample_id
miss <- setdiff(geno_ids, S$sample_id)
if (length(miss)) stop(length(miss), " genotyped samples have no environment row, e.g. ", paste(head(miss), collapse = ", "))
S <- S[match(geno_ids, sample_id)]
msg(nrow(S), " genotyped individuals matched to environment table")

# ---------------- candidate predictors
cand <- unlist(cfg$candidate_vars)
absent <- setdiff(cand, names(S))
if (length(absent)) { msg("WARNING: not in table, skipped: ", paste(absent, collapse = ", ")); cand <- intersect(cand, names(S)) }
E <- S[, ..cand]

nl <- cfg$nlcd
if (!is.null(nl$column) && nzchar(nl$column)) {
  if (identical(nl$mode, "proportion")) {
    pc <- unlist(nl$prop_cols)
    E <- cbind(E, S[, ..pc])
  } else {
    map <- unlist(nl$collapse)
    cls <- as.character(S[[nl$column]])
    grp <- unname(map[cls])
    if (anyNA(grp)) msg("WARNING: NLCD codes without a collapse entry: ", paste(unique(cls[is.na(grp)]), collapse = ", "))
    grp[is.na(grp)] <- "other"
    lv <- sort(unique(grp))
    D <- sapply(lv, function(l) as.numeric(grp == l))
    D <- matrix(D, nrow = nrow(S), dimnames = list(NULL, paste0("lc_", lv)))
    freq <- colMeans(D)
    keep <- freq >= as.numeric(nl$min_freq) & freq <= 1 - as.numeric(nl$min_freq)
    msg("NLCD groups (frequency): ", paste(sprintf("%s=%.2f", colnames(D), freq), collapse = " "))
    E <- cbind(E, as.data.table(D[, keep, drop = FALSE]))
  }
}

# ---------------- remove incomplete / near-constant
na_rows <- !complete.cases(E)
if (any(na_rows)) stop(sum(na_rows), " individuals have missing environmental values; fix or drop them in the samples table")
sdv <- sapply(E, sd)
E <- E[, names(sdv)[sdv > 1e-8], with = FALSE]

log_rows <- list()
add_log <- function(var, step, detail) log_rows[[length(log_rows) + 1]] <<- data.table(variable = var, step = step, detail = detail)

# ---------------- correlation pruning with priority
prio <- intersect(unlist(cfg$priority), names(E))
order_vars <- c(prio, setdiff(names(E), prio))
R <- abs(cor(E))
sel <- character(0)
for (v in order_vars) {
  clash <- sel[R[v, sel] >= cfg$cor_threshold]
  if (length(clash)) add_log(v, "correlation", sprintf("|r|=%.2f with %s", max(R[v, clash]), clash[which.max(R[v, clash])]))
  else sel <- c(sel, v)
}

# ---------------- VIF pruning (drop highest VIF; among ties in the
# offending set, drop the lowest-priority variable first)
vif <- function(X) { C <- cor(X); if (ncol(C) < 2) return(setNames(1, colnames(X))); diag(solve(C)) }
repeat {
  if (length(sel) < 2) break
  vv <- tryCatch(vif(as.matrix(E[, ..sel])), error = function(e) setNames(rep(Inf, length(sel)), sel))
  if (max(vv) < cfg$vif_threshold) break
  over <- names(vv)[vv >= cfg$vif_threshold]
  rank <- match(over, order_vars)
  drop <- over[which.max(rank)]              # lowest priority among offenders
  add_log(drop, "VIF", sprintf("VIF=%.2f", vv[drop]))
  sel <- setdiff(sel, drop)
}
vv <- if (length(sel) > 1) vif(as.matrix(E[, ..sel])) else setNames(1, sel)
for (v in sel) add_log(v, "selected", sprintf("VIF=%.2f", vv[v]))
msg("selected predictors: ", paste(sel, collapse = ", "))

# ---------------- write
Xs <- scale(as.matrix(E[, ..sel]))
fwrite(data.table(sample_id = S$sample_id, Xs), file.path(a$outdir, "env_selected.tsv"), sep = "\t")
geo <- unlist(cfg$geo_vars)
G <- scale(as.matrix(S[, ..geo]))
fwrite(data.table(sample_id = S$sample_id, G), file.path(a$outdir, "geo.tsv"), sep = "\t")
fwrite(rbindlist(log_rows), file.path(a$outdir, "pruning_log.tsv"), sep = "\t")
fwrite(data.table(S[, .(sample_id)], E), file.path(a$outdir, "env_all_candidates_raw.tsv"), sep = "\t")

pdf(file.path(a$outdir, "env_correlations.pdf"), width = 10, height = 9)
Rall <- cor(E)
o <- hclust(as.dist(1 - abs(Rall)))$order
image(1:ncol(Rall), 1:ncol(Rall), Rall[o, o], zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
      col = hcl.colors(41, "Blue-Red 3"), main = "Candidate predictor correlations (* = selected)")
labs <- colnames(Rall)[o]
labs <- ifelse(labs %in% sel, paste0("*", labs), labs)
axis(1, 1:ncol(Rall), labs, las = 2, cex.axis = 0.7)
axis(2, 1:ncol(Rall), labs, las = 2, cex.axis = 0.7)
if (length(sel) > 1) pairs(E[, ..sel], pch = 19, cex = 0.4, main = "Selected predictors")
dev.off()
msg("done")
