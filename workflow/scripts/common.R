# Shared helpers for the GEA pipeline scripts.
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

# --key value | --key=value | --key v1 v2 v3 (multi-value) -> named list
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(); key <- NULL
  for (a in args) {
    if (startsWith(a, "--")) {
      kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
      key <- gsub("-", "_", kv[1])
      out[[key]] <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else character(0)
    } else if (!is.null(key)) {
      out[[key]] <- c(out[[key]], a)
    }
  }
  out
}

need <- function(a, keys) {
  miss <- setdiff(keys, names(a))
  if (length(miss)) stop("missing argument(s): ", paste0("--", miss, collapse = ", "))
}

read_config <- function(path) yaml::read_yaml(path)

msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# "auto" / NULL / "" -> default ; otherwise numeric
as_auto <- function(x, default) {
  if (is.null(x) || identical(x, "") || identical(tolower(as.character(x)), "auto")) default else as.numeric(x)
}

# Read a plink2 --export A (.raw) file into an integer matrix (samples x SNPs)
# with rownames = IID and colnames = variant ID (allele suffix stripped).
read_raw <- function(path) {
  dt <- fread(path, header = TRUE, showProgress = FALSE)
  iid <- as.character(dt$IID)
  meta <- intersect(names(dt), c("#FID", "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
  dt[, (meta) := NULL]
  snp <- sub("_[^_]+$", "", names(dt))            # strip _<counted allele>
  m <- as.matrix(dt)
  storage.mode(m) <- "integer"
  dimnames(m) <- list(iid, snp)
  m
}

# Read a per-sample table keeping sample_id as character (numeric-looking IDs!)
read_samples_tab <- function(path, id = "sample_id") {
  fread(path, colClasses = setNames(list(id), "character"))
}

# Split "chrom:pos" IDs
split_id <- function(id) {
  pos <- as.integer(sub("^.*:", "", id))
  chrom <- sub(":[^:]*$", "", id)
  data.table(snp_id = id, chrom = chrom, pos = pos)
}

# q-values with a BH fallback (qvalue() fails when p-values are unusual)
safe_q <- function(p) {
  ok <- !is.na(p)
  q <- rep(NA_real_, length(p))
  if (!any(ok)) return(q)
  q[ok] <- tryCatch(qvalue::qvalue(p[ok])$qvalues,
                    error = function(e) { msg("qvalue failed (", conditionMessage(e), "); using BH"); p.adjust(p[ok], "BH") })
  q
}

# Capblancq & Forester (2021) rdadapt: Mahalanobis distance on scaled SNP
# loadings of the first K constrained axes -> chi2 p-values -> q-values
rdadapt <- function(rda_obj, K) {
  zscores <- rda_obj$CCA$v[, seq_len(K), drop = FALSE]
  resscale <- apply(zscores, 2, scale)
  resmaha <- robust::covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha) / qchisq(0.5, df = K)
  p <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  data.table(p_rda = p, q_rda = safe_q(p))
}

# Parse genes from a GFF3 *or* GTF annotation (format detected automatically).
#  GFF3 attributes look like  ID=gene1;Name=ABC1
#  GTF  attributes look like  gene_id "gene1"; gene_name "ABC1";
# If a GTF has no "gene" lines (some only list transcripts/exons), each gene's
# span is built from the min start / max end of all its features.
read_gff_genes <- function(path, feature = "gene", id_attr = "ID") {
  g <- fread(cmd = paste("gzip -cdf", shQuote(path), "| grep -v '^#'"), sep = "\t", header = FALSE,
             quote = "", col.names = c("chrom", "src", "type", "start", "end", "score",
                                       "strand", "phase", "attr"))
  g[, chrom := as.character(chrom)]
  is_gtf <- any(grepl('gene_id "', head(g$attr, 2000), fixed = TRUE))
  get_attr <- function(attr, key) {
    rx <- if (is_gtf) paste0('(^|;)[[:space:]]*', key, ' "([^"]*)"') else paste0("(^|;)", key, "=([^;]+)")
    m <- regmatches(attr, regexec(rx, attr))
    vapply(m, function(x) if (length(x)) x[3] else NA_character_, "")
  }
  if (is_gtf && identical(id_attr, "ID")) id_attr <- "gene_id"   # GFF default -> GTF equivalent
  msg("annotation format: ", if (is_gtf) "GTF" else "GFF3", "; gene ID attribute: ", id_attr)
  feat <- g[type == feature]
  if (!nrow(feat)) {
    if (!is_gtf) stop("no '", feature, "' features in ", path, "; set gene_feature in the config")
    msg("no '", feature, "' lines in GTF; building gene spans from all features per ", id_attr)
    feat <- copy(g)
  }
  feat[, gene_id := get_attr(attr, id_attr)]
  name_keys <- if (is_gtf) c("gene_name", "gene") else c("Name", "gene")
  feat[, name := NA_character_]
  for (k in name_keys) feat[is.na(name), name := get_attr(attr, k)]
  feat <- feat[!is.na(gene_id)]
  if (!nrow(feat)) stop("no gene IDs found using attribute '", id_attr, "'; check gene_id_attr in the config")
  out <- feat[, .(start = min(as.integer(start)), end = max(as.integer(end)),
                  strand = strand[1], name = name[1]), by = .(chrom, gene_id)]
  msg(nrow(out), " genes read from annotation")
  out[, .(chrom, start, end, strand, gene_id, name)]
}

# ---------------------------------------------------------------------
# Figure palette (colourblind-validated categorical order; one-hue blue ramp)
# ---------------------------------------------------------------------
PAL <- list(
  cat = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"),
  other = "#9c9a92",
  bg = "#b5b3ad",
  ink = "#0b0b0b", ink2 = "#52514e", grid = "#e4e3df",
  seq = c("#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"),
  tier = c(none = "#b5b3ad", rda_only = "#2a78d6", lfmm_only = "#1baf7a", tier1 = "#eb6834")
)
seq_ramp <- function(n) grDevices::colorRampPalette(PAL$seq)(n)
# fixed colour per predictor name so the same variable is the same colour everywhere
env_colours <- function(vars) {
  vars <- unique(vars)
  cols <- setNames(rep(PAL$other, length(vars)), vars)
  k <- min(length(vars), length(PAL$cat))
  cols[seq_len(k)] <- PAL$cat[seq_len(k)]
  cols
}
# "bio1(5);bio12(1)" -> named integer c(bio1 = 5, bio12 = 1)
parse_env_counts <- function(s) {
  s <- s[!is.na(s) & nzchar(s)]
  if (!length(s)) return(integer(0))
  tok <- unlist(strsplit(unlist(strsplit(s, "|", fixed = TRUE)), ";", fixed = TRUE))
  tok <- tok[nzchar(tok)]
  nm <- sub("\\(.*$", "", tok)
  n <- suppressWarnings(as.integer(sub("^.*\\((\\d+)\\)$", "\\1", tok)))
  n[is.na(n)] <- 1L
  tapply(n, nm, sum)
}
top_env_name <- function(s) {
  x <- parse_env_counts(s)
  if (!length(x)) NA_character_ else names(x)[which.max(x)]
}
base_theme <- function() par(family = "sans", col.axis = PAL$ink2, col.lab = PAL$ink2, fg = PAL$ink2,
                             col.main = PAL$ink, font.main = 1, las = 1, tcl = -0.25, mgp = c(2.2, 0.5, 0))
