#!/usr/bin/env Rscript
# =============================================================
# CosMx NSCLC: sigma sensitivity sweep on Lung5_Rep1
# =============================================================
#
# Refit runGeneSpaceCCA on Lung5_Rep1 at sigma in {0.02, 0.05, 0.10, 0.20}
# and check that CC1/CC2 gene weights are rank-stable across the kernel
# scale. Acceptance: cosine similarity of CC weight vectors >= 0.9 between
# adjacent sigmas, per cell type. Anchors how strongly the recovered
# biology depends on the kernel scale.
#
# Distance is computed once. Kernels for all 4 sigmas are precomputed
# in a single call (peak RSS will be roughly 4x kernel storage, ~40 GB
# for this slide; well under the box's 1.5 TB).
#
# Usage:
#   Rscript scripts/cosmx_nsclc/03_sigma_sensitivity.R [slide]
# Default slide: Lung5_Rep1.

suppressPackageStartupMessages({
  library(CoPro)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
SLIDE <- if (length(args) >= 1L) args[[1]] else "Lung5_Rep1"

DATA_ROOT <- "/home/zmiao2/spatial/cosmx_nsclc"
OUT_DIR   <- file.path(DATA_ROOT, "copro_results", paste0(SLIDE, "_sigma_sweep"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CT_OF_INTEREST <- c("tumor 5", "fibroblast", "macrophage")
SIGMAS         <- c(0.02, 0.05, 0.10, 0.20)
N_CC           <- 2
PX_TO_MM       <- 0.18 / 1000

log_step <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
}
peak_rss_gb <- function() {
  s <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
  hit <- grep("^VmHWM:", s, value = TRUE)
  if (length(hit) == 0) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", hit[1])) / 1024 / 1024
}

# ---- load data --------------------------------------------------------------
log_step(sprintf("Loading slide %s", SLIDE))
flat_dir   <- file.path(DATA_ROOT, SLIDE, paste0(SLIDE, "-Flat_files_and_images"))
expr_path  <- file.path(flat_dir, paste0(SLIDE, "_exprMat_file.csv"))
meta_path  <- file.path(flat_dir, paste0(SLIDE, "_metadata_file.csv"))
ctype_path <- file.path(DATA_ROOT, "cell_type_labels", paste0(SLIDE, "_cell_types.csv"))

expr  <- read.csv(expr_path,  check.names = FALSE)
meta  <- read.csv(meta_path,  check.names = FALSE)
ctype <- read.csv(ctype_path, check.names = FALSE)

expr <- expr[expr$cell_ID > 0, , drop = FALSE]
meta <- meta[as.integer(meta$cell_ID) > 0, , drop = FALSE]
ctype$cell_int <- as.integer(sub(".*_", "", ctype$cell_ID))

mk_uid <- function(fov, cid) paste0("fov", fov, "_c", cid)
expr$.uid  <- mk_uid(expr$fov,  expr$cell_ID)
meta$.uid  <- mk_uid(meta$fov,  meta$cell_ID)
ctype$.uid <- mk_uid(ctype$fov, ctype$cell_int)

ct_keep <- ctype[ctype$cell_type %in% CT_OF_INTEREST, , drop = FALSE]
common  <- Reduce(intersect, list(expr$.uid, meta$.uid, ct_keep$.uid))
expr    <- expr[match(common, expr$.uid), , drop = FALSE]
meta    <- meta[match(common, meta$.uid), , drop = FALSE]
ct_keep <- ct_keep[match(common, ct_keep$.uid), , drop = FALSE]

non_gene  <- c("fov", "cell_ID", ".uid")
gene_cols <- setdiff(colnames(expr), non_gene)
neg_idx <- grepl("^NegPrb",       gene_cols, ignore.case = TRUE) |
           grepl("^SystemControl", gene_cols, ignore.case = TRUE)
gene_cols <- gene_cols[!neg_idx]

mat_int <- as.matrix(expr[, gene_cols, drop = FALSE])
storage.mode(mat_int) <- "double"
rownames(mat_int) <- common

total_counts <- rowSums(mat_int)
keep_qc <- total_counts > 0
if (any(!keep_qc)) {
  mat_int      <- mat_int[keep_qc, , drop = FALSE]
  meta         <- meta[keep_qc, , drop = FALSE]
  ct_keep      <- ct_keep[keep_qc, , drop = FALSE]
  total_counts <- total_counts[keep_qc]
  common       <- common[keep_qc]
}

norm_mat <- log1p(sweep(mat_int, 1, total_counts, FUN = "/") * 10000)

loc <- data.frame(
  x = meta$CenterX_global_px * PX_TO_MM,
  y = meta$CenterY_global_px * PX_TO_MM,
  row.names = common
)
meta_df <- data.frame(
  fov          = meta$fov,
  total_counts = total_counts,
  row.names    = common
)
ct_vec <- as.character(ct_keep$cell_type)

log_step(sprintf("  %d cells; cell types: %s", length(ct_vec),
                 paste(sprintf("%s=%d", names(table(ct_vec)), as.integer(table(ct_vec))),
                       collapse = ", ")))

# ---- pipeline ---------------------------------------------------------------
log_step("Building CoProMulti")
suppressWarnings(
  obj <- newCoProMulti(
    normalizedData = norm_mat,
    locationData   = loc,
    metaData       = meta_df,
    cellTypes      = ct_vec,
    slideID        = rep(SLIDE, length(ct_vec))
  )
)
obj <- subsetData(obj, cellTypesOfInterest = CT_OF_INTEREST)

log_step("computeDistance (Euclidean2D)")
t1 <- Sys.time()
obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
log_step(sprintf("  %.1fs (peak RSS %.2f GB)",
                 as.numeric(Sys.time() - t1, units = "secs"), peak_rss_gb()))

log_step(sprintf("computeKernelMatrix for sigmas: %s",
                 paste(SIGMAS, collapse = ", ")))
t1 <- Sys.time()
obj <- computeKernelMatrix(obj, sigmaValues = SIGMAS, verbose = FALSE)
log_step(sprintf("  %.1fs (peak RSS %.2f GB)",
                 as.numeric(Sys.time() - t1, units = "secs"), peak_rss_gb()))

# Free distances now that all kernels are built (they aren't used by CCA)
obj@distances <- list()
gc(verbose = FALSE)
log_step(sprintf("  freed distances; peak RSS %.2f GB", peak_rss_gb()))

# Run runGeneSpaceCCA per sigma; geneScores keys are sigma-stamped, so we
# can keep accumulating into the same object.
for (sg in SIGMAS) {
  log_step(sprintf("runGeneSpaceCCA sigma=%g, nCC=%d", sg, N_CC))
  t1 <- Sys.time()
  obj <- runGeneSpaceCCA(obj, sigma = sg, nCC = N_CC, verbose = TRUE)
  log_step(sprintf("  %.1fs (peak RSS %.2f GB)",
                   as.numeric(Sys.time() - t1, units = "secs"), peak_rss_gb()))
}

# ---- extract per-sigma gene weights and cosine matrix -----------------------
log_step("Extracting per-sigma gene weights")

# geneScores keys: "geneScores|sigma<value>|<celltype>" (no underscore
# between 'sigma' and the value).
parse_key <- function(k) {
  parts <- strsplit(k, "\\|", perl = TRUE)[[1]]
  list(sigma = as.numeric(sub("^sigma", "", parts[2])),
       cell_type = parts[3])
}

# Build long-form weights table for all (sigma, cell_type, cc, gene)
gs_long <- do.call(rbind, lapply(names(obj@geneScores), function(k) {
  pk <- parse_key(k)
  m <- obj@geneScores[[k]]
  do.call(rbind, lapply(seq_len(N_CC), function(cc) {
    data.frame(
      sigma = pk$sigma,
      cell_type = pk$cell_type,
      cc = cc,
      gene = rownames(m),
      weight = m[, paste0("CC_", cc)],
      row.names = NULL
    )
  }))
}))

# Top-30 per (sigma, cell_type, cc) summary (handy for inspection)
top_n <- 30
gs_top <- do.call(rbind, lapply(split(gs_long,
                                       list(gs_long$sigma, gs_long$cell_type, gs_long$cc),
                                       drop = TRUE),
  function(df) {
    df <- df[order(-abs(df$weight)), , drop = FALSE]
    head(df, top_n)
  }))
write.csv(gs_top, file.path(OUT_DIR, sprintf("%s_topgene_weights_per_sigma.csv", SLIDE)),
          row.names = FALSE)

# Full long table (all genes; needed for cosine sims)
write.csv(gs_long, file.path(OUT_DIR, sprintf("%s_full_weights_per_sigma.csv", SLIDE)),
          row.names = FALSE)

# Cosine similarity of CC weight vectors between sigma pairs, per cell type.
# Sign of the eigenvector is arbitrary, so use abs(cosine).
cell_types <- sort(unique(gs_long$cell_type))
cosine_rows <- list()
for (ct in cell_types) {
  for (cc in seq_len(N_CC)) {
    # gene x sigma weight matrix
    sub <- gs_long[gs_long$cell_type == ct & gs_long$cc == cc, ]
    W <- reshape(sub[, c("gene", "sigma", "weight")],
                 idvar = "gene", timevar = "sigma", direction = "wide")
    rownames(W) <- W$gene; W$gene <- NULL
    colnames(W) <- sub("^weight\\.", "", colnames(W))
    sigma_cols <- as.numeric(colnames(W))
    ord <- order(sigma_cols)
    W <- as.matrix(W[, ord, drop = FALSE]); sigma_cols <- sigma_cols[ord]

    # all pairs cosine (abs)
    norms <- sqrt(colSums(W^2))
    M <- (t(W) %*% W) / outer(norms, norms)
    M <- abs(M)

    # also rank-based (Spearman) cosine of |weight|, since the issue
    # phrases the criterion as "rank-stable"
    R <- apply(W, 2, function(v) rank(-abs(v)))
    R <- scale(R, center = TRUE, scale = FALSE)
    rnorms <- sqrt(colSums(R^2))
    M_rank <- (t(R) %*% R) / outer(rnorms, rnorms)

    for (i in seq_along(sigma_cols)) {
      for (j in seq_along(sigma_cols)) {
        cosine_rows[[length(cosine_rows) + 1L]] <- data.frame(
          cell_type = ct, cc = cc,
          sigma_i = sigma_cols[i], sigma_j = sigma_cols[j],
          cosine_weights = M[i, j],
          cosine_ranks   = M_rank[i, j],
          row.names = NULL
        )
      }
    }
  }
}
cosine_df <- do.call(rbind, cosine_rows)
write.csv(cosine_df, file.path(OUT_DIR, sprintf("%s_sigma_cosine.csv", SLIDE)),
          row.names = FALSE)

# Adjacent-sigma cosine summary (the acceptance check)
ord_sigma <- sort(SIGMAS)
adjacent_rows <- list()
for (ct in cell_types) {
  for (cc in seq_len(N_CC)) {
    for (k in seq_len(length(ord_sigma) - 1L)) {
      mask <- cosine_df$cell_type == ct & cosine_df$cc == cc &
        cosine_df$sigma_i == ord_sigma[k] &
        cosine_df$sigma_j == ord_sigma[k + 1L]
      adjacent_rows[[length(adjacent_rows) + 1L]] <-
        cosine_df[mask, , drop = FALSE]
    }
  }
}
adjacent_df <- do.call(rbind, adjacent_rows)
adjacent_df$pass_0p9 <- adjacent_df$cosine_weights >= 0.9
write.csv(adjacent_df, file.path(OUT_DIR, sprintf("%s_sigma_adjacent.csv", SLIDE)),
          row.names = FALSE)

log_step("Adjacent-sigma cosine (acceptance: >= 0.9):")
print(adjacent_df, row.names = FALSE)
log_step(sprintf("  passes: %d / %d", sum(adjacent_df$pass_0p9), nrow(adjacent_df)))

# Save lean object
obj_lean <- obj
obj_lean@kernelMatrices <- list()
saveRDS(obj_lean, file.path(OUT_DIR, sprintf("%s_sigma_sweep.rds", SLIDE)))

log_step(sprintf("DONE. peak RSS %.2f GB. Outputs in %s",
                 peak_rss_gb(), OUT_DIR))
