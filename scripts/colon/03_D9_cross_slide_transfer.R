# =============================================================
# Colon Day 9: Cross-slide score transfer
# =============================================================
#
# Transfers learned gene weights from a reference D9 slide to
# target slides to assess co-progression consistency.
#
# Inputs:
#   - Pre-computed CoPro objects per region (from 02_run_copro_D9.R)
#   - ra_filtered_DSS9_colon_all_slides.rds
#   - meta_filtered_DSS9_colon_all_slides.rds
#
# Outputs:
#   - Transferred cell scores (RDS per target)
#   - Transferred normalized correlation values
#   - Comparison PDF
#
# Dependencies: CoPro, ggplot2

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/data"  # UPDATE THIS
OUT_DIR  <- "output/colon_D9_transfer"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_pca <- 40
cell_types <- c("Epithelial", "Fibroblast", "Immune")
sigma_choice <- c(0.01, 0.05, 0.1, 0.2)
sigma_opt <- 0.1

ref_region <- "062221_D9_m3_2_slice_3"
tar_regions <- c("062221_D9_m3_2_slice_1", "100221_D9_m5_2_slice_2")

# --- Helper: create CoPro object for a region ---
make_region_obj <- function(ra_all, meta_all, region_id) {
  keep <- meta_all$Slice_ID == region_id
  ra <- ra_all[keep, ]
  meta <- meta_all[keep, ]

  # Preprocess
  gene_frac <- colMeans(ra > 0)
  ra <- ra[, gene_frac >= 0.008]
  cell_frac <- rowSums(ra > 0)
  keep_cells <- cell_frac >= 20
  ra <- ra[keep_cells, ]
  meta <- meta[keep_cells, ]
  quant_99 <- apply(ra, 2, quantile, 0.99)
  for (g in seq_len(ncol(ra))) ra[, g] <- pmin(ra[, g], quant_99[g])
  ra <- pmin(ra, 1.0)

  location_data <- data.frame(
    x = meta$center_x / 10,
    y = meta$center_y / 10,
    row.names = rownames(ra)
  )

  obj <- newCoProSingle(
    normalizedData = ra,
    locationData = location_data,
    metaData = meta,
    cellTypes = meta$Celltype
  )
  obj <- subsetData(obj, cellTypesOfInterest = cell_types)
  obj <- computePCA(obj, nPCA = n_pca, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D")
  obj <- computeKernelMatrix(obj, sigmaValues = sigma_choice)
  return(obj)
}

# --- Load data ---
meta_all <- readRDS(file.path(DATA_DIR,
  "meta_filtered_DSS9_colon_all_slides.rds"))
ra_all   <- readRDS(file.path(DATA_DIR,
  "ra_filtered_DSS9_colon_all_slides.rds"))

# --- Reference: run full pipeline ---
ref_obj <- make_region_obj(ra_all, meta_all, ref_region)
ref_obj <- runSkrCCA(ref_obj, scalePCs = TRUE, maxIter = 500, nCC = 2)
ref_obj <- computeNormalizedCorrelation(ref_obj)
ref_obj <- computeGeneAndCellScores(ref_obj)
ref_obj <- computeRegressionGeneScores(ref_obj)

# --- Transfer to each target ---
for (tar_region in tar_regions) {
  message("Transferring to: ", tar_region)
  tar_obj <- make_region_obj(ra_all, meta_all, tar_region)

  # Transfer using regression gene weights
  tar_scores <- getTransferCellScores(
    ref_obj = ref_obj,
    tar_obj = tar_obj,
    sigma_choice = sigma_opt,
    gene_score_type = "regression"
  )

  saveRDS(tar_scores, file.path(OUT_DIR,
    paste0("transferred_scores_", tar_region, ".rds")))

  # Transfer normalized correlation
  tar_ncorr <- getTransferNormCorr(
    ref_obj = ref_obj,
    tar_obj = tar_obj,
    sigma_choice = sigma_opt,
    gene_score_type = "regression"
  )

  message("  Transferred norm. corr.: ")
  print(tar_ncorr)
}
