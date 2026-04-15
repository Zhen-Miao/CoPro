# =============================================================
# Colon Day 9: Multi-slide CoPro analysis
# =============================================================
#
# Runs CoPro on a single D9 slide. For multi-slide transfer,
# see 03_D9_cross_slide_transfer.R.
#
# Inputs:
#   - ra_filtered_DSS9_colon_all_slides.rds  (expression matrix)
#   - meta_filtered_DSS9_colon_all_slides.rds (metadata)
#
# Outputs:
#   - CoPro object (RDS)
#   - Normalized correlation PDF
#   - In situ cell score PDFs
#
# Dependencies: CoPro, ggplot2, viridis

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/data"  # UPDATE THIS
OUT_DIR  <- "output/colon_D9"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_pca <- 40
cell_types <- c("Epithelial", "Fibroblast", "Immune")
sigma_choice <- c(0.01, 0.05, 0.1, 0.2)
region_id <- "062221_D9_m3_2_slice_3"

# --- Load data ---
meta <- readRDS(file.path(DATA_DIR, "meta_filtered_DSS9_colon_all_slides.rds"))
ra   <- readRDS(file.path(DATA_DIR, "ra_filtered_DSS9_colon_all_slides.rds"))

# --- Subset to region ---
keep <- meta$Slice_ID == region_id
meta <- meta[keep, ]
ra   <- ra[keep, ]

# --- Preprocess ---
gene_frac <- colMeans(ra > 0)
ra <- ra[, gene_frac >= 0.008]

cell_frac <- rowSums(ra > 0)
keep_cells <- cell_frac >= 20
ra   <- ra[keep_cells, ]
meta <- meta[keep_cells, ]

# Cap expression at 99th percentile, then at 1.0
quant_99 <- apply(ra, 2, quantile, 0.99)
for (g in seq_len(ncol(ra))) {
  ra[, g] <- pmin(ra[, g], quant_99[g])
}
ra <- pmin(ra, 1.0)

location_data <- data.frame(
  x = meta$center_x / 10,
  y = meta$center_y / 10,
  row.names = rownames(ra)
)

# --- CoPro pipeline ---
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
obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500, nCC = 4)
obj <- computeNormalizedCorrelation(obj)
obj <- computeGeneAndCellScores(obj)
obj <- computeRegressionGeneScores(obj)

saveRDS(obj, file.path(OUT_DIR,
  paste0("colon_D9_", region_id, "_CoPro_object.rds")))

# --- Plots ---
sigma_opt <- 0.1

# Normalized correlation
ncorr <- getNormCorr(obj)
pdf(file.path(OUT_DIR, "norm_corr_D9.pdf"), width = 12, height = 8)
print(
  ggplot(ncorr, aes(x = sigmaValues, y = normalizedCorrelation)) +
    geom_point() + geom_line() +
    facet_wrap(~ ct12 + CC_index) +
    xlab("Sigma") + ylab("Normalized Correlation") +
    ggtitle("Colon D9: Normalized correlation") +
    theme_minimal()
)
dev.off()

# In situ cell scores (CC1)
cs <- getCellScoresInSitu(obj, sigmaValueChoice = sigma_opt, ccIndex = 1)
for (ct in cell_types) {
  pdf(file.path(OUT_DIR, paste0("insitu_CC1_", ct, ".pdf")),
      width = 6, height = 5)
  print(
    ggplot(cs[cs$cellType == ct, ]) +
      geom_point(aes(x = x, y = y, color = cellScores), size = 0.8) +
      scale_color_viridis_c() +
      coord_fixed() +
      ggtitle(paste("CC1 -", ct)) +
      theme_classic()
  )
  dev.off()
}
