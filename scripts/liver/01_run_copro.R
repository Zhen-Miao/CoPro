# =============================================================
# Liver: CoPro on Visium spatial transcriptomics data
# =============================================================
#
# Runs CoPro on liver Visium data to detect hepatocyte zonation
# patterns and endothelial co-progression along the
# pericentral-periportal axis.
#
# Inputs:
#   - liver_*_meta_cleaned.rds  (metadata with spatial coords)
#   - liver_*_dat_norm_cleaned.rds (normalized expression)
#
# Outputs:
#   - CoPro object (RDS)
#   - Normalized correlation and in situ plots
#
# Dependencies: CoPro, ggplot2

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/liver_data"  # UPDATE THIS
OUT_DIR  <- "output/liver"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_pca <- 40
cell_types_single <- "all_types"
sigma_choice <- c(0.02, 0.05, 0.1)

# --- Load data ---
meta <- readRDS(file.path(DATA_DIR, "liver_344_v3_meta_cleaned.rds"))
dat  <- readRDS(file.path(DATA_DIR, "liver_344_v3_dat_norm_cleaned.rds"))

# --- Spatial filtering ---
keep <- meta$x < 600 & meta$y < 500 & meta$y > 350
meta <- meta[keep, ]
dat  <- dat[keep, ]

# Prepare location data
location_data <- data.frame(
  x = meta$x,
  y = meta$y,
  row.names = rownames(dat)
)

# Create pseudo-slide division (horizontal split)
slide_id <- ifelse(meta$x < 450, "slide_left", "slide_right")

# --- CoPro pipeline (multi-slide, single cell type) ---
obj <- newCoProMulti(
  normalizedData = dat,
  locationData = location_data,
  metaData = meta,
  cellTypes = rep(cell_types_single, nrow(dat)),
  slideID = slide_id
)
obj <- subsetData(obj, cellTypesOfInterest = cell_types_single)
obj <- computePCA(obj, nPCA = n_pca, center = TRUE, scale. = TRUE)
obj <- computeDistance(obj, distType = "Euclidean2D",
                       normalizeDistance = FALSE)
obj <- computeKernelMatrix(obj, sigmaValues = sigma_choice,
                            upperQuantile = 0.85,
                            normalizeKernel = FALSE)
obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500, nCC = 4)
obj <- computeNormalizedCorrelation(obj, tol = 1e-3,
                                    calculationMode = "perSlide")
obj <- computeGeneAndCellScores(obj)
obj <- computeRegressionGeneScores(obj)

saveRDS(obj, file.path(OUT_DIR, "liver_CoPro_object.rds"))

# --- Plots ---
sigma_opt <- 0.05

# In situ cell scores
for (cc_ind in 1:2) {
  cs <- getCellScoresInSitu(obj, sigmaValueChoice = sigma_opt,
                             ccIndex = cc_ind)
  pdf(file.path(OUT_DIR, paste0("insitu_CC", cc_ind, ".pdf")),
      width = 8, height = 5)
  print(
    ggplot(cs) +
      geom_point(aes(x = x, y = y, color = cellScores), size = 0.8) +
      scale_color_viridis_c() +
      coord_fixed() +
      ggtitle(paste("Liver CC", cc_ind)) +
      theme_classic()
  )
  dev.off()
}
