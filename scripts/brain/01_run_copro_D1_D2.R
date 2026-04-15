# =============================================================
# Brain: D1/D2 neuron co-progression (MERFISH)
# =============================================================
#
# Detects co-progression between D1 and D2 GABAergic neurons
# in the striatum using brain MERFISH data.
#
# Reference: Zhang et al. Nature 624, 343-354 (2023)
#
# Inputs:
#   - Zhuang_ABCA_1.054_subset_data.rds     (expression matrix)
#   - Zhuang_ABCA_1.054_subset_metadata.rds  (metadata)
#
# Outputs:
#   - CoPro object (RDS)
#   - In situ cell scores, correlation plots, gene weight bar plots
#
# Dependencies: CoPro, ggplot2, viridis

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/brain_data"  # UPDATE THIS
OUT_DIR  <- "output/brain"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_pca <- 40
cell_types <- c("061 STR D1 Gaba", "062 STR D2 Gaba")
sigma_choice <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
sigma_opt <- 0.1

# --- Load and subset data ---
ra   <- readRDS(file.path(DATA_DIR, "Zhuang_ABCA_1.054_subset_data.rds"))
meta <- readRDS(file.path(DATA_DIR, "Zhuang_ABCA_1.054_subset_metadata.rds"))

keep <- meta$subclass %in% cell_types
ra   <- ra[keep, ]
meta <- meta[keep, ]
rownames(meta) <- meta$cell_label
meta$x <- as.numeric(meta$x)
meta$y <- as.numeric(meta$y)

location_data <- data.frame(
  x = meta$x,
  y = meta$y,
  row.names = rownames(ra)
)

# --- CoPro pipeline ---
obj <- newCoProSingle(
  normalizedData = ra,
  locationData = location_data,
  metaData = meta,
  cellTypes = meta$subclass
)
obj <- subsetData(obj, cellTypesOfInterest = cell_types)
obj <- computePCA(obj, nPCA = n_pca, center = TRUE, scale. = TRUE)
obj <- computeDistance(obj, distType = "Euclidean2D",
                       normalizeDistance = FALSE)
obj <- computeKernelMatrix(obj, sigmaValues = sigma_choice)
obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500)
obj <- computeNormalizedCorrelation(obj)
obj <- computeGeneAndCellScores(obj)
obj <- computeRegressionGeneScores(obj)

saveRDS(obj, file.path(OUT_DIR, "brain_D1_D2_CoPro_object.rds"))

# --- Plots ---
# Normalized correlation
ncorr <- getNormCorr(obj)
pdf(file.path(OUT_DIR, "norm_corr_brain.pdf"), width = 10, height = 6)
print(
  ggplot(ncorr, aes(x = sigmaValues, y = normalizedCorrelation)) +
    geom_point() + geom_line() +
    facet_wrap(~ ct12 + CC_index) +
    xlab("Sigma") + ylab("Normalized Correlation") +
    ggtitle("Brain MERFISH: D1-D2 co-progression") +
    theme_minimal()
)
dev.off()

# Cross-type correlation
df_corr <- getCorrTwoTypes(obj,
  sigmaValueChoice = sigma_opt,
  cellTypeA = cell_types[1],
  cellTypeB = cell_types[2]
)
pdf(file.path(OUT_DIR, "cross_corr_D1_D2.pdf"), width = 5, height = 5)
print(
  ggplot(df_corr) +
    geom_point(aes(x = AK, y = B), size = 0.5, alpha = 0.3) +
    xlab("D1 score (spatially smoothed)") +
    ylab("D2 score") +
    ggtitle("D1-D2 cross-type correlation") +
    theme_minimal()
)
dev.off()

# In situ
cs <- getCellScoresInSitu(obj, sigmaValueChoice = sigma_opt)
pdf(file.path(OUT_DIR, "insitu_D1_D2.pdf"), width = 7, height = 5)
print(
  ggplot(cs) +
    geom_point(aes(x = x, y = y, color = cellScores), size = 0.5) +
    scale_color_viridis_c() +
    facet_wrap(~ cellType) +
    coord_fixed() +
    ggtitle("CC1 cell scores in situ") +
    theme_minimal()
)
dev.off()
