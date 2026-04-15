# =============================================================
# Simulation: Run CoPro on simulated data
# =============================================================
#
# Runs the CoPro pipeline on each simulated dataset and computes
# Pearson correlations between CoPro cell scores and ground truth.
#
# Inputs:
#   - Parquet files from 01_generate_data.R
#
# Outputs:
#   - copro_results_[scenario].csv with correlation values
#
# Dependencies: CoPro, arrow

library(CoPro)
library(arrow)

# --- Configuration ---
DATA_DIR <- "path/to/simulation_output"  # UPDATE THIS
OUT_DIR  <- "output/simulation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_pca <- 25
sigma_choice <- c(0.05, 0.1, 0.2, 0.4, 0.8, 1, 1.5)

# --- Process each simulation run ---
# For each parquet file pair (expression + metadata):
#   1. Load expression and metadata
#   2. Create CoPro object
#   3. Run pipeline
#   4. Extract cell scores
#   5. Compute correlation with ground truth

# Example for a single configuration:
run_copro_on_sim <- function(expr_path, meta_path) {
  expr <- as.data.frame(read_parquet(expr_path))
  meta <- as.data.frame(read_parquet(meta_path))

  # Ensure matching
  rownames(expr) <- meta$cell_id
  location_data <- data.frame(
    x = meta$x, y = meta$y,
    row.names = meta$cell_id
  )

  obj <- newCoProSingle(
    normalizedData = as.matrix(expr),
    locationData = location_data,
    metaData = meta,
    cellTypes = meta$cell_type
  )
  cell_types <- unique(meta$cell_type)
  obj <- subsetData(obj, cellTypesOfInterest = cell_types)

  obj <- computePCA(obj, nPCA = n_pca, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = sigma_choice,
                              rowNormalizeKernel = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500)
  obj <- computeNormalizedCorrelation(obj)
  obj <- computeBidirCorrelation(obj, normalize_K = "row_or_col",
                                  filter_kernel = FALSE)
  obj <- computeGeneAndCellScores(obj)

  # Extract scores and compute correlation with ground truth
  ncorr <- getNormCorr(obj)
  best_sigma <- ncorr$sigmaValues[which.max(ncorr$normalizedCorrelation)]

  return(list(obj = obj, best_sigma = best_sigma))
}

message("Run CoPro on all simulation configurations.")
message("Results are saved to copro_results_[scenario].csv")
message("Columns: prop_idx, run, p_A, p_B, p_C,")
message("         CoPro_pearson_overall, CoPro_pearson_A/B/C, best_sigma")
