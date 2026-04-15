# =============================================================================
# Simulation Setting 2: Three Cell Types with Varying Proportions
# =============================================================================
#
# Generates simulated spatial transcriptomics data with three cell types across
# a ternary grid of proportion configurations (p_A, p_B, p_C). This tests
# CoPro's robustness to varying cell-type compositions.
#
# The simulation works by:
#   1. Defining a ternary grid of proportions with step = 0.10, min = 0.10
#   2. For each proportion configuration, generating multiple replicates
#   3. Using the same pipeline as Setting 1 (PCA -> weight vectors -> spatial
#      smoothing -> percentile matching -> regression ground truth), extended
#      to three cell types
#
# Inputs:
#   - A Seurat object with three cell types (annotated in meta.data$annot)
#     containing normalized expression in assays$RNA$data
#
# Outputs (parquet files per configuration x replicate):
#   - prop{i}_run{j}_expression.parquet: cells x genes
#   - prop{i}_run{j}_metadata.parquet:   coords, cell types, proportions,
#                                         ground truth
#
# Dependencies: CoPro, Seurat, irlba, arrow, fields
# =============================================================================

library(CoPro)
library(Seurat)
library(irlba)
library(arrow)

# Source helper functions (run from repo root, or adjust path)
source("simulation/sim_fun.R")

# =============================================================================
# CONFIGURATION — update these paths for your system
# =============================================================================

# Path to Seurat object with three cell types
SC_DATA_PATH <- "path/to/your_seurat_three_types.rds"

# Output directory
OUTPUT_DIR <- "output/simulation/three_types"

# Simulation parameters
nPCA             <- 25
score_sd         <- 1
n_variable_genes <- 3000
cell_types       <- c("A", "B", "C")
n_points         <- 5500
bandwidth        <- 0.4
ds_multiplier    <- 1.08   # downsample multiplier (> 1 to ensure enough cells)
seed_base        <- 33
n_runs           <- 3      # replicates per proportion configuration

# Ternary grid parameters
step     <- 0.10
min_prop <- 0.10
max_prop <- 0.80

# =============================================================================
# BUILD TERNARY PROPORTION GRID
# =============================================================================

prop_grid <- data.frame()
for (p_A in seq(min_prop, max_prop, by = step)) {
  for (p_B in seq(min_prop, max_prop, by = step)) {
    p_C <- round(1 - p_A - p_B, 2)
    if (p_C >= min_prop && p_C <= max_prop) {
      prop_grid <- rbind(prop_grid, data.frame(p_A = p_A, p_B = p_B, p_C = p_C))
    }
  }
}

message("=== Proportion Grid ===")
message(paste("Configurations:", nrow(prop_grid)))
message(paste("Runs per config:", n_runs))
message(paste("Total simulations:", nrow(prop_grid) * n_runs, "\n"))

# =============================================================================
# LOAD AND PREPROCESS SINGLE-CELL DATA
# =============================================================================

message("Loading single-cell data...")
sc_data <- readRDS(SC_DATA_PATH)
ct_names <- unique(sc_data@meta.data$annot)
stopifnot(length(ct_names) >= 3)

sc_data <- FindVariableFeatures(sc_data, nfeatures = n_variable_genes)
sc_mat <- sc_data@assays$RNA@layers$data
rownames(sc_mat) <- rownames(sc_data)
colnames(sc_mat) <- sc_data@meta.data$cell

ct_mats_full <- lapply(ct_names[1:3], function(ct) {
  sc_mat[VariableFeatures(sc_data), sc_data@meta.data$annot == ct]
})
names(ct_mats_full) <- ct_names[1:3]

for (nm in names(ct_mats_full)) {
  message(paste("  ", nm, ":", ncol(ct_mats_full[[nm]]), "cells"))
}
message("")

# =============================================================================
# Helper: compute PCA, cell scores, and ground truth for one cell type
# =============================================================================

compute_fitted_ground_truth <- function(matched_mat, matched_scores, nPCA) {
  mat_scaled <- CoPro:::center_scale_matrix_opt(as.matrix(t(matched_mat)))
  pca_matched <- prcomp_irlba(mat_scaled, center = FALSE, scale. = FALSE, n = nPCA)
  pc_mat <- scale(pca_matched$x, center = FALSE, scale = pca_matched$sdev)
  rownames(pc_mat) <- colnames(matched_mat)
  fit <- lm(matched_scores ~ pc_mat)
  fitted_scores <- fitted(fit)
  names(fitted_scores) <- colnames(matched_mat)
  return(fitted_scores)
}

# =============================================================================
# GENERATE DATA FOR EACH PROPORTION CONFIGURATION
# =============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (prop_idx in seq_len(nrow(prop_grid))) {

  label_prob <- c(prop_grid$p_A[prop_idx],
                  prop_grid$p_B[prop_idx],
                  prop_grid$p_C[prop_idx])

  message("========================================")
  message(paste("Config", prop_idx, "/", nrow(prop_grid),
                "| p_A =", label_prob[1],
                ", p_B =", label_prob[2],
                ", p_C =", label_prob[3]))

  for (run_idx in seq_len(n_runs)) {
    message(paste("  Run", run_idx, "/", n_runs))

    run_seed <- run_idx + seed_base

    # Downsample sizes proportional to label_prob
    ds_sizes <- round(n_points * label_prob * ds_multiplier)

    # Generate independent weight vectors per cell type
    weights <- lapply(1:3, function(k) {
      generate_prob_vector(length = nPCA, seed = run_seed + k * 11)
    })

    # Downsample and compute PCA + cell scores for each cell type
    set.seed(run_seed)
    ct_mats <- list()
    cell_scores <- list()

    for (k in 1:3) {
      ct_name <- ct_names[k]
      ct_full <- ct_mats_full[[ct_name]]
      ds_n <- min(ds_sizes[k], ncol(ct_full))
      ct_mats[[k]] <- ct_full[, sample(seq_len(ncol(ct_full)), ds_n)]

      ct_scaled <- CoPro:::center_scale_matrix_opt(as.matrix(t(ct_mats[[k]])))
      pca_k <- prcomp_irlba(ct_scaled, center = FALSE, scale. = FALSE, n = nPCA)
      pc_mat <- scale(pca_k$x, center = FALSE, scale = pca_k$sdev)

      scores <- as.vector(pc_mat %*% as.matrix(weights[[k]], ncol = 1))
      names(scores) <- rownames(ct_scaled)
      cell_scores[[k]] <- scores
    }

    # Simulate spatial points with shared smooth score
    sim_sp <- simulate_smooth_points(
      n_points   = n_points,
      x_range    = c(-5, 5),
      y_range    = c(-7.5, 7.5),
      bandwidth  = bandwidth,
      n_rounds   = 1,
      score_sd   = score_sd,
      seed       = run_seed,
      labels     = cell_types,
      label_prob = label_prob,
      label_name = "cell_type"
    )

    # Match each cell type to real cells
    sim_parts <- list()
    match_failed <- FALSE
    for (k in 1:3) {
      sim_parts[[k]] <- tryCatch(
        match_by_percentile(
          sim_sp[sim_sp$cell_type == cell_types[k], ],
          cell_scores[[k]], rm_outlier = TRUE
        ),
        error = function(e) NULL
      )
      if (is.null(sim_parts[[k]])) { match_failed <- TRUE; break }
    }
    if (match_failed) {
      warning(paste("Matching failed: prop_idx", prop_idx, "run", run_idx))
      next
    }

    sim_sp <- do.call(rbind, sim_parts)
    rownames(sim_sp) <- sim_sp$matched_cell

    # Ground truth = PCA-regression fitted scores
    sim_sp$fitted_score <- NA
    for (k in 1:3) {
      matched_cells_k <- sim_parts[[k]]$matched_cell
      matched_mat_k <- ct_mats[[k]][, matched_cells_k]
      fitted_k <- compute_fitted_ground_truth(
        matched_mat_k, sim_parts[[k]]$matched_score, nPCA
      )
      sim_sp$fitted_score[sim_sp$cell_type == cell_types[k]] <-
        fitted_k[sim_sp$matched_cell[sim_sp$cell_type == cell_types[k]]]
    }

    # Expression matrix
    combined_mat <- do.call(cbind, ct_mats)
    sim_exp <- as.matrix(t(combined_mat[, sim_sp$matched_cell]))

    # Metadata
    metadata <- sim_sp[, c("x", "y", "cell_type", "initial_score",
                            "smoothed_score", "matched_cell",
                            "matched_score", "fitted_score")]
    metadata$prop_idx <- prop_idx
    metadata$run      <- run_idx
    metadata$p_A      <- label_prob[1]
    metadata$p_B      <- label_prob[2]
    metadata$p_C      <- label_prob[3]

    # Save as parquet
    exp_df <- as.data.frame(sim_exp)
    exp_df$cell_id <- rownames(sim_exp)
    exp_df <- exp_df[, c("cell_id", setdiff(colnames(exp_df), "cell_id"))]

    prefix <- paste0("prop", prop_idx, "_run", run_idx)
    write_parquet(exp_df,   file.path(OUTPUT_DIR, paste0(prefix, "_expression.parquet")))
    write_parquet(metadata, file.path(OUTPUT_DIR, paste0(prefix, "_metadata.parquet")))

    message(paste("    Saved:", prefix, "| Cells:", nrow(metadata),
                  "| Genes:", ncol(sim_exp)))
  }
}

message(paste("\n=== Done ==="))
message(paste("Output:", OUTPUT_DIR))
message(paste("Generated", nrow(prop_grid) * n_runs, "datasets"))
