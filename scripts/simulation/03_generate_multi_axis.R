# =============================================================================
# Simulation Setting 3: Two Cell Types, Two Orthogonal Co-Progression Axes
# =============================================================================
#
# Generates simulated spatial transcriptomics data with two cell types and two
# independent, orthogonal co-progression axes. This tests whether CoPro can
# recover multiple latent spatial programs simultaneously (nCC >= 2).
#
# Key differences from Settings 1-2:
#   - Orthogonal PC weight vectors (via Gram-Schmidt) define two independent
#     latent axes per cell type
#   - Two independent smooth spatial fields are generated and orthogonalized
#   - 2D optimal transport matching maps simulated (score1, score2) targets
#     to real cells in 2D score space (greedy, without replacement)
#
# Inputs:
#   - A Seurat object with two cell types (annotated in meta.data$annot)
#     containing normalized expression in assays$RNA$data
#
# Outputs (parquet files per replicate):
#   - sim{i}_expression.parquet: cells x genes
#   - sim{i}_metadata.parquet:   coords, cell types, ground truth per axis
#
# Dependencies: CoPro, Seurat, irlba, arrow, fields
# =============================================================================

library(CoPro)
library(Seurat)
library(irlba)
library(arrow)

source(file.path(find.package("CoPro"), "..", "CoPro", "simulation", "sim_fun.R"))

# =============================================================================
# CONFIGURATION — update these paths for your system
# =============================================================================

# Path to Seurat object with two cell types
SC_DATA_PATH <- "path/to/your_seurat_two_types.rds"

# Output directory
OUTPUT_DIR <- "output/simulation/multi_axis"

# Simulation parameters
n_simulations    <- 5
n_axes           <- 2
nPCA             <- 25
n_variable_genes <- 3000
n_points         <- 4500
x_range          <- c(-5, 5)
y_range          <- c(-7.5, 7.5)
bandwidth        <- 0.3
score_sd         <- 1
cell_types       <- c("A", "B")
label_prob       <- c(0.4, 0.6)
down_sample_ct1  <- 1890
down_sample_ct2  <- 2835
seed_downsample  <- 458
seed_weights     <- 42
matching_method  <- "greedy_ot"
use_quantile_matching <- TRUE

# =============================================================================
# LOAD AND PREPROCESS SINGLE-CELL DATA (once)
# =============================================================================

message("Loading single-cell data...")
sc_data <- readRDS(SC_DATA_PATH)
ctypes <- unique(sc_data@meta.data$annot)
stopifnot(length(ctypes) >= 2)

sc_data <- FindVariableFeatures(sc_data, nfeatures = n_variable_genes)
hvgs <- VariableFeatures(sc_data)

sc_mat <- sc_data@assays$RNA@layers$data
rownames(sc_mat) <- rownames(sc_data)
colnames(sc_mat) <- sc_data@meta.data$cell

ct1_mat_full <- sc_mat[hvgs, sc_data@meta.data$annot == ctypes[1]]
ct2_mat_full <- sc_mat[hvgs, sc_data@meta.data$annot == ctypes[2]]

# Downsample (fixed seed across all simulations)
set.seed(seed_downsample)
ct1_mat <- ct1_mat_full[, sample(seq_len(ncol(ct1_mat_full)), down_sample_ct1)]
ct2_mat <- ct2_mat_full[, sample(seq_len(ncol(ct2_mat_full)), down_sample_ct2)]

message(paste("  CT1 (", ctypes[1], "):", ncol(ct1_mat), "cells"))
message(paste("  CT2 (", ctypes[2], "):", ncol(ct2_mat), "cells"))

# =============================================================================
# COMPUTE PCA (once)
# =============================================================================

message("Computing PCA...")
pca_results <- list()
ct_mats <- list(ct1 = ct1_mat, ct2 = ct2_mat)

for (ct_name in c("ct1", "ct2")) {
  mat_scaled <- CoPro:::center_scale_matrix_opt(as.matrix(t(ct_mats[[ct_name]])))
  pca <- prcomp_irlba(mat_scaled, center = FALSE, scale. = FALSE, n = nPCA)
  pc_mat <- scale(pca$x, center = FALSE, scale = pca$sdev)
  rownames(pc_mat) <- colnames(ct_mats[[ct_name]])
  pca_results[[ct_name]] <- list(pca = pca, pc_mat = pc_mat)
}

# Generate orthogonal weight vectors (once)
weight_list <- generate_orthogonal_weights(
  n_weights = nPCA, n_axes = n_axes, seed = seed_weights
)

# Compute cell scores per axis (once)
cell_scores <- list()
for (ct_name in c("ct1", "ct2")) {
  pc_mat <- pca_results[[ct_name]]$pc_mat
  cell_scores[[ct_name]] <- list()
  for (k in seq_len(n_axes)) {
    scores_k <- as.vector(pc_mat %*% weight_list[[k]])
    names(scores_k) <- rownames(pc_mat)
    cell_scores[[ct_name]][[k]] <- scores_k
  }
}

message("Setup complete.\n")

# =============================================================================
# GENERATE SIMULATION DATA
# =============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message(paste("Generating", n_simulations, "multi-axis simulation datasets...\n"))

for (sim_id in seq_len(n_simulations)) {
  message(paste("--- Simulation", sim_id, "/", n_simulations, "---"))

  seed_simulation <- 100 + sim_id * 111

  # --- Step 1: Generate spatial points with one smooth field ---
  sim_sp <- simulate_smooth_points(
    n_points   = n_points,
    x_range    = x_range,
    y_range    = y_range,
    bandwidth  = bandwidth,
    n_rounds   = 1,
    score_sd   = score_sd,
    seed       = seed_simulation,
    labels     = cell_types,
    label_prob = label_prob,
    label_name = "cell_type"
  )

  ct1_idx <- sim_sp$cell_type == "A"
  ct2_idx <- sim_sp$cell_type == "B"

  # --- Step 2: Generate multi-axis target scores ---
  # First axis reuses the smooth field from simulate_smooth_points.
  # Additional axes are independently smoothed, then orthogonalized.
  locations <- sim_sp[, c("x", "y")]
  dist_mat <- as.matrix(dist(locations))

  target_scores <- list()
  for (k in seq_len(n_axes)) {
    if (k == 1) {
      target_scores[[k]] <- sim_sp$smoothed_score
    } else {
      set.seed(seed_simulation + k * 1000)
      raw_scores <- rnorm(n_points, sd = score_sd)
      kernel_mat <- exp(-dist_mat^2 / (2 * bandwidth^2))
      kernel_mat <- kernel_mat / rowSums(kernel_mat)
      target_scores[[k]] <- as.vector(kernel_mat %*% raw_scores)
    }
  }
  target_scores <- orthogonalize_scores(target_scores)

  target_scores_ct1 <- lapply(target_scores, function(x) x[ct1_idx])
  target_scores_ct2 <- lapply(target_scores, function(x) x[ct2_idx])

  # --- Step 3: 2D optimal transport matching ---
  matched_results <- list()
  for (ct_name in c("ct1", "ct2")) {
    target_ct <- if (ct_name == "ct1") target_scores_ct1 else target_scores_ct2
    cell_sc   <- cell_scores[[ct_name]]
    cell_ids  <- names(cell_sc[[1]])

    matched_results[[ct_name]] <- match_2d_optimal_transport(
      target_scores_list = target_ct,
      cell_scores_list   = cell_sc,
      cell_ids           = cell_ids,
      rm_outlier         = TRUE,
      method             = matching_method,
      use_quantile       = use_quantile_matching,
      verbose            = FALSE
    )
  }

  # --- Step 4: Build matched expression matrices ---
  matched_mats <- list()
  for (ct_name in c("ct1", "ct2")) {
    matched_ids <- matched_results[[ct_name]]$ids
    mat_use <- ct_mats[[ct_name]][, matched_ids]
    colnames(mat_use) <- paste0(ct_name, "_cell_", seq_len(ncol(mat_use)))
    matched_mats[[ct_name]] <- mat_use
  }

  # --- Step 5: Recompute PCA on matched cells -> regression ground truth ---
  matched_pca <- list()
  for (ct_name in c("ct1", "ct2")) {
    mat_scaled <- CoPro:::center_scale_matrix_opt(
      as.matrix(t(matched_mats[[ct_name]]))
    )
    pca <- prcomp_irlba(mat_scaled, center = FALSE, scale. = FALSE, n = nPCA)
    pc_mat <- scale(pca$x, center = FALSE, scale = pca$sdev)
    rownames(pc_mat) <- colnames(matched_mats[[ct_name]])
    matched_pca[[ct_name]] <- list(pca = pca, pc_mat = pc_mat)

    result <- matched_results[[ct_name]]
    for (k in seq_len(n_axes)) {
      target <- result$matched_scores[, k]
      fit <- lm(target ~ pc_mat)
      matched_pca[[ct_name]][[paste0("fitted_axis_", k)]] <- fitted(fit)
    }
  }

  # --- Step 6: Assemble output ---
  exp_ct1 <- t(matched_mats[["ct1"]])
  exp_ct2 <- t(matched_mats[["ct2"]])
  sim_exp <- as.matrix(rbind(exp_ct1, exp_ct2))

  meta_ct1 <- data.frame(
    cell = rownames(exp_ct1),
    x    = sim_sp$x[ct1_idx],
    y    = sim_sp$y[ct1_idx],
    cell_type = "A",
    stringsAsFactors = FALSE
  )
  meta_ct2 <- data.frame(
    cell = rownames(exp_ct2),
    x    = sim_sp$x[ct2_idx],
    y    = sim_sp$y[ct2_idx],
    cell_type = "B",
    stringsAsFactors = FALSE
  )

  # Add fitted ground truth scores per axis
  for (k in seq_len(n_axes)) {
    col_fitted  <- paste0("fitted_score_axis_", k)
    col_matched <- paste0("matched_score_axis_", k)
    meta_ct1[[col_fitted]]  <- matched_pca[["ct1"]][[paste0("fitted_axis_", k)]]
    meta_ct2[[col_fitted]]  <- matched_pca[["ct2"]][[paste0("fitted_axis_", k)]]
    meta_ct1[[col_matched]] <- matched_results[["ct1"]]$matched_scores_raw[, k]
    meta_ct2[[col_matched]] <- matched_results[["ct2"]]$matched_scores_raw[, k]
  }

  metadata <- rbind(meta_ct1, meta_ct2)
  rownames(metadata) <- metadata$cell
  colnames(metadata)[colnames(metadata) == "cell"] <- "matched_cell"

  metadata$sim_id <- sim_id
  metadata$seed   <- seed_simulation

  # Save as parquet
  exp_df <- as.data.frame(sim_exp)
  exp_df$cell_id <- rownames(sim_exp)
  exp_df <- exp_df[, c("cell_id", setdiff(colnames(exp_df), "cell_id"))]

  prefix <- paste0("sim", sim_id)
  write_parquet(exp_df,   file.path(OUTPUT_DIR, paste0(prefix, "_expression.parquet")))
  write_parquet(metadata, file.path(OUTPUT_DIR, paste0(prefix, "_metadata.parquet")))

  message(paste("  Saved:", prefix, "| Cells:", nrow(metadata),
                "| Genes:", ncol(sim_exp)))
}

message(paste("\n=== Done ==="))
message(paste("Output:", OUTPUT_DIR))
message(paste("Generated", n_simulations, "multi-axis datasets"))
