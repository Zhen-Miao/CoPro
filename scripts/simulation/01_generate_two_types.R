# =============================================================================
# Simulation Setting 1: Two Cell Types, Single Co-Progression Axis
# =============================================================================
#
# Generates simulated spatial transcriptomics data with two cell types and a
# single co-progression axis. Produces both:
#   - Alternative hypothesis data (WITH cross-type spatial coordination)
#   - Null hypothesis data (NO cross-type coordination)
#
# The simulation works by:
#   1. Loading a real single-cell reference (for realistic gene expression)
#   2. Computing PCA on each cell type separately
#   3. Generating random PC weight vectors to define a latent cell score
#   4. Simulating 2D spatial points with smoothed scores (alternative) or
#      independent scores per cell type (null)
#   5. Matching simulated points to real cells by percentile rank
#   6. Recomputing PCA on matched cells and fitting regression to get the
#      ground truth scores (theoretical upper limit recoverable from gene
#      expression)
#
# Inputs:
#   - A Seurat object with two cell types (annotated in meta.data$annot)
#     containing normalized expression in assays$RNA$data
#
# Outputs (parquet files):
#   - alternative_round{i}_expression.parquet: cells x genes
#   - alternative_round{i}_metadata.parquet:   coords, cell types, ground truth
#   - null_round{i}_expression.parquet
#   - null_round{i}_metadata.parquet
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

# Output directory for simulation data
OUTPUT_DIR <- "output/simulation/two_types"

# Simulation parameters
nPCA             <- 25     # Number of PCs
score_sd         <- 1      # SD of initial random scores
n_variable_genes <- 3000   # Number of HVGs to use
n_rounds_sim     <- 5      # Number of simulation replicates
n_points         <- 4500   # Total cells per simulation
cell_types       <- c("A", "B")
label_prob       <- c(0.4, 0.6)

# Downsampling sizes (should be > n_points * label_prob to allow matching)
down_sample_ct1  <- 1890
down_sample_ct2  <- 2835

# =============================================================================
# LOAD AND PREPROCESS SINGLE-CELL DATA
# =============================================================================

message("Loading single-cell data...")
sc_data <- readRDS(SC_DATA_PATH)
ctypes <- unique(sc_data@meta.data$annot)
stopifnot(length(ctypes) >= 2)

sc_data <- FindVariableFeatures(sc_data, nfeatures = n_variable_genes)
sc_mat <- sc_data@assays$RNA@layers$data
rownames(sc_mat) <- rownames(sc_data)
colnames(sc_mat) <- sc_data@meta.data$cell

ct1_mat_full <- sc_mat[VariableFeatures(sc_data), sc_data@meta.data$annot == ctypes[1]]
ct2_mat_full <- sc_mat[VariableFeatures(sc_data), sc_data@meta.data$annot == ctypes[2]]

message(paste("  Cell type 1 (", ctypes[1], "):", ncol(ct1_mat_full), "cells"))
message(paste("  Cell type 2 (", ctypes[2], "):", ncol(ct2_mat_full), "cells"))

# Downsample
set.seed(458)
ct1_mat <- ct1_mat_full[, sample(seq_len(ncol(ct1_mat_full)), down_sample_ct1)]
ct2_mat <- ct2_mat_full[, sample(seq_len(ncol(ct2_mat_full)), down_sample_ct2)]

# Center, scale, PCA
ct1_scaled <- CoPro:::center_scale_matrix_opt(as.matrix(t(ct1_mat)))
ct2_scaled <- CoPro:::center_scale_matrix_opt(as.matrix(t(ct2_mat)))

pca_1 <- prcomp_irlba(ct1_scaled, center = FALSE, scale. = FALSE, n = nPCA)
pca_2 <- prcomp_irlba(ct2_scaled, center = FALSE, scale. = FALSE, n = nPCA)

PCmats1 <- scale(pca_1$x, center = FALSE, scale = pca_1$sdev)
PCmats2 <- scale(pca_2$x, center = FALSE, scale = pca_2$sdev)

# Generate random PC weight vectors -> latent cell scores
weight1 <- generate_prob_vector(length = nPCA, seed = 3)
weight2 <- generate_prob_vector(length = nPCA, seed = 33)

cell_score_1 <- as.vector(PCmats1 %*% as.matrix(weight1, ncol = 1))
cell_score_2 <- as.vector(PCmats2 %*% as.matrix(weight2, ncol = 1))
names(cell_score_1) <- rownames(ct1_scaled)
names(cell_score_2) <- rownames(ct2_scaled)

message("Data preparation complete.\n")

# =============================================================================
# Helper: recompute PCA on matched cells and fit regression for ground truth
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
# Helper: run one round of simulation and save
# =============================================================================

run_one_round <- function(round_id, condition, sim_fun, sim_seed,
                          n_smoothing_rounds) {
  message(paste("\n---", condition, "| Round", round_id, "---"))

  sim_sp <- sim_fun(
    n_points     = n_points,
    x_range      = c(-5, 5),
    y_range      = c(-7.5, 7.5),
    bandwidth    = 0.3,
    n_rounds     = n_smoothing_rounds,
    score_sd     = score_sd,
    seed         = sim_seed,
    labels       = cell_types,
    label_prob   = label_prob,
    label_name   = "cell_type"
  )

  # Match to real single-cell data
  sim_sp_1 <- match_by_percentile(
    sim_sp[sim_sp$cell_type == "A", ], cell_score_1, rm_outlier = TRUE
  )
  sim_sp_2 <- match_by_percentile(
    sim_sp[sim_sp$cell_type == "B", ], cell_score_2, rm_outlier = TRUE
  )

  sim_sp <- rbind(sim_sp_1, sim_sp_2)
  rownames(sim_sp) <- sim_sp$matched_cell

  # Matched expression matrices
  matched_ct1_mat <- ct1_mat[, sim_sp_1$matched_cell]
  matched_ct2_mat <- ct2_mat[, sim_sp_2$matched_cell]

  # Ground truth = PCA-regression fitted scores
  fitted_ct1 <- compute_fitted_ground_truth(
    matched_ct1_mat, sim_sp_1$matched_score, nPCA
  )
  fitted_ct2 <- compute_fitted_ground_truth(
    matched_ct2_mat, sim_sp_2$matched_score, nPCA
  )

  sim_sp$fitted_score <- NA
  sim_sp$fitted_score[sim_sp$cell_type == "A"] <-
    fitted_ct1[sim_sp$matched_cell[sim_sp$cell_type == "A"]]
  sim_sp$fitted_score[sim_sp$cell_type == "B"] <-
    fitted_ct2[sim_sp$matched_cell[sim_sp$cell_type == "B"]]

  # Expression matrix (genes x cells -> cells x genes)
  sim_exp <- as.matrix(t(cbind(ct1_mat, ct2_mat)[, sim_sp$matched_cell]))

  # Metadata
  metadata <- sim_sp[, c("x", "y", "cell_type", "initial_score",
                          "smoothed_score", "matched_cell",
                          "matched_score", "fitted_score")]
  metadata$round     <- round_id
  metadata$condition <- condition

  # Save as parquet
  exp_df <- as.data.frame(sim_exp)
  exp_df$cell_id <- rownames(sim_exp)
  exp_df <- exp_df[, c("cell_id", setdiff(colnames(exp_df), "cell_id"))]

  prefix <- paste0(condition, "_round", round_id)
  write_parquet(exp_df,   file.path(OUTPUT_DIR, paste0(prefix, "_expression.parquet")))
  write_parquet(metadata, file.path(OUTPUT_DIR, paste0(prefix, "_metadata.parquet")))

  message(paste("  Saved:", prefix, "| Cells:", nrow(metadata),
                "| Genes:", ncol(sim_exp)))
}

# =============================================================================
# GENERATE ALTERNATIVE (SIGNAL) AND NULL DATA
# =============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

message(paste("Generating", n_rounds_sim,
              "rounds each for alternative and null conditions...\n"))

for (i in seq_len(n_rounds_sim)) {
  # Alternative: shared smooth score -> cross-type coordination
  run_one_round(
    round_id          = i,
    condition         = "alternative",
    sim_fun           = simulate_smooth_points,
    sim_seed          = 123 + i,
    n_smoothing_rounds = 1
  )

  # Null: independent scores per cell type, no cross-type coordination
  run_one_round(
    round_id          = i,
    condition         = "null",
    sim_fun           = simulate_smooth_points_null,
    sim_seed          = 2000 + i,
    n_smoothing_rounds = 0
  )
}

message(paste("\n=== Done ==="))
message(paste("Output:", OUTPUT_DIR))
message(paste("Generated", n_rounds_sim, "rounds x 2 conditions =",
              n_rounds_sim * 2, "datasets"))
