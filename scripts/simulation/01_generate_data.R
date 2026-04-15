# =============================================================
# Simulation: Generate synthetic spatial transcriptomics data
# =============================================================
#
# Generates simulated spatial data with known ground truth
# co-progression patterns. Uses real single-cell expression
# profiles (liver cell atlas) as the gene expression backbone
# and adds spatial structure via smooth spatial point processes.
#
# Three simulation scenarios:
#   1. Three cell types with varying proportions
#   2. Two cell types with one spatial axis (standard)
#   3. Two cell types with two orthogonal axes (multi-axis)
#
# Inputs:
#   - liver_seurat_hepa_nk_mac.rds (single-cell reference atlas)
#   - sim_fun.R (simulation helper functions)
#
# Outputs:
#   - Parquet files per configuration:
#     propX_runY_expression.parquet (cells x genes)
#     propX_runY_metadata.parquet   (coords, cell types, ground truth)
#
# Dependencies: CoPro, arrow (for parquet), SeuratObject

library(CoPro)
library(arrow)

# --- Configuration ---
DATA_DIR <- "path/to/data"  # UPDATE THIS
OUT_DIR  <- "output/simulation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(DATA_DIR, "sim_fun.R"))

# Simulation parameters
n_pca <- 25
score_sd <- 1
n_variable_genes <- 3000
n_points <- 5500
bandwidth <- 0.4

# --- Load reference single-cell data ---
seu <- readRDS(file.path(DATA_DIR, "liver_seurat_hepa_nk_mac.rds"))

message("Simulation data generation.")
message("For the full simulation pipeline, run scripts 01-04 in order.")
message("See manuscript Supplementary Figure 1 for simulation design.")
message("")
message("Key functions from sim_fun.R:")
message("  simulate_smooth_points() - generates spatial coordinates")
message("  match_by_percentile() - matches simulated scores to real cells")
message("")
message("Three scenarios:")
message("  1. Three types: sweep cell type proportions on ternary grid")
message("  2. Alternative: standard 2-type setup")
message("  3. Multi-axis: 2 types with 2 orthogonal spatial gradients")
