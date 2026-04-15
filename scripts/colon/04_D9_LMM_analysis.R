# =============================================================
# Colon Day 9: Linear mixed model analysis
# =============================================================
#
# Fits linear mixed models to D9 transferred cell scores to
# assess gene-score associations across multiple slides while
# accounting for slide-level random effects.
#
# Inputs:
#   - Transferred cell scores from 03_D9_cross_slide_transfer.R
#   - Expression and metadata files
#
# Outputs:
#   - LMM results table (CSV)
#   - Supplementary gene table
#
# Dependencies: CoPro, lme4, ggplot2

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/data"  # UPDATE THIS
OUT_DIR  <- "output/colon_D9_LMM"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cell_types <- c("Epithelial", "Fibroblast", "Immune")
sigma_opt <- 0.1

# --- Load transferred scores and metadata ---
# (Load pre-computed transferred scores from all regions)
# Each file contains cell scores transferred from the reference slide
# to a target slide. The scores are used as the response variable in
# the LMM.

message("Load transferred scores and fit LMMs per gene.")
message("See manuscript methods for model specification:")
message("  gene_expression ~ cell_score + (1 | slide_id)")
message("Results in supplementary tables S6.")
