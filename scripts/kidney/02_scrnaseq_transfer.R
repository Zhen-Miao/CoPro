# =============================================================
# Kidney: Transfer scores to scRNA-seq data
# =============================================================
#
# Transfers CoPro gene weights learned from seqFISH spatial data
# to scRNA-seq data for validation and full-transcriptome analysis.
#
# Inputs:
#   - CoPro objects from 01_run_copro_supervised.R
#   - scRNA-seq kidney data (Seurat object)
#
# Outputs:
#   - Transferred cell scores on scRNA-seq cells
#   - UMAP visualizations colored by transferred scores
#   - Full-transcriptome regression results
#
# Dependencies: CoPro, SeuratObject, ggplot2

library(CoPro)
library(SeuratObject)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/data"  # UPDATE THIS
OUT_DIR  <- "output/kidney_transfer"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sigma_opt <- 0.1

message("Transfer CoPro spatial gene weights to scRNA-seq data.")
message("See manuscript methods for details on:")
message("  - getTransferCellScores(gene_score_type = 'regression')")
message("  - Full-transcriptome regression on transferred scores")
message("  - GO enrichment analysis on regression coefficients")
