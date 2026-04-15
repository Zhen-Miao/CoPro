# =============================================================
# Liver: Zonation score transfer (healthy -> mutant)
# =============================================================
#
# Transfers zonation-related gene weights from a healthy liver
# sample to a mutant sample to assess conservation of
# pericentral-periportal patterns.
#
# Inputs:
#   - CoPro object from 01_run_copro.R (healthy reference)
#   - Mutant liver data
#
# Outputs:
#   - Transferred zonation scores
#   - Rank-based zonation comparison plot
#
# Dependencies: CoPro, ggplot2

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/liver_data"  # UPDATE THIS
OUT_DIR  <- "output/liver_transfer"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sigma_opt <- 0.05

message("Transfer liver zonation scores from healthy to mutant.")
message("See manuscript Fig 7e-f for zonation transfer results.")
message("Use getTransferCellScores(gene_score_type = 'regression')")
message("for regression-based transfer.")
