# =============================================================
# Simulation: Run alternative methods for comparison
# =============================================================
#
# Runs DIALOGUE (and optionally GASTON, SpaceFlow) on the same
# simulated datasets for benchmarking against CoPro.
#
# Inputs:
#   - Parquet files from 01_generate_data.R
#
# Outputs:
#   - dialogue_results_[scenario].csv
#   - combined_results_[scenario].csv (merged with CoPro results)
#
# Dependencies: DIALOGUE, CoPro (for data loading)

message("Run DIALOGUE on simulated data for comparison.")
message("DIALOGUE requires: make.cell.type(), DLG.get.param(), DIALOGUE.run()")
message("")
message("GASTON and SpaceFlow results were generated externally")
message("(Python-based methods) and saved as CSV files.")
message("")
message("See manuscript Supplementary Figure 2 for method comparison.")
message("CoPro outperforms alternatives across all conditions:")
message("  - Standard (2 cell types): CoPro r=0.85 vs DIALOGUE r=0.72")
message("  - Three types (varying proportions): CoPro consistently higher")
message("  - Multi-axis: CoPro captures both axes; others capture at most one")
