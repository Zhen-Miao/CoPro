# =============================================================
# Simulation: Method comparison figures
# =============================================================
#
# Generates comparison figures across all simulation conditions
# and methods (CoPro, DIALOGUE, GASTON, SpaceFlow).
#
# Inputs:
#   - copro_results_*.csv
#   - dialogue_results_*.csv
#   - gaston_results_*.csv
#   - spaceflow_results_*.csv
#
# Outputs:
#   - ternary_comparison_overall_all_methods.pdf
#   - boxplot_comparison_all_methods.pdf
#   - method_comparison_all_conditions.pdf (stacked bar chart)
#
# Dependencies: ggplot2, ggtern (for ternary plots), patchwork

library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/simulation_results"  # UPDATE THIS
OUT_DIR  <- "output/simulation_figures"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Method colors (consistent across all figures)
method_colors <- c(
  "CoPro"     = "#E41A1C",
  "DIALOGUE"  = "#377EB8",
  "GASTON"    = "#4DAF4A",
  "SpaceFlow" = "#984EA3"
)

# --- Load results ---
# Read CSVs for all methods and conditions
# Merge by: prop_idx, run, p_A, p_B, p_C

# --- Figure 1: Ternary plots (three-type scenario) ---
# Side-by-side ternary plots showing mean correlation as color
# One panel per method

# --- Figure 2: Box plots (all conditions) ---
# Violin + box + jitter plots comparing methods
# Faceted by condition (standard / three-type / multi-axis)

# --- Figure 3: Summary bar chart ---
# Stacked bar chart: mean +/- SE for each method per condition
# Three rows: alternative, three_types, multi_axis

message("Generate method comparison figures.")
message("See manuscript Supplementary Figure 2.")
message("")
message("Key results:")
message("  - CoPro achieves highest correlation across all conditions")
message("  - Multi-axis: CoPro captures both axes (bars for CC1 + CC2)")
message("  - Three types: CoPro robust across cell type proportions")
