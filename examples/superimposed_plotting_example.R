# Example: Creating Superimposed Cell Type Score Plots
# 
# This script demonstrates how to create superimposed plots with different 
# color gradients for multiple cell types using the CoPro package.

library(CoPro)
library(ggplot2)
library(viridis)

# Assuming you have already run the CoPro pipeline and have:
# - br: your CoPro object with computed cell scores
# - Two cell types: "061 STR D1 Gaba" and "062 STR D2 Gaba"

# Method 1: Basic superimposed plot with default colors
# This will use viridis for the first cell type and red-pink for the second
plot1 <- plotSuperimposedCellScores(
  object = br,
  sigmaValueChoice = 0.14,  # Use your optimal sigma value
  ccIndex = 1,              # First canonical component
  cell_types = c("061 STR D1 Gaba", "062 STR D2 Gaba"),
  point_size = 0.6,
  alpha_values = list("061 STR D1 Gaba" = 0.8, "062 STR D2 Gaba" = 0.8)
)

print(plot1)

# Method 2: Custom color scales
custom_colors <- list(
  "061 STR D1 Gaba" = "viridis",
  "062 STR D2 Gaba" = "red_pink"
)

plot2 <- plotSuperimposedCellScores(
  object = br,
  sigmaValueChoice = 0.14,
  ccIndex = 1,
  cell_types = c("061 STR D1 Gaba", "062 STR D2 Gaba"),
  color_scales = custom_colors,
  point_size = 0.8,
  alpha_values = list("061 STR D1 Gaba" = 0.7, "062 STR D2 Gaba" = 0.7),
  title = "Superimposed Cell Scores: D1 vs D2 Neurons"
)

print(plot2)

# Method 3: With separate legends for each cell type
result_with_legends <- plotSuperimposedWithLegends(
  object = br,
  sigmaValueChoice = 0.14,
  ccIndex = 1,
  cell_types = c("061 STR D1 Gaba", "062 STR D2 Gaba"),
  color_scales = custom_colors,
  point_size = 0.8,
  create_legends = TRUE
)

# View the main plot
print(result_with_legends$main_plot)

# View individual legends
print(result_with_legends$legend_plots[["061 STR D1 Gaba"]])
print(result_with_legends$legend_plots[["062 STR D2 Gaba"]])

# Method 4: Side-by-side comparison (faceted plot)
comparison_plot <- plotCellTypeComparison(
  object = br,
  sigmaValueChoice = 0.14,
  ccIndex = 1,
  cell_types = c("061 STR D1 Gaba", "062 STR D2 Gaba"),
  use_viridis = TRUE,
  point_size = 0.8,
  ncol = 2
)

print(comparison_plot)

# Method 5: Advanced customization with multiple cell types
if (length(unique(br@cellTypesSub)) > 2) {
  # Get all available cell types
  all_cell_types <- unique(br@cellTypesSub)
  
  # Use different color scales for each
  multi_colors <- list(
    "viridis",
    "red_pink", 
    "blue_purple",
    "plasma",
    "inferno",
    "magma"
  )
  names(multi_colors) <- all_cell_types[1:min(length(all_cell_types), 6)]
  
  multi_plot <- plotSuperimposedCellScores(
    object = br,
    sigmaValueChoice = 0.14,
    ccIndex = 1,
    cell_types = all_cell_types[1:min(length(all_cell_types), 4)],  # Max 4 for clarity
    color_scales = multi_colors,
    point_size = 0.6,
    title = "Multi Cell Type Superimposed Scores"
  )
  
  print(multi_plot)
}

# Method 6: Creating a comprehensive visualization
# This combines the superimposed plot with the comparison plot
library(gridExtra)  # For combining plots

if (requireNamespace("gridExtra", quietly = TRUE)) {
  combined_visualization <- gridExtra::grid.arrange(
    plot2,
    comparison_plot,
    ncol = 1,
    top = "Cell Type Score Analysis: Superimposed vs Separated Views"
  )
}

# Tips for interpretation:
# 1. In the superimposed plot, areas where both cell types have high scores 
#    will show mixed colors
# 2. Use different alpha values to make overlapping regions more visible
# 3. The quantile_range parameter (default c(0.01, 0.99)) can be adjusted 
#    to enhance contrast
# 4. For publications, consider using colorblind-friendly palettes

# Example with colorblind-friendly options:
colorblind_friendly <- list(
  "061 STR D1 Gaba" = "viridis",     # Purple-blue-green
  "062 STR D2 Gaba" = "red_pink"      # Red-pink
)

final_plot <- plotSuperimposedCellScores(
  object = br,
  sigmaValueChoice = 0.14,
  ccIndex = 1,
  cell_types = c("061 STR D1 Gaba", "062 STR D2 Gaba"),
  color_scales = colorblind_friendly,
  point_size = 0.8,
  alpha_values = list("061 STR D1 Gaba" = 0.75, "062 STR D2 Gaba" = 0.75),
  title = "Cell Type Correspondence: Score Distribution Analysis",
  quantile_range = c(0.05, 0.95)  # Slightly narrower range for better contrast
)

print(final_plot)

# Save plots if needed
# ggsave("superimposed_cell_scores.png", plot = final_plot, 
#        width = 10, height = 8, dpi = 300) 