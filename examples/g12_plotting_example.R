# Example: Plotting g_12(r) pair correlation functions
# This example demonstrates how to use plotG12Functions to visualize
# cross-type pair correlation functions for spatial colocalization analysis

library(CoPro)

# Assuming you have a CoPro object with spatial data
# Replace 'your_object' with your actual CoPro object

# Basic usage - combined plot for all cell type pairs
g12_result <- plotG12Functions(your_object)

# Display the plot
print(g12_result$plot)

# View the raw data
head(g12_result$data)

# View summary statistics
print(g12_result$summary)

# Custom parameters for high-resolution data
g12_result_custom <- plotG12Functions(
  your_object,
  r_um_range = c(5, 40),        # Analyze shorter distances
  pixel_size_um = 0.325,        # High-resolution pixel size
  cell_diam_um = 8,             # Smaller cell diameter
  nsim = 99,                    # Fewer simulations for speed
  plot_type = "individual",     # Separate plots per pair
  include_confidence = TRUE,    # Include confidence bands
  confidence_level = 0.95       # 95% confidence intervals
)

# For multi-slide objects, the function automatically handles each slide
# The resulting plots will show data aggregated across slides

# Access individual components:
# - g12_result$plot: The ggplot object
# - g12_result$data: Raw g_12(r) values and confidence intervals
# - g12_result$summary: Summary statistics per cell type pair

# Save the plot
# ggsave("g12_correlation_functions.pdf", g12_result$plot, width = 10, height = 6)

# Interpretation of g_12(r) values:
# - g_12(r) = 1: Random spatial relationship at distance r
# - g_12(r) > 1: Colocalization/attraction at distance r
# - g_12(r) < 1: Segregation/repulsion at distance r

# Example with custom colors
custom_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
g12_colored <- plotG12Functions(
  your_object,
  colors = custom_colors,
  plot_type = "combined"
)

print(g12_colored$plot)
