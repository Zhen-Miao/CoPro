# Example usage of getCellScores function
# This demonstrates how to safely access cell scores from CoPro objects

# Load required libraries
library(CoPro)

# =============================================================================
# Example 1: CoProSingle object usage
# =============================================================================

# Assuming you have a CoProSingle object named 'object_single'
# (This would typically be created from your spatial transcriptomics data)

# Basic usage - get all cell scores for a specific sigma and cell type
if (exists("object_single") && is(object_single, "CoProSingle")) {
  
  # Get available sigma values and cell types
  available_sigmas <- object_single@sigmaValues
  available_cell_types <- names(object_single@cellScores[[1]])
  
  cat("Available sigma values:", paste(available_sigmas, collapse = ", "), "\n")
  cat("Available cell types:", paste(available_cell_types, collapse = ", "), "\n")
  
  # Example 1a: Get complete cell scores matrix
  if (length(available_sigmas) > 0 && length(available_cell_types) > 0) {
    scores_full <- getCellScores(object = object_single, 
                                sigma = available_sigmas[1],
                                cellType = available_cell_types[1])
    
    cat("Full scores matrix dimensions:", dim(scores_full), "\n")
    cat("First few cell IDs:", head(rownames(scores_full)), "\n")
    cat("Canonical components:", colnames(scores_full), "\n")
  }
  
  # Example 1b: Get specific canonical component
  if (length(available_sigmas) > 0 && length(available_cell_types) > 0) {
    cc1_scores <- getCellScores(object = object_single,
                               sigma = available_sigmas[1],
                               cellType = available_cell_types[1],
                               ccIndex = 1)
    
    cat("CC1 scores length:", length(cc1_scores), "\n")
    cat("First few CC1 scores:", head(cc1_scores), "\n")
  }
  
  # Example 1c: Get scores for specific cells
  if (length(available_sigmas) > 0 && length(available_cell_types) > 0) {
    all_scores <- getCellScores(object = object_single,
                               sigma = available_sigmas[1],
                               cellType = available_cell_types[1])
    
    # Get first 5 cell IDs as an example
    example_cells <- head(rownames(all_scores), 5)
    
    specific_scores <- getCellScores(object = object_single,
                                    sigma = available_sigmas[1],
                                    cellType = available_cell_types[1],
                                    cells = example_cells)
    
    cat("Specific cells scores dimensions:", dim(specific_scores), "\n")
  }
  
  # Example 1d: Combine ccIndex and cells parameters
  if (length(available_sigmas) > 0 && length(available_cell_types) > 0) {
    all_scores <- getCellScores(object = object_single,
                               sigma = available_sigmas[1],
                               cellType = available_cell_types[1])
    
    example_cells <- head(rownames(all_scores), 3)
    
    specific_cc_scores <- getCellScores(object = object_single,
                                       sigma = available_sigmas[1],
                                       cellType = available_cell_types[1],
                                       ccIndex = 1,
                                       cells = example_cells)
    
    cat("Specific cells, specific CC:", specific_cc_scores, "\n")
  }
}

# =============================================================================
# Example 2: CoProMulti object usage
# =============================================================================

# Assuming you have a CoProMulti object named 'object_multi'
if (exists("object_multi") && is(object_multi, "CoProMulti")) {
  
  # Get available sigma values, cell types, and slides
  available_sigmas <- object_multi@sigmaValues
  available_cell_types <- names(object_multi@cellScores[[1]])
  available_slides <- getSlideList(object_multi)
  
  cat("Available sigma values:", paste(available_sigmas, collapse = ", "), "\n")
  cat("Available cell types:", paste(available_cell_types, collapse = ", "), "\n")
  cat("Available slides:", paste(available_slides, collapse = ", "), "\n")
  
  # Example 2a: Get aggregated cell scores across all slides
  if (length(available_sigmas) > 0 && length(available_cell_types) > 0) {
    scores_aggregated <- getCellScores(object = object_multi,
                                      sigma = available_sigmas[1],
                                      cellType = available_cell_types[1])
    
    cat("Aggregated scores dimensions:", dim(scores_aggregated), "\n")
  }
  
  # Example 2b: Get cell scores for a specific slide
  if (length(available_sigmas) > 0 && length(available_cell_types) > 0 && 
      length(available_slides) > 0) {
    scores_slide <- getCellScores(object = object_multi,
                                 sigma = available_sigmas[1],
                                 cellType = available_cell_types[1],
                                 slide = available_slides[1])
    
    cat("Slide-specific scores dimensions:", dim(scores_slide), "\n")
  }
  
  # Example 2c: Compare slides
  if (length(available_sigmas) > 0 && length(available_cell_types) > 0 && 
      length(available_slides) >= 2) {
    
    slide1_scores <- getCellScores(object = object_multi,
                                  sigma = available_sigmas[1],
                                  cellType = available_cell_types[1],
                                  slide = available_slides[1],
                                  ccIndex = 1,
                                  verbose = FALSE)
    
    slide2_scores <- getCellScores(object = object_multi,
                                  sigma = available_sigmas[1],
                                  cellType = available_cell_types[1],
                                  slide = available_slides[2],
                                  ccIndex = 1,
                                  verbose = FALSE)
    
    cat("Slide 1 CC1 scores summary:", summary(slide1_scores), "\n")
    cat("Slide 2 CC1 scores summary:", summary(slide2_scores), "\n")
  }
}

# =============================================================================
# Example 3: Error handling demonstrations
# =============================================================================

cat("\n=== Error Handling Examples ===\n")

# These examples show the clear error messages provided by getCellScores

if (exists("object_single") && is(object_single, "CoProSingle")) {
  
  # Example 3a: Invalid sigma value
  tryCatch({
    getCellScores(object_single, sigma = 999, cellType = "TypeA")
  }, error = function(e) {
    cat("Error for invalid sigma:", e$message, "\n")
  })
  
  # Example 3b: Invalid cell type
  tryCatch({
    getCellScores(object_single, sigma = object_single@sigmaValues[1], 
                 cellType = "NonExistentType")
  }, error = function(e) {
    cat("Error for invalid cell type:", e$message, "\n")
  })
  
  # Example 3c: Invalid ccIndex
  tryCatch({
    getCellScores(object_single, sigma = object_single@sigmaValues[1],
                 cellType = names(object_single@cellScores[[1]])[1],
                 ccIndex = 999)
  }, error = function(e) {
    cat("Error for invalid ccIndex:", e$message, "\n")
  })
  
  # Example 3d: Invalid cell IDs
  tryCatch({
    getCellScores(object_single, sigma = object_single@sigmaValues[1],
                 cellType = names(object_single@cellScores[[1]])[1],
                 cells = c("fake_cell_1", "fake_cell_2"))
  }, error = function(e) {
    cat("Error for invalid cells:", e$message, "\n")
  })
}

# =============================================================================
# Example 4: Integration with existing CoPro workflows
# =============================================================================

cat("\n=== Integration Examples ===\n")

# Example 4a: Use with correlation plotting functions
if (exists("object_single") && length(object_single@sigmaValues) > 0) {
  
  sigma_val <- object_single@sigmaValues[1]
  cell_types <- names(object_single@cellScores[[1]])
  
  if (length(cell_types) >= 2) {
    # Get cell scores for correlation analysis
    typeA_scores <- getCellScores(object_single, sigma = sigma_val, 
                                 cellType = cell_types[1], ccIndex = 1, 
                                 verbose = FALSE)
    
    typeB_scores <- getCellScores(object_single, sigma = sigma_val, 
                                 cellType = cell_types[2], ccIndex = 1, 
                                 verbose = FALSE)
    
    cat("Ready for correlation analysis between", cell_types[1], "and", cell_types[2], "\n")
    cat("Type A scores range:", range(typeA_scores), "\n")
    cat("Type B scores range:", range(typeB_scores), "\n")
  }
}

# Example 4b: Use with gene expression regression
if (exists("object_single") && length(object_single@sigmaValues) > 0) {
  
  sigma_val <- object_single@sigmaValues[1]
  cell_types <- names(object_single@cellScores[[1]])
  
  if (length(cell_types) > 0) {
    # Get cell scores for regression analysis
    regression_scores <- getCellScores(object_single, sigma = sigma_val,
                                      cellType = cell_types[1], ccIndex = 1,
                                      verbose = FALSE)
    
    cat("Cell scores ready for regression analysis\n")
    cat("Number of cells:", length(regression_scores), "\n")
    cat("Score statistics:", summary(regression_scores), "\n")
  }
}

cat("\n=== getCellScores Examples Complete ===\n")

# =============================================================================
# Performance and Memory Considerations
# =============================================================================

# The getCellScores function is designed for efficiency:
# 1. No unnecessary data copying - returns references when possible
# 2. Lazy evaluation - only applies filters when requested
# 3. Clear error messages to prevent incorrect usage
# 4. Works with both single and multi-slide objects seamlessly

# For large datasets, consider:
# - Using ccIndex parameter to extract only needed components
# - Using cells parameter to extract only specific cells
# - Using slide parameter for multi-slide objects to reduce memory usage 