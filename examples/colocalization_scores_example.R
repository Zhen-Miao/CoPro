#' Example: Computing Colocalization Scores for Cell Type Pairs
#' 
#' This example demonstrates how to use the getColocScores() function to compute
#' cross-type colocalization scores using the inhomogeneous cross pair-correlation
#' function approach.

# Load required libraries
library(CoPro)

# Example workflow for colocalization analysis
example_colocalization_analysis <- function(object) {
  
  cat("=== CoPro Colocalization Analysis Example ===\n")
  
  # Step 1: Basic colocalization analysis with default parameters
  cat("Step 1: Computing colocalization scores with default parameters...\n")
  
  coloc_results_basic <- getColocScores(object, verbose = TRUE)
  
  cat("\nBasic colocalization results:\n")
  print(head(coloc_results_basic))
  
  # Step 2: High-resolution data analysis with custom parameters
  cat("\nStep 2: Computing colocalization scores for high-resolution data...\n")
  
  # Example parameters for high-resolution microscopy data
  coloc_results_hires <- getColocScores(object,
                                       r_um_range = c(5, 30),        # Smaller distance band
                                       pixel_size_um = 0.325,        # Convert pixels to microns
                                       cell_diam_um = 8,              # Smaller cell diameter
                                       nsim = 99,                     # Fewer simulations for speed
                                       r_step_um = 1,                 # Finer distance steps
                                       min_points_per_type = 15,      # Lower threshold
                                       verbose = TRUE)
  
  cat("\nHigh-resolution colocalization results:\n")
  print(head(coloc_results_hires))
  
  # Step 3: Analysis with self-type pairs included
  cat("\nStep 3: Including self-type pairs for spatial clustering analysis...\n")
  
  coloc_results_self <- getColocScores(object,
                                      include_self = TRUE,
                                      verbose = TRUE)
  
  cat("\nResults including self-type pairs:\n")
  self_pairs <- coloc_results_self[coloc_results_self$cellType1 == coloc_results_self$cellType2, ]
  print(self_pairs)
  
  # Step 4: Analysis and interpretation
  cat("\nStep 4: Analyzing colocalization patterns...\n")
  
  # Identify strongly colocalized pairs
  strong_coloc <- coloc_results_basic[coloc_results_basic$colocScore > 1, ]
  if (nrow(strong_coloc) > 0) {
    cat("Strongly colocalized pairs (score > 1):\n")
    for (i in seq_len(nrow(strong_coloc))) {
      row <- strong_coloc[i, ]
      cat(sprintf("  %s <-> %s: Score = %.3f (n1=%d, n2=%d)\n", 
                 row$cellType1, row$cellType2, row$colocScore, row$nCells1, row$nCells2))
    }
  } else {
    cat("No strongly colocalized pairs found (score > 1)\n")
  }
  
  # Identify segregated pairs
  segregated <- coloc_results_basic[coloc_results_basic$colocScore < -1, ]
  if (nrow(segregated) > 0) {
    cat("\nSegregated pairs (score < -1):\n")
    for (i in seq_len(nrow(segregated))) {
      row <- segregated[i, ]
      cat(sprintf("  %s <-> %s: Score = %.3f (n1=%d, n2=%d)\n", 
                 row$cellType1, row$cellType2, row$colocScore, row$nCells1, row$nCells2))
    }
  } else {
    cat("No strongly segregated pairs found (score < -1)\n")
  }
  
  # Summary statistics
  valid_scores <- coloc_results_basic[!is.na(coloc_results_basic$colocScore), ]
  if (nrow(valid_scores) > 0) {
    cat("\nSummary statistics:\n")
    cat(sprintf("  Mean colocalization score: %.3f\n", mean(valid_scores$colocScore)))
    cat(sprintf("  Score range: %.3f to %.3f\n", min(valid_scores$colocScore), max(valid_scores$colocScore)))
    cat(sprintf("  Number of valid pairs: %d\n", nrow(valid_scores)))
    cat(sprintf("  Number of failed pairs: %d\n", sum(is.na(coloc_results_basic$colocScore))))
  }
  
  cat("\n=== Analysis completed successfully! ===\n")
  
  # Return results for further analysis
  return(list(
    basic = coloc_results_basic,
    hires = coloc_results_hires,
    with_self = coloc_results_self
  ))
}

# Multi-slide example
example_multi_slide_colocalization <- function(multi_object) {
  
  cat("=== Multi-slide Colocalization Analysis Example ===\n")
  
  # Compute colocalization scores for multi-slide object
  coloc_results_multi <- getColocScores(multi_object, verbose = TRUE)
  
  # Analyze results by slide
  slides <- unique(coloc_results_multi$slideID)
  
  cat("\nPer-slide analysis:\n")
  for (slide in slides[1:min(3, length(slides))]) {  # Show first 3 slides
    slide_data <- coloc_results_multi[coloc_results_multi$slideID == slide, ]
    valid_data <- slide_data[!is.na(slide_data$colocScore), ]
    
    if (nrow(valid_data) > 0) {
      cat(sprintf("Slide %s: %d pairs, mean score = %.3f\n", 
                 slide, nrow(valid_data), mean(valid_data$colocScore)))
      
      # Show top colocalized pair for this slide
      top_pair <- valid_data[which.max(valid_data$colocScore), ]
      cat(sprintf("  Top pair: %s <-> %s (score = %.3f)\n", 
                 top_pair$cellType1, top_pair$cellType2, top_pair$colocScore))
    } else {
      cat(sprintf("Slide %s: No valid colocalization scores\n", slide))
    }
  }
  
  # Cross-slide comparison
  cat("\nCross-slide comparison:\n")
  if (length(slides) > 1) {
    # Compare same cell type pairs across slides
    unique_pairs <- unique(paste(coloc_results_multi$cellType1, coloc_results_multi$cellType2, sep = " vs "))
    
    for (pair in unique_pairs[1:min(3, length(unique_pairs))]) {  # Show first 3 pairs
      pair_data <- coloc_results_multi[paste(coloc_results_multi$cellType1, coloc_results_multi$cellType2, sep = " vs ") == pair, ]
      valid_pair_data <- pair_data[!is.na(pair_data$colocScore), ]
      
      if (nrow(valid_pair_data) > 1) {
        cat(sprintf("%s: %.3f ± %.3f (n=%d slides)\n", 
                   pair, mean(valid_pair_data$colocScore), sd(valid_pair_data$colocScore), nrow(valid_pair_data)))
      }
    }
  }
  
  return(coloc_results_multi)
}

# Utility function to visualize colocalization results
visualize_colocalization_results <- function(coloc_results) {
  
  cat("=== Colocalization Results Visualization ===\n")
  
  valid_results <- coloc_results[!is.na(coloc_results$colocScore), ]
  
  if (nrow(valid_results) == 0) {
    cat("No valid colocalization scores to visualize.\n")
    return(invisible(NULL))
  }
  
  # Summary by cell type pairs
  cat("Colocalization scores by cell type pair:\n")
  for (i in seq_len(nrow(valid_results))) {
    row <- valid_results[i, ]
    status <- if (row$colocScore > 1) "Strong Coloc" else if (row$colocScore > 0) "Weak Coloc" else if (row$colocScore > -1) "Weak Segregation" else "Strong Segregation"
    
    slide_info <- if ("slideID" %in% colnames(row)) paste0(" (", row$slideID, ")") else ""
    
    cat(sprintf("  %s <-> %s%s: %.3f [%s]\n", 
               row$cellType1, row$cellType2, slide_info, row$colocScore, status))
  }
  
  # Distribution summary
  cat("\nScore distribution:\n")
  cat(sprintf("  Strongly colocalized (>1): %d pairs\n", sum(valid_results$colocScore > 1)))
  cat(sprintf("  Weakly colocalized (0-1): %d pairs\n", sum(valid_results$colocScore > 0 & valid_results$colocScore <= 1)))
  cat(sprintf("  Weakly segregated (-1-0): %d pairs\n", sum(valid_results$colocScore >= -1 & valid_results$colocScore <= 0)))
  cat(sprintf("  Strongly segregated (<-1): %d pairs\n", sum(valid_results$colocScore < -1)))
  
  # Cell count analysis
  if (all(c("nCells1", "nCells2") %in% colnames(valid_results))) {
    cat("\nCell count analysis:\n")
    valid_results$min_cells <- pmin(valid_results$nCells1, valid_results$nCells2)
    valid_results$max_cells <- pmax(valid_results$nCells1, valid_results$nCells2)
    
    cat(sprintf("  Cell count range: %d - %d per type\n", 
               min(c(valid_results$nCells1, valid_results$nCells2)), 
               max(c(valid_results$nCells1, valid_results$nCells2))))
    
    # Check if colocalization correlates with cell density
    low_density <- valid_results[valid_results$min_cells < 50, ]
    high_density <- valid_results[valid_results$min_cells >= 50, ]
    
    if (nrow(low_density) > 0 && nrow(high_density) > 0) {
      cat(sprintf("  Low density pairs (<50 cells): mean score = %.3f\n", mean(low_density$colocScore)))
      cat(sprintf("  High density pairs (≥50 cells): mean score = %.3f\n", mean(high_density$colocScore)))
    }
  }
}

# Usage notes and interpretation guide
cat("=== Colocalization Analysis Usage Notes ===\n")
cat("\nKey Parameters:\n")
cat("- r_um_range: Distance band for analysis (e.g., c(10, 60) for 10-60 microns)\n")
cat("- pixel_size_um: Conversion factor from coordinates to microns\n")
cat("- cell_diam_um: Typical cell diameter for intensity smoothing\n")
cat("- nsim: Number of random simulations (more = better null, but slower)\n")
cat("- min_points_per_type: Minimum cells needed for reliable analysis\n")

cat("\nScore Interpretation:\n")
cat("- Score > 1: Strong colocalization (cells cluster together)\n")
cat("- Score 0-1: Weak colocalization\n")
cat("- Score -1-0: Weak segregation\n")
cat("- Score < -1: Strong segregation (cells avoid each other)\n")
cat("- Score ≈ 0: Random spatial arrangement\n")

cat("\nTips:\n")
cat("- Use smaller r_um_range for fine-scale interactions\n")
cat("- Adjust pixel_size_um to convert coordinates to microns\n")
cat("- Include include_self=TRUE to analyze spatial clustering within cell types\n")
cat("- For multi-slide objects, compare scores across slides for consistency\n")

cat("\nColocalization functions loaded successfully!\n")
cat("Use example_colocalization_analysis(object) to run the example.\n")
