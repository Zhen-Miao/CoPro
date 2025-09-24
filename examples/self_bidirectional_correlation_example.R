#' Example: Computing Self-Bidirectional Correlation for Multiple Cell Types
#' 
#' This example demonstrates how to use the new self-bidirectional correlation
#' functions to compute spatial autocorrelation within each cell type.

# Load the CoPro package
library(CoPro)

# Example workflow with self-bidirectional correlation
example_self_bidirectional_correlation <- function(object) {
  
  cat("=== CoPro Self-Bidirectional Correlation Example ===\n")
  
  # Step 1: Standard CoPro workflow
  cat("Step 1: Running standard CoPro workflow...\n")
  object <- computeDistance(object, distType = "Euclidean2D", verbose = TRUE)
  object <- computeKernelMatrix(object, sigmaValues = c(0.01, 0.05, 0.1), verbose = TRUE)
  object <- computePCA(object, verbose = TRUE)
  object <- runSkrCCA(object, verbose = TRUE)
  
  # Step 2: Add self-distance and self-kernel computation
  cat("\nStep 2: Computing self-distances and self-kernels...\n")
  object <- computeSelfDistance(object, verbose = TRUE)
  object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1), verbose = TRUE)
  
  # Step 3: Compute self-bidirectional correlation using native skrCCA results
  cat("\nStep 3: Computing self-bidirectional correlation (native method)...\n")
  self_bidir_native <- computeSelfBidirCorr(object, sigma_choice = 0.05, verbose = TRUE)
  
  # Display results
  cat("\nSelf-bidirectional correlation results (native method):\n")
  print(head(self_bidir_native$sigma_0.05))
  
  # Step 4: Demonstrate transfer-based self-bidirectional correlation
  cat("\nStep 4: Computing self-bidirectional correlation (transfer method)...\n")
  
  # For transfer method, we need cell scores in the right format
  # This simulates transferring scores from the same object to itself
  transfer_scores <- getTransferCellScores(ref_obj = object, 
                                          tar_obj = object, 
                                          sigma_choice = 0.05,
                                          agg_cell_type = FALSE,
                                          verbose = TRUE)
  
  # Compute self-bidirectional correlation using transferred scores
  self_bidir_transfer <- getTransferSelfBidirCorr(tar_obj = object,
                                                 transfer_cell_scores = transfer_scores,
                                                 sigma_choice = 0.05,
                                                 verbose = TRUE)
  
  cat("\nSelf-bidirectional correlation results (transfer method):\n")
  print(head(self_bidir_transfer$sigma_0.05))
  
  # Step 5: Compare with cross-type bidirectional correlation
  cat("\nStep 5: Comparing self-type vs cross-type bidirectional correlation...\n")
  
  # Compute standard cross-type bidirectional correlation
  cross_bidir <- getTransferBidirCorr(tar_obj = object,
                                     transfer_cell_scores = transfer_scores,
                                     sigma_choice = 0.05,
                                     verbose = TRUE)
  
  cat("\nCross-type bidirectional correlation results:\n")
  print(head(cross_bidir$sigma_0.05))
  
  # Step 6: Analysis and interpretation
  cat("\nStep 6: Analysis summary...\n")
  
  self_results <- self_bidir_native$sigma_0.05
  cross_results <- cross_bidir$sigma_0.05
  
  cat("Self-type correlations (spatial autocorrelation within cell types):\n")
  for (ct in unique(self_results$cellType)) {
    ct_results <- self_results[self_results$cellType == ct, ]
    avg_corr <- mean(ct_results$selfBidirCorrelation, na.rm = TRUE)
    cat(sprintf("  %s: Average self-correlation = %.3f\n", ct, avg_corr))
  }
  
  cat("\nCross-type correlations (spatial correlation between cell types):\n")
  for (i in seq_len(nrow(cross_results))) {
    row <- cross_results[i, ]
    if (row$CC_index == 1) {  # Show only first CC for brevity
      cat(sprintf("  %s <-> %s: Cross-correlation = %.3f\n", 
                 row$cellType1, row$cellType2, row$bidirCorrelation))
    }
  }
  
  cat("\n=== Analysis completed successfully! ===\n")
  
  # Return results for further analysis
  return(list(
    object = object,
    self_bidir_native = self_bidir_native,
    self_bidir_transfer = self_bidir_transfer,
    cross_bidir = cross_bidir
  ))
}

# Multi-slide example
example_multi_slide_self_bidir <- function(multi_object) {
  
  cat("=== Multi-slide Self-Bidirectional Correlation Example ===\n")
  
  # Compute self-kernels for multi-slide object
  multi_object <- computeSelfDistance(multi_object, verbose = TRUE)
  multi_object <- computeSelfKernel(multi_object, sigmaValues = c(0.01, 0.05, 0.1), verbose = TRUE)
  
  # Per-slide analysis
  cat("\nComputing per-slide self-bidirectional correlation...\n")
  self_bidir_per_slide <- computeSelfBidirCorr(multi_object, 
                                              sigma_choice = 0.05,
                                              calculationMode = "perSlide",
                                              verbose = TRUE)
  
  # Aggregate analysis
  cat("\nComputing aggregate self-bidirectional correlation...\n")
  self_bidir_aggregate <- computeSelfBidirCorr(multi_object, 
                                              sigma_choice = 0.05,
                                              calculationMode = "aggregate",
                                              verbose = TRUE)
  
  # Display results
  slides <- getSlideList(multi_object)
  cat("\nPer-slide results:\n")
  for (slide in slides[1:min(3, length(slides))]) {  # Show first 3 slides
    slide_results <- self_bidir_per_slide$sigma_0.05[
      self_bidir_per_slide$sigma_0.05$slideID == slide, ]
    cat(sprintf("Slide %s: %d correlations computed\n", slide, nrow(slide_results)))
  }
  
  cat("\nAggregate results:\n")
  print(head(self_bidir_aggregate$sigma_0.05))
  
  return(list(
    per_slide = self_bidir_per_slide,
    aggregate = self_bidir_aggregate
  ))
}

# Utility function to visualize results
visualize_self_bidir_results <- function(results) {
  
  cat("=== Self-Bidirectional Correlation Visualization ===\n")
  
  if ("sigma_0.05" %in% names(results)) {
    df <- results$sigma_0.05
    
    # Summary statistics
    if ("cellType" %in% colnames(df)) {
      # Self-correlation results
      cat("Summary of self-bidirectional correlations by cell type:\n")
      
      for (ct in unique(df$cellType)) {
        ct_data <- df[df$cellType == ct, ]
        corr_col <- if ("selfBidirCorrelation" %in% colnames(ct_data)) {
          "selfBidirCorrelation"
        } else {
          "aggregateSelfCorrelation"
        }
        
        avg_corr <- mean(ct_data[[corr_col]], na.rm = TRUE)
        sd_corr <- sd(ct_data[[corr_col]], na.rm = TRUE)
        
        cat(sprintf("  %s: Mean = %.3f, SD = %.3f, N = %d\n", 
                   ct, avg_corr, sd_corr, nrow(ct_data)))
      }
      
      # CC-specific analysis
      if ("CC_index" %in% colnames(df)) {
        cat("\nCorrelation by canonical component:\n")
        for (cc in unique(df$CC_index)) {
          cc_data <- df[df$CC_index == cc, ]
          corr_col <- if ("selfBidirCorrelation" %in% colnames(cc_data)) {
            "selfBidirCorrelation"
          } else {
            "aggregateSelfCorrelation"
          }
          avg_corr <- mean(cc_data[[corr_col]], na.rm = TRUE)
          cat(sprintf("  CC_%d: Mean correlation = %.3f\n", cc, avg_corr))
        }
      }
    }
  }
}

# Usage notes:
#
# Key differences between self-type and cross-type bidirectional correlation:
#
# 1. Self-type (spatial autocorrelation):
#    - Measures spatial clustering within the same cell type
#    - Uses self-kernel matrices (computed with computeSelfKernel)
#    - Correlates: cor(t(K) %*% A_w, A_w) and cor(A_w, K %*% A_w)
#    - Indicates how spatially organized cells of the same type are
#
# 2. Cross-type (spatial cross-correlation):
#    - Measures spatial relationships between different cell types
#    - Uses cross-kernel matrices (standard CoPro kernels)
#    - Correlates: cor(t(K) %*% A_w1, B_w2) and cor(A_w1, K %*% B_w2)
#    - Indicates how cell types spatially relate to each other
#
# Applications:
# - Self-type: Detect spatial clustering, identify organized tissue regions
# - Cross-type: Find cell type interactions, spatial niches, communication patterns
# - Combined: Comprehensive spatial analysis of tissue organization

cat("Self-bidirectional correlation functions loaded successfully!\n")
cat("Use example_self_bidirectional_correlation(object) to run the example.\n")
