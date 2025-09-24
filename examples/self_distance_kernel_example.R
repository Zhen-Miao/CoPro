#' Example: Computing Self-Distance and Self-Kernel Matrices for Multiple Cell Types
#' 
#' This example demonstrates how to use the new computeSelfDistance() and 
#' computeSelfKernel() functions to compute within-cell-type spatial relationships
#' when multiple cell types are present.

# Load the CoPro package
library(CoPro)

# Example workflow with multiple cell types
example_self_distance_kernel <- function(object) {
  
  # Assume 'object' is a CoPro object with multiple cell types
  # and has been processed through the standard workflow
  
  cat("=== CoPro Self-Distance and Self-Kernel Example ===\n")
  
  # Step 1: Compute standard cross-type distances (if not already done)
  cat("Step 1: Computing cross-type distances...\n")
  object <- computeDistance(object, 
                           distType = "Euclidean2D",
                           normalizeDistance = TRUE,
                           verbose = TRUE)
  
  # Step 2: Compute cross-type kernel matrices (if not already done)
  cat("\nStep 2: Computing cross-type kernel matrices...\n")
  sigma_values <- c(0.01, 0.05, 0.1)
  object <- computeKernelMatrix(object, 
                               sigmaValues = sigma_values,
                               verbose = TRUE)
  
  # Step 3: NEW - Compute self-distances for each cell type
  cat("\nStep 3: Computing self-distance matrices...\n")
  object <- computeSelfDistance(object,
                               distType = "Euclidean2D",
                               normalizeDistance = TRUE,
                               verbose = TRUE,
                               overwrite = FALSE)  # Add to existing distances
  
  # Step 4: NEW - Compute self-kernel matrices for each cell type
  cat("\nStep 4: Computing self-kernel matrices...\n")
  object <- computeSelfKernel(object,
                             sigmaValues = sigma_values,
                             verbose = TRUE,
                             overwrite = FALSE)  # Add to existing kernels
  
  # Step 5: Access the computed matrices
  cat("\nStep 5: Accessing self-distance and self-kernel matrices...\n")
  
  cell_types <- object@cellTypesOfInterest
  
  for (ct in cell_types) {
    cat("\nCell type:", ct, "\n")
    
    # Get self-distance matrix
    self_dist <- getSelfDistMat(object, cellType = ct)
    if (!is.null(self_dist)) {
      cat("  Self-distance matrix:", nrow(self_dist), "x", ncol(self_dist), "\n")
      cat("  Distance range:", range(self_dist[is.finite(self_dist)]), "\n")
    }
    
    # Get self-kernel matrices for different sigma values
    for (sigma in sigma_values) {
      self_kernel <- getSelfKernelMatrix(object, sigma = sigma, cellType = ct)
      if (!is.null(self_kernel)) {
        cat("  Self-kernel matrix (sigma =", sigma, "):", 
            nrow(self_kernel), "x", ncol(self_kernel), "\n")
        cat("    Kernel range:", range(self_kernel), "\n")
        cat("    Non-zero proportion:", mean(self_kernel > 1e-7), "\n")
      }
    }
  }
  
  # Step 6: Demonstrate difference between cross-type and self-type kernels
  cat("\nStep 6: Comparing cross-type vs self-type kernel structures...\n")
  
  if (length(cell_types) >= 2) {
    ct1 <- cell_types[1]
    ct2 <- cell_types[2]
    sigma <- sigma_values[1]
    
    # Cross-type kernel
    cross_kernel <- getKernelMatrix(object, sigma = sigma, 
                                   cellType1 = ct1, cellType2 = ct2)
    
    # Self-type kernels
    self_kernel1 <- getSelfKernelMatrix(object, sigma = sigma, cellType = ct1)
    self_kernel2 <- getSelfKernelMatrix(object, sigma = sigma, cellType = ct2)
    
    cat("Cross-type kernel (", ct1, "->", ct2, "):", 
        nrow(cross_kernel), "x", ncol(cross_kernel), "\n")
    cat("Self-type kernel (", ct1, "):", 
        nrow(self_kernel1), "x", ncol(self_kernel1), "\n")
    cat("Self-type kernel (", ct2, "):", 
        nrow(self_kernel2), "x", ncol(self_kernel2), "\n")
  }
  
  cat("\n=== Example completed successfully! ===\n")
  return(object)
}

# Usage notes:
# 
# The new functions extend CoPro's capabilities by allowing computation of
# within-cell-type spatial relationships even when multiple cell types are present.
# 
# Key benefits:
# 1. Standard CoPro workflow computes cross-type kernels for multiple cell types
# 2. New functions add self-type kernels for spatial autocorrelation analysis
# 3. Both kernel types can be used together for comprehensive spatial analysis
# 
# Workflow:
# 1. computeDistance() - computes cross-type distances
# 2. computeKernelMatrix() - computes cross-type kernels
# 3. computeSelfDistance() - adds self-type distances
# 4. computeSelfKernel() - adds self-type kernels
# 5. Use getKernelMatrix() for cross-type kernels
# 6. Use getSelfKernelMatrix() for self-type kernels

# For multi-slide objects, the functions work similarly:
example_multi_slide <- function(multi_object) {
  # All functions support multi-slide objects automatically
  multi_object <- computeSelfDistance(multi_object)
  multi_object <- computeSelfKernel(multi_object, sigmaValues = c(0.01, 0.05, 0.1))
  
  # Access self-kernels for specific slides
  slides <- getSlideList(multi_object)
  for (slide in slides) {
    for (ct in multi_object@cellTypesOfInterest) {
      self_kernel <- getSelfKernelMatrix(multi_object, 
                                        sigma = 0.05, 
                                        cellType = ct, 
                                        slide = slide)
      # Use the self-kernel for slide-specific analysis
    }
  }
  
  return(multi_object)
}
