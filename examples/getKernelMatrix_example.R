# Example: Using getKernelMatrix function
# 
# This example demonstrates how to use the new getKernelMatrix accessor function
# instead of manually navigating nested list structures

library(CoPro)

# Example usage with CoProSingle object
# Assume you have a CoPro object called 'copro_obj'

# BEFORE: Complex nested access 
# sigma_name <- paste("sigma", 0.1, sep = "_")
# K_ij <- copro_obj@kernelMatrices[[sigma_name]][["TypeA"]][["TypeB"]]
# if (is.null(K_ij) || length(K_ij) == 0) {
#   K_ji <- copro_obj@kernelMatrices[[sigma_name]][["TypeB"]][["TypeA"]]
#   K_ij <- t(K_ji)
# }

# AFTER: Clean accessor function
# K_ij <- getKernelMatrix(copro_obj, sigma = 0.1, 
#                        cellType1 = "TypeA", cellType2 = "TypeB")

# ============================================================================
# Function Examples
# ============================================================================

# Example 1: Basic usage for CoProSingle
get_kernel_single_example <- function(copro_single) {
  # Get kernel matrix between two different cell types
  K_AB <- getKernelMatrix(copro_single, sigma = 0.1, 
                         cellType1 = "TypeA", cellType2 = "TypeB")
  
  # Get within-cell-type kernel (diagonal)
  K_AA <- getKernelMatrix(copro_single, sigma = 0.1, 
                         cellType1 = "TypeA", cellType2 = "TypeA")
  
  return(list(K_AB = K_AB, K_AA = K_AA))
}

# Example 2: Basic usage for CoProMulti  
get_kernel_multi_example <- function(copro_multi) {
  # Get kernel matrix for specific slide
  K_AB_slide1 <- getKernelMatrix(copro_multi, sigma = 0.5, 
                                cellType1 = "TypeA", cellType2 = "TypeB",
                                slide = "slide1")
  
  # Get within-cell-type kernel for specific slide
  K_AA_slide2 <- getKernelMatrix(copro_multi, sigma = 0.5,
                                cellType1 = "TypeA", cellType2 = "TypeA", 
                                slide = "slide2")
  
  return(list(K_AB_slide1 = K_AB_slide1, K_AA_slide2 = K_AA_slide2))
}

# Example 3: Error handling and validation
get_kernel_with_error_handling <- function(copro_obj) {
  tryCatch({
    # This will provide helpful error messages if something is wrong
    K <- getKernelMatrix(copro_obj, sigma = 0.1, 
                        cellType1 = "NonExistentType", cellType2 = "TypeB")
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    return(NULL)
  })
}

# Example 4: Check what's available
check_available_kernels <- function(copro_obj) {
  # See all available kernel combinations
  available <- getAvailableKernelMatrices(copro_obj)
  print(available)
  
  # Filter by specific sigma
  available_sigma <- getAvailableKernelMatrices(copro_obj, sigma = 0.1)
  print(available_sigma)
}

# Example 5: Using in optimization functions (replaces old get_kernel_matrix)
optimize_with_new_accessor <- function(copro_obj, X_list, sigma, cell_types) {
  # Instead of passing K_list, we can get kernels on demand
  for (ct_i in cell_types) {
    for (ct_j in cell_types) {
      if (ct_i == ct_j) next
      
      # Clean access to kernel matrix
      K_ij <- getKernelMatrix(copro_obj, sigma = sigma, 
                             cellType1 = ct_i, cellType2 = ct_j,
                             verbose = FALSE)  # Turn off verbose for optimization
      
      # Use K_ij in calculations...
      X_i <- X_list[[ct_i]]
      X_j <- X_list[[ct_j]]
      # ... optimization calculations ...
    }
  }
}

# Example 6: Batch processing across slides (for CoProMulti)
process_all_slides <- function(copro_multi, sigma, cellType1, cellType2) {
  results <- list()
  
  for (slide in getSlideList(copro_multi)) {
    tryCatch({
      K <- getKernelMatrix(copro_multi, sigma = sigma,
                          cellType1 = cellType1, cellType2 = cellType2,
                          slide = slide, verbose = FALSE)
      
      # Do some processing with K
      results[[slide]] <- list(
        slide = slide,
        kernel_dims = dim(K),
        kernel_sparsity = sum(K > 1e-7) / length(K)
      )
    }, error = function(e) {
      results[[slide]] <- list(slide = slide, error = e$message)
    })
  }
  
  return(results)
}

# ============================================================================
# Benefits Summary
# ============================================================================

cat("Benefits of using getKernelMatrix():\n")
cat("1. Unified interface for single and multi-slide objects\n")
cat("2. Automatic symmetric matrix handling (K_ij or t(K_ji))\n") 
cat("3. Comprehensive error messages and validation\n")
cat("4. No need to remember complex nested structure\n")
cat("5. Consistent parameter names across the package\n")
cat("6. Easy to debug and maintain\n")
cat("7. Future-proof against internal structure changes\n") 