#' Calculate kernel matrix from distance matrix
#' The function runs the following calculation:
#' \deqn{K(x, y) = \exp\left(-\frac{\|x-y\|^2}{2 \sigma^2}\right)}. We notice
#' that the normalization factor does not affect the final results as it is
#' scale invariant, so here for easy computation we omit the scaling factor.
#'
#' @param sigma The variance parameter \eqn{\sigma}, a positive number.
#' @param dist_mat A numeric matrix representing the squared distances
#'  between cells
#' @param lower_limit A lower limit value below which the kernel value will
#'  be set to zero, default = 1e-7
#'
#' @return a matrix of the same dimensions as \code{dist_mat},
#'  containing the calculated Gaussian kernel values.
#' @noRd
kernel_from_distance <- function(
    sigma, dist_mat, lower_limit = 1e-7) {
  kernel_mat <- exp(-0.5 * (dist_mat / sigma)^2)
  kernel_mat[kernel_mat < lower_limit] <- 0
  return(kernel_mat)
}

#' Check if sigma value should be removed
#'
#' This function checks if a sigma value should be removed based on the kernel matrix.
#'
#' @param kernel_current The kernel matrix.
#' @param lowerLimit The lower limit for the kernel function.
#' @param sigma_choose The sigma value to check.
#' @param sigmaValues The vector of sigma values.
#' @param i The first cell type.
#' @param j The second cell type.
#' @param minAveCellNeighor The minimum average number of cell in the
#'  neighbor.
#' @return A logical value indicating if the sigma value should be removed.
#' @noRd
.CheckSigmaValuesToRemove <- function(kernel_current, lowerLimit,
                                      sigma_choose, sigmaValues, i, j,
                                      minAveCellNeighor) {
  n_cell1 <- nrow(kernel_current)
  n_cell2 <- ncol(kernel_current)

  minPropZero <- minAveCellNeighor * min(n_cell1, n_cell2) / (n_cell1 * n_cell2)
  if (mean(kernel_current > lowerLimit) < minPropZero) {
    warning(paste("Kernel matrix for cell types", i, "and", j,
                  "with sigma =", sigma_choose,
                  "contains too many zeros. Specifically, less than",
                  minPropZero * 100, "% total counts are above the threshold" ))
    if (length(sigmaValues) == 1) {
      stop(paste("Only one sigma value is specified,",
                 "which resulted in all Gaussian kernel being small.",
                 "Please provide a larger sigma value"))
    }else {
      warning(paste("Dropping sigma value of ",
                    sigma_choose,
                    "because all Gaussian kernel values are too small,",
                    "which will not produce meaningful results."))
      return(TRUE)
    }
  }else if (all(is.na(kernel_current))) {
    warning(paste("Kernel matrix for cell types", i, "and", j,
                  "with sigma =", sigma_choose,
                  "contains all NA."))
    if (length(sigmaValues) == 1) {
      stop(paste("Only one sigma value is specified,",
                 "which resulted in all Gaussian kernel being NA."))
    }else {
      warning(paste("Dropping sigma value of ",
                    sigma_choose,
                    "because all Gaussian kernel values are NA."))
      return(TRUE)
    }
  }

  return(FALSE)
}

# Helper function to validate kernel computation inputs
.checkInputKernel <- function(object, sigmaValues, lowerLimit, upperQuantile, 
                             normalizeKernel, minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel) {
  # Check if distance matrix exists
  if (length(object@distances) == 0) {
    stop("Please run computeDistance before computing kernel")
  }

  # Check for conflicting normalization options
  if (rowNormalizeKernel && colNormalizeKernel) {
    stop("Cannot do both row-wise and column-wise normalization.")
  }

  # Get cell types of interest
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) {
    stop("No cell types of interest specified")
  }

  # Set default sigma values if not provided
  if (length(sigmaValues) == 0) {
    warning("No Sigma specified, setting to the 5% quantile of cell distance")
    
    # Find first available distance matrix using flat structure
    dist_mat <- NULL
    if (length(cts) == 1) {
      dist_flat_name <- .createDistMatrixName(cts, cts, slide = NULL)
      if (dist_flat_name %in% names(object@distances)) {
        dist_mat <- object@distances[[dist_flat_name]]
      }
    } else {
      dist_flat_name <- .createDistMatrixName(cts[1], cts[2], slide = NULL)
      if (dist_flat_name %in% names(object@distances)) {
        dist_mat <- object@distances[[dist_flat_name]]
      }
    }
    
    if (is.null(dist_mat)) {
      stop("No distance matrices found. Please run computeDistance first.")
    }
    
    sigmaValues <- quantile(dist_mat[dist_mat > 0], 0.05)
  }

  # Validate sigma values
  if (!is.numeric(sigmaValues) || any(sigmaValues <= 0)) {
    stop("sigmaValues must be positive numeric values")
  }

  # Validate other parameters
  if (lowerLimit <= 0 || lowerLimit >= 1) {
    stop("lowerLimit must be between 0 and 1")
  }
  
  if (upperQuantile <= 0 || upperQuantile >= 1) {
    stop("upperQuantile must be between 0 and 1")
  }
  
  if (minAveCellNeighor < 1) {
    stop("minAveCellNeighor must be at least 1")
  }

  return(list(cts = cts, sigmaValues = sigmaValues))
}

# Helper function to initialize kernel matrix structure (now flat)
.initializeKernelStructure <- function(sigmaValues, cts) {
  # Initialize flat structure with informative names
  kernel_mat <- list()
  
  # Pre-allocate for all sigma-celltype pairs
  for (sigma_val in sigmaValues) {
    # Add between-celltype pairs if we have multiple cell types
    if (length(cts) > 1) {
      pair_cell_types <- combn(cts, 2)
      for (pp in seq_len(ncol(pair_cell_types))) {
        ct1 <- pair_cell_types[1, pp]
        ct2 <- pair_cell_types[2, pp]
        flat_name <- .createKernelMatrixName(sigma_val, ct1, ct2, slide = NULL)
        kernel_mat[[flat_name]] <- NULL  # Placeholder
      }
    }
    
    # Add within-celltype pairs for all cell types
    for (ct in cts) {
      flat_name <- .createKernelMatrixName(sigma_val, ct, ct, slide = NULL)
      kernel_mat[[flat_name]] <- NULL  # Placeholder
    }
  }
  
  return(kernel_mat)
}

# Helper function to process kernel matrix (clipping and normalization)
.processKernelMatrix <- function(kernel_current, lowerLimit, upperQuantile, 
                                normalizeKernel, rowNormalizeKernel, colNormalizeKernel) {
  
  # Check for empty or invalid kernel matrix
  if (nrow(kernel_current) == 0 || ncol(kernel_current) == 0) {
    stop("Cannot process empty kernel matrix")
  }
  
  # Check if there are any valid kernel values
  valid_kernel_values <- kernel_current[kernel_current >= lowerLimit & !is.na(kernel_current)]
  if (length(valid_kernel_values) == 0) {
    warning("No valid kernel values found above lowerLimit")
    kernel_current[!is.na(kernel_current)] <- 0
    return(kernel_current)
  }
  
  # Clipping large values
  upper_clip <- quantile(valid_kernel_values, upperQuantile, na.rm = TRUE)
  kernel_current[kernel_current >= upper_clip & !is.na(kernel_current)] <- upper_clip
  
  # Apply normalization
  if ((normalizeKernel && !rowNormalizeKernel) && !colNormalizeKernel) {
    # Global normalization
    rs_kernel <- rowSums(kernel_current, na.rm = TRUE)
    median_rs <- median(rs_kernel[rs_kernel > 1e-5], na.rm = TRUE)
    if (!is.na(median_rs) && median_rs > 0) {
      kernel_current <- kernel_current / median_rs
    }
    
  } else if (rowNormalizeKernel) {
    # Row-wise normalization
    rs_kernel <- rowSums(kernel_current, na.rm = TRUE)
    nz_ind <- rs_kernel > 1e-4
    if (any(nz_ind)) {
      kernel_current[nz_ind, ] <- kernel_current[nz_ind, ] / rs_kernel[nz_ind]
    }
    
  } else if (colNormalizeKernel) {
    # Column-wise normalization
    kernel_current <- t(kernel_current)
    rs_kernel <- rowSums(kernel_current, na.rm = TRUE)
    nz_ind <- rs_kernel > 1e-4
    if (any(nz_ind)) {
      kernel_current[nz_ind, ] <- kernel_current[nz_ind, ] / rs_kernel[nz_ind]
    }
    kernel_current <- t(kernel_current)
  }
  
  # Remove small values
  kernel_current[kernel_current < lowerLimit & !is.na(kernel_current)] <- 0
  
  return(kernel_current)
}

# Helper function to clean up removed sigma values (flat structure)
.cleanupSigmaValues <- function(object, kernel_mat, sigmaValuesToRemove, verbose = TRUE) {
  if (any(sigmaValuesToRemove)) {
    if (verbose) {
      cat("removing", sum(sigmaValuesToRemove), "sigma values", "\n")
    }
    
    # Get sigma values that should be removed
    sigmas_to_remove <- as.numeric(gsub("sigma_", "", names(sigmaValuesToRemove)[sigmaValuesToRemove]))
    
    # Remove flat entries that contain these sigma values
    to_remove <- sapply(names(kernel_mat), function(name) {
      parsed <- .parseKernelMatrixName(name)
      parsed$sigma %in% sigmas_to_remove
    })
    
    if (any(to_remove)) {
      kernel_mat <- kernel_mat[!to_remove]
    }
    
    object@sigmaValues <- object@sigmaValues[!sigmaValuesToRemove]
  }
  
  object@kernelMatrices <- kernel_mat
  return(object)
}

# Core dispatcher function
.computeKernelCore <- function(object, sigmaValues, lowerLimit, upperQuantile, 
                              normalizeKernel, minAveCellNeighor, rowNormalizeKernel, 
                              colNormalizeKernel, verbose) {
  
  validated_input <- .checkInputKernel(object, sigmaValues, lowerLimit, upperQuantile,
                                      normalizeKernel, minAveCellNeighor, 
                                      rowNormalizeKernel, colNormalizeKernel)
  cts <- validated_input$cts
  sigmaValues <- validated_input$sigmaValues
  
  # Save sigma values
  object@sigmaValues <- sigmaValues
  
  # Provide user feedback
  if (verbose) {
    if (length(cts) == 1) {
      cat("Computing kernel matrix for one cell type\n")
    } else {
      cat("Computing pairwise kernel matrix for", length(cts), "cell types\n")
    }
    
    if (normalizeKernel) {
      cat("normalizeKernel is set to TRUE. Kernel matrix will be",
          "normalized so that median row sums of kernel will be 1\n")
    }
  }
  
  # Determine computation approach
  if (length(cts) == 1) {
    return(.computeKernelWithin(object, cts, sigmaValues, lowerLimit, upperQuantile,
                               normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                               colNormalizeKernel, verbose))
  } else {
    return(.computeKernelPairs(object, cts, sigmaValues, lowerLimit, upperQuantile,
                              normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                              colNormalizeKernel, verbose))
  }
}

# Function for computing kernel matrices between pairs of cell types
.computeKernelPairs <- function(object, cts, sigmaValues, lowerLimit, upperQuantile,
                               normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                               colNormalizeKernel, verbose) {
  
  kernel_mat <- .initializeKernelStructure(sigmaValues, cts)
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  pair_cell_types <- combn(cts, 2)
  
  # Track sigma values to remove
  sigmaValuesToRemove <- vector(mode = "logical", length = length(sigmaValues))
  names(sigmaValuesToRemove) <- sigma_names
  
  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    sigma_choose <- sigmaValues[tt]
    
    if (verbose) {
      cat("current sigma value is", sigma_choose, "\n")
    }
    
    for (pp in seq_len(ncol(pair_cell_types))) {
      i <- pair_cell_types[1, pp]
      j <- pair_cell_types[2, pp]
      
      # Compute kernel from distance using flat structure
      dist_flat_name <- .createDistMatrixName(i, j, slide = NULL)
      if (!dist_flat_name %in% names(object@distances)) {
        warning(paste("Distance matrix not found for", i, "-", j, ". Skipping."))
        next
      }
      
      kernel_current <- kernel_from_distance(
        sigma = sigma_choose,
        dist_mat = object@distances[[dist_flat_name]],
        lower_limit = lowerLimit
      )
      
      # Check if sigma should be removed
      should_remove <- .CheckSigmaValuesToRemove(
        kernel_current = kernel_current, lowerLimit = lowerLimit,
        minAveCellNeighor = minAveCellNeighor, sigma_choose = sigma_choose,
        sigmaValues = sigmaValues, i = i, j = j
      )
      
      if (should_remove) {
        sigmaValuesToRemove[t] <- TRUE
        # Store in flat structure even if marked for removal
        flat_name <- .createKernelMatrixName(sigma_choose, i, j, slide = NULL)
        kernel_mat[[flat_name]] <- kernel_current
        next
      }
      
      # Process kernel matrix (clipping and normalization)
      kernel_current <- .processKernelMatrix(kernel_current, lowerLimit, upperQuantile,
                                           normalizeKernel, rowNormalizeKernel, colNormalizeKernel)
      
      # Store in flat structure
      flat_name <- .createKernelMatrixName(sigma_choose, i, j, slide = NULL)
      kernel_mat[[flat_name]] <- kernel_current
    }
  }
  
  # Clean up removed sigma values
  object <- .cleanupSigmaValues(object, kernel_mat, sigmaValuesToRemove, verbose)
  return(object)
}

# Function for computing kernel matrix within a single cell type
.computeKernelWithin <- function(object, cts, sigmaValues, lowerLimit, upperQuantile,
                                normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                colNormalizeKernel, verbose) {
  
  kernel_mat <- .initializeKernelStructure(sigmaValues, cts)
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  
  # Track sigma values to remove
  sigmaValuesToRemove <- vector(mode = "logical", length = length(sigmaValues))
  names(sigmaValuesToRemove) <- sigma_names
  
  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    sigma_choose <- sigmaValues[tt]
    
    if (verbose) {
      cat("current sigma value is", sigma_choose, "\n")
    }
    
    # Compute kernel from distance using flat structure
    dist_flat_name <- .createDistMatrixName(cts, cts, slide = NULL)
    if (!dist_flat_name %in% names(object@distances)) {
      warning(paste("Distance matrix not found for", cts, "-", cts, ". Skipping."))
      next
    }
    
    kernel_current <- kernel_from_distance(
      sigma = sigma_choose,
      dist_mat = object@distances[[dist_flat_name]],
      lower_limit = lowerLimit
    )
    
    # Check if sigma should be removed
    should_remove <- .CheckSigmaValuesToRemove(
      kernel_current = kernel_current, lowerLimit = lowerLimit,
      minAveCellNeighor = minAveCellNeighor, sigma_choose = sigma_choose,
      sigmaValues = sigmaValues, i = cts, j = cts
    )
    
    if (should_remove) {
      sigmaValuesToRemove[t] <- TRUE
      # Store in flat structure even if marked for removal
      flat_name <- .createKernelMatrixName(sigma_choose, cts, cts, slide = NULL)
      kernel_mat[[flat_name]] <- kernel_current
      next
    }
    
    # Process kernel matrix (clipping and normalization)
    kernel_current <- .processKernelMatrix(kernel_current, lowerLimit, upperQuantile,
                                         normalizeKernel, rowNormalizeKernel, colNormalizeKernel)
    
    # Store in flat structure
    flat_name <- .createKernelMatrixName(sigma_choose, cts, cts, slide = NULL)
    kernel_mat[[flat_name]] <- kernel_current
  }
  
  # Clean up removed sigma values
  object <- .cleanupSigmaValues(object, kernel_mat, sigmaValuesToRemove, verbose)
  return(object)
}

#' Compute Kernel Matrix for CoPro
#'
#' This method calculates the kernel matrices for pairs of cell types based on
#' their distances and a range of sigma values.
#' The formula of calculating kernel matrix is:
#' \deqn{K(x, y) = \exp\left(-\frac{\|x-y\|^2}{2 \sigma^2}\right)}
#' The matrices are adjusted by clipping the upper quantile of
#'  the values to reduce the effect of outliers. The results are stored
#'  within the object.
#'
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @param object A `CoPro` object.
#' @param sigmaValues A vector of sigma values used for kernel calculation.
#' @param lowerLimit The lower limit for the kernel function, default is 1e-7.
#' @param upperQuantile The quantile used for clipping the kernel values,
#' default is 0.85.
#' @param verbose Whether to output the progress and related information
#' @param normalizeKernel Whether to normalize the kernel matrix?
#' Default = FALSE. Note that normalization will not affect any downstream
#' analyses, it is for numerical stability and easier interpretation only.
#' @param minAveCellNeighor What is the minimum average number of cell in the
#'  neighbor? This step is to help set up the expected sparsity of the
#'  kernel matrix. If a kernel sigma value is too small, this result in too
#'  few neighbors for most cells, resulting in an overly-sparse matrix that
#'  makes the parameter estimation hard. Thus, the sigma values that results in
#'  an overly-sparse matrix will be removed for later analysis.
#' @param rowNormalizeKernel Whether the kernel matrix will be row-wise
#' normalized? Note that row or column wise normalization will result in an
#' asymmetric result in skrCCA inference.
#' @param colNormalizeKernel Whether the kernel matrix will be column-wise
#' normalized? Note that row or column wise normalization will result in an
#' asymmetric result in skrCCA inference.
#' @return The `CoPro` object with computed kernel matrices added. The kernel
#' matrices are organized into a three-layer nested list object. The first layer
#' is indexed by the sigma value, and the second and the third layers are cell
#' types
#' @rdname computeKernelMatrix
#' @aliases computeKernelMatrix,CoProSingle-method
#' @aliases computeKernelMatrix,CoProMulti-method
#' @export
#' @note To-do: Shall we include row or column normalization of the kernel?
setGeneric(
  "computeKernelMatrix",
  function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
           normalizeKernel = FALSE, minAveCellNeighor = 2,
           rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
           verbose = TRUE) standardGeneric("computeKernelMatrix"))

#' @rdname computeKernelMatrix
#' @aliases computeKernelMatrix,CoProSingle-method
#' @export
setMethod("computeKernelMatrix", "CoProSingle", 
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   verbose = TRUE) {
            .computeKernelCore(object, sigmaValues, lowerLimit, upperQuantile,
                              normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                              colNormalizeKernel, verbose)
          })

#' @rdname computeKernelMatrix
#' @aliases computeKernelMatrix,CoProMulti-method
#' @export
setMethod("computeKernelMatrix", "CoProMulti", 
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   verbose = TRUE) {
            .computeKernelCoreMulti(object, sigmaValues, lowerLimit, upperQuantile,
                                   normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                   colNormalizeKernel, verbose)
          })

# Core dispatcher for multi-slide objects
.computeKernelCoreMulti <- function(object, sigmaValues, lowerLimit, upperQuantile, 
                                   normalizeKernel, minAveCellNeighor, rowNormalizeKernel, 
                                   colNormalizeKernel, verbose) {
  
  validated_input <- .checkInputKernelMulti(object, sigmaValues, lowerLimit, upperQuantile,
                                           normalizeKernel, minAveCellNeighor, 
                                           rowNormalizeKernel, colNormalizeKernel)
  cts <- validated_input$cts
  sigmaValues <- validated_input$sigmaValues
  
  # Save sigma values
  object@sigmaValues <- sigmaValues
  
  # Provide user feedback
  if (verbose) {
    if (length(cts) == 1) {
      cat("Computing kernel matrix for one cell type across", length(getSlideList(object)), "slides\n")
    } else {
      cat("Computing pairwise kernel matrix for", length(cts), "cell types across", length(getSlideList(object)), "slides\n")
    }
    
    if (normalizeKernel) {
      cat("normalizeKernel is set to TRUE. Kernel matrix will be",
          "normalized so that median row sums of kernel will be 1\n")
    }
  }
  
  # Determine computation approach
  if (length(cts) == 1) {
    return(.computeKernelMultiWithin(object, cts, sigmaValues, lowerLimit, upperQuantile,
                                    normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                    colNormalizeKernel, verbose))
  } else {
    return(.computeKernelMultiPairs(object, cts, sigmaValues, lowerLimit, upperQuantile,
                                   normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                   colNormalizeKernel, verbose))
  }
}

# Helper function to validate kernel computation inputs for multi-slide
.checkInputKernelMulti <- function(object, sigmaValues, lowerLimit, upperQuantile, 
                                  normalizeKernel, minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel) {
  # Check if distance matrix exists
  if (length(object@distances) == 0) {
    stop("Please run computeDistance before computing kernel")
  }

  # Check for conflicting normalization options
  if (rowNormalizeKernel && colNormalizeKernel) {
    stop("Cannot do both row-wise and column-wise normalization.")
  }

  # Get cell types of interest
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) {
    stop("No cell types of interest specified")
  }

  # Check slides
  slides <- getSlideList(object)
  if (length(slides) == 0) {
    stop("No slides found in multi-slide object")
  }

  # Validate sigma values
  if (length(sigmaValues) == 0) {
    stop("sigmaValues must be provided for multi-slide kernel computation")
  }
  
  if (!is.numeric(sigmaValues) || any(sigmaValues <= 0)) {
    stop("sigmaValues must be positive numeric values")
  }

  # Validate other parameters
  if (lowerLimit <= 0 || lowerLimit >= 1) {
    stop("lowerLimit must be between 0 and 1")
  }
  
  if (upperQuantile <= 0 || upperQuantile >= 1) {
    stop("upperQuantile must be between 0 and 1")
  }
  
  if (minAveCellNeighor < 1) {
    stop("minAveCellNeighor must be at least 1")
  }

  return(list(cts = cts, sigmaValues = sigmaValues, slides = slides))
}

# Helper function to initialize kernel matrix structure for multi-slide
.initializeKernelStructureMulti <- function(sigmaValues, slides, cts) {
  # Initialize flat structure with informative names for multi-slide
  kernel_matrices_all <- list()
  
  # Pre-allocate for all sigma-slide-celltype combinations
  for (sigma_val in sigmaValues) {
    for (sID in slides) {
      # Add between-celltype pairs if we have multiple cell types
      if (length(cts) > 1) {
        pair_cell_types <- combn(cts, 2)
        for (pp in seq_len(ncol(pair_cell_types))) {
          ct1 <- pair_cell_types[1, pp]
          ct2 <- pair_cell_types[2, pp]
          flat_name <- .createKernelMatrixName(sigma_val, ct1, ct2, slide = sID)
          kernel_matrices_all[[flat_name]] <- NULL  # Placeholder
        }
      }
      
      # Add within-celltype pairs for all cell types
      for (ct in cts) {
        flat_name <- .createKernelMatrixName(sigma_val, ct, ct, slide = sID)
        kernel_matrices_all[[flat_name]] <- NULL  # Placeholder
      }
    }
  }
  
  return(kernel_matrices_all)
}

# Helper function to check kernel validity for multi-slide
.checkKernelValidityMulti <- function(kernel_current, lowerLimit, minAveCellNeighor, 
                                     sigma_val, ct_i, ct_j, sID) {
  # Check for NA values
  if (anyNA(kernel_current)) {
    warning(paste("NA values in kernel for", ct_i, "-", ct_j, "in slide", sID, "sigma =", sigma_val))
  }
  
  # Check sparsity
  n_cells1 <- nrow(kernel_current)
  n_cells2 <- ncol(kernel_current)
  min_count_nonzero <- minAveCellNeighor * min(n_cells1, n_cells2)
  
  if (sum(kernel_current > lowerLimit) < min_count_nonzero) {
    warning(paste("Kernel matrix for", ct_i, "-", ct_j, "in slide", sID,
                  "with sigma =", sigma_val, "is too sparse."))
    return(FALSE)
  }
  
  return(TRUE)
}

# Function for computing kernel matrix within a single cell type across multiple slides
.computeKernelMultiWithin <- function(object, cts, sigmaValues, lowerLimit, upperQuantile,
                                     normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                     colNormalizeKernel, verbose) {
  
  slides <- getSlideList(object)
  kernel_matrices_all <- .initializeKernelStructureMulti(sigmaValues, slides, cts)
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  
  # Track sigma values to remove
  sigmaValuesToRemove <- vector(mode = "logical", length = length(sigmaValues))
  names(sigmaValuesToRemove) <- sigma_names
  
  for (tt in seq_along(sigmaValues)) {
    sigma_val <- sigmaValues[tt]
    sigma_name <- sigma_names[tt]
    
    if (verbose) {
      cat("current sigma value is", sigma_val, "\n")
    }
    
    sigma_valid_across_slides <- TRUE
    
    for (sID in slides) {
      # Check if distance matrix exists for this slide using flat structure
      dist_flat_name <- .createDistMatrixName(cts, cts, slide = sID)
      if (!dist_flat_name %in% names(object@distances)) {
        if (verbose) message(paste("No distance matrix found for slide", sID))
        sigma_valid_across_slides <- FALSE
        next
      }
      
      dist_mat_ij <- object@distances[[dist_flat_name]]
      
      # Check if distance matrix is valid
      if (length(dist_mat_ij) == 0 || all(is.na(dist_mat_ij)) ||
          nrow(dist_mat_ij) == 0 || ncol(dist_mat_ij) == 0) {
        if (verbose) message(paste("Skipping kernel for", cts, "in slide", sID, "(invalid distance matrix)"))
        next
      }
      
      # Compute kernel from distance
      kernel_current <- kernel_from_distance(
        sigma = sigma_val,
        dist_mat = dist_mat_ij,
        lower_limit = lowerLimit
      )
      
      # Check kernel validity
      if (!.checkKernelValidityMulti(kernel_current, lowerLimit, minAveCellNeighor,
                                    sigma_val, cts, cts, sID)) {
        sigma_valid_across_slides <- FALSE
      }
      
              # Process kernel matrix (clipping and normalization)
        kernel_current <- .processKernelMatrix(kernel_current, lowerLimit, upperQuantile,
                                             normalizeKernel, rowNormalizeKernel, colNormalizeKernel)
        
        # Store in flat structure
        flat_name <- .createKernelMatrixName(sigma_val, cts, cts, slide = sID)
        kernel_matrices_all[[flat_name]] <- kernel_current
    }
    
    # Mark sigma for removal if invalid across slides
    if (!sigma_valid_across_slides) {
      sigmaValuesToRemove[sigma_name] <- TRUE
      warning(paste("Removing sigma value", sigma_val, "as it was invalid for one or more slides."))
    }
  }
  
  # Clean up removed sigma values
  object <- .cleanupSigmaValuesMulti(object, kernel_matrices_all, sigmaValuesToRemove, verbose)
  return(object)
}

# Function for computing kernel matrices between pairs of cell types across multiple slides
.computeKernelMultiPairs <- function(object, cts, sigmaValues, lowerLimit, upperQuantile,
                                    normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                    colNormalizeKernel, verbose) {
  
  slides <- getSlideList(object)
  kernel_matrices_all <- .initializeKernelStructureMulti(sigmaValues, slides, cts)
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  pair_cell_types <- combn(cts, 2)
  
  # Track sigma values to remove
  sigmaValuesToRemove <- vector(mode = "logical", length = length(sigmaValues))
  names(sigmaValuesToRemove) <- sigma_names
  
  for (tt in seq_along(sigmaValues)) {
    sigma_val <- sigmaValues[tt]
    sigma_name <- sigma_names[tt]
    
    if (verbose) {
      cat("current sigma value is", sigma_val, "\n")
    }
    
    sigma_valid_across_slides <- TRUE
    
    for (sID in slides) {
      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]
        
        # Check if specific distance matrix exists using flat structure
        dist_flat_name <- .createDistMatrixName(ct_i, ct_j, slide = sID)
        if (!dist_flat_name %in% names(object@distances)) {
          if (verbose) message(paste("No distance matrix found for", ct_i, "-", ct_j, "in slide", sID))
          sigma_valid_across_slides <- FALSE
          next
        }
        
        dist_mat_ij <- object@distances[[dist_flat_name]]
        
        # Check if distance matrix is valid
        if (length(dist_mat_ij) == 0 || all(is.na(dist_mat_ij)) ||
            nrow(dist_mat_ij) == 0 || ncol(dist_mat_ij) == 0) {
          if (verbose) message(paste("Skipping kernel for", ct_i, "-", ct_j, "in slide", sID, "(invalid distance matrix)"))
          next
        }
        
        # Compute kernel from distance
        kernel_current <- kernel_from_distance(
          sigma = sigma_val,
          dist_mat = dist_mat_ij,
          lower_limit = lowerLimit
        )
        
        # Check kernel validity
        if (!.checkKernelValidityMulti(kernel_current, lowerLimit, minAveCellNeighor,
                                      sigma_val, ct_i, ct_j, sID)) {
          sigma_valid_across_slides <- FALSE
        }
        
        # Process kernel matrix (clipping and normalization)
        kernel_current <- .processKernelMatrix(kernel_current, lowerLimit, upperQuantile,
                                             normalizeKernel, rowNormalizeKernel, colNormalizeKernel)
        
        # Store in flat structure
        flat_name <- .createKernelMatrixName(sigma_val, ct_i, ct_j, slide = sID)
        kernel_matrices_all[[flat_name]] <- kernel_current
      }
    }
    
    # Mark sigma for removal if invalid across slides
    if (!sigma_valid_across_slides) {
      sigmaValuesToRemove[sigma_name] <- TRUE
      warning(paste("Removing sigma value", sigma_val, "as it was invalid for one or more slides."))
    }
  }
  
  # Clean up removed sigma values
  object <- .cleanupSigmaValuesMulti(object, kernel_matrices_all, sigmaValuesToRemove, verbose)
  return(object)
}

# Helper function to clean up removed sigma values for multi-slide (now using flat structure)
.cleanupSigmaValuesMulti <- function(object, kernel_matrices_all, sigmaValuesToRemove, verbose = TRUE) {
  if (any(sigmaValuesToRemove)) {
    if (verbose) {
      cat("removing", sum(sigmaValuesToRemove), "sigma values", "\n")
    }
    
    # Get sigma values that should be removed
    sigmas_to_remove <- as.numeric(gsub("sigma_", "", names(sigmaValuesToRemove)[sigmaValuesToRemove]))
    
    # Remove flat entries that contain these sigma values
    to_remove <- sapply(names(kernel_matrices_all), function(name) {
      parsed <- .parseKernelMatrixName(name)
      parsed$sigma %in% sigmas_to_remove
    })
    
    if (any(to_remove)) {
      kernel_matrices_all <- kernel_matrices_all[!to_remove]
    }
    
    object@sigmaValues <- object@sigmaValues[!sigmaValuesToRemove]
  }
  
  object@kernelMatrices <- kernel_matrices_all
  return(object)
}