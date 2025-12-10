#' Helper function to get K_ij from flat kernel structure
#' @param flat_kernels Flat list of kernel matrices with names like "kernel|sigma0.1|TypeA|TypeB"
#' @param sigma Sigma value (numeric)
#' @param ct_i Name of the first cell type
#' @param ct_j Name of the second cell type
#' @param slide Slide ID (NULL for single slide)
#' @return The K_ij matrix
#' @importFrom irlba irlba
#' @noRd
get_kernel_matrix_flat <- function(flat_kernels, sigma, ct_i, ct_j, slide = NULL) {
  # Create the expected flat name using the same logic as .createKernelMatrixName
  if (is.null(slide)) {
    flat_name <- paste("kernel", paste0("sigma", sigma), ct_i, ct_j, sep = "|")
  } else {
    flat_name <- paste("kernel", paste0("sigma", sigma), slide, ct_i, ct_j, sep = "|")
  }
  
  # Try direct access
  if (flat_name %in% names(flat_kernels)) {
    return(flat_kernels[[flat_name]])
  }
  
  # Try symmetric access (swap ct_i and ct_j)
  if (is.null(slide)) {
    symmetric_name <- paste("kernel", paste0("sigma", sigma), ct_j, ct_i, sep = "|")
  } else {
    symmetric_name <- paste("kernel", paste0("sigma", sigma), slide, ct_j, ct_i, sep = "|")
  }
  
  if (symmetric_name %in% names(flat_kernels)) {
    return(t(flat_kernels[[symmetric_name]]))
  }
  
  # If we get here, kernel not found
  stop(paste("Could not find kernel matrix for pair:", ct_i, ct_j,
             "with sigma =", sigma, if (!is.null(slide)) paste("on slide", slide),
             "in flat kernel structure."))
}



#' Initialize weight vectors using SVD
#' @param X_list Named list of data matrices
#' @param cell_types Vector of cell type names
#' @return Named list of initial weight vectors
#' @importFrom irlba irlba
#' @noRd
initialize_weights_svd <- function(X_list, cell_types) {
  w_list <- setNames(vector("list", length = length(cell_types)), cell_types)
  
  for (ct in cell_types) {
    if(is.null(X_list[[ct]]) || !is.matrix(X_list[[ct]])) {
      stop(paste("Invalid or missing matrix in X_list for cell type:", ct))
    }
    svd_result <- tryCatch(irlba(X_list[[ct]], nv = 1, right_only = TRUE),
                          error = function(e) {
      stop(paste("SVD failed for cell type:", ct, "Error:", e$message))})
    if(ncol(svd_result$v) < 1) stop(paste("SVD resulted in zero singular vectors for cell type:", ct))
    w_list[[ct]] <- svd_result$v[, 1, drop = FALSE] ## orthonormal
  }
  
  return(w_list)
}

#' Check convergence of weight vectors
#' @param w_list_new New weight list
#' @param w_list_old Old weight list
#' @param cell_types Vector of cell type names
#' @return Maximum difference value
#' @noRd
check_convergence <- function(w_list_new, w_list_old, cell_types) {
  current_max_diff <- 0
  for (ct in cell_types) {
    diff_val <- max(abs(w_list_new[[ct]] - w_list_old[[ct]]))
    if(is.nan(diff_val)) diff_val <- 0
    if(diff_val > current_max_diff) {
      current_max_diff <- diff_val
    }
  }
  return(current_max_diff)
}

#' Compute update vector for cell type i in standard case
#' @param ct_i Cell type i
#' @param cell_types All cell types
#' @param X_list Data matrices
#' @param flat_kernels Flat kernel matrices
#' @param sigma Sigma value
#' @param w_list Current weights
#' @param n_features Number of features
#' @param slide Slide ID (NULL for single slide)
#' @return Update vector for ct_i
#' @noRd
compute_update_vector_standard <- function(ct_i, cell_types, X_list, flat_kernels, sigma, w_list, n_features, slide = NULL) {
  w_i_update_vec <- matrix(0, nrow = n_features, ncol = 1)
  
  for (ct_j in cell_types) {
    if (ct_i == ct_j) next
    
    # Get K_ij from flat structure
    K12 <- get_kernel_matrix_flat(flat_kernels, sigma, ct_i, ct_j, slide)
    
    # Get matrices and current weights
    X1 <- X_list[[ct_i]]
    X2 <- X_list[[ct_j]]
    w2 <- w_list[[ct_j]]
    
    # Efficient calculation: t(X1) %*% (K12 %*% (X2 %*% w2))
    v2 <- X2 %*% w2
    kv2 <- K12 %*% v2
    update_contribution <- crossprod(X1, kv2) # This is Y_ij %*% w_j
    
    # update weight vector
    w_i_update_vec <- w_i_update_vec + update_contribution
  }
  
  return(w_i_update_vec)
}

#' Compute update vector for within-cell-type case
#' @param X Data matrix
#' @param K Kernel matrix
#' @param w Current weight vector
#' @return Update vector
#' @noRd
compute_update_vector_within <- function(X, K, w) {
  # Y = t(X) %*% K %*% X for within-cell-type
  # w_new = Y %*% w = t(X) %*% K %*% X %*% w
  Xw <- X %*% w
  KXw <- K %*% Xw
  w_update <- crossprod(X, KXw)
  return(w_update)
}

#' SkrCCA optimization function for multiple groups (Single Slide) - First Component
#' Uses flat kernel structure for consistent data access
#'
#' @param X_list Named list of data matrices (cell by PC matrix)
#' @param flat_kernels Flat list of kernel matrices with names like "kernel|sigma0.1|TypeA|TypeB"
#' @param sigma Sigma value (numeric)
#' @param max_iter Maximum number of iterations
#' @param tol tolerance of accuracy
#'
#' @return Named list `w_list` containing the first weight vector component.
#' @export
optimize_bilinear <- function(X_list, flat_kernels, sigma, max_iter = 1000,
                              tol = 1e-5) {

  cell_types <- names(X_list)
  if (is.null(cell_types)) stop("Input X_list must be a named list.")
  n_mat <- length(cell_types)
  
  # Auto-detect if this is within-cell-type optimization
  is_within <- (n_mat == 1)
  
  n_features <- ncol(X_list[[cell_types[1]]])

  # Initialize w_list using SVD
  w_list <- initialize_weights_svd(X_list, cell_types)

  # Iterative refinement
  iter <- 0
  while (iter <= max_iter) {
    w_list_old <- w_list

    if (is_within) {
      # Within-cell-type optimization
      ct <- cell_types[1]
      X <- X_list[[ct]]
      K <- get_kernel_matrix_flat(flat_kernels, sigma, ct, ct, slide = NULL)
      
      w_update <- compute_update_vector_within(X, K, w_list[[ct]])
      w_list[[ct]] <- normalize_vec(w_update)
      
    } else {
      # Standard multi-cell-type optimization
      for (ct_i in cell_types) {
        w_i_update_vec <- compute_update_vector_standard(ct_i, cell_types, X_list, flat_kernels, sigma, w_list, n_features, slide = NULL)
        w_list[[ct_i]] <- normalize_vec(w_i_update_vec)
      }
    }

    # Check convergence
    current_max_diff <- check_convergence(w_list, w_list_old, cell_types)

    if (current_max_diff <= tol) {
      print(paste("Convergence reached at", iter, "iterations (Max diff =", sprintf("%.3e", current_max_diff), ")"))
      break
    }
    iter <- iter + 1
  } # end while

  if (iter > max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  # Ensure final format is list of single-column matrices
  for (ct in cell_types) {
    if (!is.matrix(w_list[[ct]]) || ncol(w_list[[ct]]) != 1) {
      w_list[[ct]] <- matrix(w_list[[ct]], ncol = 1)
    }
  }
  return(w_list)
}

#' Compute Y_resi for subsequent components
#' @param X_list Data matrices
#' @param flat_kernels Flat kernel matrices
#' @param sigma Sigma value
#' @param cell_types Cell type names
#' @param slide Slide ID (NULL for single slide)
#' @return Y_resi structure
#' @noRd
compute_Y_resi <- function(X_list, flat_kernels, sigma, cell_types, slide = NULL) {
  n_mat <- length(cell_types)
  is_within <- (n_mat == 1)
  Y_resi <- setNames(vector(mode = "list", length = n_mat), cell_types)
  
  if (is_within) {
    # For within-cell-type, we only have one Y matrix
    ct <- cell_types[1]
    Y_resi[[ct]] <- setNames(list(NULL), ct)
    X <- X_list[[ct]]
    K <- get_kernel_matrix_flat(flat_kernels, sigma, ct, ct, slide)
    Y_resi[[ct]][[ct]] <- crossprod(X, K %*% X)
  } else {
    # Standard case: initialize structure for all pairs
    for (ct_i in cell_types) {
      Y_resi[[ct_i]] <- setNames(vector(mode = "list", length = n_mat), cell_types)
    }
    
    # Compute Y matrices for all pairs
    pair_cell_types <- combn(cell_types, 2)
    for (pp in seq_len(ncol(pair_cell_types))) {
      i <- pair_cell_types[1, pp]
      j <- pair_cell_types[2, pp]
      
      K12 <- get_kernel_matrix_flat(flat_kernels, sigma, i, j, slide)
      Y_ij <- crossprod(X_list[[i]], K12 %*% X_list[[j]])
      
      Y_resi[[i]][[j]] <- Y_ij
      Y_resi[[j]][[i]] <- t(Y_ij)
    }
  }
  
  return(Y_resi)
}

#' Apply deflation to Y_resi
#' @param Y_resi Current Y_resi structure
#' @param w_list Weight list
#' @param qq Component index used for deflation
#' @param cell_types Cell type names
#' @return Updated Y_resi
#' @noRd
apply_deflation <- function(Y_resi, w_list, qq, cell_types) {
  n_mat <- length(cell_types)
  is_within <- (n_mat == 1)
  
  if (is_within) {
    ct <- cell_types[1]
    Y1 <- Y_resi[[ct]][[ct]]
    w1 <- w_list[[ct]][, qq, drop = FALSE]
    
    deflation_scalar <- (t(w1) %*% Y1 %*% w1)[1, 1]
    deflation_term <- deflation_scalar * (w1 %*% t(w1))
    Y_resi[[ct]][[ct]] <- Y1 - deflation_term
  } else {
    pair_cell_types <- combn(cell_types, 2)
    for (pp in seq_len(ncol(pair_cell_types))) {
      i <- pair_cell_types[1, pp]
      j <- pair_cell_types[2, pp]
      
      w1 <- w_list[[i]][, qq, drop = FALSE]
      w2 <- w_list[[j]][, qq, drop = FALSE]
      Y1 <- Y_resi[[i]][[j]]
      
      deflation_scalar <- (t(w1) %*% Y1 %*% w2)[1, 1]
      deflation_term <- deflation_scalar * (w1 %*% t(w2))
      Y_resi[[i]][[j]] <- Y1 - deflation_term
      Y_resi[[j]][[i]] <- t(Y_resi[[i]][[j]])
    }
  }
  
  return(Y_resi)
}

#' Initialize weights for next component using SVD
#' @param Y_resi Deflated Y matrices
#' @param cell_types Cell type names
#' @return Initial weight list for next component
#' @noRd
initialize_next_component <- function(Y_resi, cell_types) {
  n_mat <- length(cell_types)
  is_within <- (n_mat == 1)
  w_list_new <- setNames(vector("list", length = n_mat), cell_types)
  
  if (is_within) {
    ct <- cell_types[1]
    # Ensure matrix is perfectly symmetric for numerical stability
    Y_sym <- Y_resi[[ct]][[ct]]
    Y_sym <- (Y_sym + t(Y_sym)) / 2
    # For symmetric matrix, use eigen decomposition with error handling
    eigen_result <- tryCatch(
      eigen(Y_sym, symmetric = TRUE),
      error = function(e) {
        stop(paste("Eigen decomposition failed for within-cell-type:", ct, "Error:", e$message))
      }
    )
    
    # Check if eigen decomposition resulted in zero eigenvectors
    if(ncol(eigen_result$vectors) < 1) {
      stop(paste("Eigen decomposition resulted in zero eigenvectors for cell type:", ct))
    }
    w_list_new[[ct]] <- eigen_result$vectors[, 1, drop = FALSE]
  } else {
    # Use SVD on cross-cell-type Y matrices with error handling
    w_list_new[[cell_types[1]]] <- tryCatch(
      irlba(t(Y_resi[[cell_types[1]]][[cell_types[2]]]), nv = 1, right_only = TRUE)$v[, 1, drop = FALSE],
      error = function(e) {
        stop(paste("SVD failed for cell type pair initialization. Error:", e$message))
      }
    )
    for (i in cell_types[2:n_mat]) {
      w_list_new[[i]] <- tryCatch(
        irlba(Y_resi[[cell_types[1]]][[i]], nv = 1, right_only = TRUE)$v[, 1, drop = FALSE],
        error = function(e) {
          stop(paste("SVD failed for cell type:", i, "Error:", e$message))
        }
      )
    }
  }
  
  return(w_list_new)
}

#' Helper function for optimizing multiple components using precomputed Y_resi
#' @note Assumes w_list_new and Y_resi are named lists using consistent cell type names.
#' @note the structure of Y_resi is Y_resi\[[ct1\]\]\[[ct2\]\]
#' @param w_list_new Initial named list of weight vectors
#' @param Y_resi Named list of residual matrices 
#' @param n_features Number of features
#' @param max_iter Maximum iterations
#' @param tol Tolerance
#'
#' @return Updated named list of weight vectors
#' @noRd
bilinear_w_from_Y_resi <- function(w_list_new, Y_resi,
                                   n_features, max_iter, tol) {

  cell_types <- names(w_list_new)
  if (length(cell_types) == 0) stop("Input w_list_new must be a named list.")
  n_mat <- length(cell_types)
  is_within <- (n_mat == 1)

  if (length(names(Y_resi)) == 0 || !all(cell_types %in% names(Y_resi))) {
    stop("Y_resi must be a named list containing entries for all cell types in w_list_new.")
  }

  iter <- 0
  while (iter < max_iter) {
    w_list_old <- w_list_new

    if (is_within) {
      # Within-cell-type case
      ct <- cell_types[1]
      Y <- Y_resi[[ct]][[ct]]
      w <- w_list_new[[ct]]
      w_update <- Y %*% w
      w_list_new[[ct]] <- normalize_vec(w_update)
    } else {
      # Standard multi-cell-type case - use direct accumulation
      for (ct_i in cell_types) {
        # We can directly accumulate since we're just summing
        w_i_update_vec <- matrix(0, nrow = n_features, ncol = 1)
        for (ct_j in cell_types) {
          if (ct_i == ct_j) next
          w2 <- w_list_new[[ct_j]]
          Y <- Y_resi[[ct_i]][[ct_j]]
          if(is.null(Y)) stop(paste("Missing Y_resi matrix for pair:", ct_i, ct_j))
          w_i_update_vec <- w_i_update_vec + Y %*% w2
        }
        w_list_new[[ct_i]] <- normalize_vec(w_i_update_vec)
      }
    }

    # Check convergence
    current_max_diff <- check_convergence(w_list_new, w_list_old, cell_types)

    if (current_max_diff <= tol) {
      print(paste("Convergence reached at", iter, "iterations (Max diff =", sprintf("%.3e", current_max_diff), ")"))
      break
    }
    iter <- iter + 1
  } # end while

  if (iter == max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }
  return(w_list_new)
}

#' Run multi version of skrCCA to detect subsequent components (Single Slide)
#' Uses flat kernel structure for consistent data access
#'
#' @param X_list Named list of data matrices (subsetted)
#' @param flat_kernels Flat list of kernel matrices
#' @param sigma Sigma value (numeric)
#' @param w_list A named list of weights (subsetted, matrices with previous components as columns)
#' @param cellTypesOfInterest A vector specifying cell type names present in the input lists
#' @param nCC Total number of canonical vectors desired (must be >= 2)
#' @param max_iter Maximum number of iterations for helper function
#' @param tol Tolerance of accuracy for helper function
#'
#' @return A named list of weights (matrices with components 1 to nCC as columns)
#' @export
optimize_bilinear_n <- function(X_list, flat_kernels, sigma, w_list,
                                      cellTypesOfInterest,
                                      nCC = 2,
                                      max_iter = 1000,
                                      tol = 1e-5) {

  # Validate inputs based on assumption they are already subsetted
  cts <- cellTypesOfInterest
  n_mat <- length(cts)
  is_within <- (n_mat == 1)

  
  if (length(X_list) != n_mat || length(w_list) != n_mat ||
      !all(cts %in% names(X_list)) || !all(cts %in% names(w_list))) {
    stop("Input lists length or names do not match cellTypesOfInterest.")
  }
  n_features <- ncol(X_list[[cts[1]]])

  # Check input w_list dimensions
  k_start <- ncol(w_list[[cts[1]]]) ## we expect this value to be one
  if (k_start < 1) stop("Input w_list must contain at least the first component.")
  for (ct in cts) {
    if (!is.matrix(w_list[[ct]]) || ncol(w_list[[ct]]) != k_start) {
      stop(paste("w_list for", ct, "is not a matrix or has inconsistent component count."))
    }
  }
  if (nCC <= k_start) {
    stop(paste("nCC (", nCC, ") must be greater than the number of components already in w_list (", k_start, ")"))
  }

  # Initialize Y_resi structure using original data
  Y_resi <- compute_Y_resi(X_list, flat_kernels, sigma, cts, slide = NULL)

  # Loop to compute components k_start + 1 up to nCC
  for (qq in k_start:(nCC - 1)) {
    # Step 1: Apply deflation using component qq
    Y_resi <- apply_deflation(Y_resi, w_list, qq, cts)

    # Step 2: Initialize w_list_new for component qq+1
    w_list_new <- initialize_next_component(Y_resi, cts)

    # Step 3: Iterative refinement using the helper function
    w_list_qq_plus_1 <- bilinear_w_from_Y_resi(
      w_list_new = w_list_new,
      Y_resi = Y_resi,
      n_features = n_features,
      max_iter = max_iter, 
      tol = tol)

    # Step 4: Add the new component (qq+1) to w_list
    for (ct in cts) {
      w_list[[ct]] <- cbind(w_list[[ct]], w_list_qq_plus_1[[ct]])
    }
  } # end component loop qq

  return(w_list)
}

#' @importFrom parallel mclapply
#' @importFrom irlba irlba
#' @importFrom stats setNames
NULL

#' Multi-slide optimization functions for CoPro
#' 
#' These functions extend the CoPro formulation to handle multiple slides
#' where weight vectors are shared across slides.
#' 
#' The objective is:
#' \deqn{Maximize_{\{w_i,w_j\}} \sum_q w_i^T X_{i,q}^T K_{ij,q} X_{j,q} w_j}
#' subject to \eqn{||w_i|| \leq 1} for all \eqn{i}, where \eqn{q} is the slide index.
#' 
#' @name multi_slide_optimization
#' @keywords internal
NULL

# ============================================================================
# Helper Functions Specific to Multi-Slide (shared functions are in 04_optimization_function_refactored.R)
# ============================================================================

#' Validate multi-slide input data structure
#' @param X_list_all List of lists of data matrices
#' @param K_list_all List of lists of kernel matrices (deprecated, can be NULL)
#' @param expected_cell_types Expected cell types (NULL to auto-detect)
#' @param check_single_type Whether to check for single cell type
#' @return List with validated inputs and metadata
#' @noRd
validate_multi_slide_inputs <- function(X_list_all, K_list_all = NULL, 
                                       expected_cell_types = NULL,
                                       check_single_type = FALSE) {
  n_slides <- length(X_list_all)
  
  if (n_slides < 1) {
    stop("Need at least one slide.")
  }
  
  # K_list_all validation is now optional (for backward compatibility)
  if (!is.null(K_list_all) && length(K_list_all) != n_slides) {
    stop("X_list_all and K_list_all must have the same length (number of slides).")
  }
  
  # Get cell types from all slides
  cell_types_all <- lapply(X_list_all, names)
  
  # Check consistency
  if (n_slides > 1 && !all(sapply(cell_types_all[-1], function(x) identical(sort(x), sort(cell_types_all[[1]]))))) {
    stop("All slides must have the same cell types.")
  }
  
  cell_types <- cell_types_all[[1]]
  n_cell_types <- length(cell_types)
  
  # Check if single cell type when required
  if (check_single_type && n_cell_types != 1) {
    stop("This function requires exactly one cell type.")
  }
  
  # Check against expected cell types if provided
  if (!is.null(expected_cell_types) && !identical(sort(cell_types), sort(expected_cell_types))) {
    stop("Cell types in data do not match expected cell types.")
  }
  
  # Get number of features from first slide, first cell type
  n_features <- ncol(X_list_all[[1]][[cell_types[1]]])
  
  # Validate dimensions across all slides and cell types
  for (q in seq_len(n_slides)) {
    for (ct in cell_types) {
      if (!ct %in% names(X_list_all[[q]])) {
        stop(paste("Cell type", ct, "missing in slide", q))
      }
      if (!is.matrix(X_list_all[[q]][[ct]]) || ncol(X_list_all[[q]][[ct]]) != n_features) {
        stop(paste("Invalid or inconsistent feature count in slide", q, "cell type", ct))
      }
    }
  }
  
  return(list(
    n_slides = n_slides,
    cell_types = cell_types,
    n_cell_types = n_cell_types,
    n_features = n_features
  ))
}

#' Initialize weights for multi-slide using SVD on aggregated data
#' @param X_list_all List of lists of data matrices
#' @param cell_types Vector of cell type names
#' @param use_aggregation Whether to aggregate across slides for initialization
#' @return Named list of initial weight vectors
#' @noRd
initialize_weights_multi_slide <- function(X_list_all, cell_types, use_aggregation = TRUE) {
  w_list <- setNames(vector("list", length(cell_types)), cell_types)
  
  for (ct in cell_types) {
    if (use_aggregation && length(X_list_all) > 1) {
      # Stack data from all slides for better initialization
      X_stacked <- do.call(rbind, lapply(X_list_all, function(slide) slide[[ct]]))
    } else {
      # Use first slide only
      X_stacked <- X_list_all[[1]][[ct]]
    }
    
    svd_result <- tryCatch(
      irlba(X_stacked, nv = 1, right_only = TRUE),
      error = function(e) {
        stop(paste("SVD failed for cell type:", ct, "Error:", e$message))
      }
    )
    
    if (ncol(svd_result$v) < 1) {
      stop(paste("SVD resulted in zero singular vectors for cell type:", ct))
    }
    
    w_list[[ct]] <- svd_result$v[, 1, drop = FALSE]
  }
  
  return(w_list)
}

#' Compute Y matrix for a single slide
#' @param X_list Data matrices for one slide
#' @param flat_kernels Flat kernel matrices
#' @param sigma Sigma value
#' @param slide Slide ID
#' @param ct_i First cell type
#' @param ct_j Second cell type
#' @return Y_ij matrix
#' @noRd
compute_Y_slide <- function(X_list, flat_kernels, sigma, slide, ct_i, ct_j) {
  if (ct_i == ct_j) {
    # Within-cell-type case
    X <- X_list[[ct_i]]
    K <- get_kernel_matrix_flat(flat_kernels, sigma, ct_i, ct_i, slide)
    return(crossprod(X, K %*% X))
  } else {
    # Between-cell-type case
    X1 <- X_list[[ct_i]]
    X2 <- X_list[[ct_j]]
    K12 <- get_kernel_matrix_flat(flat_kernels, sigma, ct_i, ct_j, slide)
    return(crossprod(X1, K12 %*% X2))
  }
}

#' Aggregate Y matrices across slides
#' @param X_list_all List of lists of data matrices
#' @param flat_kernels Flat kernel matrices
#' @param sigma Sigma value
#' @param slides Slide IDs
#' @param cell_types Cell type names
#' @param n_cores Number of cores for parallel computation
#' @return Aggregated Y matrices structure
#' @noRd
compute_Y_multi_slide <- function(X_list_all, flat_kernels, sigma, slides, cell_types, n_cores = 1) {
  n_slides <- length(X_list_all)
  n_mat <- length(cell_types)
  is_within <- (n_mat == 1)
  
  if (is_within) {
    # Single cell type case
    ct <- cell_types[1]
    Y_list <- mclapply(seq_len(n_slides), function(q) {
      compute_Y_slide(X_list_all[[q]], flat_kernels, sigma, slides[q], ct, ct)
    }, mc.cores = n_cores)
    
    Y_sum <- Reduce("+", Y_list)
    Y_aggregate <- setNames(list(setNames(list(Y_sum), ct)), ct)
    
  } else {
    # Multiple cell types case
    Y_aggregate <- setNames(vector("list", n_mat), cell_types)
    for (ct in cell_types) {
      Y_aggregate[[ct]] <- setNames(vector("list", n_mat), cell_types)
    }
    
    # Process all pairs
    pair_indices <- combn(n_mat, 2)
    for (pp in seq_len(ncol(pair_indices))) {
      i_idx <- pair_indices[1, pp]
      j_idx <- pair_indices[2, pp]
      ct_i <- cell_types[i_idx]
      ct_j <- cell_types[j_idx]
      
      # Compute Y_ij for all slides in parallel
      Y_ij_list <- mclapply(seq_len(n_slides), function(q) {
        compute_Y_slide(X_list_all[[q]], flat_kernels, sigma, slides[q], ct_i, ct_j)
      }, mc.cores = n_cores)
      
      Y_ij_sum <- Reduce("+", Y_ij_list)
      Y_aggregate[[ct_i]][[ct_j]] <- Y_ij_sum
      Y_aggregate[[ct_j]][[ct_i]] <- t(Y_ij_sum)
    }
  }
  
  return(Y_aggregate)
}

#' Compute update vector for multi-slide optimization
#' @param ct_i Cell type to update
#' @param cell_types All cell types
#' @param X_list_all List of lists of data matrices
#' @param flat_kernels Flat kernel matrices
#' @param sigma Sigma value
#' @param slides Slide IDs
#' @param w_list Current weights
#' @param n_features Number of features
#' @param n_cores Number of cores
#' @return Update vector
#' @noRd
compute_update_vector_multi_slide <- function(ct_i, cell_types, X_list_all, flat_kernels, 
                                             sigma, slides, w_list, n_features, n_cores = 1) {
  n_slides <- length(X_list_all)
  w_i_update_vec <- matrix(0, nrow = n_features, ncol = 1)
  
  for (ct_j in cell_types) {
    if (ct_i == ct_j) next
    
    w_j <- w_list[[ct_j]]
    
    # Compute contribution from all slides in parallel
    contributions <- mclapply(seq_len(n_slides), function(q) {
      X_i <- X_list_all[[q]][[ct_i]]
      X_j <- X_list_all[[q]][[ct_j]]
      K_ij <- get_kernel_matrix_flat(flat_kernels, sigma, ct_i, ct_j, slides[q])
      
      # Efficient computation: t(X_i) %*% K_ij %*% (X_j %*% w_j)
      v_j <- X_j %*% w_j
      kv_j <- K_ij %*% v_j
      crossprod(X_i, kv_j)
    }, mc.cores = n_cores)
    
    # Sum contributions from all slides
    contribution_sum <- Reduce("+", contributions)
    w_i_update_vec <- w_i_update_vec + contribution_sum
  }
  
  return(w_i_update_vec)
}

# ============================================================================
# Main Optimization Functions
# ============================================================================

#' Multi-slide SkrCCA optimization - First Component
#' 
#' Handles both standard (multiple cell types) and within (single cell type) cases
#' Uses flat kernel structure for consistent data access
#' 
#' @param X_list_all List of lists of data matrices
#' @param flat_kernels Flat list of kernel matrices
#' @param sigma Sigma value (numeric)
#' @param slides Slide IDs
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param n_cores Number of cores for parallel computation
#' @param direct_solve For single cell type, use direct eigenvalue solution
#' @return Named list of weight vectors (first component)
#' @export
optimize_bilinear_multi_slides <- function(X_list_all, flat_kernels, sigma, slides,
                                          max_iter = 1000, tol = 1e-5,
                                          n_cores = 1, direct_solve = TRUE) {
  
  # Validate inputs
  validated <- validate_multi_slide_inputs(X_list_all, NULL)  # Don't validate K_list_all anymore
  n_slides <- validated$n_slides
  cell_types <- validated$cell_types
  n_cell_types <- validated$n_cell_types
  n_features <- validated$n_features
  
  # Check if this is within-cell-type case
  is_within <- (n_cell_types == 1)
  
  # For single slide, delegate to single-slide function
  if (n_slides == 1) {
    message("Single slide detected, using single-slide optimization")
    return(optimize_bilinear(X_list_all[[1]], flat_kernels, sigma,
                            max_iter = max_iter, tol = tol))
  }
  
  # Handle within-cell-type case with direct solution
  if (is_within && direct_solve) {
    ct <- cell_types[1]
    
    # Aggregate Y matrices across slides
    Y_aggregate <- compute_Y_multi_slide(X_list_all, flat_kernels, sigma, slides, cell_types, n_cores)
    Y_sum <- Y_aggregate[[ct]][[ct]]
    
    # Direct eigenvalue solution
    eigen_result <- tryCatch(
      eigen(Y_sum, symmetric = TRUE),
      error = function(e) {
        warning(paste("Eigen decomposition failed:", e$message, 
                     "\nFalling back to iterative method"))
        return(NULL)
      }
    )
    
    if (!is.null(eigen_result)) {
      w_list <- setNames(list(eigen_result$vectors[, 1, drop = FALSE]), ct)
      message(paste("Direct solution found, largest eigenvalue:", 
                   round(eigen_result$values[1], 6)))
      return(w_list)
    }
  }
  
  # Initialize weights
  w_list <- initialize_weights_multi_slide(X_list_all, cell_types, use_aggregation = TRUE)
  
  # Iterative optimization
  iter <- 0
  while (iter <= max_iter) {
    w_list_old <- w_list
    
    if (is_within) {
      # Within-cell-type optimization
      ct <- cell_types[1]
      
      # Compute aggregated update
      update_contributions <- mclapply(seq_len(n_slides), function(q) {
        X <- X_list_all[[q]][[ct]]
        K <- get_kernel_matrix_flat(flat_kernels, sigma, ct, ct, slides[q])
        w <- w_list[[ct]]
        
        # Compute: t(X) %*% K %*% X %*% w
        Xw <- X %*% w
        KXw <- K %*% Xw
        crossprod(X, KXw)
      }, mc.cores = n_cores)
      
      w_update <- Reduce("+", update_contributions)
      w_list[[ct]] <- normalize_vec(w_update)
      
    } else {
      # Standard multi-cell-type optimization
      for (ct_i in cell_types) {
        w_i_update_vec <- compute_update_vector_multi_slide(
          ct_i, cell_types, X_list_all, flat_kernels, sigma, slides, w_list, n_features, n_cores
        )
        w_list[[ct_i]] <- normalize_vec(w_i_update_vec)
      }
    }
    
    # Check convergence
    current_max_diff <- check_convergence(w_list, w_list_old, cell_types)
    
    if (current_max_diff <= tol) {
      message(paste("Convergence reached at", iter, "iterations (Max diff =", 
                   sprintf("%.3e", current_max_diff), ")"))
      break
    }
    iter <- iter + 1
  }
  
  if (iter > max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }
  
  # Ensure proper matrix format
  for (ct in cell_types) {
    if (!is.matrix(w_list[[ct]]) || ncol(w_list[[ct]]) != 1) {
      w_list[[ct]] <- matrix(w_list[[ct]], ncol = 1)
    }
  }
  
  return(w_list)
}

#' Multi-slide SkrCCA optimization - Multiple Components
#' 
#' Computes components 2 to nCC using deflation
#' Uses flat kernel structure for consistent data access
#' 
#' @param X_list_all List of lists of data matrices
#' @param flat_kernels Flat list of kernel matrices
#' @param sigma Sigma value (numeric)
#' @param slides Slide IDs
#' @param w_list Initial weight list with first component(s)
#' @param cellTypesOfInterest Cell types to process
#' @param nCC Total number of components desired
#' @param max_iter Maximum iterations for refinement
#' @param tol Convergence tolerance
#' @param n_cores Number of cores for parallel computation
#' @return Updated weight list with all components
#' @export
optimize_bilinear_n_multi_slides <- function(X_list_all, flat_kernels, sigma, slides, w_list,
                                            cellTypesOfInterest,
                                            nCC = 2, max_iter = 1000, 
                                            tol = 1e-5, n_cores = 1) {
  
  # Validate inputs
  validated <- validate_multi_slide_inputs(X_list_all, NULL, 
                                          expected_cell_types = cellTypesOfInterest)
  n_slides <- validated$n_slides
  cell_types <- cellTypesOfInterest
  n_cell_types <- length(cell_types)
  n_features <- validated$n_features
  is_within <- (n_cell_types == 1)
  
  # Validate w_list
  if (length(w_list) != n_cell_types || 
      !all(cell_types %in% names(w_list))) {
    stop("w_list structure does not match cellTypesOfInterest")
  }
  
  # Check existing components
  k_start <- ncol(w_list[[cell_types[1]]])
  for (ct in cell_types) {
    if (!is.matrix(w_list[[ct]]) || ncol(w_list[[ct]]) != k_start) {
      stop(paste("Inconsistent component count in w_list for cell type:", ct))
    }
  }
  
  if (nCC <= k_start) {
    stop(paste("nCC (", nCC, ") must be greater than existing components (", k_start, ")"))
  }
  
  # Initialize Y_resi for all slides
  Y_resi_all <- vector("list", n_slides)
  for (q in seq_len(n_slides)) {
    Y_resi_all[[q]] <- compute_Y_resi(X_list_all[[q]], flat_kernels, sigma, cell_types, slides[q])
  }
  
  # Compute additional components
  for (qq in k_start:(nCC - 1)) {
    
    # Step 1: Apply deflation across all slides
    Y_resi_all <- mclapply(seq_len(n_slides), function(q) {
      apply_deflation(Y_resi_all[[q]], w_list, qq, cell_types)
    }, mc.cores = n_cores)
    
    # Step 2: Aggregate deflated Y matrices
    if (is_within) {
      ct <- cell_types[1]
      Y_sum_list <- lapply(Y_resi_all, function(Y_resi) Y_resi[[ct]][[ct]])
      Y_aggregate <- setNames(list(setNames(list(Reduce("+", Y_sum_list)), ct)), ct)
    } else {
      # Initialize aggregate structure
      Y_aggregate <- setNames(vector("list", n_cell_types), cell_types)
      for (ct in cell_types) {
        Y_aggregate[[ct]] <- setNames(vector("list", n_cell_types), cell_types)
      }
      
      # Sum across slides for each pair
      for (ct_i in cell_types) {
        for (ct_j in cell_types) {
          if (ct_i == ct_j) next
          Y_ij_list <- lapply(Y_resi_all, function(Y_resi) Y_resi[[ct_i]][[ct_j]])
          Y_aggregate[[ct_i]][[ct_j]] <- Reduce("+", Y_ij_list)
        }
      }
    }
    
    # Step 3: Initialize next component
    w_list_new <- initialize_next_component(Y_aggregate, cell_types)
    
    # Step 4: Iterative refinement
    w_list_qq_plus_1 <- bilinear_w_from_Y_resi(
      w_list_new = w_list_new,
      Y_resi = Y_aggregate,
      n_features = n_features,
      max_iter = max_iter,
      tol = tol
    )
    
    # Step 5: Append new component
    for (ct in cell_types) {
      w_list[[ct]] <- cbind(w_list[[ct]], w_list_qq_plus_1[[ct]])
    }
  }
  
  return(w_list)
}


