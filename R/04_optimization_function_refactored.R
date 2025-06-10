#' Helper function to safely get K_ij or t(K_ji) from K_list
#' @note Assumes K_list is named list K_list[[ct1]][[ct2]]
#' @note Assumes if K_list[[ct_i]][[ct_j]] is NULL, K_list[[ct_j]][[ct_i]] exists.
#' @param K_list The named kernel list structure
#' @param ct_i Name of the first cell type
#' @param ct_j Name of the second cell type
#' @return The K_ij matrix
#' @importFrom irlba irlba
#' @noRd
get_kernel_matrix <- function(K_list, ct_i, ct_j) {
  K12 <- NULL
  # Try accessing K_ij directly
  if (length(K_list[[ct_i]][[ct_j]]) != 0) {
    K12 <- K_list[[ct_i]][[ct_j]]
  } else {
    # If K_ij is missing, try getting K_ji and transposing it
    K21 <- K_list[[ct_j]][[ct_i]]
    K12 <- t(K21)
  }

  # Check if we successfully obtained a matrix
  if (length(K12) == 0) {
    stop(paste("Could not find kernel matrix for pair:", ct_i, ct_j,
               "or its transpose in K_list."))
  }
  return(K12)
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
#' @param K_list Kernel matrices
#' @param w_list Current weights
#' @param n_features Number of features
#' @return Update vector for ct_i
#' @noRd
compute_update_vector_standard <- function(ct_i, cell_types, X_list, K_list, w_list, n_features) {
  w_i_update_vec <- matrix(0, nrow = n_features, ncol = 1)
  
  for (ct_j in cell_types) {
    if (ct_i == ct_j) next
    
    # Get K_ij or t(K_ji) robustly
    K12 <- get_kernel_matrix(K_list, ct_i, ct_j)
    
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
#' Assumes X_list and K_list are named lists using consistent cell type names.
#' Handles upper/lower triangle K_list. Uses w_i_left structure.
#'
#' @param X_list Named list of data matrices (cell by PC matrix)
#' @param K_list Named list of kernel matrices (potentially upper/lower triangle only)
#' @param max_iter Maximum number of iterations
#' @param tol tolerance of accuracy
#'
#' @return Named list `w_list` containing the first weight vector component.
#' @export
optimize_bilinear <- function(X_list, K_list, max_iter = 1000,
                              tol = 1e-5) {

  cell_types <- names(X_list)
  if (is.null(cell_types)) stop("Input X_list must be a named list.")
  n_mat <- length(cell_types)
  
  # Auto-detect if this is within-cell-type optimization
  is_within <- (n_mat == 1)
  
  n_features <- ncol(X_list[[cell_types[1]]])

  # Basic check on K_list names
  if (is.null(names(K_list)) || !all(cell_types %in% names(K_list))) {
    stop("K_list must be a named list containing entries for all cell types in X_list.")
  }

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
      K <- K_list[[ct]][[ct]]
      
      if (is.null(K)) {
        stop(paste("Kernel matrix not found for within-cell-type:", ct))
      }
      
      w_update <- compute_update_vector_within(X, K, w_list[[ct]])
      w_list[[ct]] <- normalize_vec(w_update)
      
    } else {
      # Standard multi-cell-type optimization
      for (ct_i in cell_types) {
        w_i_update_vec <- compute_update_vector_standard(ct_i, cell_types, X_list, K_list, w_list, n_features)
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
#' @param K_list Kernel matrices
#' @param cell_types Cell type names
#' @return Y_resi structure
#' @noRd
compute_Y_resi <- function(X_list, K_list, cell_types) {
  n_mat <- length(cell_types)
  is_within <- (n_mat == 1)
  Y_resi <- setNames(vector(mode = "list", length = n_mat), cell_types)
  
  if (is_within) {
    # For within-cell-type, we only have one Y matrix
    ct <- cell_types[1]
    Y_resi[[ct]] <- setNames(list(NULL), ct)
    X <- X_list[[ct]]
    K <- K_list[[ct]][[ct]]
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
      
      K12 <- get_kernel_matrix(K_list, i, j)
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
#' @note the structure of Y_resi is Y_resi[[ct1]][[ct2]]
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
#' Assumes X_list, K_list, w_list are named lists containing only cellTypesOfInterest.
#' Handles upper/lower triangle K_list. No internal subsetting.
#'
#' @param X_list Named list of data matrices (subsetted)
#' @param K_list Named list of kernel matrices (subsetted, potentially upper/lower triangle)
#' @param w_list A named list of weights (subsetted, matrices with previous components as columns)
#' @param cellTypesOfInterest A vector specifying cell type names present in the input lists
#' @param nCC Total number of canonical vectors desired (must be >= 2)
#' @param max_iter Maximum number of iterations for helper function
#' @param tol Tolerance of accuracy for helper function
#'
#' @return A named list of weights (matrices with components 1 to nCC as columns)
#' @export
optimize_bilinear_n <- function(X_list, K_list, w_list,
                                      cellTypesOfInterest,
                                      nCC = 2,
                                      max_iter = 1000,
                                      tol = 1e-5) {

  # Validate inputs based on assumption they are already subsetted
  cts <- cellTypesOfInterest
  n_mat <- length(cts)
  is_within <- (n_mat == 1)

  
  if (length(X_list) != n_mat || length(K_list) != n_mat || length(w_list) != n_mat ||
      !all(cts %in% names(X_list)) || !all(cts %in% names(K_list)) || !all(cts %in% names(w_list))) {
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
  Y_resi <- compute_Y_resi(X_list, K_list, cts)

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

# Backward compatibility aliases -- to be removed in the future
optimize_bilinear_multi <- optimize_bilinear
optimize_bilinear_multi_n <- optimize_bilinear_n