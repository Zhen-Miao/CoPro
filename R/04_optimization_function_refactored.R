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
    init_v <- rep(1 / ncol(X_list[[ct]]), ncol(X_list[[ct]]))
    svd_result <- tryCatch(irlba(X_list[[ct]], nv = 1, right_only = TRUE, v = init_v),
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
    # Power iteration of a symmetric eigenvalue problem has sign ambiguity:
    # if w is a fixed point so is -w, and the iterate may flip sign between
    # sweeps. Take min(||w_new - w_old||_inf, ||w_new + w_old||_inf) so a
    # pure sign flip reads as zero change rather than a spurious "change of
    # 2" that would prevent convergence forever.
    diff_fwd <- max(abs(w_list_new[[ct]] - w_list_old[[ct]]))
    diff_flip <- max(abs(w_list_new[[ct]] + w_list_old[[ct]]))
    diff_val <- min(diff_fwd, diff_flip)
    if(is.nan(diff_val)) diff_val <- 0
    if(diff_val > current_max_diff) {
      current_max_diff <- diff_val
    }
  }
  return(current_max_diff)
}

#' SkrCCA optimization function for multiple groups (Single Slide) - First Component
#' Uses flat kernel structure for consistent data access
#'
#' @param X_list Named list of data matrices (cell by PC matrix)
#' @param flat_kernels Flat list of kernel matrices with names like "kernel|sigma0.1|TypeA|TypeB"
#' @param sigma Sigma value (numeric)
#' @param max_iter Maximum number of iterations
#' @param tol tolerance of accuracy
#' @param step_size Step size for damped power iteration. Default 1 (standard
#'   power iteration). Values in (0,1) blend old and new weights for smoother
#'   convergence, which can help with many cells or many CCs.
#' @param sdev2_list Optional named list of squared standard deviations per
#'   cell type, used for weighted normalization when \code{scalePCs = FALSE}.
#'   Default \code{NULL} (unweighted).
#'
#' @return Named list `w_list` containing the first weight vector component.
#' @export
optimize_bilinear <- function(X_list, flat_kernels, sigma, max_iter = 1000,
                              tol = 1e-5, step_size = 1,
                              sdev2_list = NULL) {

  # Validate step_size
  if (!is.numeric(step_size) || length(step_size) != 1 || step_size <= 0 || step_size > 1) {
    stop("step_size must be a single numeric value in (0, 1]")
  }

  cell_types <- names(X_list)
  if (is.null(cell_types)) stop("Input X_list must be a named list.")
  feature_counts <- stats::setNames(
    vapply(X_list[cell_types], ncol, integer(1)), cell_types
  )

  # Precompute small PC-space operator matrices once, then run power iteration
  # using Y_ij %*% w_j. This avoids repeated X' K X products per iteration.
  Y_resi <- compute_Y_resi(X_list, flat_kernels, sigma, cell_types, slide = NULL)

  # The one-type problem is a symmetric Rayleigh-quotient problem, while the
  # two-type problem is an ordinary singular-vector variational problem. Both
  # have exact direct solutions and should not enter power iteration.
  if (length(cell_types) == 1L) {
    return(solve_one_type_eigen(
      Y_resi, cell_types, nCC = 1L, sdev2_list = sdev2_list
    ))
  }

  if (length(cell_types) == 2L) {
    return(solve_two_type_svd(
      Y_resi, cell_types, nCC = 1L, sdev2_list = sdev2_list
    ))
  }

  # Initialize the non-SVD cases.
  w_list <- initialize_weights_svd(X_list, cell_types)

  w_list <- bilinear_w_from_Y_resi(
    w_list_new = w_list,
    Y_resi = Y_resi,
    n_features = feature_counts,
    max_iter = max_iter,
    tol = tol,
    step_size = step_size,
    sdev2_list = sdev2_list
  )

  # Ensure final format is list of single-column matrices
  for (ct in cell_types) {
    if (!is.matrix(w_list[[ct]]) || ncol(w_list[[ct]]) != 1) {
      w_list[[ct]] <- matrix(w_list[[ct]], ncol = 1)
    }
  }
  return(w_list)
}

#' Compute PC-space operator matrices
#' @param X_list Data matrices
#' @param flat_kernels Flat kernel matrices
#' @param sigma Sigma value
#' @param cell_types Cell type names
#' @param slide Slide ID (NULL for single slide)
#' @return Y_resi structure with entries \code{Y_ij = X_i' K_ij X_j}
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
    Y_resi[[ct]][[ct]] <- compute_symmetric_Y(X, K)
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
      Y_ij <- .kernelXKY(X_list[[i]], K12, X_list[[j]])
      
      Y_resi[[i]][[j]] <- Y_ij
      Y_resi[[j]][[i]] <- t(Y_ij)
    }
  }
  
  return(Y_resi)
}

#' Form the symmetric one-type PC-space operator efficiently
#'
#' Computes \code{Y = X' K X} by first applying the (possibly triangularly
#' stored) symmetric kernel to all PC columns. This ordering uses one
#' \code{n_cells x n_PC} temporary and avoids either expanding a sparse
#' symmetric kernel or forming a much larger \code{n_PC x n_cells} product.
#' The final explicit symmetrization is a numerical no-op for an exactly
#' symmetric kernel and makes the quadratic objective well-defined if a user
#' supplied kernel has small asymmetry.
#'
#' @param X Cell-by-PC score matrix.
#' @param K Within-cell-type spatial kernel.
#' @return Symmetric base numeric matrix in PC space.
#' @noRd
compute_symmetric_Y <- function(X, K) {
  Y <- .kernelXKY(X, K, X)
  (Y + t(Y)) * 0.5
}

#' Solve the exact one-cell-type skrCCA problem by symmetric eigendecomposition
#'
#' For one cell type, skrCCA maximizes \code{w' Y w} subject to a unit-norm
#' constraint. Rayleigh--Ritz therefore selects the largest algebraic
#' eigenvalues, not the eigenvalues of largest absolute magnitude. When
#' \code{sdev2_list} is supplied, this is the generalized eigenproblem
#' \code{Y w = lambda D w}; it is solved through the symmetric transform
#' \code{D^-1/2 Y D^-1/2}.
#'
#' @param Y_resi PC-space operator structure.
#' @param cell_types The single cell-type name.
#' @param nCC Number of axes to return.
#' @param sdev2_list Optional diagonal CCA metric.
#' @return Named one-element list containing the weight matrix.
#' @noRd
solve_one_type_eigen <- function(Y_resi, cell_types, nCC = 1L,
                                 sdev2_list = NULL) {
  if (length(cell_types) != 1L) {
    stop("solve_one_type_eigen requires exactly one cell type")
  }
  if (!is.numeric(nCC) || length(nCC) != 1L || nCC < 1L ||
      nCC != as.integer(nCC)) {
    stop("nCC must be a positive integer")
  }

  ct <- cell_types[[1L]]
  Y <- as.matrix(Y_resi[[ct]][[ct]])
  if (nrow(Y) != ncol(Y)) {
    stop("The one-type PC-space operator must be square")
  }
  if (nCC > nrow(Y)) {
    stop("nCC cannot exceed the dimension of the one-type PC-space operator (",
         nrow(Y), ")")
  }

  # w'Yw depends only on the symmetric part, including when a caller supplies
  # a row- or column-normalized kernel that is not exactly symmetric.
  Y <- (Y + t(Y)) * 0.5
  inv_sqrt_d <- NULL
  if (!is.null(sdev2_list)) {
    d <- sdev2_list[[ct]]
    if (length(d) != nrow(Y) || any(!is.finite(d)) || any(d <= 0)) {
      stop("sdev2_list for ", ct,
           " must contain one finite positive value per PC")
    }
    inv_sqrt_d <- 1 / sqrt(d)
    Y <- sweep(sweep(Y, 1L, inv_sqrt_d, "*"),
               2L, inv_sqrt_d, "*")
    Y <- (Y + t(Y)) * 0.5
  }

  decomp <- eigen(Y, symmetric = TRUE)
  W <- decomp$vectors[, seq_len(nCC), drop = FALSE]
  if (!is.null(inv_sqrt_d)) {
    W <- sweep(W, 1L, inv_sqrt_d, "*")
  }
  setNames(list(W), ct)
}

#' Solve the exact two-cell-type skrCCA problem by SVD
#'
#' For two cell types, skrCCA maximizes \code{w1' Y12 w2} subject to one
#' unit-norm constraint per weight vector. This is exactly the singular-vector
#' variational problem for \code{Y12}, so one decomposition gives all axes.
#' When \code{sdev2_list} is supplied, the diagonal CCA metrics are removed
#' before the SVD and restored afterwards.
#'
#' @param Y_resi PC-space operator matrices.
#' @param cell_types The two cell-type names.
#' @param nCC Number of axes to return.
#' @param sdev2_list Optional diagonal CCA metrics.
#' @return Named list of weight matrices, one per cell type.
#' @noRd
solve_two_type_svd <- function(Y_resi, cell_types, nCC = 1L,
                               sdev2_list = NULL) {
  if (length(cell_types) != 2L) {
    stop("solve_two_type_svd requires exactly two cell types")
  }

  ct1 <- cell_types[[1L]]
  ct2 <- cell_types[[2L]]
  Y12 <- as.matrix(Y_resi[[ct1]][[ct2]])
  max_axes <- min(dim(Y12))
  if (nCC > max_axes) {
    stop("nCC cannot exceed the dimension of the two-type PC-space operator (",
         max_axes, ")")
  }

  if (!is.null(sdev2_list)) {
    inv_sqrt_d1 <- 1 / sqrt(sdev2_list[[ct1]])
    inv_sqrt_d2 <- 1 / sqrt(sdev2_list[[ct2]])
    Y12 <- sweep(sweep(Y12, 1L, inv_sqrt_d1, "*"),
                 2L, inv_sqrt_d2, "*")
  }

  decomp <- svd(Y12, nu = nCC, nv = nCC)
  W1 <- decomp$u[, seq_len(nCC), drop = FALSE]
  W2 <- decomp$v[, seq_len(nCC), drop = FALSE]

  if (!is.null(sdev2_list)) {
    W1 <- sweep(W1, 1L, inv_sqrt_d1, "*")
    W2 <- sweep(W2, 1L, inv_sqrt_d2, "*")
  }

  setNames(list(W1, W2), cell_types)
}

#' Check whether supplied first-axis weights equal the exact two-type solution
#' @noRd
matches_two_type_first_axis <- function(w_list, svd_weights, cell_types,
                                        sdev2_list = NULL,
                                        tolerance = 1e-4) {
  cosines <- vapply(cell_types, function(ct) {
    supplied <- w_list[[ct]][, 1L, drop = TRUE]
    exact <- svd_weights[[ct]][, 1L, drop = TRUE]
    if (is.null(sdev2_list)) {
      cosine <- sum(supplied * exact) /
        sqrt(sum(supplied^2) * sum(exact^2))
    } else {
      metric <- sdev2_list[[ct]]
      cosine <- sum(metric * supplied * exact) /
        sqrt(sum(metric * supplied^2) * sum(metric * exact^2))
    }
    cosine
  }, numeric(1))

  # Both blocks may flip together without changing the bilinear direction; an
  # opposite sign in only one block represents the negative singular pair and
  # must remain on the conditional sequential path.
  all(is.finite(cosines)) &&
    all(abs(abs(cosines) - 1) <= tolerance) &&
    prod(sign(cosines)) > 0
}

#' Check whether supplied first-axis weights equal the exact one-type solution
#' @noRd
matches_one_type_first_axis <- function(w_list, eigen_weights, cell_types,
                                        sdev2_list = NULL,
                                        tolerance = 1e-4) {
  ct <- cell_types[[1L]]
  supplied <- w_list[[ct]][, 1L, drop = TRUE]
  exact <- eigen_weights[[ct]][, 1L, drop = TRUE]
  if (is.null(sdev2_list)) {
    cosine <- sum(supplied * exact) /
      sqrt(sum(supplied^2) * sum(exact^2))
  } else {
    metric <- sdev2_list[[ct]]
    cosine <- sum(metric * supplied * exact) /
      sqrt(sum(metric * supplied^2) * sum(metric * exact^2))
  }
  is.finite(cosine) && abs(abs(cosine) - 1) <= tolerance
}

#' Apply deflation to Y_resi
#' @param Y_resi Current Y_resi structure
#' @param w_list Weight list
#' @param qq Component index used for deflation
#' @param cell_types Cell type names
#' @param deflation How to remove the (w1, w2) direction from Y. `"rank1"`
#'   (default) subtracts the rank-1 bilinear component
#'   `Y - (w1^T Y w2) w1 w2^T`; when `w1, w2` are singular vectors of
#'   `Y` (the ordinary two-cell-type case) this equals the full projection.
#'   `"projection"` applies the full
#'   orthogonal projection `(I - u u^T) Y (I - v v^T)` with unit `u, v` -- the
#'   Freedman-Lane / Legendre-Oksanen-ter Braak (2011) residualization. The two
#'   agree on the observed `Y` but differ on a permuted `Y` (where `w1, w2` are
#'   no longer its singular vectors), which is exactly where the conditional
#'   null lives. With `sdev2_list` (weighted / `scalePCs = FALSE` deflation),
#'   `"projection"` uses the oblique projection
#'   `(I - p1 q1^T) Y (I - q2 p2^T)` with `q = w / sqrt(w' d w)` and `p = d q`
#'   (`d = sdev^2`), the exact image of the whitened orthogonal projection under
#'   `X~ = X diag(sdev)^-1`, `w~ = diag(sdev) w`. Using it keeps every canonical
#'   axis identical to the whitened (`scalePCs = TRUE`) run, so `scalePCs` stays
#'   a pure reparametrization even for >2 cell types.
#' @return Updated Y_resi
#' @noRd
apply_deflation <- function(Y_resi, w_list, qq, cell_types, sdev2_list = NULL,
                            deflation = c("rank1", "projection")) {
  deflation <- match.arg(deflation)
  n_mat <- length(cell_types)
  is_within <- (n_mat == 1)

  ## Projection deflation. With no metric (d = NULL) this is the full
  ## orthogonal projection (I - uu^T) Y (I - vv^T) with unit u, v. With a
  ## diagonal CCA metric d = sdev^2 it is the oblique projection
  ## (I - p1 q1^T) Y (I - q2 p2^T) with q = w / sqrt(w' d w) and p = d q -- the
  ## exact image of the whitened orthogonal projection under X~ = X diag(sdev)^-1,
  ## w~ = diag(sdev) w. This makes the weighted (scalePCs = FALSE) axes equal the
  ## whitened (scalePCs = TRUE) axes, so scalePCs stays a pure reparametrization.
  proj_deflate <- function(Y1, w1, w2, d1 = NULL, d2 = NULL) {
    if (is.null(d1)) {
      q1 <- w1 / sqrt(sum(w1^2)); p1 <- q1
    } else {
      q1 <- w1 / sqrt(sum(w1 * d1 * w1)); p1 <- d1 * q1
    }
    if (is.null(d2)) {
      q2 <- w2 / sqrt(sum(w2^2)); p2 <- q2
    } else {
      q2 <- w2 / sqrt(sum(w2 * d2 * w2)); p2 <- d2 * q2
    }
    lam <- as.numeric(crossprod(q1, Y1 %*% q2))
    Y1 - p1 %*% crossprod(q1, Y1) - (Y1 %*% q2) %*% t(p2) + lam * (p1 %*% t(p2))
  }

  if (is_within) {
    ct <- cell_types[1]
    Y1 <- Y_resi[[ct]][[ct]]
    w1 <- w_list[[ct]][, qq, drop = FALSE]

    if (deflation == "projection") {
      d_ct <- if (is.null(sdev2_list)) NULL else sdev2_list[[ct]]
      Y_resi[[ct]][[ct]] <- proj_deflate(Y1, w1, w1, d_ct, d_ct)
    } else {
      deflation_scalar <- (t(w1) %*% Y1 %*% w1)[1, 1]
      if (!is.null(sdev2_list)) {
        # Weighted deflation: project out Dw direction where D = diag(sdev^2)
        Dw1 <- w1 * sdev2_list[[ct]]
        deflation_term <- deflation_scalar * (Dw1 %*% t(Dw1))
      } else {
        deflation_term <- deflation_scalar * (w1 %*% t(w1))
      }
      Y_resi[[ct]][[ct]] <- Y1 - deflation_term
    }
  } else {
    pair_cell_types <- combn(cell_types, 2)
    for (pp in seq_len(ncol(pair_cell_types))) {
      i <- pair_cell_types[1, pp]
      j <- pair_cell_types[2, pp]

      w1 <- w_list[[i]][, qq, drop = FALSE]
      w2 <- w_list[[j]][, qq, drop = FALSE]
      Y1 <- Y_resi[[i]][[j]]

      if (deflation == "projection") {
        d_i <- if (is.null(sdev2_list)) NULL else sdev2_list[[i]]
        d_j <- if (is.null(sdev2_list)) NULL else sdev2_list[[j]]
        Y_resi[[i]][[j]] <- proj_deflate(Y1, w1, w2, d_i, d_j)
      } else {
        deflation_scalar <- (t(w1) %*% Y1 %*% w2)[1, 1]
        if (!is.null(sdev2_list)) {
          Dw1 <- w1 * sdev2_list[[i]]
          Dw2 <- w2 * sdev2_list[[j]]
          deflation_term <- deflation_scalar * (Dw1 %*% t(Dw2))
        } else {
          deflation_term <- deflation_scalar * (w1 %*% t(w2))
        }
        Y_resi[[i]][[j]] <- Y1 - deflation_term
      }
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
#' @param n_features Number of features. May be a named vector with one value
#'   per cell type; retained for backward compatibility because dimensions are
#'   now derived from the weight vectors themselves.
#' @param max_iter Maximum iterations
#' @param tol Tolerance
#' @param step_size Step size for damped power iteration (default 1)
#'
#' @return Updated named list of weight vectors
#' @noRd
bilinear_w_from_Y_resi <- function(w_list_new, Y_resi,
                                   n_features, max_iter, tol,
                                   step_size = 1,
                                   sdev2_list = NULL) {

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
      sd2 <- if (!is.null(sdev2_list)) sdev2_list[[ct]] else NULL
      if (step_size < 1) {
        w_list_new[[ct]] <- normalize_vec_weighted((1 - step_size) * w_list_old[[ct]] + step_size * normalize_gradient_weighted(w_update, sd2), sd2)
      } else {
        w_list_new[[ct]] <- normalize_gradient_weighted(w_update, sd2)
      }
    } else {
      # Standard multi-cell-type case - use direct accumulation
      for (ct_i in cell_types) {
        # We can directly accumulate since we're just summing
        # Cell types can legitimately retain different PCA ranks when their
        # sample sizes differ. Allocate from the current block rather than the
        # first cell type's feature count.
        w_i_update_vec <- matrix(0, nrow = nrow(w_list_new[[ct_i]]), ncol = 1)
        for (ct_j in cell_types) {
          if (ct_i == ct_j) next
          w2 <- w_list_new[[ct_j]]
          Y <- Y_resi[[ct_i]][[ct_j]]
          if(is.null(Y)) stop(paste("Missing Y_resi matrix for pair:", ct_i, ct_j))
          w_i_update_vec <- w_i_update_vec + Y %*% w2
        }
        sd2 <- if (!is.null(sdev2_list)) sdev2_list[[ct_i]] else NULL
        if (step_size < 1) {
          w_list_new[[ct_i]] <- normalize_vec_weighted((1 - step_size) * w_list_old[[ct_i]] + step_size * normalize_gradient_weighted(w_i_update_vec, sd2), sd2)
        } else {
          w_list_new[[ct_i]] <- normalize_gradient_weighted(w_i_update_vec, sd2)
        }
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
#' @param step_size Step size for damped power iteration (default 1)
#' @param sdev2_list Optional named list of squared standard deviations per
#'   cell type for weighted normalization. Default \code{NULL}.
#'
#' @return A named list of weights (matrices with components 1 to nCC as columns)
#' @export
optimize_bilinear_n <- function(X_list, flat_kernels, sigma, w_list,
                                      cellTypesOfInterest,
                                      nCC = 2,
                                      max_iter = 1000,
                                      tol = 1e-5,
                                      step_size = 1,
                                      sdev2_list = NULL) {

  # Validate inputs based on assumption they are already subsetted
  cts <- cellTypesOfInterest
  n_mat <- length(cts)
  is_within <- (n_mat == 1)

  
  if (length(X_list) != n_mat || length(w_list) != n_mat ||
      !all(cts %in% names(X_list)) || !all(cts %in% names(w_list))) {
    stop("Input lists length or names do not match cellTypesOfInterest.")
  }
  feature_counts <- stats::setNames(
    vapply(X_list[cts], ncol, integer(1)), cts
  )

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

  # In an ordinary one-type run, one symmetric eigendecomposition returns all
  # axes. Preserve the conditional route when the caller supplied a different
  # first direction (for example, a transferred weight).
  if (n_mat == 1L) {
    eigen_weights <- solve_one_type_eigen(
      Y_resi, cts, nCC = nCC, sdev2_list = sdev2_list
    )
    if (matches_one_type_first_axis(
      w_list, eigen_weights, cts, sdev2_list = sdev2_list
    )) {
      return(eigen_weights)
    }
  }

  # The ordinary two-type run has an exact all-axis SVD. Preserve the
  # sequential path when a supplied/transferred first axis does not equal the
  # leading singular direction, because later axes are then conditional on
  # that external direction.
  if (n_mat == 2L) {
    svd_weights <- solve_two_type_svd(
      Y_resi, cts, nCC = nCC, sdev2_list = sdev2_list
    )
    if (matches_two_type_first_axis(
      w_list, svd_weights, cts, sdev2_list = sdev2_list
    )) {
      return(svd_weights)
    }
  }

  # Loop to compute components k_start + 1 up to nCC
  for (qq in k_start:(nCC - 1)) {
    # Step 1: Apply deflation using component qq
    # Rank-one subtraction is identical to projection for two-type singular
    # vectors, but not for the multi-set (>2 type) stationary equations. Full
    # projection is required there to keep later axes orthogonal within every
    # cell type. apply_deflation() implements the weighted (oblique) projection
    # for the scalePCs = FALSE case, so both metrics use projection and give the
    # same axes -- scalePCs stays a pure reparametrization.
    deflation_method <- if (n_mat == 2L) "rank1" else "projection"
    Y_resi <- apply_deflation(
      Y_resi, w_list, qq, cts, sdev2_list,
      deflation = deflation_method
    )

    if (n_mat == 1L) {
      # The residual remains a symmetric generalized eigenproblem.
      w_list_qq_plus_1 <- solve_one_type_eigen(
        Y_resi, cts, nCC = 1L, sdev2_list = sdev2_list
      )
    } else {
      # Step 2: Initialize w_list_new for component qq+1
      w_list_new <- initialize_next_component(Y_resi, cts)

      # Step 3: Iterative refinement using the helper function
      w_list_qq_plus_1 <- bilinear_w_from_Y_resi(
        w_list_new = w_list_new,
        Y_resi = Y_resi,
        n_features = feature_counts,
        max_iter = max_iter,
        tol = tol,
        step_size = step_size,
        sdev2_list = sdev2_list)
    }

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
  
  # Feature counts must be stable across slides for a given cell type, but
  # need not be identical between cell types.
  feature_counts <- stats::setNames(
    vapply(cell_types, function(ct) ncol(X_list_all[[1]][[ct]]), integer(1)),
    cell_types
  )
  
  # Validate dimensions across all slides and cell types
  for (q in seq_len(n_slides)) {
    for (ct in cell_types) {
      if (!ct %in% names(X_list_all[[q]])) {
        stop(paste("Cell type", ct, "missing in slide", q))
      }
      if (!is.matrix(X_list_all[[q]][[ct]]) ||
          ncol(X_list_all[[q]][[ct]]) != feature_counts[[ct]]) {
        stop(paste("Invalid or inconsistent feature count in slide", q, "cell type", ct))
      }
    }
  }
  
  return(list(
    n_slides = n_slides,
    cell_types = cell_types,
    n_cell_types = n_cell_types,
    n_features = feature_counts,
    feature_counts = feature_counts
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
    return(compute_symmetric_Y(X, K))
  } else {
    # Between-cell-type case
    X1 <- X_list[[ct_i]]
    X2 <- X_list[[ct_j]]
    K12 <- get_kernel_matrix_flat(flat_kernels, sigma, ct_i, ct_j, slide)
    return(.kernelXKY(X1, K12, X2))
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
#' @param step_size Step size for damped power iteration (default 1).
#'   Values in (0,1) blend old and new weights for smoother convergence.
#' @param sdev2_list Optional named list of squared standard deviations per
#'   cell type for weighted normalization. Default \code{NULL}.
#' @return Named list of weight vectors (first component)
#' @export
optimize_bilinear_multi_slides <- function(X_list_all, flat_kernels, sigma, slides,
                                          max_iter = 1000, tol = 1e-5,
                                          n_cores = 1, direct_solve = TRUE,
                                          step_size = 1,
                                          sdev2_list = NULL) {

  # Validate step_size
  if (!is.numeric(step_size) || length(step_size) != 1 || step_size <= 0 || step_size > 1) {
    stop("step_size must be a single numeric value in (0, 1]")
  }

  # Validate inputs
  validated <- validate_multi_slide_inputs(X_list_all, NULL)  # Don't validate K_list_all anymore
  n_slides <- validated$n_slides
  cell_types <- validated$cell_types
  n_cell_types <- validated$n_cell_types
  feature_counts <- validated$feature_counts
  
  # Check if this is within-cell-type case
  is_within <- (n_cell_types == 1)
  
  # For single slide, delegate to single-slide function
  if (n_slides == 1) {
    message("Single slide detected, using single-slide optimization")
    return(optimize_bilinear(X_list_all[[1]], flat_kernels, sigma,
                            max_iter = max_iter, tol = tol,
                            step_size = step_size,
                            sdev2_list = sdev2_list))
  }
  
  # Build the small PC-space operators once for all iterative first-CC updates.
  Y_aggregate <- compute_Y_multi_slide(
    X_list_all, flat_kernels, sigma, slides, cell_types, n_cores
  )

  # Stacking samples with a block-diagonal spatial kernel is algebraically
  # equivalent to summing their PC-space operators. We deliberately keep the
  # sum and never materialize that larger kernel. The resulting two-type
  # problem has the same exact SVD solution as the single-slide case.
  if (n_cell_types == 2L) {
    return(solve_two_type_svd(
      Y_aggregate, cell_types, nCC = 1L, sdev2_list = sdev2_list
    ))
  }

  # Handle within-cell-type case with direct solution
  if (is_within && direct_solve) {
    return(solve_one_type_eigen(
      Y_aggregate, cell_types, nCC = 1L, sdev2_list = sdev2_list
    ))
  }
  
  # Initialize weights
  w_list <- initialize_weights_multi_slide(X_list_all, cell_types, use_aggregation = TRUE)

  # Iterative optimization on the aggregated small operators.
  w_list <- bilinear_w_from_Y_resi(
    w_list_new = w_list,
    Y_resi = Y_aggregate,
    n_features = feature_counts,
    max_iter = max_iter,
    tol = tol,
    step_size = step_size,
    sdev2_list = sdev2_list
  )
  
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
#' @param step_size Step size for damped power iteration (default 1)
#' @param sdev2_list Optional named list of squared standard deviations per
#'   cell type for weighted normalization. Default \code{NULL}.
#' @return Updated weight list with all components
#' @export
optimize_bilinear_n_multi_slides <- function(X_list_all, flat_kernels, sigma, slides, w_list,
                                            cellTypesOfInterest,
                                            nCC = 2, max_iter = 1000,
                                            tol = 1e-5, n_cores = 1,
                                            step_size = 1,
                                            sdev2_list = NULL) {
  
  # Validate inputs
  validated <- validate_multi_slide_inputs(X_list_all, NULL, 
                                          expected_cell_types = cellTypesOfInterest)
  n_slides <- validated$n_slides
  cell_types <- cellTypesOfInterest
  n_cell_types <- length(cell_types)
  feature_counts <- validated$feature_counts
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
  
  # Shared weights make sample aggregation exact: sum the small PC-space
  # operators once, then solve and deflate that aggregate. This is algebraically
  # identical to stacking samples with a block-diagonal spatial kernel, but we
  # never materialize that larger kernel.
  Y_resi <- compute_Y_multi_slide(
    X_list_all, flat_kernels, sigma, slides, cell_types, n_cores
  )

  if (n_cell_types == 1L) {
    eigen_weights <- solve_one_type_eigen(
      Y_resi, cell_types, nCC = nCC, sdev2_list = sdev2_list
    )
    if (matches_one_type_first_axis(
      w_list, eigen_weights, cell_types, sdev2_list = sdev2_list
    )) {
      return(eigen_weights)
    }
  }

  # Obtain all ordinary two-type axes from the same exact SVD when the supplied
  # first axis is the leading singular direction. A transferred first axis
  # deliberately falls back to conditional sequential deflation.
  if (n_cell_types == 2L) {
    svd_weights <- solve_two_type_svd(
      Y_resi, cell_types, nCC = nCC, sdev2_list = sdev2_list
    )
    if (matches_two_type_first_axis(
      w_list, svd_weights, cell_types, sdev2_list = sdev2_list
    )) {
      return(svd_weights)
    }
  }
  
  # Compute additional components
  for (qq in k_start:(nCC - 1)) {

    # Full projection is needed for multi-set axes because their stationary
    # vectors are not pairwise singular vectors. apply_deflation() implements
    # the weighted (oblique) projection for scalePCs = FALSE, so both metrics
    # use projection and give the same axes.
    deflation_method <- if (n_cell_types == 2L) "rank1" else "projection"
    Y_resi <- apply_deflation(
      Y_resi, w_list, qq, cell_types, sdev2_list,
      deflation = deflation_method
    )
    
    if (n_cell_types == 1L) {
      w_list_qq_plus_1 <- solve_one_type_eigen(
        Y_resi, cell_types, nCC = 1L, sdev2_list = sdev2_list
      )
    } else {
      # Initialize and refine the next component.
      w_list_new <- initialize_next_component(Y_resi, cell_types)

      w_list_qq_plus_1 <- bilinear_w_from_Y_resi(
        w_list_new = w_list_new,
        Y_resi = Y_resi,
        n_features = feature_counts,
        max_iter = max_iter,
        tol = tol,
        step_size = step_size,
        sdev2_list = sdev2_list
      )
    }

    # Append the new component.
    for (ct in cell_types) {
      w_list[[ct]] <- cbind(w_list[[ct]], w_list_qq_plus_1[[ct]])
    }
  }
  
  return(w_list)
}
