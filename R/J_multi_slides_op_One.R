
#' SkrCCA optimization for One Cell Type, Multiple Slides - First Component
#'
#' Adapts the multi-slide logic for a single cell type.
#' Assumes X_list_all and K_list_all are lists of lists, where the inner list
#' contains data/kernel for only ONE cell type per slide.
#'
#' @importFrom parallel mclapply
#' @importFrom irlba irlba
#'
#' @param X_list_all List over slides, each element is a list containing one
#'   data matrix (cell by PC) for the single cell type.
#'   Example: X_list_all[[slide_q]][[1]]
#' @param K_list_all List over slides, each element is a list containing a list
#'   with one kernel matrix (cell by cell) for the single cell type.
#'   Example: K_list_all[[slide_q]][[1]][[1]]
#' @param max_iter Maximum number of iterations.
#' @param tol Tolerance for convergence check.
#' @param n_cores Number of cores for potential parallel computation via mclapply.
#'
#' @return A list containing one element: the first weight vector component (w1)
#'   as a single-column matrix.
#' @noRd
optimize_bilinear_multi_slides_one <- function(X_list_all, K_list_all,
                                                    max_iter = 1000, tol = 1e-5,
                                                    n_cores = 1) {

  n_slides <- length(X_list_all)
  if (n_slides < 2) {
    stop("Use optimize_bilinear_multi_One for a single slide.")
  }
  if (length(K_list_all) != n_slides) {
    stop("X_list_all and K_list_all must have the same length (N of slides).")
  }

  # Calculate update contribution from each slide q: t(Xq) %*% Kq %*% Xq %*% w1_old
  # This is equivalent to Yq %*% w1_old where Yq = t(Xq) %*% Kq %*% Xq
  update_Y_q <- mclapply(1:n_slides, function(q) {
      Xq <- X_list_all[[q]][[1]]
      Kq <- K_list_all[[q]][[1]][[1]]
      Yq <- crossprod(Xq, Kq %*% Xq) # Efficient calculation t(Xq) %*% (Kq %*% (Xq %*% w1_old))
      return(Yq)
    }, mc.cores = n_cores)

    # Aggregate contributions across slides
    sum_Y_mat <- Reduce("+", update_Y_q)

    # numerical solution by taking the right sigular vector
    w1 <- irlba(sum_Y_mat, nv = 1, maxit = max_iter)$v[, 1, drop = FALSE]

  # Return result as a list containing the single weight matrix
  w_list <- list()
  w_list[[1]] <- matrix(w1, ncol = 1) # Ensure matrix format
  return(w_list)
}

#' SkrCCA optimization for One Cell Type, Multiple Slides - First Component
#'
#' Adapts the multi-slide logic for a single cell type.
#' Assumes X_list_all and K_list_all are lists of lists, where the inner list
#' contains data/kernel for only ONE cell type per slide.
#'
#' @importFrom parallel mclapply
#' @importFrom irlba irlba
#'
#' @param X_list_all List over slides, each element is a list containing one
#'   data matrix (cell by PC) for the single cell type.
#'   Example: X_list_all[[slide_q]][[1]]
#' @param K_list_all List over slides, each element is a list containing a list
#'   with one kernel matrix (cell by cell) for the single cell type.
#'   Example: K_list_all[[slide_q]][[1]][[1]]
#' @param max_iter Maximum number of iterations.
#' @param tol Tolerance for convergence check.
#' @param n_cores Number of cores for potential parallel computation via mclapply.
#'
#' @return A list containing one element: the first weight vector component (w1)
#'   as a single-column matrix.
#' @noRd
optimize_bilinear_multi_slides_one_iter <- function(X_list_all, K_list_all,
                                         max_iter = 1000, tol = 1e-5,
                                         n_cores = 1) {

  n_slides <- length(X_list_all)
  if (n_slides < 2) {
    stop("Use optimize_bilinear_multi_One for a single slide.")
  }
  if (length(K_list_all) != n_slides) {
    stop("X_list_all and K_list_all must have the same length (N of slides).")
  }


  # --- Initialization ---
  # Initialize w1 using SVD on the first slide's data
  w1 <- tryCatch(irlba(X_list_all[[1]][[1]], nv = 1,
                       right_only = TRUE)$v[, 1, drop = FALSE],
                      error = function(e) {
                        stop("SVD failed for initialization on slide 1: ",
                             e$message)
                      })

  # --- Iterative Refinement ---
  iter <- 0
  while (iter <= max_iter) {
    w1_old <- w1

    # Calculate update contribution from each slide q: t(Xq) %*% Kq %*% Xq %*% w1_old
    # This is equivalent to Yq %*% w1_old where Yq = t(Xq) %*% Kq %*% Xq
    update_contributions <- mclapply(1:n_slides, function(q) {
      Xq <- X_list_all[[q]][[1]]
      Kq <- K_list_all[[q]][[1]][[1]]
      Yq_w1 <- crossprod(Xq, Kq %*% (Xq %*% w1_old)) # Efficient calculation t(Xq) %*% (Kq %*% (Xq %*% w1_old))
      return(Yq_w1)
    }, mc.cores = n_cores)

    # Aggregate contributions across slides [cite: 57]
    w1_update_vec <- Reduce("+", update_contributions)

    # Normalize and update w1
    w1 <- normalize_vec(w1_update_vec)

    # Check convergence
    diff_val <- max(abs(w1 - w1_old))

    if (diff_val <= tol) {
      print(paste("Convergence reached at", iter, "iterations (Max diff =", sprintf("%.3e", diff_val), ")"))
      break
    }
    iter <- iter + 1
  } # end while

  if (iter > max_iter) {
    warning("Maximum number of iterations reached without convergence.")
  }

  # Return result as a list containing the single weight matrix
  w_list <- list()
  w_list[[1]] <- matrix(w1, ncol = 1) # Ensure matrix format
  return(w_list)
}


#' Helper function for optimizing multiple components (One Cell Type, Multiple Slides)
#'
#' Used within optimize_bilinear_multi_n_slides_one. Iteratively refines a weight
#' vector using pre-computed and aggregated residual matrices.
#'
#' @param w1_new Initial guess for the new weight vector (single column matrix).
#' @param Y_resi_sum Aggregated residual matrix (Sum over slides: Sum_q[ Yq - deflation_terms ]).
#' @param n_features Number of features.
#' @param max_iter Maximum iterations for refinement.
#' @param tol Tolerance for refinement convergence.
#'
#' @return Updated weight vector (single column matrix).
#' @noRd
bilinear_w_from_Y_resi_multi_slides_one <- function(w1_new, Y_resi_sum,
                                              n_features, max_iter, tol) {

  if (!is.matrix(w1_new) || ncol(w1_new) != 1 || nrow(w1_new) != n_features) {
    stop("Input w1_new is not a single-column matrix with n_features rows.")
  }
  if (!is.matrix(Y_resi_sum) || nrow(Y_resi_sum) != n_features || ncol(Y_resi_sum) != n_features) {
    stop("Input Y_resi_sum is not a square matrix with n_features dimensions.")
  }

  w1 <- w1_new # Start with the initial guess

  # Iterative refinement based on the aggregated residual matrix
  iter <- 0
  while (iter < max_iter) {
    w1_old <- w1

    # Update using the aggregated residual matrix
    # This mirrors the update step in bilinear_w_from_Y_resi_One
    # and bilinear_w_from_Y_resi_slides, but simplified for one cell type
    update_vec <- Y_resi_sum %*% w1_old
    w1 <- normalize_vec(update_vec)

    # Check convergence
    diff_val <- max(abs(w1 - w1_old))

    if (diff_val <= tol) {
      print(paste("Helper convergence reached at", iter, "iterations."))
      break
    }
    iter <- iter + 1
  } # end while

  if (iter == max_iter) {
    warning("Fun bilinear_w_from_Y_resi_multi_slides_one reached max iterations.")
  }

  return(matrix(w1, ncol = 1)) # Ensure matrix format
}


#' Run skrCCA for One Cell Type, Multiple Slides - Subsequent Components (nCC > 1)
#'
#' Calculates components 2 to nCC using deflation.
#'
#' @importFrom stats setNames
#' @importFrom irlba irlba
#'
#' @param X_list_all List over slides of lists containing the single data matrix.
#' @param K_list_all List over slides of lists of lists containing the single kernel matrix.
#' @param w_list A list containing the previously computed weight components for the
#'   single cell type (as a matrix w_list[[1]] with components as columns).
#' @param nCC Total number of canonical vectors desired (must be >= 2).
#' @param max_iter Maximum number of iterations for the helper refinement function.
#' @param tol Tolerance of accuracy for the helper refinement function.
#'
#' @return The input w_list updated with new components appended as columns to w_list[[1]].
#' @noRd
optimize_bilinear_multi_n_slides_one <- function(X_list_all, K_list_all, w_list,
                                           nCC = 2,
                                           max_iter = 1000,
                                           tol = 1e-5) {

  n_slides <- length(X_list_all)
  if (n_slides < 1) stop("Need at least one slide.") # Allow 1 slide for consistency, though deflation is trivial then.
  if (length(K_list_all) != n_slides) {
    stop("X_list_all and K_list_all must have the same length (number of slides).")
  }

  # --- Input Validation ---
  if (length(w_list) != 1 || !is.list(w_list)) {
    stop("Input w_list should be a list containing exactly one element.")
  }
  w1_matrix <- w_list[[1]]
  if (!is.matrix(w1_matrix)) {
    stop("w_list[[1]] must be a matrix containing previous components as columns.")
  }
  k_start <- ncol(w1_matrix) # Number of components already computed
                             ## we expect this value to be always 1 based on
                             ## our current workflow
  n_features <- nrow(w1_matrix)

  if (k_start < 1) stop("Input w_list[[1]] must contain at least the first component.")
  if (nCC <= k_start) {
    stop(paste("nCC (", nCC, ") must be greater than the number of components already in w_list (", k_start, ")."))
  }

  # Verify feature consistency in X_list_all
  for (q in 1:n_slides) {
    if(length(X_list_all[[q]]) != 1 || !is.matrix(X_list_all[[q]][[1]]) ||
       ncol(X_list_all[[q]][[1]]) != n_features) {
      stop(paste("X matrix in slide", q, "is invalid or has inconsistent feature count."))
    }
    if(length(K_list_all[[q]]) != 1 || length(K_list_all[[q]][[1]]) != 1 ||
       !is.matrix(K_list_all[[q]][[1]][[1]])) {
      stop(paste("K matrix structure in slide", q, "is invalid."))
    }
  }

  # --- Deflation and Optimization Loop ---

  # Initialize list to store residual matrices Y_resi for each slide
  # Y_resi_all[[q]] will store the residual Y matrix for slide q after deflation
  Y_resi_all <- vector("list", length = n_slides)

  # Loop through components to be calculated (from k_start + 1 up to nCC)
  # qq iterates through the component index *to be used for deflation* [cite: 32]
  for (qq in k_start:(nCC - 1)) {

    # If qq == 1, we need to calculate a new component.
    # If qq <= k_start, we just need to perform deflation using the existing components.

    w1_qq <- w1_matrix[, qq, drop = FALSE] # The qq-th component weight vector

    # Step 1: Update Y_resi for each slide using component qq for deflation
    for (q in 1:n_slides) {
      Xq <- X_list_all[[q]][[1]]
      Kq <- K_list_all[[q]][[1]][[1]]

      # Get the Y matrix (original or residual from previous deflation)
      if (qq == 1) {
        # Calculate Y = t(X) %*% K %*% X for the first deflation pass
        Yq <- crossprod(Xq, Kq %*% Xq)
      } else {
        # Use the Y residual from the previous deflation step (qq-1)
        Yq <- Y_resi_all[[q]]
        if(is.null(Yq)) stop(paste("Y_resi missing for slide", q, "before deflating with component", qq))
      }

      # Apply deflation formula: Y_new = Y_old - (w^T Y_old w) * (w w^T) [cite: 48, 66]
      deflation_scalar <- (t(w1_qq) %*% Yq %*% w1_qq)[1, 1]
      deflation_term <- deflation_scalar * (w1_qq %*% t(w1_qq))
      Y_resi_all[[q]] <- Yq - deflation_term

    } # End loop over slides q for deflation

    # --- Optimization for the *next* component (qq+1) ---
    # Step 2: Aggregate the residual matrices across all slides
    # Similar aggregation as in bilinear_w_from_Y_resi_slides
    Y_resi_sum <- Reduce("+", Y_resi_all)

    # Step 3: Initialize w_list_new for component qq+1
    # Use SVD on the aggregated residual matrix, similar to[cite: 67], adapted for one type
    # Use tryCatch for robustness
    w1_new_init <- tryCatch(irlba(Y_resi_sum, nv = 1,
                                  right_only = TRUE)$v[, 1, drop = FALSE],
                            error = function(e) {
                              stop(paste("SVD for initializing component", qq + 1,
                                         "failed on aggregated Y_resi: ", e$message))
                            })

      # Step 4: Iterative refinement using the helper function
      w1_qq_plus_1 <- bilinear_w_from_Y_resi_multi_slides_one(
        w1_new = w1_new_init,
        Y_resi_sum = Y_resi_sum,
        n_features = n_features,
        max_iter = max_iter,
        tol = tol)

      # Step 5: Add the new component (qq+1) to w_list[[1]]
      # Directly modify the input w_list matrix
      w1_matrix <- cbind(w1_matrix, w1_qq_plus_1)

  } # End component loop qq

  w_list[[1]] <- w1_matrix # Update the matrix in the list

  # Return the updated list (modified w_list[[1]])
  return(w_list)
}
