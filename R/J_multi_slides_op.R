## Multi-slide
#' Multi-slide optimization function for CoPro
#'
#' @importFrom parallel mclapply
#' @importFrom irlba irlba
#'
#' @param X_list_all List of lists of data matrices (cell by PC matrix), outer list for slides, inner list for cell types
#' @param K_list_all List of lists of kernel matrices, outer list for slides, inner list for cell type pairs
#' @param max_iter Maximum number of iterations
#' @param tol Tolerance of accuracy
#' @param n_cores Number of cores for parallel computation
#'

#'
#' @return `w_list`, a list of weight vectors for each cell type
#' @noRd
optimize_bilinear_multi_slides <- function(X_list_all, K_list_all,
                                           max_iter = 1000, tol = 1e-5,
                                           n_cores = 1) {
  # Number of slides
  n_slides <- length(X_list_all)

  if(n_slides < 2){
    stop("Do not run this function if number of slides is below 2")
  }

  # Check if all slides have the same cell types
  cell_types_all <- lapply(X_list_all, names)
  if (!all(sapply(cell_types_all, function(x) all(x == cell_types_all[[1]])))) {
    stop("All slides must have the same cell types.")
  }

  # Get cell types from the first slide
  cell_types <- cell_types_all[[1]]
  n_mat <- length(cell_types)
  n_features <- ncol(X_list_all[[1]][[1]])

  # Check dimensions consistency across slides
  for (q in 1:n_slides) {
    for (i in 1:n_mat) {
      if (ncol(X_list_all[[q]][[i]]) != n_features) {
        stop(paste("Inconsistent number of features in slide", q, "cell type", i))
      }
    }
  }

  # Initialize w1 and w2 by their right singular vector from the first slide
  w_list <- rep(list(), length = n_mat)
  for (i in 1:n_mat) {
    w_list[[i]] <- irlba(X_list_all[[1]][[i]], nv = 1, right_only = 1)$v[, 1, drop = FALSE]
  }

  # Iterative refinement
  iter <- 0
  while (iter <= max_iter) {
    # After going over all j, we will update i
    for (i in 1:n_mat) {
      # Initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)
      diff_i <- vector(length = n_mat)

      for (j in 1:n_mat) {
        if (j == i) {
          next
        }

        w2 <- w_list[[j]]

        # Sum up over all slides
        Yv <- mclapply(1:n_slides, function(q) {
          if (j < i) {
            K12 <- t(K_list_all[[q]][[j]][[i]])
          } else {
            K12 <- K_list_all[[q]][[i]][[j]]
          }

          # Compute Y_ij for this slide
          X1 <- X_list_all[[q]][[i]]
          X2 <- X_list_all[[q]][[j]]
          v2 <- X2 %*% w2
          Yq_val <- crossprod(X1, K12 %*% v2)
          return(Yq_val)
        }, mc.cores = n_cores)

        ## aggregate over all slides
        Yv_sum <- Reduce("+", Yv)
        w_i_left[, j] <- Yv_sum
      }

      w_i_new <- normalize_vec(rowSums(w_i_left))
      diff_i[i] <- max(abs(w_i_new - w_list[[i]]))
      w_list[[i]] <- w_i_new
    }

    if (max(diff_i) <= tol) {
      print(paste("convergence reached at", iter, "iterations", sep = " "))
      break
    }
    iter <- iter + 1
  }

  if (iter > max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  # Convert w_list into matrix
  for (j in 1:n_mat) {
    w_list[[j]] <- matrix(w_list[[j]], ncol = 1)
  }

  # Return results
  return(w_list)
}

#' Run multi-slide version of bilinear optimization to detect the third to n_th component
#'
#' @importFrom stats setNames
#' @importFrom irlba irlba
#'
#' @param X_list_all List of lists of data matrices (cell by PC matrix), outer list for slides, inner list for cell types
#' @param K_list_all List of lists of kernel matrices, outer list for slides, inner list for cell type pairs
#' @param w_list A list of weights for each data matrix, pre-computed and to be regressed out
#' @param cellTypesOfInterest A vector specifying cell types of interest
#' @param nCC Number of canonical vectors to compute, need to be greater than two, default = 2
#' @param max_iter Maximum number of iterations
#' @param tol Tolerance of accuracy
#'
#' @return A list of weights
#' @noRd
optimize_bilinear_multi_n_slides <- function(X_list_all, K_list_all, w_list,
                                             cellTypesOfInterest,
                                             nCC = 2,
                                             max_iter = 1000,
                                             tol = 1e-5) {

  # Number of slides
  n_slides <- length(X_list_all)
  n_mat <- length(cellTypesOfInterest)
  cts <- cellTypesOfInterest
  ct1 <- cts[1]
  n_features <- ncol(X_list_all[[1]][[ct1]])

  # Check if w_list has the same length as cellTypesOfInterest
  if (n_mat != length(w_list) || n_mat != length(cellTypesOfInterest)) {
    stop(paste("the input X_list, w_list, cellTypesOfInterest",
               "are of different length!"))
  }

  # Check the dimension of the w_list
  if (length(dim(w_list[[1]])) == 0) {
    stop("the input w_list must be a matrix!")
  }

  # Get all pairs for X, K, W to obtain Y_resi
  # Let Y_resi be the same structure as K for the first slide
  Y_resi_all <- vector("list", length = n_slides)
  for (q in 1:n_slides) {
    Y_resi_all[[q]] <- setNames(vector(mode = "list", length = n_mat), cts)
    for (i in cts) {
      Y_resi_all[[q]][[i]] <- setNames(vector(mode = "list", length = n_mat), cts)
    }
  }

  # Y is between any two cell types, and each time, we update the Y_resi
  # to regress out all previous canonical vectors

  for (qq in 1:(nCC - 1)) { # qq: existing w_qq to be regressed out

    # Step 1: obtain or update Y_resi for all slides
    for (q in 1:n_slides) {
      for (pp in seq_len(ncol(combn(cts, 2)))) {
        i <- combn(cts, 2)[1, pp]
        j <- combn(cts, 2)[2, pp]
        K12 <- K_list_all[[q]][[i]][[j]]

        # Compute Y_ij
        X1 <- X_list_all[[q]][[i]]
        X2 <- X_list_all[[q]][[j]]

        w1 <- w_list[[i]][, qq, drop = FALSE]
        w2 <- w_list[[j]][, qq, drop = FALSE]

        if (qq == 1) {
          Y1 <- t(X1) %*% K12 %*% X2
        } else {
          Y1 <- Y_resi_all[[q]][[i]][[j]]
        }

        Y_resi_all[[q]][[i]][[j]] <- Y1 - ((t(w1) %*% Y1 %*% w2)[1, 1] * (w1 %*% t(w2)))

        # Fill in the matrix for easier optimization
        Y_resi_all[[q]][[j]][[i]] <- t(Y_resi_all[[q]][[i]][[j]])
      }
    }

    # Step 2: initialize w_list_new by SVD using the first slide
    w_list_new <- rep(list(), length = n_mat)

    Y_sum_12 = Reduce("+", lapply(Y_resi_all, function(y) y[[1]][[2]]))
    w_list_new[[1]] <- irlba(t(Y_sum_12), nv = 1,
                             right_only = TRUE)$v[, 1, drop = FALSE]

    for (i in 2:n_mat) {
      Y_sum_1i = Reduce("+", lapply(Y_resi_all, function(y) y[[1]][[i]]))
      w_list_new[[i]] <- irlba(Y_sum_1i, nv = 1,
                               right_only = TRUE)$v[, 1, drop = FALSE]
    }

    # Step 3: Iterative refinement
    w_list_qq <- bilinear_w_from_Y_resi_slides(w_list_new = w_list_new,
                                               Y_resi_all = Y_resi_all,
                                               n_mat = n_mat,
                                               n_features = n_features,
                                               max_iter = max_iter, tol = tol)

    # Step 4: add w_list_qq to columns in w_list
    for (ii in seq_len(n_mat)) {
      w_list[[ii]] <- cbind(w_list[[ii]], w_list_qq[[ii]])
    }
  }

  # Return results
  return(w_list)
}

#' Helper function for optimizing multiple components
#'
#' @param w_list_new Initial weight vectors
#' @param Y_resi_all Residual matrices for all slides
#' @param n_mat Number of matrices/cell types
#' @param n_features Number of features
#' @param max_iter Maximum iterations
#' @param tol Tolerance
#'
#' @return Updated weight vectors
#' @noRd
bilinear_w_from_Y_resi_slides <- function(w_list_new, Y_resi_all,
                                          n_mat, n_features, max_iter, tol) {
  # Number of slides
  n_slides <- length(Y_resi_all)

  # Iterative refinement
  iter <- 0
  while (iter < max_iter) {
    # After going over all j, we will update i
    for (i in 1:n_mat) {
      # Initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)
      diff_i <- vector(length = n_mat)

      for (j in 1:n_mat) {
        if (j == i) {
          next
        }

        ## sum up all Y_ij first, to save time
        Y_sum_ij = Reduce("+", lapply(Y_resi_all, function(y) y[[i]][[j]]))
        w2 <- w_list_new[[j]]
        w_i_left[, j] <-  Y_sum_ij %*% w2
      }

      w_i_new <- normalize_vec(rowSums(w_i_left))
      diff_i[i] <- max(abs(w_i_new - w_list_new[[i]]))
      w_list_new[[i]] <- w_i_new
    }

    if (mean(diff_i) <= tol) {
      print(paste("convergence reached at", iter, "iterations", sep = " "))
      break
    }
    iter <- iter + 1
  }

  if (iter == max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  return(w_list_new)
}

