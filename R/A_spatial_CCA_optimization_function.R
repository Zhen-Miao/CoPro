

#' centering and scaling the matrix
#' @importFrom stats sd
#' @param matrix Input matrix to be column-centered
#'
#' @return centered and scaled matrix
#' @noRd
center_scale_matrix_opt <- function(matrix) {
  # Calculate column means and standard deviations
  col_means <- colMeans(matrix)
  col_sds <- apply(matrix, 2, sd)

  # Identify columns that are not full of zeros (to avoid division by zero)
  nzc <- col_sds >= 1e-5 ## non-zero columns

  # Center and scale columns with none-zero sd column
  matrix[, nzc] <- t((t(matrix[, nzc]) - col_means[nzc]) / col_sds[nzc])

  # do not scale, but still need to center columns with zero sd
  matrix[, !nzc] <- t((t(matrix[, !nzc]) - col_means[!nzc]))

  return(matrix)
}


# Function to normalize_vec a vector to have unit norm
normalize_vec <- function(v) {
  return(v / sqrt(sum(v^2)))
}

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



#' SkrCCA optimization function for multiple groups
#'
#' @param X_list List of data matrices (cell by PC matrix)
#' @param K_list Kernel matrices between any two pairs of matrices
#' @param max_iter Maximum number of iterations
#' @param tol tolerance of accuracy
#'
#' @return `w_list`, a list of weight vectors for each cell type
#' @noRd
optimize_bilinear_multi <- function(X_list, K_list, max_iter = 1000,
                                    tol = 1e-5) {
  n_mat <- length(X_list)
  n_features <- ncol(X_list[[1]])

  # Initialize w1 and w2 by their right singular vector
  w_list <- rep(list(), length = length(X_list))
  for (i in 1:n_mat) {
    w_list[[i]] <- svd(X_list[[i]])$v[, 1, drop = FALSE]
  }

  # Iterative refinement
  iter <- 0
  while (iter <= max_iter) {
    ## after going over all j, we will update i
    for (i in 1:n_mat) {
      ## initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)
      diff_i <- vector(length = n_mat)

      for (j in 1:n_mat) {
        if (j == i) {
          next
        } else if (j < i) {
          K12 <- t(K_list[[j]][[i]])
        } else {
          K12 <- K_list[[i]][[j]]
        }

        ## compute Y_ij
        X1 <- X_list[[i]]
        X2 <- X_list[[j]]
        Y <- t(X1) %*% K12 %*% X2
        w2 <- w_list[[j]]

        w_i_left[, j] <- Y %*% w2
      }
      w_i_new <- normalize_vec(rowSums(w_i_left))
      diff_i[i] <- max(abs(w_i_new - w_list[[i]]))
      w_list[[i]] <- w_i_new
    }
    if (mean(diff_i) <= tol) {
      print(paste("convergence reached at", iter, "iterations", sep = " "))
      break
    }
    iter <- iter + 1
  }

  if (iter > max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  ## convert w_list into matrix
  for (j in 1:n_mat) {
    w_list[[j]] <- matrix(w_list[[j]], ncol = 1)
  }

  ## return results
  return(w_list)
}


bilinear_w_from_Y_resi <- function(w_list_new, Y_resi,
                                   n_mat, n_features, max_iter, tol) {
  # Iterative refinement
  iter <- 0
  while (iter < max_iter) {
    ## after going over all j, we will update i
    for (i in 1:n_mat) {
      ## initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)
      diff_i <- vector(length = n_mat)

      for (j in 1:n_mat) {
        if (j == i) {
          next
        }
        w2 <- w_list_new[[j]]
        Y <- Y_resi[[i]][[j]] ## make sure all Y_resi exist except for i == j
        w_i_left[, j] <- Y %*% w2
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

#' Run multi version of skrCCA to detect the third to n_th component
#'
#' @importFrom stats setNames
#'
#' @param X_list List of data matrices (cell by PC matrix)
#' @param K_list Kernel matrices between any two pairs of matrices
#' @param w_list A list of weights for each data matrix, pre-computed and to
#'  be regressed out
#' @param cellTypesOfInterest A vector specifying cell types of interest
#' @param nCC number of canonical vectors to compute, need to be greater than
#'  two, default = 2.
#' @param max_iter Maximum number of iterations
#' @param tol tolerance of accuracy
#'
#' @return A list of weights
#' @noRd
optimize_bilinear_multi_n <- function(X_list, K_list, w_list,
                                      cellTypesOfInterest,
                                      nCC = 2,
                                      max_iter = 1000,
                                      tol = 1e-5) {

  if (nCC < 2) {
    stop("nCC must be an integer greater than 1.")
  }
  n_mat <- length(X_list)
  n_features <- ncol(X_list[[1]])

  ## check if w_list has the same length of x_list
  if (n_mat != length(w_list) || n_mat != length(cellTypesOfInterest) ) {
    stop(paste("the input X_list, w_list, cellTypesOfInterest",
               "are of different length!"))
  }

  ## check the dimension of the w_list
  if (length(dim(w_list[[1]])) == 0) {
    stop("the input w_list must be a matrix after Zhen updated the function!")
  }

  ## get all pairs for X, K, W to obtain Y_resi
  ## let Y_resi be the same structure as K
  Y_resi <- setNames(vector(mode = "list", length = n_mat), cellTypesOfInterest)
  for (i in cellTypesOfInterest) {
    Y_resi[[i]] <- setNames(vector(mode = "list", length = n_mat),
                            cellTypesOfInterest)
  }

  ## computations are between all pairs
  pair_cell_types <- combn(cellTypesOfInterest, 2)

  ## Y is between any two cell types, and each time, we update the Y_resi
  ## to regress out all previous canonical vectors

  for (qq in 1:(nCC - 1)) { ## qq: existing w_qq to be regressed out

    ## step 1: obtain or update Y_resi
    for (pp in seq_len(ncol(pair_cell_types))) {
      i <- pair_cell_types[1, pp]
      j <- pair_cell_types[2, pp]
      K12 <- K_list[[i]][[j]]

      ## compute Y_ij
      X1 <- X_list[[i]]
      X2 <- X_list[[j]]

      w1 <- w_list[[i]][, qq, drop = FALSE]
      w2 <- w_list[[j]][, qq, drop = FALSE]

      if (qq == 1) {
        Y1 <- t(X1) %*% K12 %*% X2
      }else {
        Y1 <- Y_resi[[i]][[j]]
      }

      Y_resi[[i]][[j]] <- Y1 - ((t(w1) %*% Y1 %*% w2)[1, 1] * (w1 %*% t(w2)))

      ## fill in the matrix, this will make the optimization step easier
      Y_resi[[j]][[i]] <- t(Y_resi[[i]][[j]])
    }


    ## step 2: initialize w_list_new by SVD
    w_list_new <- rep(list(), length = n_mat)
    w_list_new[[1]] <- svd(t(Y_resi[[1]][[2]]))$v[, 1, drop = FALSE]
    for (i in 2:n_mat) {
      w_list_new[[i]] <- svd(Y_resi[[1]][[i]])$v[, 1, drop = FALSE]
    }

    ## step 3: Iterative refinement
    w_list_qq <- bilinear_w_from_Y_resi(w_list_new = w_list_new,
                                        Y_resi = Y_resi, n_mat = n_mat,
                                        n_features = n_features,
                                        max_iter = max_iter, tol = tol)

    ## step 4: add w_list_qq to columns in w_list
    for (ii in seq_len(n_mat)) {
      w_list[[ii]] <- cbind(w_list[[ii]], w_list_qq[[ii]])
    }
  }

  ## return results
  return(w_list)
}
