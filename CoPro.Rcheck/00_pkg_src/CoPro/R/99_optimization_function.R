
# #' Helper function to safely get K_ij or t(K_ji) from K_list
# #' Assumes K_list is named list K_list[[ct1]][[ct2]]
# #' Assumes if K_list[[ct_i]][[ct_j]] is NULL, K_list[[ct_j]][[ct_i]] exists.
# #' @param K_list The named kernel list structure
# #' @param ct_i Name of the first cell type
# #' @param ct_j Name of the second cell type
# #' @return The K_ij matrix
# #' @noRd
# get_kernel_matrix <- function(K_list, ct_i, ct_j) {
#   K12 <- NULL
#   # Try accessing K_ij directly
#   if (length(K_list[[ct_i]][[ct_j]]) != 0) {
#     K12 <- K_list[[ct_i]][[ct_j]]
#   } else {
#     # If K_ij is missing, try getting K_ji and transposing it
#     K21 <- K_list[[ct_j]][[ct_i]]
#     K12 <- t(K21)
#   }

#   # Check if we successfully obtained a matrix
#   if (length(K12) == 0) {
#     stop(paste("Could not find kernel matrix for pair:", ct_i, ct_j,
#                "or its transpose in K_list."))
#   }
#   return(K12)
# }





# #' SkrCCA optimization function for multiple groups (Single Slide) - First Component
# #' Assumes X_list and K_list are named lists using consistent cell type names.
# #' Handles upper/lower triangle K_list. Uses w_i_left structure.
# #'
# #' @param X_list Named list of data matrices (cell by PC matrix)
# #' @param K_list Named list of kernel matrices (potentially upper/lower triangle only)
# #' @param max_iter Maximum number of iterations
# #' @param tol tolerance of accuracy
# #'
# #' @return Named list `w_list` containing the first weight vector component.
# #' @noRd
# optimize_bilinear <- function(X_list, K_list, max_iter = 1000,
#                                     tol = 1e-5) {

#   cell_types <- names(X_list)
#   if (is.null(cell_types)) stop("Input X_list must be a named list.")
#   n_mat <- length(cell_types)
#   if (n_mat < 2) stop("Need at least two cell types (entries in X_list).")
#   n_features <- ncol(X_list[[cell_types[1]]])

#   # Basic check on K_list names
#   if (is.null(names(K_list)) || !all(cell_types %in% names(K_list))) {
#     stop("K_list must be a named list containing entries for all cell types in X_list.")
#   }

#   # Initialize w_list (named list) using SVD
#   w_list <- setNames(vector("list", length = n_mat), cell_types)
#   for (ct in cell_types) {
#     if(is.null(X_list[[ct]]) || !is.matrix(X_list[[ct]]) ||
#        ncol(X_list[[ct]]) != n_features) {
#       stop(paste("Invalid or inconsistent matrix in X_list for cell type:", ct))
#     }
#     svd_result <- tryCatch(irlba(X_list[[ct]], nv = 1, right_only = TRUE),
#                            error = function(e) {
#       stop(paste("SVD failed for cell type:", ct, "Error:", e$message))})
#     if(ncol(svd_result$v) < 1) stop(paste("SVD resulted in zero singular vectors for cell type:", ct))
#     w_list[[ct]] <- svd_result$v[, 1, drop = FALSE] ## orthonormal
#   }

#   # Iterative refinement
#   iter <- 0
#   while (iter <= max_iter) {
#     w_list_old <- w_list
#     # max_diff_iter <- 0 # Track max change in this iteration

#     for (ct_i in cell_types) {
#       # Initialize the w_i_left matrix to store contributions column-wise
#       # Columns correspond to ct_j (excluding ct_i)
#       w_i_update_vec <- matrix(0, nrow = n_features, ncol = 1)

#       for (ct_j in cell_types) {
#         if (ct_i == ct_j) {
#           next
#         }

#         # Get K_ij or t(K_ji) robustly
#         K12 <- get_kernel_matrix(K_list, ct_i, ct_j)

#         # Get matrices and current weights
#         X1 <- X_list[[ct_i]]
#         X2 <- X_list[[ct_j]]
#         w2 <- w_list[[ct_j]] # Use current w_j from this iteration's w_list

#         # Efficient calculation: t(X1) %*% (K12 %*% (X2 %*% w2))
#         v2 <- X2 %*% w2
#         kv2 <- K12 %*% v2
#         update_contribution <- crossprod(X1, kv2) # This is Y_ij %*% w_j

#         # update weight vector
#         w_i_update_vec <- w_i_update_vec + update_contribution

#       } # end loop ct_j

#       # Normalize the summed update vector
#       w_i_new <- normalize_vec(w_i_update_vec)

#       # Update w_list for the next inner iteration
#       w_list[[ct_i]] <- w_i_new

#     } # end loop ct_i

#     # Check convergence after updating all w_i
#     current_max_diff <- 0
#     for (ct in cell_types) {
#       diff_val <- max(abs(w_list[[ct]] - w_list_old[[ct]]))
#       if(is.nan(diff_val)) diff_val <- 0
#       if(diff_val > current_max_diff) {
#         current_max_diff <- diff_val
#       }
#     }

#     if (current_max_diff <= tol) {
#       print(paste("Convergence reached at", iter, "iterations (Max diff =", sprintf("%.3e", current_max_diff), ")"))
#       break
#     }
#     iter <- iter + 1
#   } # end while

#   if (iter > max_iter) {
#     warning("Maximum number of iterations reached without convergence")
#   }

#   # Ensure final format is list of single-column matrices
#   for (ct in cell_types) {
#     if (!is.matrix(w_list[[ct]]) || ncol(w_list[[ct]]) != 1) {
#       w_list[[ct]] <- matrix(w_list[[ct]], ncol = 1)
#     }
#   }
#   return(w_list)
# }


# #' Helper function for optimizing multiple components using precomputed Y_resi
# #' Assumes w_list_new and Y_resi are named lists using consistent cell type names.
# #'
# #' @param w_list_new Initial named list of weight vectors
# #' @param Y_resi Named list of residual matrices (Y_resi[[ct1]][[ct2]])
# #' @param n_features Number of features
# #' @param max_iter Maximum iterations
# #' @param tol Tolerance
# #'
# #' @return Updated named list of weight vectors
# #' @noRd
# bilinear_w_from_Y_resi <- function(w_list_new, Y_resi,
#                                    n_features, max_iter, tol) {

#   cell_types <- names(w_list_new)
#   if (length(cell_types) == 0) stop("Input w_list_new must be a named list.")
#   n_mat <- length(cell_types)

#   if (length(names(Y_resi)) == 0 || !all(cell_types %in% names(Y_resi))) {
#     stop("Y_resi must be a named list containing entries for all cell types in w_list_new.")
#   }

#   iter <- 0
#   while (iter < max_iter) {
#     w_list_old <- w_list_new
#     # max_diff_iter <- 0

#     for (ct_i in cell_types) {
#       # Use w_i_left structure for accumulation
#       w_i_left <- matrix(0, nrow = n_features, ncol = n_mat)
#       colnames(w_i_left) <- cell_types

#       for (ct_j in cell_types) {
#         if (ct_i == ct_j) {
#           next
#         }
#         w2 <- w_list_new[[ct_j]]
#         # Y_resi should be fully populated (Y_ij and Y_ji = t(Y_ij))
#         # by the calling function
#         Y <- Y_resi[[ct_i]][[ct_j]]
#         if(is.null(Y)) stop(paste("Missing Y_resi matrix for pair:", ct_i, ct_j))

#         update_contribution <- Y %*% w2

#         j_idx <- match(ct_j, cell_types)
#         w_i_left[, j_idx] <- update_contribution
#       } # end ct_j loop

#       # i_idx <- match(ct_i, cell_types)
#       # w_i_update_vec <- rowSums(w_i_left[, -i_idx, drop = FALSE])
#       w_i_update_vec <- rowSums(w_i_left[, , drop = FALSE])
#       w_i_new <- normalize_vec(w_i_update_vec)
#       w_list_new[[ct_i]] <- w_i_new
#     } # end ct_i loop

#     current_max_diff <- 0
#     for(ct in cell_types){
#       diff_val <- max(abs(w_list_new[[ct]] - w_list_old[[ct]]))
#       if(is.nan(diff_val)) diff_val <- 0
#       if(diff_val > current_max_diff){
#         current_max_diff <- diff_val
#       }
#     }

#     if (current_max_diff <= tol) {
#       print(paste("Convergence reached at", iter, "iterations (Max diff =", sprintf("%.3e", current_max_diff), ")"))
#       break
#     }
#     iter <- iter + 1
#   } # end while

#   if (iter == max_iter) {
#     warning("Maximum number of iterations reached without convergence")
#   }
#   return(w_list_new)
# }


# #' Run multi version of skrCCA to detect subsequent components (Single Slide)
# #' Assumes X_list, K_list, w_list are named lists containing only cellTypesOfInterest.
# #' Handles upper/lower triangle K_list. No internal subsetting.
# #'
# #' @param X_list Named list of data matrices (subsetted)
# #' @param K_list Named list of kernel matrices (subsetted, potentially upper/lower triangle)
# #' @param w_list A named list of weights (subsetted, matrices with previous components as columns)
# #' @param cellTypesOfInterest A vector specifying cell type names present in the input lists
# #' @param nCC Total number of canonical vectors desired (must be >= 2)
# #' @param max_iter Maximum number of iterations for helper function
# #' @param tol Tolerance of accuracy for helper function
# #'
# #' @return A named list of weights (matrices with components 1 to nCC as columns)
# #' @noRd

# optimize_bilinear_n <- function(X_list, K_list, w_list,
#                                       cellTypesOfInterest,
#                                       nCC = 2,
#                                       max_iter = 1000,
#                                       tol = 1e-5) {

#   # Validate inputs based on assumption they are already subsetted
#   cts <- cellTypesOfInterest
#   n_mat <- length(cts)
#   if (n_mat < 2) stop("Need at least two cell types.")
#   if (length(X_list) != n_mat || length(K_list) != n_mat || length(w_list) != n_mat ||
#       !all(cts %in% names(X_list)) || !all(cts %in% names(K_list)) || !all(cts %in% names(w_list))) {
#     stop("Input lists length or names do not match cellTypesOfInterest.")
#   }
#   n_features <- ncol(X_list[[cts[1]]])

#   # Check input w_list dimensions
#   k_start <- ncol(w_list[[cts[1]]]) ## we expect this value to be one
#   ## based on the structure of our function
#   if (k_start < 1) stop("Input w_list must contain at least the first component.")
#   for (ct in cts) {
#     if (!is.matrix(w_list[[ct]]) || ncol(w_list[[ct]]) != k_start) {
#       stop(paste("w_list for", ct, "is not a matrix or has inconsistent component count."))
#     }
#   }
#   if (nCC <= k_start) {
#     stop(paste("nCC (", nCC, ") must be greater than the number of components already in w_list (", k_start, ")"))
#   }


#   # Initialize Y_resi structure using names
#   Y_resi <- setNames(vector(mode = "list", length = n_mat), cts)
#   for (ct_i in cts) {
#     Y_resi[[ct_i]] <- setNames(vector(mode = "list", length = n_mat), cts)
#   }

#   # Generate pairs using names
#   pair_cell_types <- combn(cts, 2)

#   # Loop to compute components k_start + 1 up to nCC
#   # qq iterates through the component index *to be used for deflation*
#   for (qq in k_start:(nCC - 1)) {

#     ## qq == 1 means we have the weight vector of 1st component,
#     ## so the residual is calculated based on the 1st weight

#     # Step 1: Update Y_resi using component qq for deflation
#     for (pp in seq_len(ncol(pair_cell_types))) {
#       i <- pair_cell_types[1, pp] # Name ct1
#       j <- pair_cell_types[2, pp] # Name ct2

#       # Get X matrices
#       X1 <- X_list[[i]]
#       X2 <- X_list[[j]]

#       # Get the qq-th weight vectors
#       ## initially, we should have the first component ready (qq == 1)
#       w1 <- w_list[[i]][, qq, drop = FALSE]
#       w2 <- w_list[[j]][, qq, drop = FALSE]

#       # Get the cross-product matrix (use original if first deflation, else use deflated)
#       if (qq == 1) {
#         # Calculate Y from original data for the first deflation pass
#         # Use robust kernel access
#         K12 <- get_kernel_matrix(K_list, i, j)
#         Y1 <- t(X1) %*% K12 %*% X2
#       } else {
#         # Use the Y deflated by previous component
#         ## say that qq == 2, then Y_resi is already calculated for 1st component
#         ## Then, we start from there to further regress out the residual
#         Y1 <- Y_resi[[i]][[j]]
#         if (is.null(Y1)) stop(paste("Y_resi missing for pair", i, j, "before deflating with component", qq))
#       }

#       # Apply deflation formula
#       deflation_scalar <- (t(w1) %*% Y1 %*% w2)[1, 1]
#       deflation_term <- deflation_scalar * (w1 %*% t(w2))
#       Y_resi[[i]][[j]] <- Y1 - deflation_term

#       # Update the transpose entry
#       Y_resi[[j]][[i]] <- t(Y_resi[[i]][[j]])
#     } # end deflation pair loop pp


#     # Step 2: Initialize w_list_new for component qq+1 (using names)
#     # Using SVD initialization
#     w_list_new <- setNames(vector("list", length = n_mat), cts)
#     w_list_new[[cts[1]]] <- irlba(t(Y_resi[[cts[1]]][[cts[2]]]),
#                                   nv = 1, right_only = TRUE)$v[, 1, drop = FALSE]
#     for (i in cts[2:n_mat]) {
#       w_list_new[[i]] <- irlba(Y_resi[[cts[1]]][[i]],
#                                nv = 1, right_only = TRUE)$v[, 1, drop = FALSE]
#     }

#     # Step 3: Iterative refinement using the helper function
#     w_list_qq_plus_1 <- bilinear_w_from_Y_resi(
#       w_list_new = w_list_new,
#       Y_resi = Y_resi,
#       n_features = n_features,
#       max_iter = max_iter, tol = tol)

#     # Step 4: Add the new component (qq+1) to w_list using names
#     for (ct in cts) {
#       # Directly modify the input w_list
#       w_list[[ct]] <- cbind(w_list[[ct]], w_list_qq_plus_1[[ct]])
#     }
#   } # end component loop qq

#   # Return the updated list (modified in place)
#   return(w_list)
# }


# #' SkrCCA optimization function for multiple groups
# #'
# #' @importFrom irlba irlba
# #' @param X_list List of data matrices (cell by PC matrix)
# #' @param K_list Kernel matrices between any two pairs of matrices
# #' @param max_iter Maximum number of iterations
# #' @param tol tolerance of accuracy
# #'
# #' @return `w_list`, a list of weight vectors for each cell type
# #' @noRd
# optimize_bilinear_One <- function(X_list, K_list, max_iter = 1000,
#                                     tol = 1e-5) {
#   n_mat <- 1
#   n_features <- ncol(X_list[[1]])

#   # Initialize w1 by its right singular vector
#   w_list <- list()
#   w_list[[1]] <- irlba::irlba(X_list[[1]], nv = 1,
#                               right_only = TRUE)$v[, 1, drop = FALSE]

#   # Iterative refinement
#   iter <- 0
#   diff_i <- vector(length = n_mat)
#   while (iter <= max_iter) {
#     ## after going over all j, we will update i
#       ## initialize the w_i_left matrix
#       w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)

#       K11 <- K_list[[1]][[1]]
#       ## compute Y_ij
#       X1 <- X_list[[1]]
#       Y <- t(X1) %*% K11 %*% X1
#       w2 <- w_list[[1]]

#       w_i_left[, 1] <- Y %*% w2

#       w_i_new <- normalize_vec(rowSums(w_i_left))
#       diff_i[1] <- max(abs(w_i_new - w_list[[1]]))
#       w_list[[1]] <- w_i_new

#     if (max(diff_i) <= tol) {
#       print(paste("convergence reached at", iter, "iterations"))
#       break
#     }
#     iter <- iter + 1
#   }

#   if (iter > max_iter) {
#     warning("Maximum number of iterations reached without convergence")
#   }

#   ## convert w_list into matrix
#   w_list[[1]] <- matrix(w_list[[1]], ncol = 1)

#   ## return results
#   return(w_list)
# }




# bilinear_w_from_Y_resi_One <- function(w_list_new, Y_resi,
#                                     n_features, max_iter, tol) {
#   # Iterative refinement
#   iter <- 0
#   while (iter < max_iter) {
#     ## after going over all j, we will update i

#       ## initialize the w_i_left matrix
#       w_i_left <- matrix(data = 0, ncol = 1, nrow = n_features)
#       diff_i <- vector(length = 1)

#       w2 <- w_list_new[[1]]
#       Y <- Y_resi[[1]][[1]] ## make sure all Y_resi exist except for i == j
#       w_i_left[, 1] <- Y %*% w2

#       w_i_new <- normalize_vec(rowSums(w_i_left))
#       diff_i[1] <- max(abs(w_i_new - w_list_new[[1]]))
#       w_list_new[[1]] <- w_i_new

#     if (mean(diff_i) <= tol) {
#       print(paste("convergence reached at", iter, "iterations", sep = " "))
#       break
#     }
#     iter <- iter + 1
#   }

#   if (iter == max_iter) {
#     warning("Maximum number of iterations reached without convergence")
#   }

#   return(w_list_new)
# }



# #' Run multi version of skrCCA to detect the third to n_th component
# #'
# #' @importFrom stats setNames
# #'
# #' @param X_list List of data matrices (cell by PC matrix)
# #' @param K_list Kernel matrices between any two pairs of matrices
# #' @param w_list A list of weights for each data matrix, pre-computed and to
# #'  be regressed out
# #' @param cellTypesOfInterest A vector specifying cell types of interest
# #' @param nCC number of canonical vectors to compute, need to be greater than
# #'  two, default = 2.
# #' @param max_iter Maximum number of iterations
# #' @param tol tolerance of accuracy
# #'
# #' @return A list of weights
# #' @noRd
# optimize_bilinear_n_One <- function(X_list, K_list, w_list,
#                                       cellTypesOfInterest,
#                                       nCC = 2,
#                                       max_iter = 1000,
#                                       tol = 1e-5) {

#   n_mat <- 1
#   n_features <- ncol(X_list[[1]])
#   cts <- cellTypesOfInterest

#   ## check if w_list has the same length of x_list
#   if (n_mat != length(w_list) || n_mat != length(cellTypesOfInterest) ) {
#     stop(paste("the input X_list, w_list, cellTypesOfInterest",
#                "are of different length!"))
#   }

#   ## check the dimension of the w_list
#   if (length(dim(w_list[[1]])) == 0) {
#     stop("the input w_list must be a matrix!")
#   }

#   ## get all pairs for X, K, W to obtain Y_resi
#   ## let Y_resi be the same structure as K
#   Y_resi <- setNames(vector(mode = "list", length = n_mat), cts)
#   Y_resi[[cts]] <- setNames(vector(mode = "list", length = n_mat), cts)

#   ## Y is between any two cell types, and each time, we update the Y_resi
#   ## to regress out all previous canonical vectors

#   for (qq in 1:(nCC - 1)) { ## qq: existing w_qq to be regressed out

#     ## step 1: obtain or update Y_resi
#       K11 <- K_list[[cts]][[cts]]

#       ## compute Y_ij
#       X1 <- X_list[[cts]]
#       w1 <- w_list[[1]][, qq, drop = FALSE]

#       if (qq == 1) {
#         Y1 <- t(X1) %*% K11 %*% X1
#       }else {
#         Y1 <- Y_resi[[cts]][[cts]]
#       }

#       Y_resi[[cts]][[cts]] <- Y1 - ((t(w1) %*% Y1 %*% w1)[1, 1] * (w1 %*% t(w1)))

#     ## step 2: initialize w_list_new by SVD
#     w_list_new <- rep(list(), length = n_mat)
#     w_list_new[[1]] <- irlba::irlba(t(Y_resi[[cts]][[cts]]),
#                                     nv = 1, right_only = TRUE)$v[, 1, drop = FALSE]

#     ## step 3: Iterative refinement
#     w_list_qq <- bilinear_w_from_Y_resi_One(w_list_new = w_list_new,
#                                         Y_resi = Y_resi,
#                                         n_features = n_features,
#                                         max_iter = max_iter, tol = tol)

#     ## step 4: add w_list_qq to columns in w_list
#     w_list[[1]] <- cbind(w_list[[1]], w_list_qq[[1]])

#   }

#   ## return results
#   return(w_list)
# }
