

## function from distance to kernel matrix
kernel_from_distance <- function(
    sigma_square, dist_mat, lower_limit = 0.001
){
  kernel_mat <- exp(-1 * dist_mat^2 / (2 * sigma_square)) / (2 * pi * sigma_square)
  kernel_mat[kernel_mat < lower_limit] <- 0
  return(kernel_mat)
}

## function for centering and scaling the matrix 
center_scale_matrix_opt <- function(matrix) {
  # Calculate column means and standard deviations
  col_means <- colMeans(matrix)
  col_sds <- apply(matrix, 2, sd)
  
  # Identify columns that are not full of zeros (to avoid division by zero)
  non_zero_cols <- col_sds >= 1e-5
  
  # Center and scale columns with none-zero sd column
  matrix[, non_zero_cols] <- t((t(matrix[, non_zero_cols]) - col_means[non_zero_cols]) / col_sds[non_zero_cols])
  
  # do not scale, but still need to center columns with zero sd 
  matrix[, !non_zero_cols] <- t((t(matrix[, !non_zero_cols]) - col_means[!non_zero_cols]))
  
  return(matrix)
}


# Function to normalize_vec a vector to have unit norm
normalize_vec <- function(v) {
  return(v / sqrt(sum(v^2)))
}


# Multi-group optimization 
optimize_bilinear_multi <- function(X_list, K_list, max_iter = 1000,
                              tol = 1e-5) {
  
  n_mat <- length(X_list)
  n_features <- ncol(X_list[[1]])
  
  ## computations are between all pairs
  ij_pairs <- combn(n_mat, 2)
  
  ## check all pairwise kernel matrices exist 
  for(t in ncol(ij_pairs)){
    i <- ij_pairs[1,t]
    j <- ij_pairs[2,t]
    
    if(is.null(K_list[[i]][[j]])){
      stop(paste('Kernel matrix between', i, 'and', j, 'not available', 
                 sep = ' ' ))
    }
  }

  
  # Initialize w1 and w2 by their right singular vector
  w_list = rep(list(), length = length(X_list))
  for(i in 1:n_mat){
    w_list[[i]] <- svd(X_list[[i]])$v[,1,drop = FALSE]
  }
  
  # Iterative refinement
  iter <- 0
  while (iter <= max_iter) {

    ## after going over all j, we will update i
    for(i in 1:n_mat){
      ## initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)
      diff_i <- vector(length = n_mat)
      
      for (j in 1:n_mat){
        if (j == i){
          next
        } else if (j < i){
          K12 <- t(K_list[[j]][[i]])
        } else {
          K12 <- K_list[[i]][[j]]
        }
  
        ## compute Y_ij
        X1 <- X_list[[i]]
        X2 <- X_list[[j]]
        Y = t(X1) %*% K12 %*% X2 
        
        w1 <- w_list[[i]]
        w2 <- w_list[[j]]
        
        w_i_left[,j] <- Y %*% w2
      }
      w_i_new <- normalize_vec(rowSums(w_i_left))
      diff_i[i] <- max(abs(w_i_new - w_list[[i]]))
      w_list[[i]] <- w_i_new
    }
    if (mean(diff_i) <= tol){
      print(paste('convergence reached at', iter, 'iterations', sep = ' '))
      break
    }
    iter <- iter + 1
  }
  
  if (iter > max_iter){
    warning("Maximum number of iterations reached without convergence")
  }
  
  ## return results 
  return(w_list)
}


## multi-group optimization to identify the functions beyond the first component
## This is like the PCA where we can compute second PC and so on 
optimize_bilinear_multi_w <- function(X_list, K_list, w_list, max_iter = 1000,
                                    tol = 1e-5) {
  # Function to normalize_vec a vector to have unit norm
  normalize_vec <- function(v) {
    return(v / sqrt(sum(v^2)))
  }
  
  n_mat <- length(X_list)
  n_features <- ncol(X_list[[1]])
  
  ## n_mat has to be greater than 2
  if(n_mat <= 2){
    stop(paste('Make sure number of matrices to be greater than 2', 
    'to run the multi version'))
  }
  
  ## check all pairwise kernel matrices exist 
  for(i in 1:(n_mat-1)){
    for(j in 1:n_mat){
      if(i >= j){
        next
      }else if(is.null(K_list[[i]][[j]])){
        stop(paste('Kernel matrix between', i, 'and', j, 'not available', 
                   sep = ' ' ))
      }
    }
  }
  
  ## check if w_list has the same length of x_list 
  if(length(X_list) != length(w_list)){
    stop('the input X_list and w_list are of different length!')
  }
  
  ## get all pairs for X, K, W to obtain Y_resi
  ## let Y_resi be the same structure as K
  # pair_types = combn(n_mat, 2)
  # Y_resi <- rep(list(), length = ncol(pair_types))
  Y_resi <- rep(list(), length = n_mat)
  for(i in 1:n_mat){
    Y_resi <- rep(list(), length = n_mat)
  }
  
  ## 
  for(i in 1:n_mat){
    for (j in 1:n_mat){
      if (j == i){
        next
      } else if (j < i){
        K12 <- t(K_list[[j]][[i]])
      } else {
        K12 <- K_list[[i]][[j]]
      }
      
      ## compute Y_ij
      X1 <- X_list[[i]]
      X2 <- X_list[[j]]
      Y <- t(X1) %*% K12 %*% X2 
      w1 <- matrix(w_list[[i]], ncol = 1)
      w2 <- matrix(w_list[[j]], ncol = 1)
      
      Y_resi[[i]][[j]] <- Y - ( ( t(w1) %*% Y %*% w2 )[1,1] * (w1 %*% t(w2) ) )
    }
  }
    
  
  
  
  ## initialize w_list_new
  w_list_new <- rep(list(), length = n_mat)
  w_list_new[[1]] <- svd(t(Y_resi[[1]][[2]]))$v[,1,drop = FALSE]
  for(i in 2:n_mat){
    w_list_new[[i]] <- svd(Y_resi[[1]][[i]])$v[,1,drop = FALSE]
  }
  
  # Iterative refinement
  iter <- 0
  while (iter <= max_iter) {
    
    ## after going over all j, we will update i
    for(i in 1:n_mat){
      ## initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)
      diff_i <- vector(length = n_mat)
      
      for (j in 1:n_mat){
        if (j == i){
          next
        } 
        # w1 <- w_list[[i]]
        w2 <- w_list[[j]]
        Y <- Y_resi[[i]][[j]]
        w_i_left[,j] <- Y %*% w2
      }
      w_i_new <- normalize_vec(rowSums(w_i_left))
      diff_i[i] <- max(abs(w_i_new - w_list[[i]]))
      w_list[[i]] <- w_i_new
    }
    if (mean(diff_i) <= tol){
      print(paste('convergence reached at', iter, 'iterations', sep = ' '))
      break
    }
    iter <- iter + 1
  }
  
  if (iter == max_iter){
    warning("Maximum number of iterations reached without convergence")
  }
  
  ## return results 
  return(w_list)
  
}

