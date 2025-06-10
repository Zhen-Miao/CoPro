
## get PC matrices
.getAllPCMats <- function(allPCs, scalePCs) {

  if (length(allPCs) == 0) {
    stop("PCA results do not exist, run computePCA() first.")
  }

  PCmats <- setNames(
    vector("list", length = length(allPCs)),
    names(allPCs)
  )

  ## optionally, scale the PCs before running CCA
  if (scalePCs) {
    for (i in names(allPCs)) {
      pca_A_sd <- allPCs[[i]]$sdev
      PCmats[[i]] <- scale(allPCs[[i]]$x,
                           center = FALSE,
                           scale = pca_A_sd
      )
    }
  } else {
    for (i in names(allPCs)) {
      PCmats[[i]] <- allPCs[[i]]$x
    }
  }
  return(PCmats)
}


#' centering and scaling the matrix
#' @importFrom stats sd
#' @importFrom matrixStats colSds
#' @param matrix Input matrix to be column-centered
#'
#' @return centered and scaled matrix
#' @noRd
center_scale_matrix_opt <- function(input_matrix,
                                    zero_sd_threshold = 1e-3,
                                    nz_propion_threshold = 0.01) {
  # Calculate column standard deviations
  col_means <- colMeans(input_matrix)
  col_sds <- apply(input_matrix, 2, sd)

  # non-zero proportion
  col_nz <- colSums(input_matrix != 0) / nrow(input_matrix)

  # Identify columns that are not full of zeros (to avoid division by zero)

  zero_sd_cols <- which(col_sds < zero_sd_threshold |
                        col_nz < nz_propion_threshold)

  ## do not scale if the sd is too small, or if the proportion of non-zero
  # values is too low
  col_sds_safe <- col_sds
  if (length(zero_sd_cols) > 0) {
    col_sds_safe[zero_sd_cols] <- 1.0
  }

  scaled_matrix <- scale(input_matrix, center = col_means, scale = col_sds_safe)


  return(scaled_matrix)
}



# Function to normalize a vector to have unit norm
normalize_vec <- function(v) {
  norm_v <- sqrt(sum(v^2))
  # Check if norm is near zero to avoid division by zero
  if (norm_v < 1e-7) {
    # If close to zero, return itself, trigger warning
    warning("vector is close to a zero vector, check potential issue")
    return(v)
  } else {
    return(v / norm_v)
  }
}
