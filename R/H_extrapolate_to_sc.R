
#' quantile normalize data
#'
#' @param A The reference cell by feature matrix
#' @param B The cell by feature matrix we want to normalize to A
#' @param save_Sparse whether to save as a sparse matrix
#' @param ties_method If there are ties, which method to choose
#'
#' @returns The normalized B matrix
#' @export
#'
#' @importFrom stats approx
quantile_normalize <- function(A, B, save_Sparse = FALSE,
                               ties_method = "min") {
  # Number of features (columns)
  n_features <- ncol(A)
  if(n_features != ncol(B)){
    stop("The number of features (columns) must match between A and B")
  }

  ## check ties_method
  if(!(ties_method %in% c("average", "first", "last", "random", "max", "min"))){
    stop(paste("ties_methods must be one of the following:",
               "average, first, last, random, max, min."))
  }

  ## if sparse matrix, convert to dense matrix
  if(inherits(A, "Matrix")){
    A <- as.matrix(A)
  }

  if(inherits(B, "Matrix")){
    B <- as.matrix(B)
  }

  # Initialize the output matrix for B with the same dimensions
  B_normalized <- matrix(nrow = nrow(B), ncol = n_features)
  rownames(B_normalized) <- rownames(B)
  colnames(B_normalized) <- colnames(B)

  # Obtain the quantiles from A based on B's distribution
  prob_A <- (0:(nrow(A)-1))/(nrow(A)-1)

  ## sort all A values
  sorted_A <- apply(A, 2, sort)

  ## rank B values
  rank_B_1 <- apply(B, 2, rank, ties.method = ties_method) - 1
  n_row_B_1 = nrow(B) - 1

  # Loop through each feature
  for (i in 1:n_features) {
    prob_B <- rank_B_1[, i] / n_row_B_1

    quantiles_A <- approx(x = prob_A, y = sorted_A[, i],
                          xout = prob_B, method = "linear")$y

    # Apply the calculated quantiles
    B_normalized[, i] <- quantiles_A
  }

  if(save_Sparse){
    B_normalized = Matrix::Matrix(data = B_normalized, sparse = TRUE)
  }

  return(B_normalized)
}

# # Example usage:
# # Example data matrices
# A <- matrix(c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50), ncol = 2)
# B <- matrix(c(2, 3, 4, 5, 6, 7, 15, 30, 45, 60, 75, 80), ncol = 2)
#
# # Run the quantile normalization
# B_new <- quantile_normalize(A, B)
#
# # Print the original and normalized matrices
# print("Original B:")
# print(B)
# print("Normalized B:")
# print(B_new)
#
# B_sd <- apply(B, MARGIN = 2, sd)
# B_new_sd <- apply(B_new, MARGIN = 2, sd)
# A_sd  <- apply(A, MARGIN = 2, sd)
#
# B_mean <- apply(B, MARGIN = 2, mean)
# B_new_mean <- apply(B_new, MARGIN = 2, mean)
# A_mean  <- apply(A, MARGIN = 2, mean)
