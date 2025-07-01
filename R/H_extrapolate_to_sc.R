
#' quantile normalize data
#'
#' @param A The reference cell-by-feature matrix Values in this matrix define
#'   the target distribution for each feature.
#' @param B The cell-by-feature matrix we want to normalize to A
#' @param save_Sparse whether to save as a sparse matrix
#' @param ties_method If there are ties, which method to choose
#' @param verbose whether to print progress updates (default TRUE)
#'
#' @returns The normalized B matrix
#' @export
#'
#' @importFrom stats approx
quantile_normalize <- function(A, B, save_Sparse = FALSE,
                               ties_method = "min", verbose = TRUE) {

  # --- Input Validation ---
  if (!is.matrix(A) && !inherits(A, "Matrix")) {
    stop("Input A must be a matrix or a Matrix object.")
  }
  if (!is.matrix(B) && !inherits(B, "Matrix")) {
    stop("Input B must be a matrix or a Matrix object.")
  }

  # Number of features (columns)
  n_features <- ncol(A)
  if(n_features != ncol(B)){
    stop("The number of features (columns) must match between A and B")
  }

  ## number of rows
  if (nrow(A) < 2 || nrow(B) < 2) {
    stop("Matrix A and B must have at least 2 rows to define quantiles.")
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

  # Obtain the quantiles from A based on B's distribution
  prob_A <- (0:(nrow(A)-1))/(nrow(A)-1)

  ## sort all A values
  A <- apply(A, 2, sort)

  # Pre-calculate ranks for B (vectorized operation)
  rank_B <- apply(B, 2, rank, ties.method = ties_method)
  prob_B <- (rank_B - 1) / (nrow(B) - 1)

  # Use vapply for better performance and type safety
  B_normalized <- vapply(seq_len(n_features), function(i) {
    if (verbose && i %% max(1, n_features %/% 10) == 0) {
      message(sprintf("Processing feature %d/%d", i, n_features))
    }

    # Use approx to interpolate quantiles
    result <- approx(x = prob_A, y = A[, i],
                     xout = prob_B[, i], method = "linear", rule = 2)
    result$y
  }, numeric(nrow(B)))

  if(!is.null(rownames(B))){
    rownames(B_normalized) <- rownames(B)
  }

  if(!is.null(colnames(B))){
    colnames(B_normalized) <- colnames(B)
  }


  # # Loop through each feature
  # for (i in 1:n_features) {
  #   prob_B <- rank_B_1[, i] / n_row_B_1
  #
  #   quantiles_A <- approx(x = prob_A, y = A[, i],
  #                         xout = prob_B, method = "linear")$y
  #
  #   # Apply the calculated quantiles
  #   B_normalized[, i] <- quantiles_A
  # }

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

transfer_scores <- function(mat_A, mat_B, gs_ct, verbose = TRUE){
  B_qn <- quantile_normalize(A = mat_A, B = mat_B, verbose = verbose)

  ## transform B with mean and sd from A
  A_mean  <- apply(mat_A, MARGIN = 2, mean)
  A_sd  <- apply(mat_A, MARGIN = 2, sd)
  B_qn_cs <- scale(B_qn, center = A_mean, scale = A_sd)

  ## get cell weight by marix multiplication
  B_cell_score <- B_qn_cs %*% gs_ct

  return(B_cell_score)
}

#' Get cell score by transferring gene weights from another slide
#' By default, quantile normalization is used to ensure distribution match
#' @param ref_obj Reference object (where the gene weights will be obtained)
#' @param tar_obj Target object (where the cell scores will be obtained)
#' @param sigma_choice Sigma value to be used
#' @param verbose verbose
#'
#' @returns cell scores as a matrix
#' @export
getTransferCellScores <- function(ref_obj, tar_obj, sigma_choice,
                            verbose = TRUE){
  ## check object
  if (!(is(ref_obj, "CoProm") || is(ref_obj, "CoPro"))) {
    stop("ref_obj must be a CoPro or CoProm object")
  }
  if (!(is(tar_obj, "CoProm") || is(tar_obj, "CoPro"))) {
    stop("tar_obj must be a CoPro or CoProm object")
  }

  ## make sure gene names match
  if( length(ref_obj@geneList) != length(tar_obj@geneList) ||
      any(ref_obj@geneList != tar_obj@geneList)){
    warning("gene names mismatch between the two object, use the intersect")
    genes_sel <- intersect(ref_obj@geneList, tar_obj@geneList)
  }else{
    genes_sel <- ref_obj@geneList
  }

  ## get gene weights for reference object
  gs <- ref_obj@geneScores
  sigma_name <- paste0("sigma_",sigma_choice)
  if(length(gs) == 0){
    stop(paste("geneScores slot in reference object does not exist.",
               "Run `computeGeneAndCellScores()` first"))
  }
  if(!(sigma_name %in% names(gs))){
    stop("sigma value does not exist in geneScores slot.")
  }else{
    gs <- gs[[sigma_name]]
  }

  ## check number of cell types
  cts <- ref_obj@cellTypesOfInterest
  cts_tar <- tar_obj@cellTypesOfInterest
  if(!all(cts == cts_tar)){
    stop("cell types mismatch between the two objects")
  }

  B_cs <- rep(list(), length = length(cts))
  names(B_cs) = cts

  for(ct in cts){
    ## get data matrices and run quantile normalization
    mat_A = ref_obj@normalizedDataSub[ref_obj@cellTypesSub == ct, genes_sel]
    mat_B = tar_obj@normalizedDataSub[tar_obj@cellTypesSub == ct, genes_sel]
    if(verbose) cat("running quantile normalization\n")
    B_cs[[ct]] <- transfer_scores(mat_A = mat_A, mat_B = mat_B,
                                  gs_ct = gs[[ct]][genes_sel,],
                                  verbose = verbose)

  }
  names(B_cs) <- NULL

  B_cs_all <- do.call(rbind, B_cs)

  ## re-order the B_cs by cell names
  B_cs_all <- B_cs_all[rownames(tar_obj@normalizedDataSub),]
  return(B_cs_all)

}


