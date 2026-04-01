
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



utils::globalVariables(c("getKernelMatrix", "getSlideList"))

#' Transfer cell scores between matrices using gene weights
#'
#' Given reference (`mat_A`) and target (`mat_B`) cell-by-gene matrices with the
#' same set and order of genes (columns), this function optionally quantile
#' normalizes the target to match the reference per feature, standardizes the
#' target using the reference mean and standard deviation, filters small-magnitude
#' gene weights, and computes cell scores via matrix multiplication
#' \eqn{B_{standardized} \\times W}.
#'
#' - If `use_quantile_normalization = TRUE`, `mat_B` is quantile-normalized to
#'   the distribution of `mat_A` for each gene.
#' - Columns of the normalized `mat_B` are then centered and scaled using the
#'   column means and standard deviations computed from `mat_A` (with a small
#'   safeguard for near-zero standard deviations).
#' - Gene weights in `gs_ct` are thresholded by absolute value per column; values
#'   with `abs(weight) < gs_weight_threshold` are set to 0.
#' - Final cell scores are obtained as the matrix product of the standardized
#'   `mat_B` with the filtered weight matrix `gs_ct`.
#'
#' @param mat_A Numeric matrix of shape cells-by-genes (reference). Column order
#'   must correspond to genes in `gs_ct`.
#' @param mat_B Numeric matrix of shape cells-by-genes (target). Must have the
#'   same genes (columns) and order as `mat_A`.
#' @param gs_ct Numeric matrix of gene weights of shape genes-by-K, where K is
#'   the number of signatures or cell types. Column names, if present, are used
#'   in verbose messages.
#' @param use_quantile_normalization Logical; whether to quantile-normalize
#'   `mat_B` to `mat_A` per feature before standardization (default `TRUE`).
#' @param gs_weight_threshold Numeric; absolute threshold used to zero out small
#'   gene weights in `gs_ct` (default `0.005`).
#' @param verbose Logical; whether to print progress messages (default `TRUE`).
#'
#' @returns A numeric matrix of cell scores with shape (nrow(`mat_B`) x ncol(`gs_ct`)).
#'
#' @details This function assumes that `mat_A`, `mat_B`, and `gs_ct` are aligned
#' on the same set and ordering of genes (columns of `mat_A`/`mat_B`, rows of
#' `gs_ct`). Callers are responsible for subsetting/reordering if needed.
#'
#' @seealso [getTransferCellScores]
#'
#' @keywords internal
transfer_scores <- function(mat_A, mat_B, gs_ct,
                            use_quantile_normalization = TRUE,
                            gs_weight_threshold = 0.005,
                            verbose = TRUE){
  ## quantile normalize if needed
  if(use_quantile_normalization){
    B_qn <- quantile_normalize(A = mat_A, B = mat_B, verbose = verbose)
  }else{
    B_qn <- mat_B
  }

  ## transform B with mean and sd from A
  A_mean  <- apply(mat_A, MARGIN = 2, mean)
  A_sd  <- apply(mat_A, MARGIN = 2, sd)
  # guard against zero (or near-zero) sd to avoid division by zero
  A_sd_safe <- A_sd
  A_sd_safe[is.na(A_sd_safe) | A_sd_safe < 1e-8] <- 1.0
  B_qn_cs <- scale(B_qn, center = A_mean, scale = A_sd_safe)

  # fitler gs_ct by column. Look at each column, and get the elements that are greater than the threshold
  # replace smaller elements with 0
  gs_ct_filt <- gs_ct
  for(i in seq_len(ncol(gs_ct_filt))){
    ind <- abs(gs_ct_filt[, i, drop = TRUE]) < gs_weight_threshold
    if(verbose){
      sum_ind <- sum(!ind)
      cat(paste("retaining", sum_ind, "genes for", colnames(gs_ct_filt)[i],
       "with threshold", gs_weight_threshold, "\n"))
    }
    gs_ct_filt[,i][ind] <- 0
  }

  ## get cell weight by marix multiplication
  B_cell_score <- B_qn_cs %*% gs_ct_filt

  return(B_cell_score)
}

#' Select genes by intersection of gene names
#' @param ref_geneList A vector of gene names
#' @param tar_geneList A vector of gene names
#' @return A vector of gene names
#' @keywords internal
#' @noRd
.trans_gene_sel <- function(ref_geneList, tar_geneList){
  # if gene names mismatch, use the intersect
  if( length(ref_geneList) != length(tar_geneList) ||
      any(ref_geneList != tar_geneList)){
    warning("gene names mismatch between the two object, use the intersect")
    genes_sel <- intersect(ref_geneList, tar_geneList)
  }else{
    genes_sel <- ref_geneList
  }

  # select genes by weight 

  return(genes_sel)
}

#' Get gene score names
#' @param sigma_choice Sigma value to be used
#' @param cts Cell types
#' @param slide Slide name
#' @return A vector of gene score names
#' @keywords internal
#' @noRd
.get_gs_names <- function(sigma_choice, cts, slide = NULL){
  gs_names <- vector("character", length = length(cts))
  for(i in seq_along(cts)){
    gs_names[i] <- .createGeneScoresName(sigma_choice, cts[i], slide = slide)
  }
  return(gs_names)
}

#' Get cell score by transferring gene weights from another slide
#' By default, quantile normalization is used to ensure distribution match
#' @param ref_obj Reference object (where the gene weights will be obtained)
#' @param tar_obj Target object (where the cell scores will be obtained)
#' @param sigma_choice Sigma value to be used
#' @param use_quantile_normalization Logical; apply quantile normalization of target to reference distribution (default TRUE)
#' @param agg_cell_type Logical; if TRUE, returns a single matrix aggregated across cell types (default FALSE)
#' @param gs_weight_threshold Numeric; absolute gene-weight threshold for filtering prior to transfer (default 0.005)
#' @param sigma_choice_tar Numeric; sigma value for target object. If NULL (default), uses sigma_choice. Not recommended for general use.
#' @param verbose verbose
#'
#' @returns cell scores as a matrix
#' @export
getTransferCellScores <- function(ref_obj, tar_obj, sigma_choice,
                            use_quantile_normalization = TRUE,
                            agg_cell_type = FALSE,
                            gs_weight_threshold = 0.005,
                            sigma_choice_tar = NULL, verbose = TRUE){
  ## check object
  if (!(is(ref_obj, "CoProMulti") || is(ref_obj, "CoProSingle"))) {
    stop("ref_obj must be a CoProSingle or CoProMulti object")
  }
  if (!(is(tar_obj, "CoProMulti") || is(tar_obj, "CoProSingle"))) {
    stop("tar_obj must be a CoProSingle or CoProMulti object")
  }
  
  ## handle sigma_choice_tar parameter
  if (is.null(sigma_choice_tar)) {
    sigma_choice_tar <- sigma_choice
  } else {
    warning("Using different sigma values for reference and target objects is not recommended and is intended for development use only.")
  }

  ## make sure gene names match
  genes_sel <- .trans_gene_sel(ref_geneList = ref_obj@geneList,
   tar_geneList = tar_obj@geneList)
  if(length(genes_sel) == 0){
    stop("No overlapping genes between reference and target objects after matching")
  }


  ## check number of cell types
  cts <- ref_obj@cellTypesOfInterest
  cts_tar <- tar_obj@cellTypesOfInterest
  if(!all(cts == cts_tar)){
    stop("cell types mismatch between the two objects")
  }

  ## get gene weights for reference object
  gs <- ref_obj@geneScores
  if(length(gs) == 0){
    stop(paste("geneScores slot in reference object does not exist.",
               "Run `computeGeneAndCellScores()` first"))
  }
  # subset geneScores by sigma_choice and cell types
  gs_names <- .get_gs_names(sigma_choice, cts, slide = NULL)
  missing_gs <- setdiff(gs_names, names(gs))
  if(length(missing_gs) > 0){
    stop(paste0("Missing gene score entries for: ", paste(missing_gs, collapse = ", "),
                ". Ensure computeGeneAndCellScores() has been run for these."))
  }
  gs <- gs[gs_names]

  # initialize the list to store cell scores
  B_cs <- rep(list(), length = length(cts))
  names(B_cs) = cts

  for(ct in cts){
    # get data matrices
    mat_A = ref_obj@normalizedDataSub[ref_obj@cellTypesSub == ct, genes_sel]
    mat_B = tar_obj@normalizedDataSub[tar_obj@cellTypesSub == ct, genes_sel]
    if(verbose) cat("transferring gene scores for cell type", ct, "\n")
    # slide = NULL for both single and multi slide objects
    gs_name_ct <- .createGeneScoresName(sigma_choice, ct, slide = NULL)
    B_cs[[ct]] <- transfer_scores(mat_A = mat_A, mat_B = mat_B,
                                  gs_ct = gs[[gs_name_ct]][genes_sel,],
                                  use_quantile_normalization = use_quantile_normalization,
                                  gs_weight_threshold = gs_weight_threshold,
                                  verbose = verbose)

  }

  if(agg_cell_type){
    names(B_cs) <- NULL
    B_cs_all <- do.call(rbind, B_cs)
    ## re-order the B_cs by cell names
    B_cs_all <- B_cs_all[rownames(tar_obj@normalizedDataSub),]
  }else{
    B_cs_all <- B_cs
  }

  return(B_cs_all)

}




#' Compute Normalized Correlation from Transferred Cell Scores
#'
#' Given a target object and a list of transferred cell scores (per cell type),
#' compute the normalized correlation for each pair of cell types using the
#' provided kernel matrices at a selected sigma value.
#'
#' This function mirrors the correlation calculation used in
#' `computeNormalizedCorrelation()` but operates on precomputed cell scores
#' (e.g., obtained from `getTransferCellScores(agg_cell_type = FALSE)`).
#'
#' The normalized correlation is computed as:
#'   numerator = t(A_w1) %*% K %*% B_w2
#'   denominator = sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * ||K||_2
#' where A_w1 and B_w2 are the cell score vectors (for the same CC index)
#' for cell types A and B respectively, K is the kernel matrix between the two
#' cell types at the chosen sigma, and ||K||_2 is the spectral norm of K.
#'
#' @param tar_obj A `CoProSingle` or `CoProMulti` object containing kernel matrices
#'   and metadata needed for alignment.
#' @param transfer_cell_scores A named list of matrices, with one entry per cell type
#'   (names must be the cell type names). Each matrix should be cells-by-CCs, where
#'   rows are cell IDs and columns are `CC_1`, `CC_2`, ..., as returned by
#'   `getTransferCellScores(agg_cell_type = FALSE)`.
#' @param sigma_choice Numeric scalar specifying the sigma value of the kernel to use.
#' @param tol Numeric tolerance passed to the truncated SVD for spectral norm
#'   estimation (default 1e-4).
#' @param calculationMode For `CoProMulti` objects only, either "perSlide" or
#'   "aggregate". Ignored for `CoProSingle`. Default "perSlide" if `tar_obj`
#'   is multi-slide.
#' @param sigma_choice_tar Numeric; sigma value for target object kernel matrices. If NULL (default), uses sigma_choice. Not recommended for general use.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list with one element named `paste0("sigma_", sigma_choice)`, whose
#'   value is a data.frame of results. For single-slide objects, the data.frame has
#'   columns `sigmaValue`, `cellType1`, `cellType2`, `CC_index`, `normalizedCorrelation`.
#'   For multi-slide objects in `perSlide` mode, the data.frame additionally includes
#'   `slideID`. For `aggregate` mode, the correlation column is named
#'   `aggregateCorrelation`.
#'
#' @examples
#' # Assuming `tar_obj` is prepared and `trans_scores` was computed with
#' # getTransferCellScores(..., agg_cell_type = FALSE)
#' # res <- getTransferNormCorr(tar_obj, trans_scores, sigma_choice = 2.0)
#'
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @export
getTransferNormCorr <- function(tar_obj,
                                transfer_cell_scores,
                                sigma_choice,
                                tol = 1e-4,
                                calculationMode = NULL,
                                sigma_choice_tar = NULL,
                                verbose = TRUE) {
  # --- Input validation ---
  if (!(is(tar_obj, "CoProMulti") || is(tar_obj, "CoProSingle"))) {
    stop("tar_obj must be a CoProSingle or CoProMulti object")
  }
  if (!is.list(transfer_cell_scores) || length(transfer_cell_scores) == 0) {
    stop("transfer_cell_scores must be a non-empty named list of matrices")
  }
  if (is.null(names(transfer_cell_scores)) || any(names(transfer_cell_scores) == "")) {
    stop("transfer_cell_scores must be a named list with cell type names")
  }
  if (!is.numeric(sigma_choice) || length(sigma_choice) != 1 || is.na(sigma_choice) || sigma_choice <= 0) {
    stop("sigma_choice must be a positive numeric scalar")
  }
  
  ## handle sigma_choice_tar parameter
  if (is.null(sigma_choice_tar)) {
    sigma_choice_tar <- sigma_choice
  } else {
    warning("Using different sigma values for reference and target objects is not recommended and is intended for development use only.")
  }

  # Determine cell types to use (from the provided scores)
  cts <- names(transfer_cell_scores)
  if (length(cts) == 0) stop("No cell types found in transfer_cell_scores")

  # Infer number of CCs from the first matrix
  first_mat <- transfer_cell_scores[[cts[1]]]
  if (!is.matrix(first_mat)) stop("Each entry of transfer_cell_scores must be a numeric matrix")
  nCC <- ncol(first_mat)
  if (is.null(nCC) || nCC < 1) stop("transferred cell score matrices must have >=1 columns (CCs)")

  # Validate structure of matrices
  for (ct in cts) {
    mat <- transfer_cell_scores[[ct]]
    if (!is.matrix(mat)) stop(paste0("transfer_cell_scores[[", ct, "]] is not a matrix"))
    if (ncol(mat) != nCC) stop("All transferred score matrices must have the same number of CC columns")
    if (is.null(rownames(mat))) stop(paste0("Row names (cell IDs) are required for cell type ", ct))
  }

  # Decide calculation mode for multi-slide objects
  is_multi <- is(tar_obj, "CoProMulti")
  if (is_multi) {
    if (is.null(calculationMode)) calculationMode <- "perSlide"
    if (!calculationMode %in% c("perSlide", "aggregate")) {
      stop("calculationMode must be either 'perSlide' or 'aggregate' for CoProMulti")
    }
  }

  # Build pairwise cell type combinations
  if (length(cts) == 1) {
    pair_cell_types <- matrix(c(cts, cts), nrow = 2, ncol = 1)
  } else {
    pair_cell_types <- combn(cts, 2)
  }

  sigma_name <- paste0("sigma_", sigma_choice)
  sigma_name_tar <- paste0("sigma_", sigma_choice_tar)

  # Helper: align a score vector to matrix dimension names if available
  .align_scores <- function(scores_mat, target_names) {
    if (!is.null(target_names) && !is.null(rownames(scores_mat))) {
      idx <- match(target_names, rownames(scores_mat))
      if (any(is.na(idx))) {
        stop("Mismatch between kernel dimension names and score row names.")
      }
      return(scores_mat[idx, , drop = FALSE])
    }
    return(scores_mat)
  }

  if (!is_multi) {
    # Precompute spectral norms using existing helper
    norm_K12 <- .computeSpecNorm(
      object = tar_obj, tol = tol, cts = cts,
      scalePCs = if (length(tar_obj@scalePCs) == 0) FALSE else tar_obj@scalePCs,
      sigmaValues = sigma_choice_tar,
      nCC = if (length(tar_obj@nCC) == 0) nCC else tar_obj@nCC,
      pair_cell_types = pair_cell_types
    )

    # --- Single slide object ---
    df <- data.frame(
      sigmaValue = numeric(),
      cellType1 = character(), cellType2 = character(),
      CC_index = integer(), normalizedCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (pp in seq_len(ncol(pair_cell_types))) {
      ct_i <- pair_cell_types[1, pp]
      ct_j <- pair_cell_types[2, pp]

      # Kernel and its spectral norm
      K_ij <- getKernelMatrix(tar_obj, sigma = sigma_choice_tar,
                              cellType1 = ct_i, cellType2 = ct_j,
                              verbose = FALSE)
      norm_K_ij <- norm_K12[[sigma_name_tar]][[ct_i]][[ct_j]]
      if (is.na(norm_K_ij) || norm_K_ij < 1e-9) next

      # Align transferred scores if kernel has dimnames
      A_scores <- transfer_cell_scores[[ct_i]]
      B_scores <- transfer_cell_scores[[ct_j]]
      A_scores <- .align_scores(A_scores, rownames(K_ij))
      B_scores <- .align_scores(B_scores, colnames(K_ij))

      for (cc in seq_len(nCC)) {
        A_w1 <- A_scores[, cc, drop = FALSE]
        B_w2 <- B_scores[, cc, drop = FALSE]
        num <- as.numeric(t(A_w1) %*% K_ij %*% B_w2)
        denom <- sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K_ij
        nc_val <- ifelse(is.na(denom) || abs(denom) < 1e-9, 0, num / denom)
        df <- rbind(df, data.frame(
          sigmaValue = sigma_choice,
          cellType1 = ct_i, cellType2 = ct_j,
          CC_index = cc, normalizedCorrelation = as.numeric(nc_val),
          stringsAsFactors = FALSE
        ))
      }
    }

    return(setNames(list(df), sigma_name))
  }

  # --- Multi-slide object ---
  slides <- getSlideList(tar_obj)
  # Precompute spectral norms using existing multi-slide helper
  norm_K_all <- .computeSpecNormMulti(
    object = tar_obj, tol = tol, cts = cts, slides = slides,
    sigmas_run = sigma_name_tar, nCC = if (length(tar_obj@nCC) == 0) nCC else tar_obj@nCC,
    pair_cell_types = pair_cell_types
  )

  if (calculationMode == "perSlide") {
    df_all <- data.frame(
      sigmaValue = numeric(),
      slideID = character(),
      cellType1 = character(), cellType2 = character(),
      CC_index = integer(), normalizedCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (sID in slides) {
      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]

        K_ij <- tryCatch({
          getKernelMatrix(tar_obj, sigma = sigma_choice_tar,
                          cellType1 = ct_i, cellType2 = ct_j,
                          slide = sID, verbose = FALSE)
        }, error = function(e) NULL)
        if (is.null(K_ij)) next
        norm_K_ij <- norm_K_all[[sigma_name_tar]][[sID]][[ct_i]][[ct_j]]
        if (is.na(norm_K_ij) || norm_K_ij < 1e-9) next

        # Extract slide-specific rows by cell IDs
        cells_i <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_i)
        cells_j <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_j)
        if (length(cells_i) == 0 || length(cells_j) == 0) next

        A_scores <- transfer_cell_scores[[ct_i]][cells_i, , drop = FALSE]
        B_scores <- transfer_cell_scores[[ct_j]][cells_j, , drop = FALSE]
        A_scores <- .align_scores(A_scores, rownames(K_ij))
        B_scores <- .align_scores(B_scores, colnames(K_ij))

        for (cc in seq_len(nCC)) {
          A_w1 <- A_scores[, cc, drop = FALSE]
          B_w2 <- B_scores[, cc, drop = FALSE]
          num <- as.numeric(t(A_w1) %*% K_ij %*% B_w2)
          denom <- sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K_ij
          nc_val <- ifelse(is.na(denom) || abs(denom) < 1e-9, 0, num / denom)
          df_all <- rbind(df_all, data.frame(
            sigmaValue = sigma_choice,
            slideID = sID,
            cellType1 = ct_i, cellType2 = ct_j,
            CC_index = cc, normalizedCorrelation = as.numeric(nc_val),
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    return(setNames(list(df_all), sigma_name))

  } else { # aggregate
    df_agg <- data.frame(
      sigmaValue = numeric(),
      cellType1 = character(), cellType2 = character(),
      CC_index = integer(), aggregateCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (pp in seq_len(ncol(pair_cell_types))) {
      ct_i <- pair_cell_types[1, pp]
      ct_j <- pair_cell_types[2, pp]

      for (cc in seq_len(nCC)) {
        total_numerator <- 0
        total_norm_sum_i <- 0
        total_norm_sum_j <- 0
        total_K_norm <- 0
        valid_slides <- 0

        for (sID in slides) {
          K_ij <- tryCatch({
            getKernelMatrix(tar_obj, sigma = sigma_choice_tar,
                            cellType1 = ct_i, cellType2 = ct_j,
                            slide = sID, verbose = FALSE)
          }, error = function(e) NULL)
          if (is.null(K_ij)) next
          norm_K_ij <- norm_K_all[[sigma_name_tar]][[sID]][[ct_i]][[ct_j]]
          if (is.na(norm_K_ij) || norm_K_ij < 1e-9) next

          cells_i <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_i)
          cells_j <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_j)
          if (length(cells_i) == 0 || length(cells_j) == 0) next

          A_scores <- transfer_cell_scores[[ct_i]][cells_i, , drop = FALSE]
          B_scores <- transfer_cell_scores[[ct_j]][cells_j, , drop = FALSE]
          A_scores <- .align_scores(A_scores, rownames(K_ij))
          B_scores <- .align_scores(B_scores, colnames(K_ij))

          A_w1 <- A_scores[, cc, drop = FALSE]
          B_w2 <- B_scores[, cc, drop = FALSE]
          total_numerator <- total_numerator + as.numeric(t(A_w1) %*% K_ij %*% B_w2)
          total_norm_sum_i <- total_norm_sum_i + sum(A_w1^2)
          total_norm_sum_j <- total_norm_sum_j + sum(B_w2^2)
          total_K_norm <- total_K_norm + norm_K_ij
          valid_slides <- valid_slides + 1
        }

        if (valid_slides > 0) {
          avg_K_norm <- total_K_norm / valid_slides
          denom <- sqrt(total_norm_sum_i) * sqrt(total_norm_sum_j) * avg_K_norm
          agg_val <- ifelse(is.na(denom) || abs(denom) < 1e-9, 0, total_numerator / denom)
          df_agg <- rbind(df_agg, data.frame(
            sigmaValue = sigma_choice,
            cellType1 = ct_i, cellType2 = ct_j,
            CC_index = cc, aggregateCorrelation = as.numeric(agg_val),
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    return(setNames(list(df_agg), sigma_name))
  }
}

#' Compute Bidirectional Correlation from Transferred Cell Scores
#'
#' Given a target object and a list of transferred cell scores (per cell type),
#' compute the bidirectional correlation for each pair of cell types using the
#' provided kernel matrices at a selected sigma value.
#'
#' This function mirrors the correlation calculation used in
#' `computeBidirCorrelation()` but operates on precomputed cell scores
#' (e.g., obtained from `getTransferCellScores(agg_cell_type = FALSE)`).
#'
#' The bidirectional correlation is computed as the mean of two correlations:
#' cor(t(K) %*% A_w1, B_w2) and cor(A_w1, K %*% B_w2), where A_w1 and B_w2
#' are the cell score vectors (for the same CC index) for cell types A and B
#' respectively, and K is the kernel matrix between the two cell types at the
#' chosen sigma.
#'
#' @param tar_obj A `CoProSingle` or `CoProMulti` object containing kernel matrices
#'   and metadata needed for alignment.
#' @param transfer_cell_scores A named list of matrices, with one entry per cell type
#'   (names must be the cell type names). Each matrix should be cells-by-CCs, where
#'   rows are cell IDs and columns are `CC_1`, `CC_2`, ..., as returned by
#'   `getTransferCellScores(agg_cell_type = FALSE)`.
#' @param sigma_choice Numeric scalar specifying the sigma value of the kernel to use.
#' @param calculationMode For `CoProMulti` objects only, either "perSlide" or
#'   "aggregate". Ignored for `CoProSingle`. Default "perSlide" if `tar_obj`
#'   is multi-slide.
#' @param normalize_K Character; method for normalizing the kernel matrix, one of 
#'   "row_or_col", "sinkhorn_knopp", or "none". Default "row_or_col".
#' @param filter_kernel Logical; whether to filter the kernel matrix. Default TRUE.
#' @param K_row_sum_cutoff Numeric; cutoff for row sums when normalizing kernel matrix.
#'   Default 5e-3.
#' @param K_col_sum_cutoff Numeric; cutoff for column sums when normalizing kernel matrix.
#'   Default 5e-3.
#' @param sigma_choice_tar Numeric; sigma value for target object kernel matrices. If NULL (default), uses sigma_choice. Not recommended for general use.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list with one element named `paste0("sigma_", sigma_choice)`, whose
#'   value is a data.frame of results. For single-slide objects, the data.frame has
#'   columns `sigmaValue`, `cellType1`, `cellType2`, `CC_index`, `bidirCorrelation`.
#'   For multi-slide objects in `perSlide` mode, the data.frame additionally includes
#'   `slideID`. For `aggregate` mode, the correlation column is named
#'   `aggregateCorrelation`.
#'
#' @examples
#' # Assuming `tar_obj` is prepared and `trans_scores` was computed with
#' # getTransferCellScores(..., agg_cell_type = FALSE)
#' # res <- getTransferBidirCorr(tar_obj, trans_scores, sigma_choice = 2.0)
#'
#' @importFrom utils combn
#' @export
getTransferBidirCorr <- function(tar_obj,
                                 transfer_cell_scores,
                                 sigma_choice,
                                 calculationMode = NULL,
                                 normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
                                 filter_kernel = TRUE,
                                 K_row_sum_cutoff = 5e-3,
                                 K_col_sum_cutoff = 5e-3,
                                 sigma_choice_tar = NULL,
                                 verbose = TRUE) {
  normalize_K <- match.arg(normalize_K)
  # --- Input validation ---
  if (!(is(tar_obj, "CoProMulti") || is(tar_obj, "CoProSingle"))) {
    stop("tar_obj must be a CoProSingle or CoProMulti object")
  }
  if (!is.list(transfer_cell_scores) || length(transfer_cell_scores) == 0) {
    stop("transfer_cell_scores must be a non-empty named list of matrices")
  }
  if (is.null(names(transfer_cell_scores)) || any(names(transfer_cell_scores) == "")) {
    stop("transfer_cell_scores must be a named list with cell type names")
  }
  if (!is.numeric(sigma_choice) || length(sigma_choice) != 1 || is.na(sigma_choice) || sigma_choice <= 0) {
    stop("sigma_choice must be a positive numeric scalar")
  }
  
  ## handle sigma_choice_tar parameter
  if (is.null(sigma_choice_tar)) {
    sigma_choice_tar <- sigma_choice
  } else {
    warning("Using different sigma values for reference and target objects is not recommended and is intended for development use only.")
  }

  # Determine cell types to use (from the provided scores)
  cts <- names(transfer_cell_scores)
  if (length(cts) == 0) stop("No cell types found in transfer_cell_scores")

  # Infer number of CCs from the first matrix
  first_mat <- transfer_cell_scores[[cts[1]]]
  if (!is.matrix(first_mat)) stop("Each entry of transfer_cell_scores must be a numeric matrix")
  nCC <- ncol(first_mat)
  if (is.null(nCC) || nCC < 1) stop("transferred cell score matrices must have >=1 columns (CCs)")

  # Validate structure of matrices
  for (ct in cts) {
    mat <- transfer_cell_scores[[ct]]
    if (!is.matrix(mat)) stop(paste0("transfer_cell_scores[[", ct, "]] is not a matrix"))
    if (ncol(mat) != nCC) stop("All transferred score matrices must have the same number of CC columns")
    if (is.null(rownames(mat))) stop(paste0("Row names (cell IDs) are required for cell type ", ct))
  }

  # Decide calculation mode for multi-slide objects
  is_multi <- is(tar_obj, "CoProMulti")
  if (is_multi) {
    if (is.null(calculationMode)) calculationMode <- "perSlide"
    if (!calculationMode %in% c("perSlide", "aggregate")) {
      stop("calculationMode must be either 'perSlide' or 'aggregate' for CoProMulti")
    }
  }

  # Build pairwise cell type combinations
  if (length(cts) == 1) {
    pair_cell_types <- matrix(c(cts, cts), nrow = 2, ncol = 1)
  } else {
    pair_cell_types <- combn(cts, 2)
  }

  sigma_name <- paste0("sigma_", sigma_choice)
  sigma_name_tar <- paste0("sigma_", sigma_choice_tar)

  # Helper: align a score vector to matrix dimension names if available
  .align_scores <- function(scores_mat, target_names) {
    if (!is.null(target_names) && !is.null(rownames(scores_mat))) {
      idx <- match(target_names, rownames(scores_mat))
      if (any(is.na(idx))) {
        stop("Mismatch between kernel dimension names and score row names.")
      }
      return(scores_mat[idx, , drop = FALSE])
    }
    return(scores_mat)
  }

  if (!is_multi) {
    # --- Single slide object ---
    df <- data.frame(
      sigmaValue = numeric(),
      cellType1 = character(), cellType2 = character(),
      CC_index = integer(), bidirCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (pp in seq_len(ncol(pair_cell_types))) {
      ct_i <- pair_cell_types[1, pp]
      ct_j <- pair_cell_types[2, pp]

      # Kernel
      K_ij <- getKernelMatrix(tar_obj, sigma = sigma_choice_tar,
                              cellType1 = ct_i, cellType2 = ct_j,
                              verbose = FALSE)

      # Align transferred scores if kernel has dimnames
      A_scores <- transfer_cell_scores[[ct_i]]
      B_scores <- transfer_cell_scores[[ct_j]]
      A_scores <- .align_scores(A_scores, rownames(K_ij))
      B_scores <- .align_scores(B_scores, colnames(K_ij))

      for (cc in seq_len(nCC)) {
        A_w1 <- A_scores[, cc, drop = FALSE]
        B_w2 <- B_scores[, cc, drop = FALSE]
        
        bidir_val <- .computeSpatialCrossCorrelation(A_w1, B_w2, K_ij, normalize_K = normalize_K, filter_kernel = filter_kernel, K_row_sum_cutoff = K_row_sum_cutoff, K_col_sum_cutoff = K_col_sum_cutoff)
        
        df <- rbind(df, data.frame(
          sigmaValue = sigma_choice,
          cellType1 = ct_i, cellType2 = ct_j,
          CC_index = cc, bidirCorrelation = as.numeric(bidir_val),
          stringsAsFactors = FALSE
        ))
      }
    }

    return(setNames(list(df), sigma_name))
  }

  # --- Multi-slide object ---
  slides <- getSlideList(tar_obj)

  if (calculationMode == "perSlide") {
    df_all <- data.frame(
      sigmaValue = numeric(),
      slideID = character(),
      cellType1 = character(), cellType2 = character(),
      CC_index = integer(), bidirCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (sID in slides) {
      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]

        K_ij <- tryCatch({
          getKernelMatrix(tar_obj, sigma = sigma_choice_tar,
                          cellType1 = ct_i, cellType2 = ct_j,
                          slide = sID, verbose = FALSE)
        }, error = function(e) NULL)
        if (is.null(K_ij)) next

        # Extract slide-specific rows by cell IDs
        cells_i <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_i)
        cells_j <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_j)
        if (length(cells_i) == 0 || length(cells_j) == 0) next

        A_scores <- transfer_cell_scores[[ct_i]][cells_i, , drop = FALSE]
        B_scores <- transfer_cell_scores[[ct_j]][cells_j, , drop = FALSE]
        A_scores <- .align_scores(A_scores, rownames(K_ij))
        B_scores <- .align_scores(B_scores, colnames(K_ij))

        for (cc in seq_len(nCC)) {
          A_w1 <- A_scores[, cc, drop = FALSE]
          B_w2 <- B_scores[, cc, drop = FALSE]
          
          bidir_val <- .computeSpatialCrossCorrelation(A_w1, B_w2, K_ij, normalize_K = normalize_K, filter_kernel = filter_kernel, K_row_sum_cutoff = K_row_sum_cutoff, K_col_sum_cutoff = K_col_sum_cutoff)
          
          df_all <- rbind(df_all, data.frame(
            sigmaValue = sigma_choice,
            slideID = sID,
            cellType1 = ct_i, cellType2 = ct_j,
            CC_index = cc, bidirCorrelation = as.numeric(bidir_val),
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    return(setNames(list(df_all), sigma_name))

  } else { # aggregate
    df_agg <- data.frame(
      sigmaValue = numeric(),
      cellType1 = character(), cellType2 = character(),
      CC_index = integer(), aggregateCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (pp in seq_len(ncol(pair_cell_types))) {
      ct_i <- pair_cell_types[1, pp]
      ct_j <- pair_cell_types[2, pp]

      for (cc in seq_len(nCC)) {
        sum_corr <- 0
        valid_slides <- 0

        for (sID in slides) {
          K_ij <- tryCatch({
            getKernelMatrix(tar_obj, sigma = sigma_choice_tar,
                            cellType1 = ct_i, cellType2 = ct_j,
                            slide = sID, verbose = FALSE)
          }, error = function(e) NULL)
          if (is.null(K_ij)) next

          cells_i <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_i)
          cells_j <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct_j)
          if (length(cells_i) == 0 || length(cells_j) == 0) next

          A_scores <- transfer_cell_scores[[ct_i]][cells_i, , drop = FALSE]
          B_scores <- transfer_cell_scores[[ct_j]][cells_j, , drop = FALSE]
          A_scores <- .align_scores(A_scores, rownames(K_ij))
          B_scores <- .align_scores(B_scores, colnames(K_ij))

          A_w1 <- A_scores[, cc, drop = FALSE]
          B_w2 <- B_scores[, cc, drop = FALSE]
          
          bidir_val <- .computeSpatialCrossCorrelation(A_w1, B_w2, K_ij, normalize_K = normalize_K, filter_kernel = filter_kernel, K_row_sum_cutoff = K_row_sum_cutoff, K_col_sum_cutoff = K_col_sum_cutoff)

          if (is.finite(bidir_val) && !is.na(bidir_val)) {
            sum_corr <- sum_corr + as.numeric(bidir_val)
            valid_slides <- valid_slides + 1
          }
        }

        if (valid_slides > 0) {
          df_agg <- rbind(df_agg, data.frame(
            sigmaValue = sigma_choice,
            cellType1 = ct_i, cellType2 = ct_j,
            CC_index = cc, aggregateCorrelation = sum_corr / valid_slides,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    return(setNames(list(df_agg), sigma_name))
  }
}
