#' Compute Self-Bidirectional Correlation for Multiple Cell Types
#'
#' This file contains functions to compute within-cell-type (self) bidirectional
#' correlation using the self-kernel matrices computed by computeSelfKernel().
#' This extends the CoPro package's capability to analyze spatial autocorrelation
#' patterns within individual cell types.
#'
#' @name self_bidirectional_correlation
#' @keywords internal
NULL

#' Compute Self-Bidirectional Correlation from Transferred Cell Scores
#'
#' Given transferred cell scores for each cell type, compute the within-cell-type bidirectional
#' correlation using self-kernel matrices. This function computes spatial autocorrelation
#' within each cell type using transferred scores, complementing the cross-type analysis 
#' provided by `getTransferBidirCorr()`.
#'
#' The self-bidirectional correlation is computed as the mean of two correlations:
#' cor(t(K) %*% A_w, A_w) and cor(A_w, K %*% A_w), where A_w is the transferred cell score
#' vector for a cell type and K is the self-kernel matrix for that cell type.
#'
#' @param tar_obj A `CoProSingle` or `CoProMulti` object containing self-kernel matrices
#'   and metadata needed for alignment. Must have self-kernel matrices computed using
#'   `computeSelfKernel()`.
#' @param transfer_cell_scores A named list of matrices, with one entry per cell type
#'   (names must be the cell type names). Each matrix should be cells-by-CCs, where
#'   rows are cell IDs and columns are `CC_1`, `CC_2`, ..., as returned by
#'   `getTransferCellScores(agg_cell_type = FALSE)`.
#' @param sigma_choice Numeric scalar specifying the sigma value of the self-kernel to use.
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
#' @param sigma_choice_tar Numeric; sigma value for target object self-kernel matrices. 
#'   If NULL (default), uses sigma_choice. Not recommended for general use.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list with one element named `paste0("sigma_", sigma_choice)`, whose
#'   value is a data.frame of results. For single-slide objects, the data.frame has
#'   columns `sigmaValue`, `cellType`, `CC_index`, `selfBidirCorrelation`.
#'   For multi-slide objects in `perSlide` mode, the data.frame additionally includes
#'   `slideID`. For `aggregate` mode, the correlation column is named
#'   `aggregateSelfCorrelation`.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a CoPro object with multiple cell types
#' # First compute standard workflow
#' object <- computeDistance(object)
#' object <- computeKernelMatrix(object, sigmaValues = c(0.01, 0.05, 0.1))
#' 
#' # Add self-distances and self-kernels
#' object <- computeSelfDistance(object)
#' object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))
#' 
#' # Compute transferred cell scores
#' trans_scores <- getTransferCellScores(ref_obj, tar_obj, sigma_choice = 0.05, 
#'                                      agg_cell_type = FALSE)
#' 
#' # Compute self-bidirectional correlation from transferred scores
#' self_bidir <- getTransferSelfBidirCorr(tar_obj, trans_scores, sigma_choice = 0.05)
#' }
#'
#' @importFrom utils combn
#' @export
getTransferSelfBidirCorr <- function(tar_obj,
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
  
  # Handle sigma_choice_tar parameter
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

  # Check that self-kernel matrices exist
  missing_self_kernels <- character(0)
  is_multi <- is(tar_obj, "CoProMulti")
  
  if (is_multi) {
    slides <- getSlideList(tar_obj)
    for (ct in cts) {
      for (sID in slides) {
        self_kernel <- tryCatch({
          getSelfKernelMatrix(tar_obj, sigma = sigma_choice_tar, cellType = ct, 
                             slide = sID, verbose = FALSE)
        }, error = function(e) NULL)
        if (is.null(self_kernel)) {
          missing_self_kernels <- c(missing_self_kernels, paste(ct, sID, sep = "@"))
        }
      }
    }
  } else {
    for (ct in cts) {
      self_kernel <- tryCatch({
        getSelfKernelMatrix(tar_obj, sigma = sigma_choice_tar, cellType = ct, verbose = FALSE)
      }, error = function(e) NULL)
      if (is.null(self_kernel)) {
        missing_self_kernels <- c(missing_self_kernels, ct)
      }
    }
  }
  
  if (length(missing_self_kernels) > 0) {
    stop(paste("Self-kernel matrices missing for:", 
               paste(missing_self_kernels, collapse = ", "), 
               ". Run computeSelfKernel() first."))
  }

  # Decide calculation mode for multi-slide objects
  if (is_multi) {
    if (is.null(calculationMode)) calculationMode <- "perSlide"
    if (!calculationMode %in% c("perSlide", "aggregate")) {
      stop("calculationMode must be either 'perSlide' or 'aggregate' for CoProMulti")
    }
  }

  sigma_name <- paste0("sigma_", sigma_choice)

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
      cellType = character(),
      CC_index = integer(), 
      selfBidirCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (ct in cts) {
      if (verbose) cat("Computing self-bidirectional correlation for cell type:", ct, "\n")
      
      # Get self-kernel matrix
      K_self <- getSelfKernelMatrix(tar_obj, sigma = sigma_choice_tar, 
                                   cellType = ct, verbose = FALSE)
      
      # Align transferred scores if kernel has dimnames
      A_scores <- transfer_cell_scores[[ct]]
      A_scores <- .align_scores(A_scores, rownames(K_self))

      for (cc in seq_len(nCC)) {
        A_w <- A_scores[, cc, drop = FALSE]
        
        # Compute self-bidirectional correlation
        # This is the spatial autocorrelation within the cell type
        self_bidir_val <- .computeSpatialSelfCorrelation(A_w, K_self, 
                                                        normalize_K = normalize_K, 
                                                        filter_kernel = filter_kernel, 
                                                        K_row_sum_cutoff = K_row_sum_cutoff, 
                                                        K_col_sum_cutoff = K_col_sum_cutoff)
        
        df <- rbind(df, data.frame(
          sigmaValue = sigma_choice,
          cellType = ct,
          CC_index = cc, 
          selfBidirCorrelation = as.numeric(self_bidir_val),
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
      cellType = character(),
      CC_index = integer(), 
      selfBidirCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (sID in slides) {
      if (verbose) cat("Processing slide:", sID, "\n")
      
      for (ct in cts) {
        if (verbose) cat("  Cell type:", ct, "\n")
        
        # Get self-kernel matrix for this slide and cell type
        K_self <- tryCatch({
          getSelfKernelMatrix(tar_obj, sigma = sigma_choice_tar, cellType = ct, 
                             slide = sID, verbose = FALSE)
        }, error = function(e) NULL)
        if (is.null(K_self)) next

        # Extract slide-specific rows by cell IDs
        cells_ct <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct)
        if (length(cells_ct) == 0) next

        A_scores <- transfer_cell_scores[[ct]][cells_ct, , drop = FALSE]
        A_scores <- .align_scores(A_scores, rownames(K_self))

        for (cc in seq_len(nCC)) {
          A_w <- A_scores[, cc, drop = FALSE]
          
          self_bidir_val <- .computeSpatialSelfCorrelation(A_w, K_self, 
                                                          normalize_K = normalize_K, 
                                                          filter_kernel = filter_kernel,
                                                          K_row_sum_cutoff = K_row_sum_cutoff, 
                                                          K_col_sum_cutoff = K_col_sum_cutoff)
          
          df_all <- rbind(df_all, data.frame(
            sigmaValue = sigma_choice,
            slideID = sID,
            cellType = ct,
            CC_index = cc, 
            selfBidirCorrelation = as.numeric(self_bidir_val),
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    return(setNames(list(df_all), sigma_name))

  } else { # aggregate
    df_agg <- data.frame(
      sigmaValue = numeric(),
      cellType = character(),
      CC_index = integer(), 
      aggregateSelfCorrelation = numeric(),
      stringsAsFactors = FALSE
    )

    for (ct in cts) {
      if (verbose) cat("Computing aggregate self-correlation for cell type:", ct, "\n")
      
      for (cc in seq_len(nCC)) {
        sum_corr <- 0
        valid_slides <- 0

        for (sID in slides) {
          K_self <- tryCatch({
            getSelfKernelMatrix(tar_obj, sigma = sigma_choice_tar, cellType = ct, 
                               slide = sID, verbose = FALSE)
          }, error = function(e) NULL)
          if (is.null(K_self)) next

          cells_ct <- .getSlideCellTypeIDs(tar_obj, slide = sID, cellType = ct)
          if (length(cells_ct) == 0) next

          A_scores <- transfer_cell_scores[[ct]][cells_ct, , drop = FALSE]
          A_scores <- .align_scores(A_scores, rownames(K_self))

          A_w <- A_scores[, cc, drop = FALSE]
          
          self_bidir_val <- .computeSpatialSelfCorrelation(A_w, K_self, 
                                                          normalize_K = normalize_K, 
                                                          filter_kernel = filter_kernel,
                                                          K_row_sum_cutoff = K_row_sum_cutoff, 
                                                          K_col_sum_cutoff = K_col_sum_cutoff)

          if (is.finite(self_bidir_val) && !is.na(self_bidir_val)) {
            sum_corr <- sum_corr + as.numeric(self_bidir_val)
            valid_slides <- valid_slides + 1
          }
        }

        if (valid_slides > 0) {
          df_agg <- rbind(df_agg, data.frame(
            sigmaValue = sigma_choice,
            cellType = ct,
            CC_index = cc, 
            aggregateSelfCorrelation = sum_corr / valid_slides,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    return(setNames(list(df_agg), sigma_name))
  }
}

#' Compute Spatial Self-Correlation (Within Cell Type)
#'
#' This is the core function that computes bidirectional correlation within a single
#' cell type using its self-kernel matrix. It computes the mean of two correlations:
#' cor(t(K) %*% A_w, A_w) and cor(A_w, K %*% A_w), where A_w is the cell score
#' vector and K is the self-kernel matrix.
#'
#' @param A_w Cell score vector/matrix (cells x 1)
#' @param K_self Self-kernel matrix (cells x cells) for the same cell type
#' @param normalize_K Character; method for normalizing the kernel matrix
#' @param filter_kernel Logical; whether to filter the kernel matrix
#' @param K_row_sum_cutoff Numeric; cutoff for row sums when filtering
#' @param K_col_sum_cutoff Numeric; cutoff for column sums when filtering
#'
#' @return A single numeric value: the mean of the two self-correlations
#' @keywords internal
#' @noRd
.computeSpatialSelfCorrelation <- function(A_w, K_self, 
                                          normalize_K = c("row_or_col", "sinkhorn_knopp", "none"), 
                                          filter_kernel = TRUE,
                                          K_row_sum_cutoff = 5e-3, 
                                          K_col_sum_cutoff = 5e-3) {

  normalize_K <- match.arg(normalize_K)
  
  # Ensure vector is a column matrix
  if (!is.matrix(A_w)) A_w <- matrix(A_w, ncol = 1)
  
  # For self-correlation, we need specialized filtering to maintain square matrix
  if (filter_kernel) {
    filtered_result <- .filter_self_kernel_matrix(K_self, A_w, K_row_sum_cutoff, K_col_sum_cutoff)
    K_self <- filtered_result$K
    A_w <- filtered_result$A_w
  }

  if (normalize_K == "row_or_col") {
    # Optimized row/column normalization
    K_row_sum <- rowSums(K_self)
    K_col_sum <- colSums(K_self)

    # More efficient normalization using vectorized operations
    K_row_norm <- K_self / K_row_sum  # Broadcasting division
    K_col_norm <- sweep(K_self, 2, K_col_sum, "/")  # More efficient than t(t(K) / K_col_sum)

    # Use crossprod for more efficient matrix multiplication
    KA <- crossprod(K_row_norm, A_w)  # Equivalent to t(K_row_norm) %*% A_w
    KB <- K_col_norm %*% A_w

  } else if (normalize_K == "none") {
    # Compute kernel-weighted vectors - use crossprod for efficiency
    KA <- crossprod(K_self, A_w)  # Equivalent to t(K_self) %*% A_w but faster
    KB <- K_self %*% A_w
    
  } else if (normalize_K == "sinkhorn_knopp") {
    K_self <- sinkhorn_knopp(K_self)  # Use the optimized version
    KA <- crossprod(K_self, A_w)  # More efficient transpose multiplication
    KB <- K_self %*% A_w
  }

  # Optimized correlation computation
  # For self-correlation: cor(t(K) %*% A_w, A_w) and cor(A_w, K %*% A_w)
  if (length(KA) == 1 && length(A_w) == 1) {
    # Single values case
    cor1 <- 1.0
  } else {
    # Use more efficient correlation calculation
    cor1 <- stats::cor(as.vector(KA), as.vector(A_w))
  }
  
  if (length(A_w) == 1 && length(KB) == 1) {
    # Single values case  
    cor2 <- 1.0
  } else {
    cor2 <- stats::cor(as.vector(A_w), as.vector(KB))
  }
  
  return((cor1 + cor2) * 0.5)  # Multiplication is faster than division
}

#' Filter Self-Kernel Matrix (Square Matrix)
#'
#' For self-correlation, we need to filter both rows and columns consistently
#' to maintain a square matrix. This function ensures that the same cells are
#' kept for both dimensions.
#'
#' @param K_self Square kernel matrix (cells x cells)
#' @param A_w Cell score vector/matrix (cells x 1)
#' @param K_row_sum_cutoff Cutoff for row sums
#' @param K_col_sum_cutoff Cutoff for column sums
#'
#' @return List with filtered K and A_w
#' @keywords internal
#' @noRd
.filter_self_kernel_matrix <- function(K_self, A_w, K_row_sum_cutoff = 5e-3, K_col_sum_cutoff = 5e-3) {
  
  if (!is.matrix(A_w)) A_w <- matrix(A_w, ncol = 1)
  
  # For self-kernel, we need to keep the same cells for both rows and columns
  # Use the more restrictive of row and column filtering
  K_row_sum <- rowSums(K_self)
  K_col_sum <- colSums(K_self)
  
  # Keep cells that pass both row and column thresholds
  cells_keep <- (K_row_sum > K_row_sum_cutoff) & (K_col_sum > K_col_sum_cutoff)
  
  if (any(!cells_keep)) {
    n_original <- length(cells_keep)
    n_kept <- sum(cells_keep)
    
    # Filter both dimensions consistently
    K_self <- K_self[cells_keep, cells_keep, drop = FALSE]
    A_w <- A_w[cells_keep, , drop = FALSE]
    
    cat("After row filtering: kept", n_kept, "of", n_original, "rows\n")
    cat("After column filtering: kept", n_kept, "of", n_original, "columns\n")
  }
  
  return(list(K = K_self, A_w = A_w))
}

#' Compute Self-Bidirectional Correlation using skrCCA Results
#'
#' This function computes self-bidirectional correlation directly from a CoPro object
#' that has skrCCA results, using the object's own cell scores rather than transferred
#' scores. This is useful for computing spatial autocorrelation patterns within each
#' cell type using the object's native skrCCA results.
#'
#' @param object A `CoProSingle` or `CoProMulti` object with skrCCA results and
#'   self-kernel matrices computed using `computeSelfKernel()`.
#' @param sigma_choice Numeric scalar specifying the sigma value to use.
#' @param calculationMode For `CoProMulti` objects only, either "perSlide" or
#'   "aggregate". Default "perSlide".
#' @param normalize_K Character; method for normalizing the kernel matrix, one of 
#'   "row_or_col", "sinkhorn_knopp", or "none". Default "row_or_col".
#' @param filter_kernel Logical; whether to filter the kernel matrix. Default TRUE.
#' @param K_row_sum_cutoff Numeric; cutoff for row sums when normalizing kernel matrix.
#'   Default 5e-3.
#' @param K_col_sum_cutoff Numeric; cutoff for column sums when normalizing kernel matrix.
#'   Default 5e-3.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list with one element named `paste0("sigma_", sigma_choice)`, whose
#'   value is a data.frame of results with the same structure as `getTransferSelfBidirCorr()`.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a CoPro object with skrCCA results and self-kernels
#' object <- runSkrCCA(object)
#' object <- computeSelfDistance(object)
#' object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))
#' 
#' # Compute self-bidirectional correlation using native skrCCA results
#' self_bidir <- computeSelfBidirCorr(object, sigma_choice = 0.05)
#' }
#'
#' @export
computeSelfBidirCorr <- function(object,
                                sigma_choice,
                                calculationMode = "perSlide",
                                normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
                                filter_kernel = TRUE,
                                K_row_sum_cutoff = 5e-3,
                                K_col_sum_cutoff = 5e-3,
                                verbose = TRUE) {
  
  normalize_K <- match.arg(normalize_K)
  
  # Input validation
  if (!(is(object, "CoProMulti") || is(object, "CoProSingle"))) {
    stop("object must be a CoProSingle or CoProMulti object")
  }
  
  # Check for skrCCA results
  if (length(object@skrCCAOut) == 0) {
    stop("No skrCCA results found. Run runSkrCCA() first.")
  }
  
  # Check for PCA results
  if (length(object@pcaGlobal) == 0) {
    stop("No PCA results found. Run computePCA() first.")
  }
  
  # Get cell types and other parameters
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) {
    stop("No cell types of interest specified")
  }
  
  sigma_name <- paste0("sigma_", sigma_choice)
  if (!sigma_name %in% names(object@skrCCAOut)) {
    stop(paste("Sigma value", sigma_choice, "not found in skrCCA results"))
  }
  
  nCC <- object@nCC
  if (length(nCC) == 0 || nCC < 1) {
    stop("Invalid number of canonical components")
  }
  
  scalePCs <- if (length(object@scalePCs) == 0) FALSE else object@scalePCs
  
  # Get PCA matrices
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
  
  # Get skrCCA weights
  W_list <- object@skrCCAOut[[sigma_name]]
  
  # Convert cell scores from PCA * weights for each cell type
  cell_scores_list <- list()
  for (ct in cts) {
    if (!ct %in% names(PCmats)) {
      stop(paste("PCA results missing for cell type:", ct))
    }
    if (!ct %in% names(W_list)) {
      stop(paste("skrCCA weights missing for cell type:", ct))
    }
    
    # Compute cell scores: PCA_matrix %*% skrCCA_weights
    cell_scores_list[[ct]] <- PCmats[[ct]] %*% W_list[[ct]]
    
    # Set row names using cell names from the object (PCA doesn't preserve cell names)
    ct_indices <- which(object@cellTypesSub == ct)
    if (length(ct_indices) != nrow(cell_scores_list[[ct]])) {
      stop(paste("Dimension mismatch for cell type", ct, 
                ": expected", length(ct_indices), "cells but got", nrow(cell_scores_list[[ct]])))
    }
    
    # Get cell names from the object
    cell_names <- rownames(object@normalizedDataSub)[ct_indices]
    if (is.null(cell_names)) {
      # Fallback: create row names from cell indices
      cell_names <- paste0(ct, "_cell_", ct_indices)
    }
    
    # Set row names to match cell IDs
    rownames(cell_scores_list[[ct]]) <- cell_names
    
    # Set column names for CCs
    colnames(cell_scores_list[[ct]]) <- paste0("CC_", seq_len(ncol(cell_scores_list[[ct]])))
    
    if (verbose) {
      cat("Cell scores for", ct, ": dimensions =", dim(cell_scores_list[[ct]]), 
          ", row names =", length(rownames(cell_scores_list[[ct]])), "\n")
    }
  }
  
  # Use the main function with computed cell scores
  return(getTransferSelfBidirCorr(tar_obj = object,
                                 transfer_cell_scores = cell_scores_list,
                                 sigma_choice = sigma_choice,
                                 calculationMode = calculationMode,
                                 normalize_K = normalize_K,
                                 filter_kernel = filter_kernel,
                                 K_row_sum_cutoff = K_row_sum_cutoff,
                                 K_col_sum_cutoff = K_col_sum_cutoff,
                                 sigma_choice_tar = sigma_choice,
                                 verbose = verbose))
}
