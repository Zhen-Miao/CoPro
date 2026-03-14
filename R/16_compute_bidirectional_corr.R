#' Sinkhorn-Knopp algorithm for matrix scaling (optimized)
#'
#' @param A Non-negative input matrix
#' @param tol Convergence tolerance
#' @param max_iter Maximum number of iterations
#' @return Scaled matrix
#' @noRd
sinkhorn_knopp <- function(A, tol = 1e-8, max_iter = 1000) {
  if (any(A < 0, na.rm = TRUE)) {
    stop("Input matrix A must be non-negative.")
  }
  
  n <- nrow(A)
  m <- ncol(A)
  eps <- .Machine$double.eps
  
  u <- rep(1.0, n)
  v <- rep(1.0, m)
  
  # Pre-allocate
  Av <- numeric(n)
  Atu <- numeric(m)
  
  for (i in seq_len(max_iter)) {
    u_prev <- u
    v_prev <- v
    
    # Row normalization
    Av <- A %*% v
    u <- ifelse(Av > eps, 1.0 / Av, 0.0)
    
    # Column normalization
    Atu <- crossprod(A, u)
    v <- ifelse(Atu > eps, 1.0 / Atu, 0.0)
    
    # Simple convergence check (like fast1)
    if (max(abs(u - u_prev), abs(v - v_prev)) < tol) {
      break
    }
  }
  
  # Use fast1's scaling method since it's proven faster
  scaled_matrix <- (u %*% t(rep(1, m))) * A * (rep(1, n) %*% t(v))
  
  return(scaled_matrix)
}



#' Vectorized correlation calculation for all CC indices
#' @param A_W1_all Matrix of projected scores for cell type 1 (cells x nCC)
#' @param B_W2_all Matrix of projected scores for cell type 2 (cells x nCC)  
#' @param K Kernel matrix (already filtered)
#' @param normalize_K Normalization method
#' @return Vector of correlation values for each CC
#' @noRd
.computeAllCCCorrelations <- function(A_W1_all, B_W2_all, K, normalize_K) {
  nCC <- ncol(A_W1_all)
  corr_values <- numeric(nCC)
  
  # Apply kernel normalization once based on method
  if (normalize_K == "row_or_col") {
    # Optimized row/column normalization (guard against zero sums)
    K_row_sum <- rowSums(K)
    K_col_sum <- colSums(K)
    K_row_sum[K_row_sum == 0] <- 1  # avoid division by zero
    K_col_sum[K_col_sum == 0] <- 1  # avoid division by zero
    K_row_norm <- K / K_row_sum
    K_col_norm <- sweep(K, 2, K_col_sum, "/")

    # Compute kernel-weighted matrices for all CC at once
    KA_all <- crossprod(K_row_norm, A_W1_all)  # nCC columns
    KB_all <- K_col_norm %*% B_W2_all         # nCC columns
    
  } else if (normalize_K == "sinkhorn_knopp") {
    K_norm <- sinkhorn_knopp(K)
    KA_all <- crossprod(K_norm, A_W1_all)
    KB_all <- K_norm %*% B_W2_all
    
  } else { # "none"
    KA_all <- crossprod(K, A_W1_all)
    KB_all <- K %*% B_W2_all
  }
  
  # Ultra-fast vectorized correlation calculation
  # Compute correlations for all CC at once using matrix operations
  if (nrow(KA_all) > 1 && nrow(B_W2_all) > 1) {
    # Vectorized correlation using efficient covariance calculation
    cor1_all <- numeric(nCC)
    cor2_all <- numeric(nCC)
    
    # Center the matrices once
    KA_centered <- scale(KA_all, center = TRUE, scale = FALSE)
    B_centered <- scale(B_W2_all, center = TRUE, scale = FALSE)
    A_centered <- scale(A_W1_all, center = TRUE, scale = FALSE)
    KB_centered <- scale(KB_all, center = TRUE, scale = FALSE)
    
    # Compute correlations efficiently
    for (cc in seq_len(nCC)) {
      # Faster correlation using centered data
      ka <- KA_centered[, cc]
      b2 <- B_centered[, cc]
      a1 <- A_centered[, cc]
      kb <- KB_centered[, cc]
      
      # Correlation = sum(x*y) / sqrt(sum(x^2) * sum(y^2))
      denom1 <- sqrt(sum(ka^2) * sum(b2^2))
      denom2 <- sqrt(sum(a1^2) * sum(kb^2))
      cor1_all[cc] <- if (denom1 < 1e-12) 0 else sum(ka * b2) / denom1
      cor2_all[cc] <- if (denom2 < 1e-12) 0 else sum(a1 * kb) / denom2
    }
    
    corr_values <- (cor1_all + cor2_all) * 0.5
    
  } else {
    # Fallback for edge cases
    corr_values[] <- 1.0
  }
  
  return(corr_values)
}

# Optimized filtering function
.filter_kernel_matrix <- function(K, A_w1, B_w2, filter_kernel = TRUE,
 K_row_sum_cutoff = 5e-3, K_col_sum_cutoff = 5e-3){

  if(filter_kernel){
    # Vectorized row filtering
    K_row_sum <- rowSums(K)
    row_keep <- K_row_sum > K_row_sum_cutoff
    
    if (any(!row_keep)) {
        K <- K[row_keep, , drop = FALSE]
        A_w1 <- A_w1[row_keep, , drop = FALSE]
        if (length(row_keep) > 0) {
            cat("After row filtering: kept", sum(row_keep), "of", length(row_keep), "rows\n")
        }
    }
    
    # Vectorized column filtering
    K_col_sum <- colSums(K)
    col_keep <- K_col_sum > K_col_sum_cutoff
    
    if (any(!col_keep)) {
        K <- K[, col_keep, drop = FALSE]
        B_w2 <- B_w2[col_keep, , drop = FALSE]
        if (length(col_keep) > 0) {
            cat("After column filtering: kept", sum(col_keep), "of", length(col_keep), "columns\n")
        }
    }

    return(list(K = K, A_w1 = A_w1, B_w2 = B_w2))
  }

  # When filter_kernel=FALSE, return inputs unchanged
  return(list(K = K, A_w1 = A_w1, B_w2 = B_w2))
}

#' Mean bidirectional spatial cross-correlation (internal helper)
#'
#' Computes the average of two correlations involving the kernel-weighted
#' vectors: cor(t(K) %*% A_w1, B_w2) and cor(A_w1, K %*% B_w2).
#'
#' @param A_w1 Numeric column vector (or matrix coerced to column) of scores for cell type 1
#' @param B_w2 Numeric column vector (or matrix coerced to column) of scores for cell type 2
#' @param K Numeric kernel matrix relating cells of type 1 and type 2
#' @param normalize_K whether to normalize the kernel matrix, one of "none", "sinkhorn_knopp", "row_or_col" 
#' @param K_row_sum_cutoff Numeric; cutoff for row sums when normalizing kernel matrix
#' @param K_col_sum_cutoff Numeric; cutoff for column sums when normalizing kernel matrix
#' @return A single numeric value: the mean of the two correlations
#' @keywords internal
#' @noRd
.computeSpatialCrossCorrelation <- function(A_w1, B_w2, K, 
normalize_K = c("row_or_col", "sinkhorn_knopp", "none"), filter_kernel = TRUE,
 K_row_sum_cutoff = 5e-2, K_col_sum_cutoff = 5e-2) {

  normalize_K <- match.arg(normalize_K)
  # Ensure vectors are column matrices
  if (!is.matrix(A_w1)) A_w1 <- matrix(A_w1, ncol = 1)
  if (!is.matrix(B_w2)) B_w2 <- matrix(B_w2, ncol = 1)
  
     # row and column filtering
   if(filter_kernel){
     filtered_result <- .filter_kernel_matrix(K, A_w1, B_w2, TRUE, K_row_sum_cutoff, K_col_sum_cutoff)
     K <- filtered_result$K
     A_w1 <- filtered_result$A_w1
     B_w2 <- filtered_result$B_w2
   }

  if(normalize_K == "row_or_col"){
    # Optimized row/column normalization (guard against zero sums)
    K_row_sum <- rowSums(K)
    K_col_sum <- colSums(K)
    K_row_sum[K_row_sum == 0] <- 1  # avoid division by zero
    K_col_sum[K_col_sum == 0] <- 1  # avoid division by zero

    # More efficient normalization using vectorized operations
    K_row_norm <- K / K_row_sum  # Broadcasting division
    K_col_norm <- sweep(K, 2, K_col_sum, "/")  # More efficient than t(t(K) / K_col_sum)

    # Use crossprod for more efficient matrix multiplication
    KA <- crossprod(K_row_norm, A_w1)  # Equivalent to t(K_row_norm) %*% A_w1
    KB <- K_col_norm %*% B_w2

  }else if(normalize_K == "none"){
    # Compute kernel-weighted vectors - use crossprod for efficiency
    KA <- crossprod(K, A_w1)  # Equivalent to t(K) %*% A_w1 but faster
    KB <- K %*% B_w2
  }else if(normalize_K == "sinkhorn_knopp"){
    K <- sinkhorn_knopp(K)  # Use the optimized version
    KA <- crossprod(K, A_w1)  # More efficient transpose multiplication
    KB <- K %*% B_w2
  }

  # Optimized correlation computation using covariance
  # For vectors, correlation can be computed more efficiently
  if (length(KA) == 1 && length(B_w2) == 1) {
    # Single values case
    cor1 <- 1.0
  } else {
    # Use more efficient correlation calculation
    cor1 <- stats::cor(as.vector(KA), as.vector(B_w2))
  }
  
  if (length(A_w1) == 1 && length(KB) == 1) {
    # Single values case  
    cor2 <- 1.0
  } else {
    cor2 <- stats::cor(as.vector(A_w1), as.vector(KB))
  }
  
  return((cor1 + cor2) * 0.5)  # Multiplication is faster than division
}

#' Compute Bidirectional Correlation
#'
#' Computes the mean of two correlations: cor(t(A_w1) %*% K, B_w2) and
#' cor(A_w1, K %*% B_w2), where A_w1 and B_w2 are cell scores derived from
#' PCA matrices multiplied by skrCCA weights for the corresponding cell types,
#' and K is the kernel matrix between the two cell types.
#'
#' The results are stored in the slot `bidirCorrelation`.
#'
#' For multi-slide (`CoProMulti`) objects, results can be computed per slide
#' (default, returns a data frame with a `slideID` column) or aggregated across
#' slides (returns a data frame with an `aggregateCorrelation` column).
#'
#' @name computeBidirCorrelation
#' @param object A `CoPro` or `CoProMulti` object containing skrCCA results,
#'  PCA results, and kernel matrices.
#' @param calculationMode (CoProMulti only) "perSlide" or "aggregate".
#'  Ignored for single-slide objects.
#' @param normalize_K whether to normalize the kernel matrix, one of "none", "sinkhorn_knopp", "row_or_col" 
#' @param filter_kernel whether to filter the kernel matrix, default is TRUE
#' @param K_row_sum_cutoff Numeric; cutoff for row sums when normalizing kernel matrix.
#'  Default 5e-3.
#' @param K_col_sum_cutoff Numeric; cutoff for column sums when normalizing kernel matrix.
#'  Default 5e-3.
#' @return The input object with `@bidirCorrelation` populated (a list keyed by
#'  sigma names). For single-slide objects, each entry is a data frame with
#'  columns `sigmaValues`, `cellType1`, `cellType2`, `CC_index`, `bidirCorrelation`.
#'  For multi-slide objects, the columns depend on `calculationMode`.
#' @examples
#' # Assuming `obj` is a prepared CoProSingle with PCA, kernels, and skrCCA:
#' # obj <- computeBidirCorrelation(obj)
#'
#' # For CoProMulti per-slide results:
#' # objm <- computeBidirCorrelation(objm, calculationMode = "perSlide")
#'
#' # For CoProMulti aggregate results:
#' # objm <- computeBidirCorrelation(objm, calculationMode = "aggregate")
#' @export
setGeneric(
  "computeBidirCorrelation",
  function(object, calculationMode = "perSlide", 
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
   K_row_sum_cutoff = 5e-2, K_col_sum_cutoff = 5e-2) standardGeneric("computeBidirCorrelation")
)

.checkInputBidirCorr <- function(object) {
  if (length(object@skrCCAOut) == 0) {
    stop("CCA results are not available. Please run CCA first.")
  }
  if (length(object@kernelMatrices) == 0) {
    stop("Kernel matrices are not available. Please compute the kernel matrices first.")
  }
  if (length(object@pcaGlobal) == 0) {
    stop("PCA results missing. Please run computePCA first.")
  }

  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    warning("no cell type of interest specified, using all cell types to run the analysis")
    cts <- unique(object@cellTypesSub)
  }

  if (length(object@scalePCs) == 0) stop("object@scalePCs not specified")
  scalePCs <- object@scalePCs

  if (length(object@sigmaValues) == 0) stop("`sigmaValues` is empty, please specify")
  sigmaValues <- object@sigmaValues
  nCC <- object@nCC

  return(list(cts = cts, scalePCs = scalePCs, sigmaValues = sigmaValues, nCC = nCC))
}

.computeBidirCorrCore <- function(object, cts, scalePCs, sigmaValues, nCC, normalize_K = c("row_or_col", "sinkhorn_knopp", "none"), filter_kernel = TRUE, K_row_sum_cutoff = 5e-3, K_col_sum_cutoff = 5e-3) {
  normalize_K <- match.arg(normalize_K)
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

  if (length(cts) == 1) {
    pair_cell_types <- matrix(c(cts, cts), nrow = 2, ncol = 1)
  } else {
    pair_cell_types <- utils::combn(cts, 2)
  }

  correlation_value <- vector("list", length = length(sigmaValues))
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  names(correlation_value) <- sigma_names
  n_pairs <- ncol(pair_cell_types)

  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    sigma_val <- sigmaValues[tt]
    
    # Pre-allocate results as vectors (much faster than dataframe indexing)
    n_total <- n_pairs * nCC
    result_corr <- numeric(n_total)
    result_idx <- 0L
    
    # Cache kernel matrices for this sigma
    kernel_cache <- list()
    
    for (pp in seq_len(n_pairs)) {
      cellType1 <- pair_cell_types[1, pp]
      cellType2 <- pair_cell_types[2, pp]
      
      # Get cached kernel matrix
      cache_key <- paste(cellType1, cellType2, sep = ":")
      if (!cache_key %in% names(kernel_cache)) {
        kernel_cache[[cache_key]] <- getKernelMatrix(object, sigma = sigma_val, 
                                                    cellType1 = cellType1, 
                                                    cellType2 = cellType2, 
                                                    verbose = FALSE)
      }
      K_orig <- kernel_cache[[cache_key]]
      
      # Get PCA matrices
      A <- PCmats[[cellType1]]
      B <- PCmats[[cellType2]]
      
      # Get all weight vectors for this pair
      W1_all <- object@skrCCAOut[[t]][[cellType1]]
      W2_all <- object@skrCCAOut[[t]][[cellType2]]
      
      # Compute all projections at once
      A_W1_all <- A %*% W1_all  # Matrix: cells x nCC
      B_W2_all <- B %*% W2_all  # Matrix: cells x nCC
      
      # MAJOR OPTIMIZATION: Filter kernel once per pair, not per CC
      if (filter_kernel) {
        # Compute row and column masks from the original kernel consistently
        row_keep <- rowSums(K_orig) > K_row_sum_cutoff
        K_row_filtered <- K_orig[row_keep, , drop = FALSE]
        col_keep <- colSums(K_row_filtered) > K_col_sum_cutoff

        K <- K_row_filtered[, col_keep, drop = FALSE]
        A_W1_all <- A_W1_all[row_keep, , drop = FALSE]
        B_W2_all <- B_W2_all[col_keep, , drop = FALSE]
      } else {
        K <- K_orig
      }
      
      # Vectorized correlation calculation for all CC at once
      corr_values <- .computeAllCCCorrelations(A_W1_all, B_W2_all, K, normalize_K)
      
      # Store results efficiently
      start_idx <- result_idx + 1L
      end_idx <- result_idx + nCC
      result_corr[start_idx:end_idx] <- corr_values
      result_idx <- end_idx
    }
    
    # Build dataframe once at the end (much faster)
    correlation_value[[t]] <- data.frame(
      sigmaValues = rep(sigma_val, n_total),
      cellType1 = rep(pair_cell_types[1, ], times = nCC),
      cellType2 = rep(pair_cell_types[2, ], times = nCC),
      CC_index = rep(seq_len(nCC), each = n_pairs),
      bidirCorrelation = result_corr,
      stringsAsFactors = FALSE
    )
  }

  object@bidirCorrelation <- correlation_value
  return(object)
}

#' @rdname computeBidirCorrelation
#' @aliases computeBidirCorrelation,CoPro-method
#' @importFrom utils combn
#' @export
setMethod(
  "computeBidirCorrelation", "CoPro",
  function(object, calculationMode = "perSlide", 
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
   K_row_sum_cutoff = 5e-2, K_col_sum_cutoff = 5e-2) { # calculationMode ignored for single-slide
    normalize_K <- match.arg(normalize_K)
    input_check <- .checkInputBidirCorr(object)
    cts <- input_check$cts
    scalePCs <- input_check$scalePCs
    sigmaValues <- input_check$sigmaValues
    nCC <- input_check$nCC

    object <- .computeBidirCorrCore(object, cts = cts, scalePCs = scalePCs,
                                    sigmaValues = sigmaValues, nCC = nCC, normalize_K = normalize_K,
                                    filter_kernel = filter_kernel,
                                    K_row_sum_cutoff = K_row_sum_cutoff, K_col_sum_cutoff = K_col_sum_cutoff)
    return(object)
  }
)

.checkInputBidirCorrMulti <- function(object) {
  if (length(object@skrCCAOut) == 0) stop("skrCCA results missing. Run runSkrCCAMulti.")
  if (length(object@pcaResults) == 0) stop("PCA results missing. Please run computePCA first.")
  if (length(object@pcaGlobal) == 0) stop("PCA global results missing. Please run computePCA first.")
  if (length(object@kernelMatrices) == 0) stop("Kernel matrices missing.")
  cts <- object@cellTypesOfInterest
  if (length(cts) < 1) stop("Need at least one cell type.")
  slides <- getSlideList(object)
  sigmas_run <- names(object@skrCCAOut)
  if (length(sigmas_run) == 0) stop("No skrCCA results found.")
  nCC <- object@nCC
  if (is.null(nCC) || length(nCC) == 0) {
    first_sigma <- sigmas_run[1]
    first_ct <- cts[1]
    if (!is.null(object@skrCCAOut[[first_sigma]][[first_ct]])) {
      nCC <- ncol(object@skrCCAOut[[first_sigma]][[first_ct]])
    } else {
      stop("Cannot infer nCC from skrCCA results")
    }
  }
  return(list(cts = cts, slides = slides, sigmas_run = sigmas_run, nCC = nCC))
}

.computeBidirCorrCoreMulti <- function(object, cts, slides, sigmas_run, nCC, calculationMode = "perSlide", normalize_K = c("row_or_col", "sinkhorn_knopp", "none"), filter_kernel = TRUE, K_row_sum_cutoff = 5e-3, K_col_sum_cutoff = 5e-3) {
  normalize_K <- match.arg(normalize_K)
  if (length(cts) == 1) {
    pair_cell_types <- matrix(c(cts, cts), nrow = 2, ncol = 1)
  } else {
    pair_cell_types <- utils::combn(cts, 2)
  }

  correlation_results <- setNames(vector("list", length = length(sigmas_run)), sigmas_run)

  for (sig_name in sigmas_run) {
    W_list_sigma <- object@skrCCAOut[[sig_name]]

    if (calculationMode == "perSlide") {
      correlation_per_slide <- setNames(vector("list", length = length(slides)), slides)
      for (sID in slides) {
        df_slide <- data.frame(
          sigmaValue = character(),
          slideID = character(),
          cellType1 = character(), cellType2 = character(),
          CC_index = integer(), bidirCorrelation = numeric(),
          stringsAsFactors = FALSE
        )

        X_list_slide <- object@pcaResults[[sID]]
        if (is.null(X_list_slide)) next

        for (pp in seq_len(ncol(pair_cell_types))) {
          ct_i <- pair_cell_types[1, pp]
          ct_j <- pair_cell_types[2, pp]

          X_i <- X_list_slide[[ct_i]]
          X_j <- X_list_slide[[ct_j]]
          K_ij <- tryCatch({
            getKernelMatrix(object,
                            sigma = as.numeric(gsub("sigma_", "", sig_name)),
                            cellType1 = ct_i,
                            cellType2 = ct_j,
                            slide = sID,
                            verbose = FALSE)
          }, error = function(e) NULL)

          if (is.null(X_i) || is.null(X_j) || is.null(K_ij) || nrow(X_i) == 0 || nrow(X_j) == 0) next

          for (cc in seq_len(nCC)) {
            w_i <- W_list_sigma[[ct_i]][, cc, drop = FALSE]
            w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]

            Xiw <- X_i %*% w_i
            Xjw <- X_j %*% w_j

            corr_val <- .computeSpatialCrossCorrelation(Xiw, Xjw, K_ij, normalize_K = normalize_K, filter_kernel = filter_kernel, K_row_sum_cutoff = K_row_sum_cutoff, K_col_sum_cutoff = K_col_sum_cutoff)

            df_slide <- rbind(df_slide, data.frame(
              sigmaValue = as.numeric(gsub("sigma_", "", sig_name)),
              slideID = sID,
              cellType1 = ct_i, cellType2 = ct_j,
              CC_index = cc, bidirCorrelation = as.numeric(corr_val),
              stringsAsFactors = FALSE
            ))
          }
        }
        correlation_per_slide[[sID]] <- df_slide
      }
      correlation_results[[sig_name]] <- do.call(rbind, correlation_per_slide)

    } else { # aggregate mode
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
            X_list_slide <- object@pcaResults[[sID]]
            if (is.null(X_list_slide)) next
            X_i <- X_list_slide[[ct_i]]
            X_j <- X_list_slide[[ct_j]]
            K_ij <- tryCatch({
              getKernelMatrix(object,
                              sigma = as.numeric(gsub("sigma_", "", sig_name)),
                              cellType1 = ct_i,
                              cellType2 = ct_j,
                              slide = sID,
                              verbose = FALSE)
            }, error = function(e) NULL)
            if (is.null(X_i) || is.null(X_j) || is.null(K_ij) || nrow(X_i) == 0 || nrow(X_j) == 0) next

            w_i <- W_list_sigma[[ct_i]][, cc, drop = FALSE]
            w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]

            Xiw <- X_i %*% w_i
            Xjw <- X_j %*% w_j

            cat("current cell type pair: ", ct_i, " and ", ct_j, "\n")

            corr_val <- .computeSpatialCrossCorrelation(Xiw, Xjw, K_ij, normalize_K = normalize_K, filter_kernel = filter_kernel, K_row_sum_cutoff = K_row_sum_cutoff, K_col_sum_cutoff = K_col_sum_cutoff)

            if (is.finite(corr_val) && !is.na(corr_val)) {
              sum_corr <- sum_corr + as.numeric(corr_val)
              valid_slides <- valid_slides + 1
            }
          }
          if (valid_slides > 0) {
            df_agg <- rbind(df_agg, data.frame(
              sigmaValue = as.numeric(gsub("sigma_", "", sig_name)),
              cellType1 = ct_i, cellType2 = ct_j,
              CC_index = cc, aggregateCorrelation = sum_corr / valid_slides,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      correlation_results[[sig_name]] <- df_agg
    }
  }

  object@bidirCorrelation <- correlation_results
  return(object)
}

#' @rdname computeBidirCorrelation
#' @aliases computeBidirCorrelation,CoProMulti-method
#' @importFrom utils combn
#' @export
setMethod("computeBidirCorrelation", "CoProMulti", function(object, calculationMode = "perSlide", 
normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
filter_kernel = TRUE,
K_row_sum_cutoff = 5e-2, K_col_sum_cutoff = 5e-2) {
  normalize_K <- match.arg(normalize_K)
  if (!calculationMode %in% c("perSlide", "aggregate")) {
    stop("calculationMode must be either 'perSlide' or 'aggregate'")
  }

  input_check <- .checkInputBidirCorrMulti(object)
  cts <- input_check$cts
  slides <- input_check$slides
  sigmas_run <- input_check$sigmas_run
  nCC <- input_check$nCC

  object <- .computeBidirCorrCoreMulti(object, cts = cts, slides = slides,
                                       sigmas_run = sigmas_run, nCC = nCC,
                                       calculationMode = calculationMode, normalize_K = normalize_K,
                                       filter_kernel = filter_kernel,
                                       K_row_sum_cutoff = K_row_sum_cutoff, K_col_sum_cutoff = K_col_sum_cutoff)
  return(object)
})

#' Ensure object has bidirCorrelation slot
#'
#' Safely upgrades a `CoPro`/`CoProSingle`/`CoProMulti` object created with an
#' older package version (without the `bidirCorrelation` slot) to the current
#' class definition by recreating an instance of the same class and copying
#' over all available slots. If `bidirCorrelation` already exists, the object is
#' returned unchanged.
#'
#' @param object A `CoProSingle`, `CoProMulti`, or `CoProm` object
#' @return The same object class with a valid `bidirCorrelation` slot
#' @examples
#' # Upgrade legacy object to include bidirCorrelation slot
#' # obj <- ensureBidirCorrelationSlot(obj)
#' @export
#' @importFrom methods new slot getSlots is
ensureBidirCorrelationSlot <- function(object) {
  if (!methods::is(object, "CoPro")) {
    stop("ensureBidirCorrelationSlot expects a CoPro-derived object")
  }

  # If slot is already accessible, return as-is
  has_bidir <- tryCatch({ methods::slot(object, "bidirCorrelation"); TRUE },
                        error = function(e) FALSE)
  if (has_bidir) return(object)

  # Recreate a fresh instance of the same concrete class and copy slots
  cls <- class(object)[1]
  new_obj <- methods::new(cls)
  all_slots <- names(methods::getSlots(cls))

  for (s in all_slots) {
    if (s == "bidirCorrelation") next
    val <- tryCatch(methods::slot(object, s), error = function(e) NULL)
    if (!is.null(val)) methods::slot(new_obj, s) <- val
  }

  # Initialize new slot(s)
  methods::slot(new_obj, "bidirCorrelation") <- list()
  if ("bidirCorrelationPermu" %in% all_slots) {
    permu_val <- tryCatch(methods::slot(object, "bidirCorrelationPermu"), error = function(e) NULL)
    if (is.null(permu_val)) methods::slot(new_obj, "bidirCorrelationPermu") <- list()
  }

  return(new_obj)
}

