#' Compute Normalized Correlation (approximation)
#'
#' This method calculates the normalized correlation between pairs of cell types
#' based on CCA weights and the respective kernel matrix. It uses the
#' whitened-Frobenius norm ||R_x^(1/2) K_c R_y^(1/2)||_F of the kernel for
#' normalization (R_x, R_y = matched-sigma within-type kernels).
#'
#' @param object A `CoPro` or `CoProMulti` object containing CCA results and kernel matrices.
#' @param tol tolerance for approximate SVD calculation
#' @param calculationMode (for CoProMulti only) either "perSlide" or "aggregate",
#'   for single slide analysis, it is ignored, with default value "perSlide".
#' @return The object with the normalized correlation value
#' between any pair of cell types
#' added as a new slot, `normalizedCorrelation`.
#' @family scores-and-correlation
#' @seealso [runSkrCCA()], [computeBidirCorrelation()],
#'   [computeGeneAndCellScores()]
#' @export
#'
setGeneric(
  "computeNormalizedCorrelation",
  function(object, tol = 1e-4, calculationMode = "perSlide") standardGeneric("computeNormalizedCorrelation")
)


.checkInputNormCorr <- function(object) {
  ## Check for required components
  if (length(object@skrCCAOut) == 0) {
    stop("CCA results are not available. Please run CCA first.")
  }
  if (length(object@kernelMatrices) == 0) {
    stop(paste(
      "Kernel matrices are not available.",
      "Please compute the kernel matrices first."
    ))
  }
  if (length(object@pcaGlobal) == 0) {
    stop("PCA results missing. Please run computePCA first.")
  }
  
  ## choose cell types
  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    warning(paste(
      "no cell type of interest specified,",
      "using all cell types to run the analysis"
    ))
    cts <- unique(object@cellTypesSub)
  }

  ## check whether the PCs are being scaled prior to CCA
  if (length(object@scalePCs) == 0) {
    stop("object@scalePCs not specified")
  }
  scalePCs <- object@scalePCs

  ## check sigmaValues
  if (length(object@sigmaValues) == 0) {
    stop("`sigmaValues` is empty, please specify")
  }

  sigmaValues <- object@sigmaValues
  nCC <- object@nCC

  return(list(cts = cts, scalePCs = scalePCs, sigmaValues = sigmaValues, nCC = nCC))
}

#' Whitened-Frobenius null SD of the bilinear statistic a' K b
#'
#' Internal. Returns \eqn{\|R_x^{1/2} K_c R_y^{1/2}\|_F}, the distribution-free
#' null standard deviation of \eqn{T = a' K b}, where \eqn{K_c} is the
#' double-centered cross-kernel and \eqn{R_x, R_y} are the within-type
#' correlation operators. This replaces the spectral norm \eqn{\|K\|_2}, which
#' is scale-blind and rails bandwidth selection to the noise floor. With
#' \code{Rx}/\code{Ry} omitted it degrades to the un-whitened \eqn{\|K_c\|_F}
#' (i.e. \eqn{R_x = R_y = I}).
#'
#' @param K Cross-type kernel matrix (double-centered internally).
#' @param Rx,Ry Within-type kernels used as correlation operators; symmetrized
#'   defensively (a no-op when the kernel is un-normalized). `NULL` to skip
#'   whitening.
#' @return Scalar whitened-Frobenius norm.
#' @keywords internal
.whitenedFrobNorm <- function(K, Rx = NULL, Ry = NULL) {
  ## Sparse fixed-radius kernels can be very large. Double-centering a sparse
  ## matrix explicitly makes it dense, and coercing the matched within-type
  ## kernels to base matrices compounds that memory cost. Use the equivalent
  ## sparse-plus-low-rank expansion whenever all three kernels are sparse.
  ##
  ## K_c = K + U V', with
  ##   U = [-rowMeans(K), 1]
  ##   V = [1, mean(K) - colMeans(K)].
  ## For symmetric Rx/Ry,
  ##   <Rx K_c Ry, K_c>
  ## = <Rx K Ry, K> + 2 <Rx K Ry, U V'>
  ##   + sum((U' Rx U) * (V' Ry V)).
  ## This is algebraically identical to materializing K_c, but retains sparse
  ## matrix multiplication for the expensive term.
  all_sparse <- inherits(K, "sparseMatrix") &&
    (is.null(Rx) || inherits(Rx, "sparseMatrix")) &&
    (is.null(Ry) || inherits(Ry, "sparseMatrix"))

  if (all_sparse) {
    nr <- nrow(K)
    nc <- ncol(K)
    rmean <- as.numeric(rowMeans(K))
    cmean <- as.numeric(colMeans(K))
    grand_mean <- mean(cmean)

    if (is.null(Rx) || is.null(Ry)) {
      ## ||H_r K H_c||_F^2 without forming the dense centered matrix.
      norm2 <- sum(K * K) - nc * sum(rmean^2) - nr * sum(cmean^2) +
        nr * nc * grand_mean^2
      return(sqrt(max(as.numeric(norm2), 0)))
    }

    Rx <- (Rx + t(Rx)) / 2
    Ry <- (Ry + t(Ry)) / 2
    M <- (Rx %*% K) %*% Ry

    U <- cbind(-rmean, rep.int(1, nr))
    V <- cbind(rep.int(1, nc), grand_mean - cmean)
    base_term <- as.numeric(sum(M * K))
    cross_term <- as.numeric(sum(U * (M %*% V)))
    rank_term <- as.numeric(sum(crossprod(U, Rx %*% U) *
                                  crossprod(V, Ry %*% V)))
    return(sqrt(max(base_term + 2 * cross_term + rank_term, 0)))
  }

  K <- as.matrix(K)
  ## double-center: a' K b = a' K_c b for centered scores; Var_0(T) uses K_c
  Kc <- K - rowMeans(K) - rep(colMeans(K), each = nrow(K)) + mean(K)
  if (is.null(Rx) || is.null(Ry)) {
    return(sqrt(sum(Kc * Kc)))                 # ||K_c||_F  (R_x = R_y = I)
  }
  ## symmetrize within-type kernels into correlation operators (no-op when the
  ## kernel is un-normalized, i.e. already symmetric PSD with unit diagonal)
  Rx <- as.matrix(Rx); Ry <- as.matrix(Ry)
  Rx <- (Rx + t(Rx)) / 2
  Ry <- (Ry + t(Ry)) / 2
  ## ||R_x^{1/2} K_c R_y^{1/2}||_F^2 = tr(R_x K_c R_y K_c^T) = sum((R_x K_c R_y) * K_c)
  M <- (Rx %*% Kc) %*% Ry
  sqrt(max(sum(M * Kc), 0))
}

.kernelNormalizerKey <- function(sigma, cellType1, cellType2, slide = NULL) {
  paste(
    format(sigma, scientific = FALSE, trim = TRUE),
    if (is.null(slide)) "single" else slide,
    cellType1, cellType2, sep = "|"
  )
}

.kernelMatrixSignature <- function(x) {
  if (is.null(x)) return("NULL")
  values <- if (inherits(x, "sparseMatrix")) x@x else as.numeric(x)
  paste(
    nrow(x), ncol(x), length(values),
    format(sum(values), digits = 16),
    format(sum(values^2), digits = 16),
    sep = ":"
  )
}

.kernelNormalizerSignature <- function(K, Rx = NULL, Ry = NULL) {
  paste(
    .kernelMatrixSignature(K),
    .kernelMatrixSignature(Rx),
    .kernelMatrixSignature(Ry),
    sep = "||"
  )
}

.cacheKernelNormalizer <- function(cache, key, K, Rx, Ry, value) {
  cache[[key]] <- list(
    value = as.numeric(value),
    signature = .kernelNormalizerSignature(K, Rx, Ry)
  )
  cache
}

.readKernelNormalizer <- function(cache, key, K, Rx, Ry) {
  entry <- cache[[key]]
  if (is.null(entry) || !is.finite(entry$value)) return(NULL)
  if (!identical(entry$signature, .kernelNormalizerSignature(K, Rx, Ry))) {
    return(NULL)
  }
  entry$value
}

.computeCrossKernelNorm <- function(object, tol = 1e-4, cts, scalePCs,
 sigmaValues, nCC, pair_cell_types) {
  ## Whitened-Frobenius normalizer ||R_x^{1/2} K_c R_y^{1/2}||_F per (sigma,
  ## pair), with R_x, R_y the matched-sigma within-type kernels. Replaces the
  ## old spectral norm ||K||_2.
  message("Calculating whitened-Frobenius normalizers, this may take a while.")

  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  norm_K12 <- setNames(vector(mode = "list", length = length(sigma_names)),
                       sigma_names)
  normalizer_cache <- attr(object, "kernelNormalizerCache", exact = TRUE)
  if (is.null(normalizer_cache)) normalizer_cache <- list()

  for (t in sigma_names) {
    norm_K12[[t]] <- setNames(vector(mode = "list", length = length(cts)),
                              cts)
    for (i in cts) {
      norm_K12[[t]][[i]] <- setNames(vector(mode = "list",
                                            length = length(cts)), cts)
    }
  }

  for (t in sigma_names) {
    sigma_val <- as.numeric(gsub("sigma_", "", t))
    for (pp in seq_len(ncol(pair_cell_types))) {
      cellType1 <- pair_cell_types[1, pp]
      cellType2 <- pair_cell_types[2, pp]
      K <- getKernelMatrix(object, sigma = sigma_val,
                           cellType1 = cellType1, cellType2 = cellType2,
                           verbose = FALSE)
      ## matched-sigma within-type kernels serve as the whitening operators;
      ## if unavailable, .whitenedFrobNorm falls back to ||K_c||_F
      Rx <- tryCatch(getKernelMatrix(object, sigma = sigma_val,
                       cellType1 = cellType1, cellType2 = cellType1,
                       verbose = FALSE), error = function(e) NULL)
      Ry <- tryCatch(getKernelMatrix(object, sigma = sigma_val,
                       cellType1 = cellType2, cellType2 = cellType2,
                       verbose = FALSE), error = function(e) NULL)
      cache_key <- .kernelNormalizerKey(
        sigma_val, cellType1, cellType2, slide = NULL
      )
      nrm <- .readKernelNormalizer(normalizer_cache, cache_key, K, Rx, Ry)
      if (is.null(nrm)) {
        nrm <- .whitenedFrobNorm(K, Rx, Ry)
        normalizer_cache <- .cacheKernelNormalizer(
          normalizer_cache, cache_key, K, Rx, Ry, nrm
        )
      }
      norm_K12[[t]][[cellType1]][[cellType2]] <- nrm
      # Frobenius norm is transpose-invariant: store the mirror for symmetry
      if (cellType1 != cellType2) {
        norm_K12[[t]][[cellType2]][[cellType1]] <- nrm
      }
    }
  }

  message("Finished calculating whitened-Frobenius normalizers.")

  attr(norm_K12, "kernelNormalizerCache") <- normalizer_cache

  return(norm_K12)
}

.computeNormCorrCore <- function(object, tol = 1e-4, cts, scalePCs, sigmaValues, nCC) {
  
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

  # Check if there are at least 2 cell types for pairwise analysis
  if (length(cts) == 1) {
    pair_cell_types <- matrix(c(cts, cts), nrow = 2, ncol = 1)
  } else {
    pair_cell_types <- combn(cts, 2)
  }
  

  correlation_value <- vector("list", length = length(sigmaValues))
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  names(correlation_value) <- sigma_names

  norm_K12 <- .computeCrossKernelNorm(object, tol = tol, cts = cts,
   scalePCs = scalePCs, sigmaValues = sigmaValues, nCC = nCC,
   pair_cell_types = pair_cell_types)
  attr(object, "kernelNormalizerCache") <-
    attr(norm_K12, "kernelNormalizerCache", exact = TRUE)

  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    correlation_value[[t]] <- data.frame(
      sigmaValues = sigmaValues[tt],
      cellType1 = rep(pair_cell_types[1, ], times = nCC),
      cellType2 = rep(pair_cell_types[2, ], times = nCC),
      CC_index = rep(x = 1:nCC, each = ncol(pair_cell_types)),
      normalizedCorrelation = numeric(length = ncol(pair_cell_types) * nCC),
      stringsAsFactors = FALSE
    )
    
    sigma_val <- as.numeric(gsub("sigma_", "", t))

    for (pp in seq_len(ncol(pair_cell_types))) {
      cellType1 <- pair_cell_types[1, pp]
      cellType2 <- pair_cell_types[2, pp]

      A <- PCmats[[cellType1]]
      B <- PCmats[[cellType2]]

      K <- getKernelMatrix(object, sigma = sigma_val,
                           cellType1 = cellType1, cellType2 = cellType2,
                           verbose = FALSE)
      norm_K12_sel <- norm_K12[[t]][[cellType1]][[cellType2]]
      if (is.null(norm_K12_sel) || !is.finite(norm_K12_sel)) {
        warning(paste("Normalizer unavailable for", cellType1, "-", cellType2, "at", t, "- skipping"))
        next
      }

      # Check that skrCCA results exist for this pair
      if (is.null(object@skrCCAOut[[t]][[cellType1]]) ||
          is.null(object@skrCCAOut[[t]][[cellType2]])) {
        warning(paste("skrCCA results missing for", cellType1, "-", cellType2, "at", t, "- skipping"))
        next
      }

      for (cc_index in seq_len(nCC)) {
        w_1 <- object@skrCCAOut[[t]][[cellType1]][, cc_index, drop = FALSE]
        w_2 <- object@skrCCAOut[[t]][[cellType2]][, cc_index, drop = FALSE]

        A_w1 <- A %*% w_1
        B_w2 <- B %*% w_2

        numerator <- (t(A_w1) %*% K %*% B_w2)[1, 1]
        denominator <- sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K12_sel

        correlation_value[[t]]$"normalizedCorrelation"[
          pp + (cc_index - 1) * ncol(pair_cell_types)] <-
          ifelse(abs(denominator) < 1e-9, 0, numerator / denominator)
      }
    }
  }

  ## Store the result in the object
  object@normalizedCorrelation <- correlation_value

  ## obtain the sigma value with the highest
  ## normalized correlation
  ncorr <- do.call(rbind, correlation_value)
  ncorr$ct12 <- paste(ncorr$cellType1, ncorr$cellType2, sep = "-")

  # Calculate the mean of column 2 for each unique value in column 1
  ## only for cc_index == 1
  meanCorr <- tapply(
    ncorr$"normalizedCorrelation"[ncorr$"CC_index" == 1],
    ncorr$"sigmaValues"[ncorr$"CC_index" == 1], mean, na.rm = TRUE
  )

  # Find the value of column 1 with the highest mean in column 2
  if (length(meanCorr) > 0 && any(!is.na(meanCorr))) {
    sigmaValueChoice <- as.numeric(names(which.max(meanCorr)))
    object@sigmaValueChoice <- sigmaValueChoice
  } else {
    warning("Could not determine optimal sigma: All correlations were NA.")
    object@sigmaValueChoice <- numeric()
  }

  ## Return the modified object
  return(object)
}

#' @rdname computeNormalizedCorrelation
#' @aliases computeNormalizedCorrelation,CoPro-method
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @export
setMethod(
  "computeNormalizedCorrelation", "CoPro",
  function(object, tol = 1e-4) {
    input_check <- .checkInputNormCorr(object)
    cts <- input_check$cts
    scalePCs <- input_check$scalePCs
    sigmaValues <- input_check$sigmaValues
    nCC <- input_check$nCC

    object <- .computeNormCorrCore(object, tol = tol, cts = cts, 
                                   scalePCs = scalePCs, sigmaValues = sigmaValues, nCC = nCC)
    return(object)
  }
)

.checkInputNormCorrMulti <- function(object) {
  
  # --- Input Checks ---
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
  if(is.null(nCC) || length(nCC)==0) {
    # Try to infer nCC from the first available result
    first_sigma <- sigmas_run[1]
    first_ct <- cts[1]
    if (!is.null(object@skrCCAOut[[first_sigma]][[first_ct]])) {
      nCC <- ncol(object@skrCCAOut[[first_sigma]][[first_ct]])
    } else {
      stop("Cannot infer nCC from skrCCA results")
    }
  }

  if (length(object@scalePCs) == 0) stop("object@scalePCs not specified")
  scalePCs <- object@scalePCs

  return(list(cts = cts, slides = slides, sigmas_run = sigmas_run, nCC = nCC, scalePCs = scalePCs))
}

.computeCrossKernelNormMulti <- function(object, tol = 1e-4, cts, slides, sigmas_run, nCC, pair_cell_types) {

  # --- Precompute whitened-Frobenius normalizers (Per Slide, Per Sigma) ---
  message("Calculating whitened-Frobenius normalizers (can take time)...")
  norm_K_all <- setNames(vector("list", length = length(sigmas_run)), sigmas_run)
  normalizer_cache <- attr(object, "kernelNormalizerCache", exact = TRUE)
  if (is.null(normalizer_cache)) normalizer_cache <- list()

  for (sig_name in sigmas_run) {
    sigma_val <- as.numeric(gsub("sigma_", "", sig_name))
    norm_K_sigma <- setNames(vector("list", length = length(slides)), slides)

    for (sID in slides) {
      # Create empty (ct_i, ct_j) list structure for this slide
      norm_K_slide <- setNames(vector("list", length = length(cts)), cts)
      for (ct_i in cts) {
        norm_K_slide[[ct_i]] <- setNames(vector("list", length = length(cts)), cts)
      }

      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]

        # Retrieve kernel matrix via accessor (works with flat storage)
        K <- tryCatch({
          getKernelMatrix(object, sigma = sigma_val, cellType1 = ct_i, cellType2 = ct_j,
                          slide = sID, verbose = FALSE)
        }, error = function(e) NULL)

        if (!is.null(K) && nrow(K) > 0 && ncol(K) > 0) {
          # matched-sigma within-type kernels for this slide are the whitening
          # operators; fall back to ||K_c||_F if either is unavailable
          Rx <- tryCatch(getKernelMatrix(object, sigma = sigma_val, cellType1 = ct_i,
                          cellType2 = ct_i, slide = sID, verbose = FALSE),
                         error = function(e) NULL)
          Ry <- tryCatch(getKernelMatrix(object, sigma = sigma_val, cellType1 = ct_j,
                          cellType2 = ct_j, slide = sID, verbose = FALSE),
                         error = function(e) NULL)
          cache_key <- .kernelNormalizerKey(
            sigma_val, ct_i, ct_j, slide = sID
          )
          nrm <- .readKernelNormalizer(
            normalizer_cache, cache_key, K, Rx, Ry
          )
          if (is.null(nrm)) nrm <- tryCatch(.whitenedFrobNorm(K, Rx, Ry), error = function(e) {
            warning(paste("Normalizer failed for K[", sID, ",", ct_i, ",", ct_j, "]:", e$message))
            NA
          })
          if (is.finite(nrm)) {
            normalizer_cache <- .cacheKernelNormalizer(
              normalizer_cache, cache_key, K, Rx, Ry, nrm
            )
          }
          norm_K_slide[[ct_i]][[ct_j]] <- nrm
          norm_K_slide[[ct_j]][[ct_i]] <- nrm # symmetry (Frobenius transpose-invariant)
        } else {
          norm_K_slide[[ct_i]][[ct_j]] <- NA
          norm_K_slide[[ct_j]][[ct_i]] <- NA
        }
      }
      norm_K_sigma[[sID]] <- norm_K_slide
    }
    norm_K_all[[sig_name]] <- norm_K_sigma
  }

  message("Finished calculating whitened-Frobenius normalizers.")
  attr(norm_K_all, "kernelNormalizerCache") <- normalizer_cache
  return(norm_K_all)
}

.computeNormCorrCoreMulti <- function(object, tol = 1e-4, cts, slides, sigmas_run, nCC, scalePCs, calculationMode = "perSlide") {

  if (length(cts) == 1) {
    pair_cell_types <- matrix(c(cts, cts), nrow = 2, ncol = 1)
  } else {
    pair_cell_types <- combn(cts, 2)
  }

  # Scale per-slide PCA matrices to match optimization (whitening)
  X_scaled <- .preparePCMatrices(
    pc_data = object@pcaResults,
    pca_global = object@pcaGlobal,
    scalePCs = scalePCs,
    slides = slides,
    cts = cts
  )

  norm_K_all <- .computeCrossKernelNormMulti(object, tol = tol, cts = cts,
   slides = slides, sigmas_run = sigmas_run, nCC = nCC, pair_cell_types = pair_cell_types)
  attr(object, "kernelNormalizerCache") <-
    attr(norm_K_all, "kernelNormalizerCache", exact = TRUE)
  
  # --- Calculate Correlation ---
  correlation_results <- setNames(vector("list", length = length(sigmas_run)), sigmas_run)

  for (sig_name in sigmas_run) {
    W_list_sigma <- object@skrCCAOut[[sig_name]] # Shared weights
    sigma_val <- as.numeric(gsub("sigma_", "", sig_name))

    if (calculationMode == "perSlide") {
      correlation_per_slide <- setNames(vector("list", length=length(slides)), slides)
      for(sID in slides) {
        row_buffer <- vector("list", ncol(pair_cell_types) * nCC)
        row_idx <- 1
        X_list_slide <- X_scaled[[sID]]
        norm_K_slide <- norm_K_all[[sig_name]][[sID]]

        if (is.null(X_list_slide)) {
          correlation_per_slide[[sID]] <- data.frame(
            sigmaValue = numeric(),
            slideID = character(),
            cellType1 = character(), cellType2 = character(),
            CC_index = integer(), normalizedCorrelation = numeric(),
            stringsAsFactors = FALSE
          )
          next  # Skip if PCA data missing
        }

        for (pp in seq_len(ncol(pair_cell_types))) {
          ct_i <- pair_cell_types[1, pp]
          ct_j <- pair_cell_types[2, pp]

          X_i <- X_list_slide[[ct_i]]
          X_j <- X_list_slide[[ct_j]]
          K_ij <- tryCatch({
            getKernelMatrix(object,
                           sigma = sigma_val,
                           cellType1 = ct_i,
                           cellType2 = ct_j,
                           slide = sID,
                           verbose = FALSE)
          }, error = function(e) NULL)
          norm_K_ij <- norm_K_slide[[ct_i]][[ct_j]]

          if(is.null(X_i) || is.null(X_j) || is.null(K_ij) ||
             is.na(norm_K_ij) || norm_K_ij < 1e-9 ||
             nrow(X_i)==0 || nrow(X_j)==0) next # Skip if data missing/invalid

          for (cc in 1:nCC) {
            w_i <- W_list_sigma[[ct_i]][, cc, drop = FALSE]
            w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]

            Xiw <- X_i %*% w_i
            Xjw <- X_j %*% w_j

            numerator <- (t(Xiw) %*% K_ij %*% Xjw)[1, 1]
            denom_norm <- sqrt(sum(Xiw^2)) * sqrt(sum(Xjw^2)) * norm_K_ij

            norm_corr_val <- ifelse(abs(denom_norm) < 1e-9, 0, numerator / denom_norm)
            row_buffer[[row_idx]] <- list(
              sigmaValue = sigma_val,
              slideID = sID,
              cellType1 = ct_i,
              cellType2 = ct_j,
              CC_index = cc,
              normalizedCorrelation = as.numeric(norm_corr_val)
            )
            row_idx <- row_idx + 1
          } # end CC loop
        } # end pair loop
        if (row_idx > 1) {
          df_slide <- do.call(rbind.data.frame, c(row_buffer[seq_len(row_idx - 1)], stringsAsFactors = FALSE))
        } else {
          df_slide <- data.frame(
            sigmaValue = numeric(),
            slideID = character(),
            cellType1 = character(), cellType2 = character(),
            CC_index = integer(), normalizedCorrelation = numeric(),
            stringsAsFactors = FALSE
          )
        }
        correlation_per_slide[[sID]] <- df_slide
      } # end slide loop
      correlation_results[[sig_name]] <- do.call(rbind, correlation_per_slide)

    } else { # calculationMode == "aggregate"
      # For aggregate mode, calculate correlations across all slides for each cell type pair
      row_buffer <- vector("list", ncol(pair_cell_types) * nCC)
      row_idx <- 1
      
      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]

        valid_slide_data <- vector("list", length(slides))
        valid_slide_idx <- 1
        for(sID in slides) {
          X_list_slide <- X_scaled[[sID]]
          norm_K_slide <- norm_K_all[[sig_name]][[sID]]
          
          if (is.null(X_list_slide)) next
          
          X_i <- X_list_slide[[ct_i]]
          X_j <- X_list_slide[[ct_j]]
          K_ij <- tryCatch({
            getKernelMatrix(object,
                           sigma = sigma_val,
                           cellType1 = ct_i,
                           cellType2 = ct_j,
                           slide = sID,
                           verbose = FALSE)
          }, error = function(e) NULL)
          norm_K_ij <- norm_K_slide[[ct_i]][[ct_j]]
          
          if(is.null(X_i) || is.null(X_j) || is.null(K_ij) ||
             is.na(norm_K_ij) || norm_K_ij < 1e-9 ||
             nrow(X_i)==0 || nrow(X_j)==0) next

          valid_slide_data[[valid_slide_idx]] <- list(
            X_i = X_i,
            X_j = X_j,
            K_ij = K_ij,
            norm_K_ij = norm_K_ij
          )
          valid_slide_idx <- valid_slide_idx + 1
        }

        valid_slides_count <- valid_slide_idx - 1
        if (valid_slides_count == 0) next
        
        for (cc in 1:nCC) {
          # Aggregate numerator and denominator components across slides
          total_numerator <- 0
          total_norm_sum_i <- 0
          total_norm_sum_j <- 0
          total_K_norm <- 0
          w_i <- W_list_sigma[[ct_i]][, cc, drop = FALSE]
          w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]

          for(slide_idx in seq_len(valid_slides_count)) {
            slide_data <- valid_slide_data[[slide_idx]]
            Xiw <- slide_data$X_i %*% w_i
            Xjw <- slide_data$X_j %*% w_j
            
            total_numerator <- total_numerator + (t(Xiw) %*% slide_data$K_ij %*% Xjw)[1, 1]
            total_norm_sum_i <- total_norm_sum_i + sum(Xiw^2)
            total_norm_sum_j <- total_norm_sum_j + sum(Xjw^2)
            total_K_norm <- total_K_norm + slide_data$norm_K_ij
          }
          
          if (valid_slides_count > 0) {
            avg_K_norm <- total_K_norm / valid_slides_count
            denom_norm <- sqrt(total_norm_sum_i) * sqrt(total_norm_sum_j) * avg_K_norm
            
            agg_corr_val <- ifelse(abs(denom_norm) < 1e-9, 0, total_numerator / denom_norm)
            row_buffer[[row_idx]] <- list(
              sigmaValue = sigma_val,
              cellType1 = ct_i,
              cellType2 = ct_j,
              CC_index = cc,
              aggregateCorrelation = as.numeric(agg_corr_val)
            )
            row_idx <- row_idx + 1
          }
        }
      }
      if (row_idx > 1) {
        df_agg <- do.call(rbind.data.frame, c(row_buffer[seq_len(row_idx - 1)], stringsAsFactors = FALSE))
      } else {
        df_agg <- data.frame(
          sigmaValue = numeric(),
          cellType1 = character(), cellType2 = character(),
          CC_index = integer(), aggregateCorrelation = numeric(),
          stringsAsFactors = FALSE
        )
      }
      correlation_results[[sig_name]] <- df_agg
    }
  } # End sigma loop

  object@normalizedCorrelation <- correlation_results

  # --- Choose Optimal Sigma (based on aggregate or mean per-slide CC1 correlation) ---
  if(length(correlation_results) > 0) {
    first_sigma_res <- correlation_results[[1]]
    if(!is.null(first_sigma_res) && nrow(first_sigma_res) > 0) {
      corr_col_name <- ifelse("aggregateCorrelation" %in% names(first_sigma_res), "aggregateCorrelation", "normalizedCorrelation")
      all_corrs <- do.call(rbind, correlation_results)
      cc1_corrs <- all_corrs[all_corrs$CC_index == 1, ]

      if(nrow(cc1_corrs) > 0) {
        # Average correlation across pairs and potentially slides for CC1 for each sigma
        mean_corr_per_sigma <- tapply(cc1_corrs[[corr_col_name]], cc1_corrs$sigmaValue, mean, na.rm=TRUE)
        if(any(!is.na(mean_corr_per_sigma))){
          object@sigmaValueChoice <- as.numeric(names(which.max(mean_corr_per_sigma)))
        } else {
          warning("Could not determine optimal sigma: All correlations were NA.")
          object@sigmaValueChoice <- numeric() # Empty
        }

      } else {
        warning("Could not determine optimal sigma: No CC1 correlation values found.")
        object@sigmaValueChoice <- numeric() # Empty
      }
    } else {
      warning("Could not determine optimal sigma: Correlation results are empty.")
      object@sigmaValueChoice <- numeric() # Empty
    }
  } else {
    object@sigmaValueChoice <- numeric() # Empty if no results
  }

  return(object)
}

#' @rdname computeNormalizedCorrelation
#' @aliases computeNormalizedCorrelation,CoProMulti-method
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @export
setMethod("computeNormalizedCorrelation", "CoProMulti", function(
    object, tol = 1e-4, calculationMode = "perSlide") {

  # Validate calculationMode
  if (!calculationMode %in% c("perSlide", "aggregate")) {
    stop("calculationMode must be either 'perSlide' or 'aggregate'")
  }

  input_check <- .checkInputNormCorrMulti(object)
  cts <- input_check$cts
  slides <- input_check$slides
  sigmas_run <- input_check$sigmas_run
  nCC <- input_check$nCC
  scalePCs <- input_check$scalePCs

  object <- .computeNormCorrCoreMulti(object, tol = tol, cts = cts, slides = slides,
                                      sigmas_run = sigmas_run, nCC = nCC, scalePCs = scalePCs,
                                      calculationMode = calculationMode)
  return(object)
})
