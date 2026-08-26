# The \describe{} items below continue at a two-space indent, not four. A
# paragraph break followed by a four-space indent reads as a markdown code
# block, so roxygen wraps the rest of the item in \preformatted{} and escapes
# the brace that was meant to close \item{}{}. The .Rd then has unbalanced
# braces and R CMD check reports a WARNING from both "Rd files" and "can be
# installed".

#' Compute Normalized Correlation (approximation)
#'
#' This method calculates the normalized correlation between pairs of cell types
#' based on CCA weights and the respective kernel matrix. It divides the
#' bilinear statistic \eqn{T = a' K b} by the whitened-Frobenius norm
#' \eqn{\|R_x^{1/2} K_c R_y^{1/2}\|_F}, which is the null standard deviation of
#' \eqn{T} when the score vectors carry within-type covariance proportional to
#' \eqn{R_x} and \eqn{R_y}.
#'
#' @section Choosing the normalizer:
#' The whitening operators decide how the criterion behaves across the sigma
#' grid, and therefore which bandwidth `sigmaValueChoice` ends up at.
#'
#' \describe{
#'   \item{`"legacy"` (default)}{Use the historical selection rule: use the
#'     matched-sigma within-type kernels when the
#'     object happens to contain them, and \eqn{R = I} otherwise. Because
#'     `computeKernelMatrix()` builds only cross-type kernels, this is normally
#'     \eqn{\|K_c\|_F}; but it silently becomes the whitened norm on an object
#'     that has been through `computeSelfKernel()`. The whitening copy always
#'     has its unit diagonal restored so that it is a valid correlation
#'     operator. Which path applied is reported by [getNormalizerInfo()].}
#'   \item{`"unwhitened"`}{Force \eqn{R_x = R_y = I}, i.e. \eqn{\|K_c\|_F}. This
#'     assumes spatially independent scores, so it under-counts noise at large
#'     sigma and biases the selected bandwidth upward.}
#'   \item{`"kernel"`}{Force \eqn{R = K(\sigma)}, the matched-sigma self-kernel,
#'     with the unit diagonal restored. Errors if the self-kernels are absent
#'     rather than falling back. This over-counts noise at large sigma and
#'     biases the selected bandwidth downward. Provided as a diagnostic, not
#'     for analysis.}
#'   \item{`"variogram"`}{Estimate one within-type autocorrelation range per
#'   cell type from the feature-averaged spatial autocorrelation of the PC
#'   scores, and use \eqn{R = \exp(-d^2/2\ell^2)}. The range is a property of
#'   the score field rather than of the bandwidth, so it is fitted once and
#'   reused across the grid. It is fitted from all PC columns, never from a
#'   fitted canonical score, which would leak the association under test into
#'   the denominator.
#'
#'   **This flattens the null when the scores share one correlation length,
#'   and should be treated as a diagnostic otherwise.** On a real targeted
#'   panel the principal components within a cell type generally do *not*
#'   share one: measured per-component ranges spanned 0.3-4.1 cell-spacings on
#'   a colon MERFISH dataset, where the flat feature-average returned a
#'   sub-spacing range (so \eqn{R \approx I}, reproducing `"unwhitened"`)
#'   while pooling over the leading five components gave a range 6x longer.
#'   The selected bandwidth is steeply sensitive to that choice -- sweeping
#'   the assumed range over 0.4-5 cell-spacings moved it from 3.4 to 0.3
#'   spacings on the same data. Check `getNormalizerInfo()$ranges` against
#'   what you believe about the tissue, supply `normalizerControl$range`
#'   if you have a better estimate, and prefer a permutation null when the
#'   selected bandwidth matters.}
#' }
#'
#' @section Multi-slide aggregation:
#' In `calculationMode = "aggregate"`, let \eqn{T_s} be the numerator and
#' \eqn{d_s} its per-slide denominator. `aggregateDenominator = "sum"` returns
#' \eqn{\sum_s T_s / \sum_s d_s}, the \eqn{d_s}-weighted mean of the per-slide
#' normalized correlations. `aggregateDenominator = "rss"` instead returns
#' \eqn{\sum_s T_s / \sqrt{\sum_s d_s^2}}, the null-standardized summed
#' statistic when slides are independent. The latter is not an aggregate
#' correlation and can grow with the number of concordant slides.
#'
#' The numerator is not centered again here. Its correlation interpretation
#' relies on the PC scores already being centered; objects fitted with
#' `computePCA(center = FALSE)` retain score-mean coupling in the numerator.
#' Sigma selection averages the available pair/slide rows with `na.rm = TRUE`,
#' so users should inspect missingness when valid-pair coverage varies by sigma.
#'
#' @param object A `CoPro` or `CoProMulti` object containing CCA results and kernel matrices.
#' @param tol tolerance for approximate SVD calculation
#' @param calculationMode (for CoProMulti only) either "perSlide" or "aggregate",
#'   for single slide analysis, it is ignored, with default value "perSlide".
#' @param normalizer Which whitening operators to use in the denominator; one of
#'   `"legacy"`, `"unwhitened"`, `"kernel"`, `"variogram"`. See details.
#' @param normalizerControl A named list tuning the `"variogram"` normalizer.
#'   `distType`, `xDistScale`, `yDistScale`, and `zDistScale` default to `NULL`
#'   and inherit the geometry recorded by the distance/kernel builder. Explicit
#'   values must match that record.
#'   `range` accepts a named numeric vector of per-cell-type ranges, in
#'   normalized distance units, to skip estimation. Remaining entries
#'   (`maxCells`, `nBins`, `maxLagQuantile`, `minCorrelation`, `minBins`,
#'   `lowerLimit`) tune the fit and the operator's truncation.
#' @param aggregateDenominator For multi-slide aggregate mode, either `"sum"`
#'   (default; a denominator-weighted aggregate correlation) or `"rss"` (a
#'   null-standardized summed statistic under independent slides). Ignored for
#'   single-slide and per-slide calculations.
#' @return The object with the normalized correlation value
#' between any pair of cell types
#' added as a new slot, `normalizedCorrelation`. The resolved normalizer is
#' attached as an attribute and can be read back with [getNormalizerInfo()]. In
#' aggregate mode the selected `aggregateDenominator` is also attached to that
#' result list.
#' @family scores-and-correlation
#' @seealso [runSkrCCA()], [computeBidirCorrelation()],
#'   [computeGeneAndCellScores()], [getNormalizerInfo()]
#' @export
#'
setGeneric(
  "computeNormalizedCorrelation",
  function(object, tol = 1e-4, calculationMode = "perSlide",
           normalizer = c("legacy", "unwhitened", "kernel", "variogram"),
           normalizerControl = list(),
           aggregateDenominator = c("sum", "rss")) {
    standardGeneric("computeNormalizedCorrelation")
  }
)

#' Combine per-slide normalized-correlation components
#' @noRd
.aggregateNormCorr <- function(numerators, slideDenominators,
                               denominator = c("sum", "rss")) {
  denominator <- match.arg(denominator)
  if (length(numerators) != length(slideDenominators)) {
    stop("numerators and slideDenominators must have the same length")
  }
  if (length(numerators) == 0L) return(NA_real_)
  if (any(!is.finite(numerators)) || any(!is.finite(slideDenominators)) ||
      any(slideDenominators < 0)) {
    return(NA_real_)
  }
  aggregate_denominator <- switch(
    denominator,
    sum = sum(slideDenominators),
    rss = sqrt(sum(slideDenominators^2))
  )
  if (abs(aggregate_denominator) < 1e-9) return(0)
  sum(numerators) / aggregate_denominator
}


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
#' Sum of squares of a sparse matrix's represented entries
#'
#' `sum(K * K)` allocates a whole second sparse matrix -- 12 bytes per nonzero,
#' several GB on a large kernel -- to produce one scalar. Read the value slot
#' instead. A `dsCMatrix` stores one triangle, so each stored off-diagonal
#' entry represents two.
#' @noRd
.sparseSumSquares <- function(K) {
  squares <- sum(K@x^2)
  if (!inherits(K, "symmetricMatrix")) return(squares)
  diagonal <- Matrix::diag(K)
  2 * squares - sum(diagonal^2)
}

#' `<Rx K Ry, K>` without materializing `Rx K Ry`
#'
#' Both sparse products fill in heavily: on a 40k-cell kernel with ~30
#' neighbours per cell, `Rx %*% K` grows nnz 7x and `(Rx K) %*% Ry` 11x, which
#' at 200k cells and ~40 neighbours extrapolates to roughly 1.5 GB for a
#' quantity that reduces to a scalar.
#'
#' With `Rx` and `Ry` symmetric (the caller symmetrizes them),
#' `sum((Rx K Ry) * K) = tr(K' Rx K Ry) = sum((Rx K) * (K Ry))`, and that form
#' splits over column blocks of `K`: block `J` needs only `Rx K[, J]` and
#' `K Ry[, J]`. Peak memory becomes proportional to the block rather than to
#' the whole filled-in product.
#'
#' Blocking turns out to be faster as well as smaller, because the intermediates
#' stay in cache instead of streaming a multi-hundred-MB product through memory.
#' Measured on a 60k-cell kernel (nnz 1.56M, 11x fill unblocked):
#'
#' ```
#'   unblocked            2.43 s   821 MB
#'   block_nnz = 2e5      1.72 s   296 MB      8 blocks
#'   block_nnz = 5e4      1.58 s   299 MB     32 blocks
#' ```
#'
#' `block_nnz` is a budget on nonzeros of `K` per block, not on the product, so
#' with the ~7x fill of one operator it caps each intermediate near 1.4M
#' nonzeros (~17 MB) regardless of how large `K` is. Setting it high enough to
#' yield a single block is worse than not blocking at all -- both one-sided
#' products are then live at once -- so the default stays at the measured knee.
#'
#' @param K Cross-cell-type kernel.
#' @param Rx,Ry Symmetric within-type whitening operators.
#' @param block_nnz Target nonzeros of `K` per block.
#' @return Scalar `<Rx K Ry, K>`.
#' @noRd
.sparseWhitenedInner <- function(K, Rx, Ry, block_nnz = 2e5) {
  n_columns <- ncol(K)
  if (n_columns == 0L) return(0)
  per_column <- max(1, length(K@x) / n_columns)
  block <- max(1L, min(n_columns, as.integer(block_nnz / per_column)))

  total <- 0
  start <- 1L
  while (start <= n_columns) {
    columns <- seq.int(start, min(n_columns, start + block - 1L))
    total <- total + sum(
      (Rx %*% K[, columns, drop = FALSE]) *
        (K %*% Ry[, columns, drop = FALSE])
    )
    start <- start + block
  }
  as.numeric(total)
}

.whitenedFrobNorm <- function(K, Rx = NULL, Ry = NULL) {
  if (.isFloat32SparseKernel(K)) {
    if (!is.null(Rx) || !is.null(Ry)) {
      K_double <- asDoubleSparseMatrix(K)
      Rx_double <- if (identical(Rx, K)) {
        K_double
      } else if (.isFloat32SparseKernel(Rx)) {
        asDoubleSparseMatrix(Rx)
      } else {
        Rx
      }
      Ry_double <- if (identical(Ry, K)) {
        K_double
      } else if (identical(Ry, Rx)) {
        Rx_double
      } else if (.isFloat32SparseKernel(Ry)) {
        asDoubleSparseMatrix(Ry)
      } else {
        Ry
      }
      return(.whitenedFrobNorm(K_double, Rx_double, Ry_double))
    }
    nr <- nrow(K)
    nc <- ncol(K)
    sums <- .float32KernelSums(K)
    rmean <- sums$rowSums / nc
    cmean <- sums$colSums / nr
    grand_mean <- sum(sums$rowSums) / (as.numeric(nr) * nc)
    norm2 <- sums$sumSquares -
      nc * sum(rmean^2) -
      nr * sum(cmean^2) +
      as.numeric(nr) * nc * grand_mean^2
    return(sqrt(max(as.numeric(norm2), 0)))
  }

  ## K is double here, but a within-type whitening operator may still be an
  ## encoded float32 kernel -- e.g. a cross kernel built in double while the
  ## self-kernels supplied as Rx/Ry were built with computeSparseKernelFloat32,
  ## or vice versa on the same object. Decode any float32 operator to double so
  ## the formulas below apply; whitening needs only the values, so the temporary
  ## double copy is exact. (.isFloat32SparseKernel(NULL) is FALSE, so NULL Rx/Ry
  ## pass through untouched.)
  if (.isFloat32SparseKernel(Rx)) Rx <- asDoubleSparseMatrix(Rx)
  if (.isFloat32SparseKernel(Ry)) Ry <- asDoubleSparseMatrix(Ry)

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
      norm2 <- .sparseSumSquares(K) - nc * sum(rmean^2) - nr * sum(cmean^2) +
        nr * nc * grand_mean^2
      return(sqrt(max(as.numeric(norm2), 0)))
    }

    Rx <- (Rx + t(Rx)) / 2
    Ry <- (Ry + t(Ry)) / 2

    U <- cbind(-rmean, rep.int(1, nr))
    V <- cbind(rep.int(1, nc), grand_mean - cmean)
    ## <Rx K Ry, K>, streamed over column blocks; see .sparseWhitenedInner().
    base_term <- .sparseWhitenedInner(K, Rx, Ry)
    ## <Rx K Ry, U V'> = <U, (Rx K Ry) V>. V has two columns, so applying the
    ## three operators to it right-to-left never widens past an n x 2 dense
    ## block -- M itself is not needed here.
    cross_term <- as.numeric(sum(U * (Rx %*% (K %*% (Ry %*% V)))))
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

#' Compute an encoded float32 kernel's content signature from its values
#'
#' Walks every nonzero, so callers should prefer the copy
#' `.newFloat32SparseKernel()` stores at construction.
#' @noRd
.computeFloat32KernelSignature <- function(x) {
  sums <- .float32KernelSums(x)
  paste(
    nrow(x), ncol(x), .float32KernelNnz(x),
    format(sum(sums$rowSums), digits = 16),
    format(sums$sumSquares, digits = 16),
    sep = ":"
  )
}

.kernelNormalizerKey <- function(sigma, cellType1, cellType2, slide = NULL) {
  paste(
    .formatSigmaValue(sigma),
    if (is.null(slide)) "single" else slide,
    cellType1, cellType2, sep = "|"
  )
}

#' Content signature used to invalidate a cached normalizer
#'
#' `.readKernelNormalizer()` and `.cacheKernelNormalizer()` each signature `K`,
#' `Rx` and `Ry`, so validating one cache entry used to walk every nonzero six
#' times over. For an encoded float32 kernel the signature is a pure function
#' of the stored values, and those are immutable once built, so
#' `.newFloat32SparseKernel()` computes it once at construction and this reads
#' it back. The full scan remains for kernels without a stored signature:
#' double `sparseMatrix` kernels, and float32 kernels inside objects saved
#' before the field existed.
#' @noRd
.kernelMatrixSignature <- function(x) {
  if (is.null(x)) return("NULL")
  if (.isFloat32SparseKernel(x)) {
    stored <- x$signature
    if (!is.null(stored)) {
      # A transposed view represents a different matrix, and t() only flips the
      # flag, so distinguish the two orientations. Values, and therefore every
      # summary in the signature, are identical either way.
      return(if (isTRUE(x$transposed)) paste0(stored, ":t") else stored)
    }
    return(.computeFloat32KernelSignature(x))
  }
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
 sigmaValues, nCC, pair_cell_types, resolver) {
  ## Whitened-Frobenius normalizer ||R_x^{1/2} K_c R_y^{1/2}||_F per (sigma,
  ## pair). `resolver` supplies R_x and R_y; see .makeWhiteningResolver().
  message("Calculating whitened-Frobenius normalizers, this may take a while.")
  whitened_pairs <- 0L
  total_pairs <- 0L

  sigma_names <- .sigmaName(sigmaValues)
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
                           verbose = FALSE, materialize = FALSE)
      Rx <- resolver$get(sigma_val, cellType1)
      Ry <- resolver$get(sigma_val, cellType2)
      total_pairs <- total_pairs + 1L
      if (!is.null(Rx) && !is.null(Ry)) {
        whitened_pairs <- whitened_pairs + 1L
      }
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

  description <- .describeNormalizer(resolver, whitened_pairs, total_pairs)
  message("Normalizer: ", description)
  if (resolver$mode == "variogram" && length(resolver$ranges) > 0) {
    message("  fitted autocorrelation ranges: ",
            paste(names(resolver$ranges),
                  format(resolver$ranges, digits = 3),
                  sep = " = ", collapse = ", "))
  }
  message("Finished calculating whitened-Frobenius normalizers.")

  attr(norm_K12, "kernelNormalizerCache") <- normalizer_cache
  attr(norm_K12, "normalizerInfo") <- list(
    mode = resolver$mode, description = description, ranges = resolver$ranges
  )

  return(norm_K12)
}

.computeNormCorrCore <- function(object, tol = 1e-4, cts, scalePCs, sigmaValues,
                                 nCC, normalizer = "legacy",
                                 normalizerControl = list()) {

  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

  # Check if there are at least 2 cell types for pairwise analysis
  if (length(cts) == 1) {
    pair_cell_types <- matrix(c(cts, cts), nrow = 2, ncol = 1)
  } else {
    pair_cell_types <- combn(cts, 2)
  }
  

  correlation_value <- vector("list", length = length(sigmaValues))
  sigma_names <- .sigmaName(sigmaValues)
  names(correlation_value) <- sigma_names

  resolver <- .makeWhiteningResolver(
    object, normalizer, normalizerControl, scoreMats = PCmats, cts = cts
  )
  norm_K12 <- .computeCrossKernelNorm(object, tol = tol, cts = cts,
   scalePCs = scalePCs, sigmaValues = sigmaValues, nCC = nCC,
   pair_cell_types = pair_cell_types, resolver = resolver)
  attr(object, "kernelNormalizerCache") <-
    attr(norm_K12, "kernelNormalizerCache", exact = TRUE)
  normalizer_info <- attr(norm_K12, "normalizerInfo", exact = TRUE)

  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    correlation_value[[t]] <- data.frame(
      sigmaValues = sigmaValues[tt],
      cellType1 = rep(pair_cell_types[1, ], times = nCC),
      cellType2 = rep(pair_cell_types[2, ], times = nCC),
      CC_index = rep(x = 1:nCC, each = ncol(pair_cell_types)),
      normalizedCorrelation = rep(NA_real_, ncol(pair_cell_types) * nCC),
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
                           verbose = FALSE, materialize = FALSE)
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

        numerator <- .kernelXKY(A_w1, K, B_w2)[1, 1]
        denominator <- sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K12_sel

        correlation_value[[t]]$"normalizedCorrelation"[
          pp + (cc_index - 1) * ncol(pair_cell_types)] <-
          ifelse(abs(denominator) < 1e-9, 0, numerator / denominator)
      }
    }
  }

  ## Store the result in the object, tagged with the denominator that produced
  ## it -- the numbers are not comparable across normalizers.
  attr(correlation_value, "normalizer") <- normalizer_info
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
#' @export
setMethod(
  "computeNormalizedCorrelation", "CoPro",
  function(object, tol = 1e-4, calculationMode = "perSlide",
           normalizer = c("legacy", "unwhitened", "kernel", "variogram"),
           normalizerControl = list(),
           aggregateDenominator = c("sum", "rss")) {
    normalizer <- match.arg(normalizer)
    aggregateDenominator <- match.arg(aggregateDenominator)
    input_check <- .checkInputNormCorr(object)
    cts <- input_check$cts
    scalePCs <- input_check$scalePCs
    sigmaValues <- input_check$sigmaValues
    nCC <- input_check$nCC

    object <- .computeNormCorrCore(object, tol = tol, cts = cts,
                                   scalePCs = scalePCs, sigmaValues = sigmaValues, nCC = nCC,
                                   normalizer = normalizer,
                                   normalizerControl = normalizerControl)
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

.computeCrossKernelNormMulti <- function(object, tol = 1e-4, cts, slides, sigmas_run, nCC, pair_cell_types, resolver) {

  # --- Precompute whitened-Frobenius normalizers (Per Slide, Per Sigma) ---
  message("Calculating whitened-Frobenius normalizers (can take time)...")
  whitened_pairs <- 0L
  total_pairs <- 0L
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
                          slide = sID, verbose = FALSE,
                          materialize = FALSE)
        }, error = function(e) NULL)

        if (!is.null(K) && nrow(K) > 0 && ncol(K) > 0) {
          Rx <- resolver$get(sigma_val, ct_i, slide = sID)
          Ry <- resolver$get(sigma_val, ct_j, slide = sID)
          total_pairs <- total_pairs + 1L
          if (!is.null(Rx) && !is.null(Ry)) {
            whitened_pairs <- whitened_pairs + 1L
          }
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

  description <- .describeNormalizer(resolver, whitened_pairs, total_pairs)
  message("Normalizer: ", description)
  message("Finished calculating whitened-Frobenius normalizers.")
  attr(norm_K_all, "kernelNormalizerCache") <- normalizer_cache
  attr(norm_K_all, "normalizerInfo") <- list(
    mode = resolver$mode, description = description, ranges = resolver$ranges
  )
  return(norm_K_all)
}

.computeNormCorrCoreMulti <- function(object, tol = 1e-4, cts, slides, sigmas_run, nCC, scalePCs, calculationMode = "perSlide",
                                      normalizer = "legacy",
                                      normalizerControl = list(),
                                      aggregateDenominator = "sum") {

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

  resolver <- .makeWhiteningResolver(
    object, normalizer, normalizerControl, scoreMats = X_scaled, cts = cts,
    slides = slides
  )
  norm_K_all <- .computeCrossKernelNormMulti(object, tol = tol, cts = cts,
   slides = slides, sigmas_run = sigmas_run, nCC = nCC, pair_cell_types = pair_cell_types,
   resolver = resolver)
  attr(object, "kernelNormalizerCache") <-
    attr(norm_K_all, "kernelNormalizerCache", exact = TRUE)
  normalizer_info <- attr(norm_K_all, "normalizerInfo", exact = TRUE)

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
                           verbose = FALSE,
                           materialize = FALSE)
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

            numerator <- .kernelXKY(Xiw, K_ij, Xjw)[1, 1]
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
                           verbose = FALSE,
                           materialize = FALSE)
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
          slide_numerators <- numeric(valid_slides_count)
          slide_denominators <- numeric(valid_slides_count)
          w_i <- W_list_sigma[[ct_i]][, cc, drop = FALSE]
          w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]

          for(slide_idx in seq_len(valid_slides_count)) {
            slide_data <- valid_slide_data[[slide_idx]]
            Xiw <- slide_data$X_i %*% w_i
            Xjw <- slide_data$X_j %*% w_j
            
            slide_numerators[slide_idx] <-
              .kernelXKY(Xiw, slide_data$K_ij, Xjw)[1, 1]
            slide_denominators[slide_idx] <-
              sqrt(sum(Xiw^2)) * sqrt(sum(Xjw^2)) *
              slide_data$norm_K_ij
          }
          
          if (valid_slides_count > 0) {
            agg_corr_val <- .aggregateNormCorr(
              slide_numerators, slide_denominators,
              denominator = aggregateDenominator
            )
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

  attr(correlation_results, "normalizer") <- normalizer_info
  if (calculationMode == "aggregate") {
    attr(correlation_results, "aggregateDenominator") <- aggregateDenominator
  }
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
#' @export
setMethod("computeNormalizedCorrelation", "CoProMulti", function(
    object, tol = 1e-4, calculationMode = "perSlide",
    normalizer = c("legacy", "unwhitened", "kernel", "variogram"),
    normalizerControl = list(),
    aggregateDenominator = c("sum", "rss")) {

  normalizer <- match.arg(normalizer)
  aggregateDenominator <- match.arg(aggregateDenominator)

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
                                      calculationMode = calculationMode,
                                      normalizer = normalizer,
                                      normalizerControl = normalizerControl,
                                      aggregateDenominator = aggregateDenominator)
  return(object)
})
