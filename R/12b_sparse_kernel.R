# =============================================================================
# Sparse, exact Gaussian-kernel construction via fixed-radius neighbor search
# -----------------------------------------------------------------------------
# A fused replacement for computeDistance() + computeKernelMatrix() that goes
# directly from spatial coordinates to a sparse `dgCMatrix` kernel, never
# materializing a dense n x n distance or kernel matrix.
#
# Exactness: the Gaussian kernel exp(-0.5 (d/sigma)^2) is below `lowerLimit`
# precisely when d > sigma * sqrt(-2 log(lowerLimit)). The dense path already
# zeroes those entries, so enumerating only pairs within that radius (via the
# grid search in 11b_sparse_neighbors.R) and applying the identical transforms
# reproduces the dense result up to floating-point rounding, while cost and
# memory scale with the number of near pairs instead of n^2.
#
# Every transform in the dense pipeline is mirrored here, in the same order:
#   coordinates -> distance -> zero handling -> low-percentile flooring
#   -> global distance scaling -> Gaussian kernel -> upper-quantile clip
#   -> optional normalization -> lower-limit drop -> sigma sparsity check.
# =============================================================================

#' Radius beyond which a Gaussian kernel falls below `lowerLimit`
#' @noRd
.kernelSupportMultiplier <- function(lowerLimit) {
  sqrt(-2 * log(lowerLimit))
}

#' Largest number of cells in any single cell type of interest (subset).
#' Used by method = "auto" to decide dense vs sparse.
#' @noRd
.maxCellTypeCount <- function(object) {
  cts <- object@cellTypesOfInterest
  sub <- object@cellTypesSub
  if (length(cts) == 0 || length(sub) == 0) return(0L)
  tt <- table(sub[sub %in% cts])
  if (length(tt) == 0) return(0L)
  as.integer(max(tt))
}

#' Dispatch computeKernelMatrix() to the dense or sparse path and optionally
#' drop the @distances slot afterward.
#' @noRd
.computeKernelDispatch <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                   normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                   colNormalizeKernel, verbose, method, dropDistances,
                                   autoThreshold, distType, xDistScale, yDistScale, zDistScale,
                                   normalizeDistance, normalizeTarget, truncateLowDist,
                                   is_multi) {
  # Infer distType from coordinates when not supplied (sparse path only uses it).
  if (is.null(distType)) {
    distType <- if ("z" %in% tolower(colnames(object@locationDataSub))) {
      "Euclidean3D"
    } else {
      "Euclidean2D"
    }
  }

  if (method == "auto") {
    n_max <- .maxCellTypeCount(object)
    method <- if (n_max > autoThreshold) "sparse" else "dense"
    if (verbose) {
      message(sprintf(
        "computeKernelMatrix: method='auto' -> '%s' (largest cell type = %d cells, threshold = %d).",
        method, n_max, autoThreshold))
    }
  }

  if (method == "dense") {
    if (is_multi) {
      object <- .computeKernelCoreMulti(object, sigmaValues, lowerLimit, upperQuantile,
                                        normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                        colNormalizeKernel, verbose)
    } else {
      object <- .computeKernelCore(object, sigmaValues, lowerLimit, upperQuantile,
                                   normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                   colNormalizeKernel, verbose)
    }
  } else {  # sparse
    if (is_multi) {
      object <- .computeSparseKernelCoreMulti(object, sigmaValues, lowerLimit, upperQuantile,
                                              normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                              colNormalizeKernel, verbose,
                                              distType, xDistScale, yDistScale, zDistScale,
                                              normalizeDistance, normalizeTarget, truncateLowDist)
    } else {
      object <- .computeSparseKernelCore(object, sigmaValues, lowerLimit, upperQuantile,
                                         normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                         colNormalizeKernel, verbose,
                                         distType, xDistScale, yDistScale, zDistScale,
                                         normalizeDistance, normalizeTarget, truncateLowDist)
    }
  }

  if (dropDistances && length(object@distances) > 0) {
    if (verbose) message("Cleared @distances after kernel computation (dropDistances = TRUE).")
    object@distances <- list()
  }
  object
}

#' Build, for one cell-type block, the cached neighbor triplets with distances
#' already zero-handled, floored, and globally scaled (kernel-ready for any
#' sigma). Returns NULL if no near pairs exist.
#'
#' @param A,B coordinate matrices (B = NULL for the within-type case).
#' @param percentile this block's low-distance percentile (original units); used
#'   for `truncateLowDist` flooring. May be NA when not needed.
#' @param scaling_factor global distance scaling factor (1 when normalizeDistance
#'   is FALSE).
#' @param max_sigma largest sigma value (sets the search radius).
#' @param lowerLimit kernel floor.
#' @param truncateLowDist whether to floor small distances.
#' @return list(i, j, dscaled, n_i, n_j) or NULL.
#' @noRd
.buildBlockTriplets <- function(A, B, percentile, scaling_factor, max_sigma,
                                lowerLimit, truncateLowDist) {
  within <- is.null(B)
  n_i <- nrow(A)
  n_j <- if (within) n_i else nrow(B)

  mult <- .kernelSupportMultiplier(lowerLimit)
  # search radius in ORIGINAL units; slightly inflated so floating-point
  # boundary pairs are captured and then settled by the exact >= lowerLimit
  # filter when kernelizing.
  r_orig <- mult * max_sigma / scaling_factor * (1 + 1e-6)

  tri <- .frnnGrid(A, B, r_orig)
  if (length(tri$i) == 0) return(NULL)

  d <- tri$d
  # zero distances -> smallest non-zero distance (mirrors .processDistanceMatrix)
  if (any(d == 0)) {
    nz <- d[d > 0]
    if (length(nz) > 0) {
      warning(paste("Zero distances detected, replacing with",
                    "the smallest non-zero distances, please",
                    "consider checking the location of cells",
                    "for potential errors"))
      d[d == 0] <- min(nz)
    }
  }
  # floor small distances
  if (truncateLowDist && !is.na(percentile)) {
    d[d < percentile] <- percentile
  }
  d_scaled <- d * scaling_factor

  list(i = tri$i, j = tri$j, dscaled = d_scaled, n_i = n_i, n_j = n_j)
}

#' Assemble the raw sparse Gaussian kernel for one block at one sigma.
#' @noRd
.sparseKernelFromTriplets <- function(bt, sigma, lowerLimit) {
  k <- exp(-0.5 * (bt$dscaled / sigma)^2)
  keep <- k >= lowerLimit
  Matrix::sparseMatrix(
    i = bt$i[keep], j = bt$j[keep], x = k[keep],
    dims = c(bt$n_i, bt$n_j)
  )
}

#' Sparse analogue of .CheckSigmaValuesToRemove (single-slide, proportion-based)
#' @noRd
.checkSparseSigmaRemove <- function(K, lowerLimit, sigma_choose, sigmaValues,
                                    i, j, minAveCellNeighor) {
  n1 <- nrow(K); n2 <- ncol(K)
  minPropZero <- minAveCellNeighor * min(n1, n2) / (n1 * n2)
  prop_above <- sum(K@x > lowerLimit) / (n1 * n2)
  if (prop_above < minPropZero) {
    warning(paste("Kernel matrix for cell types", i, "and", j,
                  "with sigma =", sigma_choose,
                  "contains too many zeros. Specifically, less than",
                  minPropZero * 100, "% total counts are above the threshold"))
    if (length(sigmaValues) == 1) {
      stop(paste("Only one sigma value is specified,",
                 "which resulted in all Gaussian kernel being small.",
                 "Please provide a larger sigma value"))
    } else {
      warning(paste("Dropping sigma value of ", sigma_choose,
                    "because all Gaussian kernel values are too small,",
                    "which will not produce meaningful results."))
      return(TRUE)
    }
  }
  FALSE
}

#' Sparse analogue of .checkKernelValidityMulti (count-based)
#' @noRd
.checkSparseKernelValidityMulti <- function(K, lowerLimit, minAveCellNeighor,
                                            sigma_val, ct_i, ct_j, sID) {
  n1 <- nrow(K); n2 <- ncol(K)
  if (sum(K@x > lowerLimit) < minAveCellNeighor * min(n1, n2)) {
    warning(paste("Kernel matrix for", ct_i, "-", ct_j, "in slide", sID,
                  "with sigma =", sigma_val, "is too sparse."))
    return(FALSE)
  }
  TRUE
}

#' Sparse analogue of .processKernelMatrix: upper-quantile clip, optional
#' normalization, lower-limit drop. Operates on the @x slot so it never densifies.
#' @noRd
.processSparseKernelMatrix <- function(K, lowerLimit, upperQuantile,
                                       normalizeKernel, rowNormalizeKernel,
                                       colNormalizeKernel) {
  if (nrow(K) == 0 || ncol(K) == 0) {
    stop("Cannot process empty kernel matrix")
  }
  K <- as(K, "CsparseMatrix")
  valid <- K@x  # all stored values are >= lowerLimit, no NA
  if (length(valid) == 0) {
    warning("No valid kernel values found above lowerLimit")
    return(K)
  }

  # Clip large values (matches dense: clip computed from valid values)
  upper_clip <- as.numeric(stats::quantile(valid, upperQuantile, na.rm = TRUE))
  K@x[K@x >= upper_clip] <- upper_clip

  if (normalizeKernel && !rowNormalizeKernel && !colNormalizeKernel) {
    rs_kernel <- Matrix::rowSums(K)
    median_rs <- stats::median(rs_kernel[rs_kernel > 1e-5], na.rm = TRUE)
    if (!is.na(median_rs) && median_rs > 0) {
      K@x <- K@x / median_rs
    }
  } else if (rowNormalizeKernel) {
    rs_kernel <- Matrix::rowSums(K)
    scl <- ifelse(rs_kernel > 1e-4, 1 / rs_kernel, 1)
    K <- Matrix::Diagonal(x = scl) %*% K
  } else if (colNormalizeKernel) {
    cs_kernel <- Matrix::colSums(K)
    scl <- ifelse(cs_kernel > 1e-4, 1 / cs_kernel, 1)
    K <- K %*% Matrix::Diagonal(x = scl)
  }

  # Remove small values (mirrors final dense lower-limit zeroing)
  K <- as(K, "CsparseMatrix")
  K@x[K@x < lowerLimit] <- 0
  Matrix::drop0(K)
}

#' Validate inputs for the sparse kernel path (does NOT require @distances).
#' @noRd
.checkInputSparseKernel <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                    minAveCellNeighor, rowNormalizeKernel,
                                    colNormalizeKernel, distType) {
  if (rowNormalizeKernel && colNormalizeKernel) {
    stop("Cannot do both row-wise and column-wise normalization.")
  }
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) stop("No cell types of interest specified")

  if (distType == "Morphology-Aware") {
    stop(paste("method = 'sparse' supports 'Euclidean2D' / 'Euclidean3D' only.",
               "Use method = 'dense' for Morphology-Aware distances."))
  }
  if (length(sigmaValues) == 0) {
    stop(paste("sigmaValues must be provided for the sparse kernel method",
               "(no distance matrix is built to derive a default from)."))
  }
  if (!is.numeric(sigmaValues) || any(sigmaValues <= 0)) {
    stop("sigmaValues must be positive numeric values")
  }
  if (lowerLimit <= 0 || lowerLimit >= 1) stop("lowerLimit must be between 0 and 1")
  if (upperQuantile <= 0 || upperQuantile >= 1) stop("upperQuantile must be between 0 and 1")
  if (minAveCellNeighor < 1) stop("minAveCellNeighor must be at least 1")
  invisible(cts)
}

#' Per-cell-type-pair quantile probability matching the dense path
#' @noRd
.pairPercentileProb <- function(n_i, n_j) min(1e-3, 2 / max(n_i, n_j))

# -----------------------------------------------------------------------------
# Single-slide core
# -----------------------------------------------------------------------------
.computeSparseKernelCore <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                     normalizeKernel, minAveCellNeighor,
                                     rowNormalizeKernel, colNormalizeKernel, verbose,
                                     distType, xDistScale, yDistScale, zDistScale,
                                     normalizeDistance, normalizeTarget, truncateLowDist) {

  cts <- .checkInputSparseKernel(object, sigmaValues, lowerLimit, upperQuantile,
                                 minAveCellNeighor, rowNormalizeKernel,
                                 colNormalizeKernel, distType)
  object@sigmaValues <- sigmaValues
  max_sigma <- max(sigmaValues)

  # coordinates per cell type (computed once)
  ct_coords <- stats::setNames(
    lapply(cts, function(ct) .getCoordinateMatrix(object, ct, distType,
                                                  xDistScale, yDistScale, zDistScale)),
    cts
  )

  within_only <- length(cts) == 1
  if (within_only) {
    blocks <- list(list(i = cts, j = cts, within = TRUE))
  } else {
    pct <- utils::combn(cts, 2)
    blocks <- lapply(seq_len(ncol(pct)),
                     function(k) list(i = pct[1, k], j = pct[2, k], within = FALSE))
  }

  need_pct <- truncateLowDist || normalizeDistance

  if (verbose) {
    cat(sprintf("Computing sparse kernel for %d cell type(s) [%s]\n",
                length(cts), if (within_only) "within" else "pairwise"))
  }

  # PASS 1: per-block low-distance percentile + global scaling factor
  pctls <- rep(NA_real_, length(blocks))
  if (need_pct) {
    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      A <- ct_coords[[blk$i]]
      B <- if (blk$within) NULL else ct_coords[[blk$j]]
      p <- if (blk$within) 1e-4 else .pairPercentileProb(nrow(A), nrow(B))
      pctls[b] <- .lowPercentileBlock(A, B, p)$percentile
    }
  }
  scaling_factor <- if (normalizeDistance) normalizeTarget / min(pctls) else 1
  if (normalizeDistance && verbose) {
    message(sprintf("Distance normalization scaling factor: %g", scaling_factor))
  }

  # PASS 2a: cache kernel-ready triplets per block
  block_tri <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    A <- ct_coords[[blk$i]]
    B <- if (blk$within) NULL else ct_coords[[blk$j]]
    block_tri[[b]] <- .buildBlockTriplets(A, B, pctls[b], scaling_factor,
                                          max_sigma, lowerLimit, truncateLowDist)
  }

  # PASS 2b: sigma-outer (matches dense ordering for sigma-removal semantics)
  kernel_mat <- list()
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  sigmaValuesToRemove <- stats::setNames(logical(length(sigmaValues)), sigma_names)

  for (tt in seq_along(sigmaValues)) {
    sigma_choose <- sigmaValues[tt]
    t <- sigma_names[tt]
    if (verbose) cat("current sigma value is", sigma_choose, "\n")

    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      bt <- block_tri[[b]]
      flat_name <- .createKernelMatrixName(sigma_choose, blk$i, blk$j, slide = NULL)

      if (is.null(bt)) {
        # no near pairs at all: empty kernel
        Kraw <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                     dims = c(nrow(ct_coords[[blk$i]]),
                                              if (blk$within) nrow(ct_coords[[blk$i]])
                                              else nrow(ct_coords[[blk$j]])))
      } else {
        Kraw <- .sparseKernelFromTriplets(bt, sigma_choose, lowerLimit)
      }

      should_remove <- .checkSparseSigmaRemove(
        Kraw, lowerLimit, sigma_choose, sigmaValues, blk$i, blk$j, minAveCellNeighor
      )
      if (should_remove) {
        sigmaValuesToRemove[t] <- TRUE
        kernel_mat[[flat_name]] <- Kraw
        next
      }

      kernel_mat[[flat_name]] <- .processSparseKernelMatrix(
        Kraw, lowerLimit, upperQuantile, normalizeKernel,
        rowNormalizeKernel, colNormalizeKernel
      )
    }
  }

  object <- .cleanupSigmaValues(object, kernel_mat, sigmaValuesToRemove, verbose)
  object
}

# -----------------------------------------------------------------------------
# Multi-slide core
# -----------------------------------------------------------------------------
.computeSparseKernelCoreMulti <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                          normalizeKernel, minAveCellNeighor,
                                          rowNormalizeKernel, colNormalizeKernel, verbose,
                                          distType, xDistScale, yDistScale, zDistScale,
                                          normalizeDistance, normalizeTarget, truncateLowDist) {

  cts <- .checkInputSparseKernel(object, sigmaValues, lowerLimit, upperQuantile,
                                 minAveCellNeighor, rowNormalizeKernel,
                                 colNormalizeKernel, distType)
  slides <- getSlideList(object)
  if (length(slides) == 0) stop("No slides found in multi-slide object")
  object@sigmaValues <- sigmaValues
  max_sigma <- max(sigmaValues)
  within_only <- length(cts) == 1
  need_pct <- truncateLowDist || normalizeDistance

  if (verbose) {
    cat(sprintf("Computing sparse kernel for %d cell type(s) across %d slides [%s]\n",
                length(cts), length(slides), if (within_only) "within" else "pairwise"))
  }

  # Enumerate valid blocks across slides (mirrors dense cell-count skips)
  blocks <- list()
  for (sID in slides) {
    if (within_only) {
      cnt <- .countSlideCellType(object, slide = sID, cellType = cts)
      if (cnt <= 5) next
      blocks[[length(blocks) + 1L]] <- list(slide = sID, i = cts, j = cts, within = TRUE)
    } else {
      pct <- utils::combn(cts, 2)
      for (pp in seq_len(ncol(pct))) {
        ci <- pct[1, pp]; cj <- pct[2, pp]
        if (.countSlideCellType(object, sID, ci) <= 5 ||
            .countSlideCellType(object, sID, cj) <= 5) next
        blocks[[length(blocks) + 1L]] <- list(slide = sID, i = ci, j = cj, within = FALSE)
      }
    }
  }
  if (length(blocks) == 0) stop("No slide/cell-type blocks with enough cells to compute kernels.")

  # coordinate cache keyed by slide|celltype
  coord_key <- function(sID, ct) paste(sID, ct, sep = "|")
  coord_cache <- new.env(parent = emptyenv())
  get_coords <- function(sID, ct) {
    key <- coord_key(sID, ct)
    if (is.null(coord_cache[[key]])) {
      coord_cache[[key]] <- .getCoordinateMatrix(object, ct, distType,
                                                 xDistScale, yDistScale, zDistScale,
                                                 slideID = sID)
    }
    coord_cache[[key]]
  }

  # PASS 1: per-block percentiles + GLOBAL min across all slides/pairs
  pctls <- rep(NA_real_, length(blocks))
  if (need_pct) {
    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      A <- get_coords(blk$slide, blk$i)
      B <- if (blk$within) NULL else get_coords(blk$slide, blk$j)
      p <- if (blk$within) 1e-4 else .pairPercentileProb(nrow(A), nrow(B))
      pctls[b] <- .lowPercentileBlock(A, B, p)$percentile
    }
  }
  global_min_pct <- if (need_pct) min(pctls, na.rm = TRUE) else NA_real_
  scaling_factor <- if (normalizeDistance) normalizeTarget / global_min_pct else 1
  if (normalizeDistance && verbose) {
    message(sprintf("Global distance scaling factor: %g", scaling_factor))
  }

  # PASS 2a: cache triplets per block
  block_tri <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    A <- get_coords(blk$slide, blk$i)
    B <- if (blk$within) NULL else get_coords(blk$slide, blk$j)
    block_tri[[b]] <- .buildBlockTriplets(A, B, pctls[b], scaling_factor,
                                          max_sigma, lowerLimit, truncateLowDist)
  }

  # PASS 2b: sigma-outer; drop a sigma only if NO valid kernel across all blocks
  kernel_mat <- list()
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  sigmaValuesToRemove <- stats::setNames(logical(length(sigmaValues)), sigma_names)

  for (tt in seq_along(sigmaValues)) {
    sigma_val <- sigmaValues[tt]
    sigma_name <- sigma_names[tt]
    if (verbose) cat("current sigma value is", sigma_val, "\n")
    sigma_has_valid_kernel <- FALSE

    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      bt <- block_tri[[b]]
      if (is.null(bt)) next
      Kraw <- .sparseKernelFromTriplets(bt, sigma_val, lowerLimit)
      if (!.checkSparseKernelValidityMulti(Kraw, lowerLimit, minAveCellNeighor,
                                           sigma_val, blk$i, blk$j, blk$slide)) {
        next
      }
      Kp <- .processSparseKernelMatrix(Kraw, lowerLimit, upperQuantile,
                                       normalizeKernel, rowNormalizeKernel,
                                       colNormalizeKernel)
      flat_name <- .createKernelMatrixName(sigma_val, blk$i, blk$j, slide = blk$slide)
      kernel_mat[[flat_name]] <- Kp
      sigma_has_valid_kernel <- TRUE
    }

    if (!sigma_has_valid_kernel) {
      sigmaValuesToRemove[sigma_name] <- TRUE
      warning(paste("Removing sigma value", sigma_val,
                    "as no valid kernels were generated across slides."))
    }
  }

  object <- .cleanupSigmaValuesMulti(object, kernel_mat, sigmaValuesToRemove, verbose)
  object
}

#' Compute sparse Gaussian kernels directly from coordinates
#'
#' A fused, memory-efficient alternative to [computeDistance()] +
#' [computeKernelMatrix()] for large datasets. It builds, for every cell-type
#' pair (and within-type), a sparse `dgCMatrix` Gaussian kernel using a
#' fixed-radius neighbor search, never forming a dense `n x n` matrix. Results
#' are numerically equivalent to the dense path (every pair beyond the kernel's
#' support radius is zero anyway). Distances are not stored.
#'
#' @inheritParams computeKernelMatrix
#' @param distType "Euclidean2D" or "Euclidean3D" (Morphology-Aware is not
#'   supported by the sparse path).
#' @param xDistScale,yDistScale,zDistScale per-axis coordinate scales.
#' @param normalizeDistance,normalizeTarget,truncateLowDist distance-processing
#'   options, matching [computeDistance()].
#' @return The `CoPro` object with sparse kernel matrices in `@kernelMatrices`.
#' @family spatial-pipeline
#' @seealso [computeKernelMatrix()], [computeDistance()]
#' @export
#' @rdname computeSparseKernel
setGeneric(
  "computeSparseKernel",
  function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
           normalizeKernel = FALSE, minAveCellNeighor = 2,
           rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
           distType = c("Euclidean2D", "Euclidean3D"),
           xDistScale = 1, yDistScale = 1, zDistScale = 1,
           normalizeDistance = TRUE, normalizeTarget = 0.01,
           truncateLowDist = TRUE, verbose = TRUE) standardGeneric("computeSparseKernel")
)

#' @rdname computeSparseKernel
#' @export
setMethod("computeSparseKernel", "CoProSingle",
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   distType = c("Euclidean2D", "Euclidean3D"),
                   xDistScale = 1, yDistScale = 1, zDistScale = 1,
                   normalizeDistance = TRUE, normalizeTarget = 0.01,
                   truncateLowDist = TRUE, verbose = TRUE) {
            distType <- match.arg(distType)
            .computeSparseKernelCore(object, sigmaValues, lowerLimit, upperQuantile,
                                     normalizeKernel, minAveCellNeighor,
                                     rowNormalizeKernel, colNormalizeKernel, verbose,
                                     distType, xDistScale, yDistScale, zDistScale,
                                     normalizeDistance, normalizeTarget, truncateLowDist)
          })

#' @rdname computeSparseKernel
#' @export
setMethod("computeSparseKernel", "CoProMulti",
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   distType = c("Euclidean2D", "Euclidean3D"),
                   xDistScale = 1, yDistScale = 1, zDistScale = 1,
                   normalizeDistance = TRUE, normalizeTarget = 0.01,
                   truncateLowDist = TRUE, verbose = TRUE) {
            distType <- match.arg(distType)
            .computeSparseKernelCoreMulti(object, sigmaValues, lowerLimit, upperQuantile,
                                          normalizeKernel, minAveCellNeighor,
                                          rowNormalizeKernel, colNormalizeKernel, verbose,
                                          distType, xDistScale, yDistScale, zDistScale,
                                          normalizeDistance, normalizeTarget, truncateLowDist)
          })
