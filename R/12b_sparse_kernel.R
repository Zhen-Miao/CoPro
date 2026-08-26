# =============================================================================
# Sparse, exact Gaussian-kernel construction via fixed-radius neighbor search
# -----------------------------------------------------------------------------
# A fused replacement for computeDistance() + computeKernelMatrix() that goes
# directly from spatial coordinates to a sparse kernel, never materializing a
# dense n x n distance or kernel matrix. Symmetric within-type kernels retain
# one triangle in a `dsCMatrix`; cross-type kernels use a `dgCMatrix`.
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

#' Largest number of cells in any kernel block dimension.
#'
#' Multi-slide kernels are built slide by slide, so use the largest per-slide,
#' per-cell-type count rather than the total count across all slides. This keeps
#' `method = "auto"` on the fast dense route for many small slides while still
#' protecting genuinely large blocks.
#' @noRd
.maxCellTypeCount <- function(object) {
  cts <- object@cellTypesOfInterest
  sub <- object@cellTypesSub
  if (length(cts) == 0 || length(sub) == 0) return(0L)

  if (isMultiSlide(object)) {
    slides <- getSlideList(object)
    if (length(slides) == 0) return(0L)
    counts <- vapply(
      slides,
      function(sID) max(vapply(
        cts,
        function(ct) .countSlideCellType(object, slide = sID, cellType = ct),
        numeric(1)
      )),
      numeric(1)
    )
    return(as.integer(max(counts)))
  }

  tt <- table(sub[sub %in% cts])
  if (length(tt) == 0) return(0L)
  as.integer(max(tt))
}

#' Number of entries the dense distance path would materialize.
#'
#' Counts every cross-type block (or the within-type block for a one-type
#' analysis), including all slides. Returned as a double to avoid integer
#' overflow. This catches workloads with many medium-sized blocks whose total
#' memory cost is large even though no single cell type crosses the threshold.
#' @noRd
.denseKernelEntryCount <- function(object) {
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0 || length(object@cellTypesSub) == 0) return(0)

  slides <- if (isMultiSlide(object)) getSlideList(object) else NULL
  slide_keys <- if (length(slides) > 0) slides else NA_character_

  sum(vapply(slide_keys, function(sID) {
    slide <- if (is.na(sID)) NULL else sID
    counts <- vapply(
      cts,
      function(ct) .countSlideCellType(object, slide = slide, cellType = ct),
      numeric(1)
    )
    if (length(counts) == 1) {
      counts[1]^2
    } else {
      pair_idx <- utils::combn(seq_along(counts), 2)
      sum(counts[pair_idx[1, ]] * counts[pair_idx[2, ]])
    }
  }, numeric(1)))
}

#' Number of entries required by dense multitype self-kernels
#' @noRd
.denseSelfKernelEntryCount <- function(object) {
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0 || length(object@cellTypesSub) == 0) return(0)

  slides <- if (isMultiSlide(object)) getSlideList(object) else NULL
  slide_keys <- if (length(slides) > 0) slides else NA_character_
  sum(vapply(slide_keys, function(sID) {
    slide <- if (is.na(sID)) NULL else sID
    counts <- vapply(
      cts,
      function(ct) .countSlideCellType(object, slide = slide, cellType = ct),
      numeric(1)
    )
    sum(counts^2)
  }, numeric(1)))
}


#' Whether every dense self-distance block required by computeSelfKernel exists
#' @noRd
.hasSelfDistanceMatrices <- function(object, is_multi) {
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) return(FALSE)
  if (!is_multi) {
    required <- vapply(cts, function(ct) {
      .createDistMatrixName(ct, ct, slide = NULL)
    }, character(1))
  } else {
    slides <- getSlideList(object)
    required <- unlist(lapply(slides, function(sID) {
      vapply(cts, function(ct) {
        .createDistMatrixName(ct, ct, slide = sID)
      }, character(1))
    }), use.names = FALSE)
  }
  length(required) > 0 && all(required %in% names(object@distances))
}


#' Dispatch multitype self-kernel construction
#' @noRd
.computeSelfKernelDispatch <- function(
    object, sigmaValues, lowerLimit, upperQuantile,
    normalizeKernel, minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel,
    method, autoThreshold, distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist,
    verbose, overwrite, is_multi) {
  if (!is.null(distType)) {
    distType <- match.arg(
      distType,
      c("Euclidean2D", "Euclidean3D", "Morphology-Aware")
    )
  }
  geometry <- .resolveDistanceGeometry(
    object,
    requested = list(distType = distType, xDistScale = xDistScale,
                     yDistScale = yDistScale, zDistScale = zDistScale,
                     normalizeDistance = normalizeDistance,
                     normalizeMethod = normalizeMethod,
                     normalizeTarget = normalizeTarget,
                     truncateLowDist = truncateLowDist),
    what = "computeSelfKernel", verbose = verbose
  )
  distType <- geometry$distType
  xDistScale <- geometry$xDistScale
  yDistScale <- geometry$yDistScale
  zDistScale <- geometry$zDistScale
  normalizeDistance <- geometry$normalizeDistance
  normalizeMethod <- geometry$normalizeMethod
  normalizeTarget <- geometry$normalizeTarget
  truncateLowDist <- geometry$truncateLowDist

  if (method == "auto") {
    n_max <- .maxCellTypeCount(object)
    dense_entries <- .denseSelfKernelEntryCount(object)
    has_dense_inputs <- .hasSelfDistanceMatrices(object, is_multi)
    use_sparse <- !has_dense_inputs || n_max >= autoThreshold ||
      dense_entries >= as.numeric(autoThreshold)^2
    method <- if (use_sparse) "float32" else "dense"
    if (verbose) {
      message(sprintf(
        paste0("computeSelfKernel: method='auto' -> '%s' ",
               "(largest block = %d cells, estimated dense entries = %.3g, ",
               "self-distances available = %s)."),
        method, n_max, dense_entries, has_dense_inputs
      ))
    }
  }

  if (method == "dense") {
    if (is_multi) {
      return(.computeSelfKernelCoreMulti(
        object, sigmaValues, lowerLimit, upperQuantile, normalizeKernel,
        minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel,
        verbose, overwrite
      ))
    }
    return(.computeSelfKernelCore(
      object, sigmaValues, lowerLimit, upperQuantile, normalizeKernel,
      minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel,
      verbose, overwrite
    ))
  }

  # Sparse routes only: see the matching note in .computeKernelDispatch().
  # `normalizeDistance = "inherit"` says where the scaling factor comes from,
  # not which coordinates were used, so the record keeps the resulting scaling.
  stamped <- geometry
  stamped$normalizeDistance <- .recordedNormalizeDistance(
    geometry$normalizeDistance, .getDistanceGeometry(object)
  )
  object@distanceGeometry <- stamped

  if (method == "float32") {
    return(.computeSparseKernelFloat32Core(
      object, sigmaValues, lowerLimit, upperQuantile, normalizeKernel,
      minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel,
      distType, xDistScale, yDistScale, zDistScale,
      normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist,
      overwrite = overwrite, verbose = verbose, is_multi = is_multi,
      nThreads = NULL, self_only = TRUE
    ))
  }

  .computeSparseSelfKernelCore(
    object, sigmaValues, lowerLimit, upperQuantile, normalizeKernel,
    minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel,
    distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist,
    verbose, overwrite, is_multi
  )
}

#' Choose a sparse storage format, and say so when sparse is a poor fit
#'
#' Called once `method = "auto"` has ruled out the dense path on size. Two
#' things are decided here:
#'
#' 1. Which sparse representation. `float32` is the default: it stores 8 bytes
#'    per entry against `dgCMatrix`'s 12, streams one block at a time instead of
#'    caching every block's neighbor triplets, and its `X' K X` operator agrees
#'    with the float64 path to about 1e-5 relative -- far inside the tolerance
#'    that matters for a CCA objective. `method = "sparse"` still selects
#'    float64 explicitly, for exactness checks.
#' 2. Whether sparse makes sense at all. A fixed-radius kernel is only sparse
#'    while its support radius stays well below the tissue scale; the radius
#'    grows linearly in sigma, so density grows as sigma^d and saturates. Past
#'    saturation "sparse" storage costs more than dense and buys nothing, and
#'    the user should hear about it rather than silently pay for it.
#' @noRd
.resolveSparseStorage <- function(object, sigmaValues, lowerLimit, geometry,
                                  autoThreshold, denseThreshold, is_multi,
                                  verbose) {
  cts <- if (length(object@cellTypesOfInterest) != 0L) {
    object@cellTypesOfInterest
  } else {
    unique(object@cellTypesSub)
  }

  probe <- tryCatch(
    .kernelStorageProbe(object, cts, is_multi, geometry, sigmaValues, lowerLimit),
    error = function(e) NULL
  )
  if (is.null(probe) || !is.finite(probe$maxDensity)) {
    if (verbose) {
      message("computeKernelMatrix: method='auto' -> 'float32' (density probe unavailable).")
    }
    return("float32")
  }

  bytes <- probe$nnz * 8          # float32 value + int32 column index
  dense_bytes <- probe$denseEntries * 8
  if (verbose) {
    message(sprintf(
      paste0("computeKernelMatrix: method='auto' -> 'float32' ",
             "(predicted kernel density %.1f%% at sigma = %g, ~%.2f GB per ",
             "sigma against ~%.2f GB dense)."),
      100 * probe$maxDensity, max(sigmaValues),
      bytes / 1e9, dense_bytes / 1e9))
  }

  if (probe$maxDensity >= denseThreshold) {
    # density ~ sigma^d, so invert for the sigma that would come back under it
    sigma_ok <- max(sigmaValues) *
      (denseThreshold / probe$maxDensity)^(1 / probe$dim)
    warning(sprintf(
      paste0(
        "The kernel at sigma = %g is predicted to be %.0f%% dense, so a ",
        "fixed-radius sparse kernel saves little or no memory here (sparse ",
        "storage costs more than dense past ~67%% density).\n",
        "  Sigma below about %.3g would keep it under %.0f%%. Use ",
        "detectSigmaRange() to pick sigma from the data, or build one sigma ",
        "at a time to bound peak memory."
      ),
      max(sigmaValues), 100 * probe$maxDensity, sigma_ok,
      100 * denseThreshold
    ), call. = FALSE)
  }
  "float32"
}

#' Dispatch computeKernelMatrix() to the dense or sparse path and optionally
#' drop the @distances slot afterward.
#' @noRd
.computeKernelDispatch <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                   normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                   colNormalizeKernel, verbose, method, dropDistances,
                                   autoThreshold, denseThreshold, distType, xDistScale, yDistScale, zDistScale,
                                   normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist,
                                   is_multi) {
  # Resolve the coordinate geometry once, before choosing a path. Any geometry
  # already recorded on the object wins over this function's own defaults --
  # otherwise the sparse path would rebuild coordinates from unscaled x,y (and
  # flip to 3-D whenever a z column exists) regardless of what computeDistance()
  # was told. Contradictory arguments are an error, not a silent preference.
  geometry <- .resolveDistanceGeometry(
    object,
    requested = list(distType = distType, xDistScale = xDistScale,
                     yDistScale = yDistScale, zDistScale = zDistScale,
                     normalizeDistance = normalizeDistance,
                     normalizeMethod = normalizeMethod,
                     normalizeTarget = normalizeTarget,
                     truncateLowDist = truncateLowDist),
    what = "computeKernelMatrix", verbose = verbose
  )
  distType <- geometry$distType
  xDistScale <- geometry$xDistScale
  yDistScale <- geometry$yDistScale
  zDistScale <- geometry$zDistScale
  normalizeDistance <- geometry$normalizeDistance
  normalizeMethod <- geometry$normalizeMethod
  normalizeTarget <- geometry$normalizeTarget
  truncateLowDist <- geometry$truncateLowDist

  if (method == "auto") {
    n_max <- .maxCellTypeCount(object)
    dense_entries <- .denseKernelEntryCount(object)
    sparse_by_block <- n_max >= autoThreshold
    sparse_by_total <- dense_entries >= as.numeric(autoThreshold)^2
    # The dense route reads distance matrices someone else built. With none on
    # the object it cannot run at all, so size is irrelevant: take a route that
    # builds its own coordinates. This is what lets the documented workflow --
    # detectSigmaRange() then computeKernelMatrix() -- work without a separate
    # computeDistance() call, on data of any size.
    has_distances <- length(object@distances) > 0
    if (!has_distances) {
      sparse_by_block <- TRUE
    }
    if (!sparse_by_block && !sparse_by_total) {
      # Small enough that dense is both affordable and fastest; skip the probe.
      method <- "dense"
      if (verbose) {
        message(sprintf(
          paste0("computeKernelMatrix: method='auto' -> 'dense' ",
                 "(largest block dimension = %d cells, estimated dense entries ",
                 "= %.3g, threshold = %d cells / %.3g entries)."),
          n_max, dense_entries, autoThreshold, as.numeric(autoThreshold)^2))
      }
    } else {
      method <- .resolveSparseStorage(
        object, sigmaValues, lowerLimit, geometry, autoThreshold,
        denseThreshold, is_multi, verbose
      )
    }
  }

  # Stamp the resolved geometry before dispatching, but only on the routes that
  # build coordinates here. The dense path reads distances someone else built,
  # so on a legacy object with no record there is nothing authoritative to
  # write and a guess would be worse than a blank.
  if (method %in% c("sparse", "float32")) {
    object@distanceGeometry <- geometry
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
  } else if (method == "float32") {
    object <- .computeSparseKernelFloat32Core(
      object, sigmaValues, lowerLimit, upperQuantile, normalizeKernel,
      minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel,
      distType, xDistScale, yDistScale, zDistScale,
      normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist,
      overwrite = TRUE, verbose = verbose, is_multi = is_multi
    )
  } else {  # sparse (float64)
    if (is_multi) {
      object <- .computeSparseKernelCoreMulti(object, sigmaValues, lowerLimit, upperQuantile,
                                              normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                              colNormalizeKernel, verbose,
                                              distType, xDistScale, yDistScale, zDistScale,
                                              normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist)
    } else {
      object <- .computeSparseKernelCore(object, sigmaValues, lowerLimit, upperQuantile,
                                         normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                         colNormalizeKernel, verbose,
                                         distType, xDistScale, yDistScale, zDistScale,
                                         normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist)
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
    warning(paste("Zero distances detected; applying dense-compatible",
                  "nearest-distance handling, please",
                  "consider checking the location of cells",
                  "for potential errors"))
    nz <- d[d > 0]
    if (length(nz) > 0) {
      d[d == 0] <- min(nz)
    } else {
      # Every nonzero pair is beyond the largest kernel's retained support.
      # The block-global minimum used by the dense path is therefore also
      # beyond support, so the coincident pairs must be dropped rather than
      # retained at distance zero with weight one.
      keep <- d != 0
      tri$i <- tri$i[keep]
      tri$j <- tri$j[keep]
      d <- d[keep]
      if (length(d) == 0L) return(NULL)
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
#' @param symmetric If \code{TRUE}, retain only the upper triangle and return a
#'   symmetric compressed-column matrix. Within-type triplets contain both
#'   directions, so this nearly halves persistent sparse-kernel storage.
#' @noRd
.sparseKernelFromTriplets <- function(bt, sigma, lowerLimit,
                                      symmetric = FALSE) {
  k <- exp(-0.5 * (bt$dscaled / sigma)^2)
  keep <- k >= lowerLimit
  if (symmetric) {
    if (bt$n_i != bt$n_j) {
      stop("Symmetric sparse kernel storage requires a square block")
    }
    keep <- keep & bt$i <= bt$j
  }
  Matrix::sparseMatrix(
    i = bt$i[keep], j = bt$j[keep], x = k[keep],
    dims = c(bt$n_i, bt$n_j), symmetric = symmetric
  )
}

#' Count represented sparse entries above a threshold
#'
#' A \code{dsCMatrix} stores one triangle, while kernel-validity checks are
#' defined on the represented full matrix. Self-kernels have a zero diagonal,
#' so every stored value represents two off-diagonal entries.
#' @noRd
.countSparseKernelAbove <- function(K, threshold) {
  stored <- sum(K@x > threshold)
  if (inherits(K, "symmetricMatrix")) {
    diagonal <- Matrix::diag(K)
    diagonal_count <- sum(diagonal > threshold)
    return(2 * stored - diagonal_count)
  }
  stored
}

#' Sparse analogue of .CheckSigmaValuesToRemove (single-slide, proportion-based)
#' @noRd
.checkSparseSigmaRemove <- function(K, lowerLimit, sigma_choose, sigmaValues,
                                    i, j, minAveCellNeighor) {
  # as.numeric() avoids integer overflow of n1 * n2 past ~46k cells (NA would
  # make the proportion comparison silently incorrect on large data).
  n1 <- as.numeric(nrow(K)); n2 <- as.numeric(ncol(K))
  minPropZero <- minAveCellNeighor * min(n1, n2) / (n1 * n2)
  prop_above <- .countSparseKernelAbove(K, lowerLimit) / (n1 * n2)
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
  if (.countSparseKernelAbove(K, lowerLimit) <
      minAveCellNeighor * min(n1, n2)) {
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
  # The dense/full sparse path contains both copies of every off-diagonal value,
  # so a symmetric kernel's quantile is defined on twice what it stores.
  # .type7QuantileRepeated() reproduces that exactly from the stored triangle
  # rather than materializing rep(valid, each = 2L), which would double a vector
  # that reaches hundreds of millions of entries on large data.
  repetitions <- if (inherits(K, "symmetricMatrix")) 2L else 1L
  if (length(valid) == 0) {
    warning("No valid kernel values found above lowerLimit")
    return(K)
  }

  # Clip large values (matches dense: clip computed from valid values)
  upper_clip <- as.numeric(
    .type7QuantileRepeated(valid, upperQuantile, repetitions)
  )
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

#' Resolve geometry for a direct sparse-kernel entry point
#'
#' Direct sparse APIs use `NULL` defaults so omission means inheritance, just as
#' it does for `computeKernelMatrix()`. This helper keeps their public methods
#' aligned and leaves the low-level builders with fully resolved scalar values.
#' @noRd
.resolveDirectSparseGeometry <- function(
    object, distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist,
    what, verbose) {
  if (!is.null(distType)) {
    distType <- match.arg(distType, c("Euclidean2D", "Euclidean3D"))
  }
  .resolveDistanceGeometry(
    object,
    requested = list(
      distType = distType,
      xDistScale = xDistScale,
      yDistScale = yDistScale,
      zDistScale = zDistScale,
      normalizeDistance = normalizeDistance,
      normalizeMethod = normalizeMethod,
      normalizeTarget = normalizeTarget,
      truncateLowDist = truncateLowDist
    ),
    what = what, verbose = verbose
  )
}

#' Record the geometry a sparse kernel core built its coordinates on
#'
#' Every public route resolves its arguments before reaching a core. Resolve
#' once more at this invariant boundary so an internal caller cannot create a
#' mixed-coordinate object either; contradictions are errors.
#' @noRd
.recordSparseKernelGeometry <- function(object, distType, xDistScale, yDistScale,
                                        zDistScale, normalizeDistance,
                                        normalizeMethod, normalizeTarget,
                                        truncateLowDist, what) {
  requested <- list(
    distType = distType, xDistScale = xDistScale, yDistScale = yDistScale,
    zDistScale = zDistScale, normalizeDistance = normalizeDistance,
    normalizeMethod = normalizeMethod,
    normalizeTarget = normalizeTarget, truncateLowDist = truncateLowDist
  )
  recorded <- .getDistanceGeometry(object)
  geometry <- .resolveDistanceGeometry(
    object, requested = requested, what = what, verbose = FALSE
  )
  # A record describes the scaling the matrices ended up on, so the self-kernel
  # paths' "inherit" instruction resolves to whatever it inherited from.
  geometry$normalizeDistance <- .recordedNormalizeDistance(
    geometry$normalizeDistance, recorded
  )
  object@distanceGeometry <- geometry
  object
}

# -----------------------------------------------------------------------------
# Single-slide core
# -----------------------------------------------------------------------------
.computeSparseKernelCore <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                     normalizeKernel, minAveCellNeighor,
                                     rowNormalizeKernel, colNormalizeKernel, verbose,
                                     distType, xDistScale, yDistScale, zDistScale,
                                     normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist) {

  cts <- .checkInputSparseKernel(object, sigmaValues, lowerLimit, upperQuantile,
                                 minAveCellNeighor, rowNormalizeKernel,
                                 colNormalizeKernel, distType)
  object <- .recordSparseKernelGeometry(
    object, distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist, "computeSparseKernel"
  )
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

  # The low percentile floors small distances; the scaling reference sets the
  # unit. They were one number before 1.2.0, which is why one dense block could
  # rescale the whole object. They are now computed separately.
  scaling_mode <- .normalizeDistanceMode(normalizeDistance)
  derive_own <- identical(scaling_mode, "own")
  need_pct <- truncateLowDist ||
    (derive_own && identical(normalizeMethod, "percentile"))
  need_spacing <- derive_own && identical(normalizeMethod, "spacing")

  if (verbose) {
    cat(sprintf("Computing sparse kernel for %d cell type(s) [%s]\n",
                length(cts), if (within_only) "within" else "pairwise"))
  }

  # PASS 1: per-block reference distances + global scaling factor
  pctls <- rep(NA_real_, length(blocks))
  spacings <- rep(NA_real_, length(blocks))
  if (need_pct || need_spacing) {
    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      A <- ct_coords[[blk$i]]
      B <- if (blk$within) NULL else ct_coords[[blk$j]]
      if (need_pct) {
        p <- if (blk$within) 1e-4 else .pairPercentileProb(nrow(A), nrow(B))
        pctls[b] <- .lowPercentileBlock(A, B, p)$percentile
      }
      if (need_spacing) {
        spacings[b] <- .blockNearestSpacing(
          A, if (blk$within) A else B, within = blk$within
        )
      }
    }
  }
  scaling_factor <- .selfScaleFactor(
    object, normalizeDistance,
    blockValues = if (need_spacing) spacings else pctls,
    normalizeMethod = normalizeMethod, normalizeTarget = normalizeTarget,
    distType = distType, xDistScale = xDistScale, yDistScale = yDistScale,
    zDistScale = zDistScale, what = "computeSparseKernel", verbose = verbose
  )
  if (!identical(scaling_mode, "none") && verbose) {
    message(sprintf("Distance normalization scaling factor: %g", scaling_factor))
  }
  if (!identical(scaling_mode, "none") &&
      length(object@distanceScaleFactor) == 0L) {
    object@distanceScaleFactor <- scaling_factor
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
  sigma_names <- .sigmaName(sigmaValues)
  sigmaValuesToRemove <- stats::setNames(logical(length(sigmaValues)), sigma_names)

  for (tt in seq_along(sigmaValues)) {
    sigma_choose <- sigmaValues[tt]
    t <- sigma_names[tt]
    if (verbose) cat("current sigma value is", sigma_choose, "\n")

    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      bt <- block_tri[[b]]
      flat_name <- .createKernelMatrixName(sigma_choose, blk$i, blk$j, slide = NULL)
      symmetric_storage <- blk$within &&
        !rowNormalizeKernel && !colNormalizeKernel

      if (is.null(bt)) {
        # no near pairs at all: empty kernel
        Kraw <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                     dims = c(nrow(ct_coords[[blk$i]]),
                                              if (blk$within) nrow(ct_coords[[blk$i]])
                                              else nrow(ct_coords[[blk$j]])),
                                     symmetric = symmetric_storage)
      } else {
        Kraw <- .sparseKernelFromTriplets(
          bt, sigma_choose, lowerLimit, symmetric = symmetric_storage
        )
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
                                          normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist) {

  cts <- .checkInputSparseKernel(object, sigmaValues, lowerLimit, upperQuantile,
                                 minAveCellNeighor, rowNormalizeKernel,
                                 colNormalizeKernel, distType)
  slides <- getSlideList(object)
  if (length(slides) == 0) stop("No slides found in multi-slide object")
  object <- .recordSparseKernelGeometry(
    object, distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist, "computeSparseKernel"
  )
  object@sigmaValues <- sigmaValues
  max_sigma <- max(sigmaValues)
  within_only <- length(cts) == 1
  scaling_mode <- .normalizeDistanceMode(normalizeDistance)
  derive_own <- identical(scaling_mode, "own")
  need_pct <- truncateLowDist ||
    (derive_own && identical(normalizeMethod, "percentile"))
  need_spacing <- derive_own && identical(normalizeMethod, "spacing")

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

  # PASS 1: per-block reference distances, combined into one global scale.
  # "spacing" takes the median across slide/pair blocks, so a single dense
  # block cannot set the unit for every other slide; "percentile" keeps the
  # pre-1.2.0 minimum.
  pctls <- rep(NA_real_, length(blocks))
  spacings <- rep(NA_real_, length(blocks))
  if (need_pct || need_spacing) {
    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      A <- get_coords(blk$slide, blk$i)
      B <- if (blk$within) NULL else get_coords(blk$slide, blk$j)
      if (need_pct) {
        p <- if (blk$within) 1e-4 else .pairPercentileProb(nrow(A), nrow(B))
        pctls[b] <- .lowPercentileBlock(A, B, p)$percentile
      }
      if (need_spacing) {
        spacings[b] <- .blockNearestSpacing(
          A, if (blk$within) A else B, within = blk$within
        )
      }
    }
  }
  scaling_factor <- .selfScaleFactor(
    object, normalizeDistance,
    blockValues = if (need_spacing) spacings else pctls,
    normalizeMethod = normalizeMethod, normalizeTarget = normalizeTarget,
    distType = distType, xDistScale = xDistScale, yDistScale = yDistScale,
    zDistScale = zDistScale, what = "computeSparseKernel", verbose = verbose
  )
  if (!identical(scaling_mode, "none") && verbose) {
    message(sprintf("Global distance scaling factor: %g", scaling_factor))
  }
  if (!identical(scaling_mode, "none") &&
      length(object@distanceScaleFactor) == 0L) {
    object@distanceScaleFactor <- scaling_factor
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
  sigma_names <- .sigmaName(sigmaValues)
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
      symmetric_storage <- blk$within &&
        !rowNormalizeKernel && !colNormalizeKernel
      Kraw <- .sparseKernelFromTriplets(
        bt, sigma_val, lowerLimit, symmetric = symmetric_storage
      )
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


#' Build multitype within-cell-type kernels directly from coordinates
#'
#' This is the sparse counterpart of `computeSelfDistance()` followed by
#' `computeSelfKernel()`. All self blocks share one distance-normalization
#' factor, matching the dense self-kernel workflow, while existing cross-type
#' kernels are retained unless `overwrite = TRUE`.
#' @noRd
.computeSparseSelfKernelCore <- function(
    object, sigmaValues, lowerLimit, upperQuantile,
    normalizeKernel, minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel,
    distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist,
    verbose, overwrite, is_multi) {
  cts <- .checkInputSparseKernel(
    object, sigmaValues, lowerLimit, upperQuantile, minAveCellNeighor,
    rowNormalizeKernel, colNormalizeKernel, distType
  )
  if (length(cts) == 1L) {
    warning("Only one cell type detected. Use computeKernelMatrix() instead ",
            "for single-cell-type kernels.")
    return(object)
  }
  object <- .recordSparseKernelGeometry(
    object, distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeMethod, normalizeTarget, truncateLowDist, "computeSelfKernel"
  )

  slides <- if (is_multi) getSlideList(object) else NULL
  if (is_multi && length(slides) == 0L) {
    stop("No slides found in multi-slide object")
  }

  coord_cache <- new.env(parent = emptyenv())
  coord_key <- function(slide, ct) {
    paste(if (is.null(slide)) "single" else slide, ct, sep = "|")
  }
  get_coords <- function(slide, ct) {
    key <- coord_key(slide, ct)
    if (is.null(coord_cache[[key]])) {
      coord_cache[[key]] <- .getCoordinateMatrix(
        object, ct, distType, xDistScale, yDistScale, zDistScale,
        slideID = slide
      )
    }
    coord_cache[[key]]
  }

  blocks <- list()
  slide_keys <- if (is_multi) slides else list(NULL)
  for (slide in slide_keys) {
    for (ct in cts) {
      coords <- get_coords(slide, ct)
      if (nrow(coords) <= 5L) {
        if (verbose) {
          message(sprintf(
            "Skipping self-kernel for %s%s: only %d cells.",
            ct,
            if (is.null(slide)) "" else paste0(" in slide ", slide),
            nrow(coords)
          ))
        }
        next
      }
      blocks[[length(blocks) + 1L]] <- list(
        slide = slide, i = ct, j = ct, within = TRUE
      )
    }
  }
  if (length(blocks) == 0L) {
    stop("No cell-type blocks have enough cells to compute self-kernels.")
  }

  if (verbose) {
    message(sprintf(
      "Computing sparse self-kernels for %d cell types%s.",
      length(cts),
      if (is_multi) paste0(" across ", length(slides), " slides") else ""
    ))
  }

  # normalizeDistance may be "inherit", so branch on the mode rather than on
  # the argument itself.
  scaling_mode <- .normalizeDistanceMode(normalizeDistance)
  derive_own <- identical(scaling_mode, "own")
  need_pct <- truncateLowDist ||
    (derive_own && identical(normalizeMethod, "percentile"))
  need_spacing <- derive_own && identical(normalizeMethod, "spacing")
  pctls <- rep(NA_real_, length(blocks))
  spacings <- rep(NA_real_, length(blocks))
  if (need_pct || need_spacing) {
    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      A <- get_coords(blk$slide, blk$i)
      if (need_pct) {
        pctls[b] <- .lowPercentileBlock(A, NULL, 1e-4)$percentile
      }
      if (need_spacing) {
        spacings[b] <- .blockNearestSpacing(A, A, within = TRUE)
      }
    }
  }
  # .selfScaleFactor() raises the diagnostic itself when no block yielded a
  # usable reference, rather than quietly falling back to 1 and skipping the
  # normalization the caller asked for.
  scaling_factor <- .selfScaleFactor(
    object, normalizeDistance,
    blockValues = if (need_spacing) spacings else pctls,
    normalizeMethod = normalizeMethod, normalizeTarget = normalizeTarget,
    distType = distType, xDistScale = xDistScale, yDistScale = yDistScale,
    zDistScale = zDistScale, what = "computeSelfKernel", verbose = verbose
  )
  if (!identical(scaling_mode, "none") &&
      (!is.finite(scaling_factor) || scaling_factor <= 0)) {
    stop("Cannot normalize self-kernel distances: no valid low-distance ",
         "reference was found.")
  }
  if (!identical(scaling_mode, "none")) {
    ## Only record a factor when there is not already one. Overwriting the
    ## value computeDistance() derived from the cross-type distances would
    ## leave the object claiming a scale its cross-kernels were not built at,
    ## and everything downstream that reads the slot -- the variogram
    ## normalizer in particular -- would build operators in the wrong units.
    if (length(object@distanceScaleFactor) == 0L) {
      object@distanceScaleFactor <- scaling_factor
    }
    if (verbose) message(sprintf("Self-distance scaling factor: %g", scaling_factor))
  }

  max_sigma <- max(sigmaValues)
  block_tri <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    block_tri[[b]] <- .buildBlockTriplets(
      get_coords(blk$slide, blk$i), NULL, pctls[b], scaling_factor,
      max_sigma, lowerLimit, truncateLowDist
    )
  }

  kernel_matrices <- if (overwrite || length(object@kernelMatrices) == 0L) {
    list()
  } else {
    object@kernelMatrices
  }
  sigma_names <- .sigmaName(sigmaValues)
  sigma_invalid <- stats::setNames(logical(length(sigmaValues)), sigma_names)

  for (tt in seq_along(sigmaValues)) {
    sigma <- sigmaValues[tt]
    sigma_name <- sigma_names[tt]
    if (verbose) message("Processing self-kernel sigma: ", sigma)

    for (b in seq_along(blocks)) {
      blk <- blocks[[b]]
      coords <- get_coords(blk$slide, blk$i)
      bt <- block_tri[[b]]
      symmetric_storage <- !rowNormalizeKernel && !colNormalizeKernel
      Kraw <- if (is.null(bt)) {
        Matrix::sparseMatrix(
          i = integer(0), j = integer(0), x = numeric(0),
          dims = c(nrow(coords), nrow(coords)),
          symmetric = symmetric_storage
        )
      } else {
        .sparseKernelFromTriplets(
          bt, sigma, lowerLimit, symmetric = symmetric_storage
        )
      }

      valid <- if (is_multi) {
        .checkSparseKernelValidityMulti(
          Kraw, lowerLimit, minAveCellNeighor, sigma,
          blk$i, blk$j, blk$slide
        )
      } else {
        !.checkSparseSigmaRemove(
          Kraw, lowerLimit, sigma, sigmaValues,
          blk$i, blk$j, minAveCellNeighor
        )
      }
      if (!valid) sigma_invalid[sigma_name] <- TRUE

      flat_name <- .createKernelMatrixName(
        sigma, blk$i, blk$j, slide = blk$slide
      )
      kernel_matrices[[flat_name]] <- .processSparseKernelMatrix(
        Kraw, lowerLimit, upperQuantile, normalizeKernel,
        rowNormalizeKernel, colNormalizeKernel
      )
    }
  }

  if (any(sigma_invalid)) {
    invalid_sigmas <- sigmaValues[sigma_invalid]
    is_invalid_self <- vapply(names(kernel_matrices), function(name) {
      parsed <- .parseKernelMatrixName(name)
      parsed$sigma %in% invalid_sigmas &&
        identical(parsed$cellType1, parsed$cellType2)
    }, logical(1))
    kernel_matrices <- kernel_matrices[!is_invalid_self]
    if (verbose) {
      message("Removed ", sum(sigma_invalid),
              " invalid sigma value(s) from self-kernels.")
    }
  }

  surviving <- sigmaValues[!sigma_invalid]
  object <- .pruneSigmaValues(
    object, surviving,
    invalid = sigmaValues[sigma_invalid]
  )
  object@kernelMatrices <- kernel_matrices
  object
}

#' Compute sparse Gaussian kernels directly from coordinates
#'
#' A fused, memory-efficient alternative to [computeDistance()] +
#' [computeKernelMatrix()] for large datasets. It builds, for every cell-type
#' pair (and within-type), a sparse Gaussian kernel using a fixed-radius
#' neighbor search, never forming a dense `n x n` matrix. Within-type kernels
#' are stored with one triangle as symmetric `dsCMatrix` objects; cross-type
#' kernels, and kernels made asymmetric by row/column normalization, use
#' `dgCMatrix`. Results are numerically equivalent to the dense path (every pair
#' beyond the kernel's support radius is zero anyway). Distances are not stored.
#'
#' @inheritParams computeKernelMatrix
#' @param distType "Euclidean2D" or "Euclidean3D" (Morphology-Aware is not
#'   supported by the sparse path). `NULL` inherits the recorded geometry.
#' @param xDistScale,yDistScale,zDistScale Per-axis coordinate scales. `NULL`
#'   inherits the recorded geometry.
#' @param normalizeDistance,normalizeTarget,truncateLowDist distance-processing
#'   options, matching [computeDistance()]. `NULL` inherits the recorded
#'   geometry.
#' @param normalizeMethod How the reference distance is estimated when
#'   `normalizeDistance = TRUE`. `"global"` uses the median
#'   nearest-neighbor distance over all cells of interest, ignoring their type
#'   labels, so the unit is a property of the tissue rather than of whichever
#'   blocks this call builds. `"spacing"` measures each cell-type block and
#'   takes the median across blocks. `"percentile"` reproduces the pre-1.2.0
#'   behavior: the minimum, across blocks, of a low quantile of pairwise
#'   distances.
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
           distType = NULL,
           xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
           normalizeDistance = NULL, normalizeMethod = NULL, normalizeTarget = NULL,
           truncateLowDist = NULL, verbose = TRUE) standardGeneric("computeSparseKernel")
)

#' @rdname computeSparseKernel
#' @export
setMethod("computeSparseKernel", "CoProSingle",
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   distType = NULL,
                   xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
                   normalizeDistance = NULL, normalizeMethod = NULL, normalizeTarget = NULL,
                   truncateLowDist = NULL, verbose = TRUE) {
            geometry <- .resolveDirectSparseGeometry(
              object, distType, xDistScale, yDistScale, zDistScale,
              normalizeDistance, normalizeMethod, normalizeTarget,
              truncateLowDist, "computeSparseKernel", verbose
            )
            object <- .invalidateCoProState(object, "kernel")
            .computeSparseKernelCore(object, sigmaValues, lowerLimit, upperQuantile,
                                     normalizeKernel, minAveCellNeighor,
                                     rowNormalizeKernel, colNormalizeKernel, verbose,
                                     geometry$distType, geometry$xDistScale,
                                     geometry$yDistScale, geometry$zDistScale,
                                     geometry$normalizeDistance,
                                     geometry$normalizeMethod,
                                     geometry$normalizeTarget,
                                     geometry$truncateLowDist)
          })

#' @rdname computeSparseKernel
#' @export
setMethod("computeSparseKernel", "CoProMulti",
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   distType = NULL,
                   xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
                   normalizeDistance = NULL, normalizeMethod = NULL, normalizeTarget = NULL,
                   truncateLowDist = NULL, verbose = TRUE) {
            geometry <- .resolveDirectSparseGeometry(
              object, distType, xDistScale, yDistScale, zDistScale,
              normalizeDistance, normalizeMethod, normalizeTarget,
              truncateLowDist, "computeSparseKernel", verbose
            )
            object <- .invalidateCoProState(object, "kernel")
            .computeSparseKernelCoreMulti(object, sigmaValues, lowerLimit, upperQuantile,
                                          normalizeKernel, minAveCellNeighor,
                                          rowNormalizeKernel, colNormalizeKernel, verbose,
                                          geometry$distType, geometry$xDistScale,
                                          geometry$yDistScale, geometry$zDistScale,
                                          geometry$normalizeDistance,
                                          geometry$normalizeMethod,
                                          geometry$normalizeTarget,
                                          geometry$truncateLowDist)
          })
