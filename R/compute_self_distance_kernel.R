#' Compute Self-Distance and Self-Kernel for Multiple Cell Types
#'
#' These functions extend the CoPro package to compute within-cell-type (self) 
#' distance and kernel matrices even when multiple cell types are present.
#' This is useful for analyses that need both cross-type and within-type 
#' spatial relationships.
#'
#' @name self_distance_kernel
#' @keywords internal
NULL

#' Compute Self-Distance Matrices for Multiple Cell Types
#'
#' This function computes within-cell-type distance matrices for each cell type
#' when multiple cell types are present. Unlike the standard computeDistance
#' which only computes cross-type distances for multiple cell types, this
#' function computes self-distances (within each cell type).
#'
#' `computeSelfDistance()` is normally additive: it adds within-type blocks to
#' an object whose cross-type blocks [computeDistance()] already built. Its
#' distance arguments therefore default to `NULL`, meaning "use whatever the
#' existing distances were built on" ([getDistanceGeometry()]); supplying a
#' value that contradicts that record is an error rather than a silent
#' second basis. When nothing is recorded -- a self-distance-only workflow --
#' the package defaults apply.
#'
#' When `normalizeDistance = TRUE` and the object already carries a scale
#' factor, that factor is reused rather than re-derived from the within-type
#' blocks, so cross-type and within-type distances stay on one unit. Pass
#' `overwrite = TRUE` to discard the existing blocks and re-derive.
#'
#' @param object A `CoPro` object with multiple cell types
#' @param distType Type of distance to compute: `"Euclidean2D"`,
#'  `"Euclidean3D"`, or `"Morphology-Aware"`. `NULL` (default) inherits the
#'  recorded geometry, falling back to the coordinate columns present.
#' @param xDistScale Scale for x distance. `NULL` (default) inherits the
#'  recorded geometry, falling back to 1.
#' @param yDistScale Scale for y distance. `NULL` (default) inherits the
#'  recorded geometry, falling back to 1.
#' @param zDistScale Scale for z distance. `NULL` (default) inherits the
#'  recorded geometry, falling back to 1.
#' @param normalizeDistance Whether to rescale distances so that the reference
#'  cell-cell distance becomes `normalizeTarget`. `NULL` (default) inherits the
#'  recorded geometry, falling back to `FALSE` as of CoPro 1.2.0.
#' @param normalizeMethod How the reference distance is estimated when
#'  `normalizeDistance = TRUE` and no scale factor is already pinned on the
#'  object. `"global"` uses the median nearest-neighbor distance over all cells,
#'  ignoring type labels, which gives the same unit here as in
#'  [computeDistance()]; `"spacing"` measures the within-type blocks instead;
#'  `"percentile"` reproduces the pre-1.2.0 behavior (the minimum, across
#'  blocks, of a low quantile of pairwise distances). `NULL` (default) inherits
#'  the recorded geometry, falling back to `"global"`.
#' @param normalizeTarget Numeric scalar. The value the reference distance is
#'  rescaled to. `NULL` (default) inherits the recorded geometry, falling back
#'  to 0.01.
#' @param truncateLowDist Whether to truncate small distances so that cells
#'  that are nearly overlapping do not have super small distances. `NULL`
#'  (default) inherits the recorded geometry, falling back to `TRUE`.
#' @param verbose Whether to print info about the computation progress
#' @param overwrite Whether to overwrite existing distance matrices. If FALSE,
#'  will add self-distance matrices to existing cross-type distances. Default = FALSE
#'
#' @return `CoPro` object with self-distance matrices added to the distances slot
#' @seealso [computeDistance()], [getDistanceGeometry()]
#' @export
#' @rdname computeSelfDistance
#' @aliases computeSelfDistance,CoProSingle-method
#' @aliases computeSelfDistance,CoProMulti-method
#'
#' @examples
#' \dontrun{
#' # Assume you have a CoPro object with multiple cell types
#' # First compute cross-type distances
#' object <- computeDistance(object)
#' 
#' # Then add self-distances
#' object <- computeSelfDistance(object)
#' 
#' # Now you have both cross-type and self-type distance matrices
#' }
setGeneric(
  "computeSelfDistance",
  function(object, distType = NULL,
           xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
           normalizeDistance = NULL, normalizeMethod = NULL,
           normalizeTarget = NULL, truncateLowDist = NULL,
           verbose = TRUE, overwrite = FALSE) standardGeneric("computeSelfDistance")
)

#' @rdname computeSelfDistance
#' @export
setMethod("computeSelfDistance", "CoProSingle",
          function(object, distType = NULL,
                   xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
                   normalizeDistance = NULL, normalizeMethod = NULL,
                   normalizeTarget = NULL, truncateLowDist = NULL,
                   verbose = TRUE, overwrite = FALSE) {
            .computeSelfDistanceCore(object, distType, xDistScale, yDistScale, zDistScale,
                                    normalizeDistance, normalizeMethod, normalizeTarget,
                                    truncateLowDist, verbose, overwrite)
          })

#' @rdname computeSelfDistance
#' @export
setMethod("computeSelfDistance", "CoProMulti",
          function(object, distType = NULL,
                   xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
                   normalizeDistance = NULL, normalizeMethod = NULL,
                   normalizeTarget = NULL, truncateLowDist = NULL,
                   verbose = TRUE, overwrite = FALSE) {
            .computeSelfDistanceCoreMulti(object, distType, xDistScale, yDistScale, zDistScale,
                                         normalizeDistance, normalizeMethod, normalizeTarget,
                                         truncateLowDist, verbose, overwrite)
          })

#' Resolve the geometry a self-distance step should build on
#'
#' `computeSelfDistance()` runs on top of whatever [computeDistance()] left
#' behind, so its geometry arguments are `NULL`-defaulted and inherited rather
#' than re-defaulted. Wraps [.resolveDistanceGeometry()] so both the single-
#' and multi-slide cores enter with the same resolved values.
#' @noRd
.resolveSelfDistanceGeometry <- function(object, distType, xDistScale, yDistScale,
                                         zDistScale, normalizeDistance,
                                         normalizeMethod, normalizeTarget,
                                         truncateLowDist, verbose,
                                         rebuild = FALSE) {
  geometry <- .resolveDistanceGeometry(
    object,
    requested = list(
      distType = distType, xDistScale = xDistScale, yDistScale = yDistScale,
      zDistScale = zDistScale, normalizeDistance = normalizeDistance,
      normalizeMethod = normalizeMethod, normalizeTarget = normalizeTarget,
      truncateLowDist = truncateLowDist
    ),
    what = "computeSelfDistance", verbose = verbose, rebuild = rebuild
  )
  if (isTRUE(geometry$normalizeDistance) &&
      geometry$normalizeMethod %in% c("global", "spacing") &&
      identical(geometry$distType, "Morphology-Aware")) {
    warning(sprintf(paste("normalizeMethod = '%s' is not defined for",
                          "Morphology-Aware distances; using 'percentile'."),
                    geometry$normalizeMethod))
    geometry$normalizeMethod <- "percentile"
  }
  geometry
}

#' Core function for computing self-distances (single slide)
#' @noRd
.computeSelfDistanceCore <- function(object, distType, xDistScale, yDistScale, zDistScale,
                                    normalizeDistance, normalizeMethod, normalizeTarget,
                                    truncateLowDist, verbose, overwrite) {
  # Resolve against the record before clearing it: overwrite means "re-derive
  # the scale", not "stop normalizing", so the recorded intent still supplies
  # the defaults for arguments the caller did not give.
  geometry <- .resolveSelfDistanceGeometry(
    object, distType, xDistScale, yDistScale, zDistScale, normalizeDistance,
    normalizeMethod, normalizeTarget, truncateLowDist, verbose,
    rebuild = overwrite
  )
  if (overwrite) object <- .clearDistanceRecord(object)
  distType <- geometry$distType
  xDistScale <- geometry$xDistScale
  yDistScale <- geometry$yDistScale
  zDistScale <- geometry$zDistScale
  normalizeDistance <- geometry$normalizeDistance
  normalizeMethod <- geometry$normalizeMethod
  normalizeTarget <- geometry$normalizeTarget
  truncateLowDist <- geometry$truncateLowDist

  # Validate inputs
  cts <- .checkInputDistance(object, distType, xDistScale, yDistScale, zDistScale)

  if (length(cts) == 1) {
    warning("Only one cell type detected. Use computeDistance() instead for single cell type.")
    return(object)
  }

  # Initialize or preserve existing distances
  if (overwrite || length(object@distances) == 0) {
    distances <- list()
  } else {
    distances <- object@distances
  }
  
  if (verbose) {
    cat("Computing self-distance matrices for", length(cts), "cell types\n")
  }
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance && verbose) {
    message(sprintf(
      paste0("normalizeDistance is TRUE, so self-distances will be rescaled ",
             "so that the reference cell-cell distance becomes %g."),
      normalizeTarget
    ))
  }

  # Compute self-distances for each cell type
  all_percentiles <- numeric(length(cts))
  names(all_percentiles) <- cts
  
  for (i in seq_along(cts)) {
    ct <- cts[i]
    
    if (verbose) {
      cat("Computing self-distance for cell type:", ct, "\n")
    }
    
    # Get coordinate matrix for this cell type
    mat <- .getCoordinateMatrix(object, ct, distType, xDistScale, yDistScale, zDistScale)
    
    # Check if we have enough cells
    if (nrow(mat) <= 5) {
      warning(paste("Insufficient cells for cell type", ct, "(", nrow(mat), "cells). Skipping."))
      next
    }
    
    # Compute self-distance matrix
    distances_ii <- fields::rdist(mat)
    
    # Process distance matrix (set diagonal to Inf)
    processed <- .processDistanceMatrix(distances_ii, truncateLowDist, 
                                       percentile_choice = 1e-4, set_diag_inf = TRUE)
    distances_ii <- processed$distances
    all_percentiles[i] <- processed$percentile
    
    # Save using flat structure
    flat_name <- .createDistMatrixName(ct, ct, slide = NULL)
    distances[[flat_name]] <- distances_ii
    
    if (verbose) {
      cat("Cell type:", ct, "- Cells:", nrow(mat), "\n")
      print(quantile(distances_ii[is.finite(distances_ii)], na.rm = TRUE))
    }
  }
  
  # Apply normalization if requested
  scaling_factor <- 1
  normalized <- FALSE
  if (normalizeDistance) {
    # "global" measures the cells rather than these blocks, so it needs no
    # per-block reference at all.
    reference_values <- if (identical(normalizeMethod, "spacing")) {
      vapply(cts, function(ct) {
        coords_ct <- .getCoordinateMatrix(object, ct, distType,
                                          xDistScale, yDistScale, zDistScale)
        .blockNearestSpacing(coords_ct, coords_ct, within = TRUE)
      }, numeric(1))
    } else if (identical(normalizeMethod, "percentile")) {
      all_percentiles
    } else {
      numeric(0)
    }
    valid_percentiles <- reference_values[!is.na(reference_values) & is.finite(reference_values)]
    if (identical(normalizeMethod, "global") || length(valid_percentiles) > 0) {
      scaling_factor <- .normalizationScaleFactor(
        object, blockValues = valid_percentiles,
        normalizeMethod = normalizeMethod, normalizeTarget = normalizeTarget,
        distType = distType, xDistScale = xDistScale, yDistScale = yDistScale,
        zDistScale = zDistScale, what = "computeSelfDistance", verbose = verbose
      )
      if (verbose) {
        message(sprintf("Self-distance scaling factor: %g", scaling_factor))
      }

      for (ct in cts) {
        flat_name <- .createDistMatrixName(ct, ct, slide = NULL)
        if (flat_name %in% names(distances) && !is.null(distances[[flat_name]])) {
          distances[[flat_name]] <- distances[[flat_name]] * scaling_factor
          normalized <- TRUE
        }
      }
    }
  }

  object@distances <- distances
  # Pin the factor so a later step normalizes the blocks it adds the same way,
  # and so .recoverDistanceScaleFactor() can map back to raw coordinates. The
  # record has to agree with the slot: see .recordNormalizationOutcome().
  object@distanceScaleFactor <- scaling_factor
  object@distanceGeometry <- .recordNormalizationOutcome(
    geometry, normalized, "computeSelfDistance"
  )
  return(object)
}

#' Core function for computing self-distances (multi-slide)
#' @noRd
.computeSelfDistanceCoreMulti <- function(object, distType, xDistScale, yDistScale, zDistScale,
                                         normalizeDistance, normalizeMethod, normalizeTarget,
                                         truncateLowDist, verbose, overwrite) {
  # See .computeSelfDistanceCore(): resolve first, clear second.
  geometry <- .resolveSelfDistanceGeometry(
    object, distType, xDistScale, yDistScale, zDistScale, normalizeDistance,
    normalizeMethod, normalizeTarget, truncateLowDist, verbose,
    rebuild = overwrite
  )
  if (overwrite) object <- .clearDistanceRecord(object)
  distType <- geometry$distType
  xDistScale <- geometry$xDistScale
  yDistScale <- geometry$yDistScale
  zDistScale <- geometry$zDistScale
  normalizeDistance <- geometry$normalizeDistance
  normalizeMethod <- geometry$normalizeMethod
  normalizeTarget <- geometry$normalizeTarget
  truncateLowDist <- geometry$truncateLowDist

  # Validate inputs
  cts <- .checkInputDistance(object, distType, xDistScale, yDistScale, zDistScale)
  
  if (length(cts) == 1) {
    warning("Only one cell type detected. Use computeDistance() instead for single cell type.")
    return(object)
  }
  
  slides <- getSlideList(object)
  
  # Initialize or preserve existing distances
  if (overwrite || length(object@distances) == 0) {
    distances_all <- list()
  } else {
    distances_all <- object@distances
  }
  
  if (verbose) {
    cat("Computing self-distance matrices for", length(cts), "cell types across", length(slides), "slides\n")
  }
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance && verbose) {
    message(sprintf(
      paste0("normalizeDistance is TRUE, so self-distances will be rescaled ",
             "across all slides so that the reference cell-cell distance ",
             "becomes %g."),
      normalizeTarget
    ))
  }

  global_min_percentile <- Inf
  self_spacings <- numeric(0)
  
  # Compute self-distances for each cell type across all slides
  for (ct in cts) {
    if (verbose) {
      cat("Computing self-distance for cell type:", ct, "\n")
    }
    
    for (sID in slides) {
      # Check if there are enough cells for this cell type in this slide
      slide_ct_count <- .countSlideCellType(object, slide = sID, cellType = ct)
      if (slide_ct_count <= 5) {
        if (verbose) message(paste("Skipping", ct, "in slide", sID, "- insufficient cells (", slide_ct_count, "cells)"))
        next
      }
      
      # Get coordinate matrix for this slide and cell type
      mat <- .getCoordinateMatrix(object, ct, distType, xDistScale, yDistScale, zDistScale, slideID = sID)
      
      # Compute self-distance matrix
      distances_ii <- fields::rdist(mat)
      
      # Process distance matrix (set diagonal to Inf)
      processed <- .processDistanceMatrix(distances_ii, truncateLowDist, 
                                         percentile_choice = 1e-4, set_diag_inf = TRUE)
      distances_ii <- processed$distances
      dist_percentile <- processed$percentile
      
      if (!is.na(dist_percentile) && is.finite(dist_percentile)) {
        global_min_percentile <- min(global_min_percentile, dist_percentile, na.rm = TRUE)
        if (normalizeDistance && identical(normalizeMethod, "spacing")) {
          self_spacings <- c(self_spacings, .blockSpacingForBlock(
            object, distType, xDistScale, yDistScale, zDistScale, sID, ct, ct
          ))
        }
      }
      
      # Save using flat structure
      flat_name <- .createDistMatrixName(ct, ct, slide = sID)
      distances_all[[flat_name]] <- distances_ii
      
      if (verbose) {
        cat("Slide:", sID, ", Cell type:", ct, "- Cells:", nrow(mat), "\n")
        print(quantile(distances_ii[is.finite(distances_ii)], na.rm = TRUE))
      }
    }
  }
  
  # Apply normalization if requested
  scaling_factor <- 1
  normalized <- FALSE
  if (normalizeDistance && is.finite(global_min_percentile)) {
    scaling_factor <- .normalizationScaleFactor(
      object,
      blockValues = if (identical(normalizeMethod, "spacing")) self_spacings
                    else global_min_percentile,
      normalizeMethod = normalizeMethod, normalizeTarget = normalizeTarget,
      distType = distType, xDistScale = xDistScale, yDistScale = yDistScale,
      zDistScale = zDistScale, what = "computeSelfDistance", verbose = verbose
    )
    if (verbose) {
      message(sprintf("Global self-distance scaling factor: %g", scaling_factor))
    }

    for (ct in cts) {
      for (sID in slides) {
        flat_name <- .createDistMatrixName(ct, ct, slide = sID)
        if (flat_name %in% names(distances_all) && !is.null(distances_all[[flat_name]])) {
          distances_all[[flat_name]] <- distances_all[[flat_name]] * scaling_factor
          normalized <- TRUE
        }
      }
    }
  }

  object@distances <- distances_all
  object@distanceScaleFactor <- scaling_factor
  object@distanceGeometry <- .recordNormalizationOutcome(
    geometry, normalized, "computeSelfDistance"
  )
  return(object)
}

#' Compute Self-Kernel Matrices for Multiple Cell Types
#'
#' This function computes within-cell-type kernel matrices for each cell type
#' when multiple cell types are present. The default `method = "auto"` uses
#' existing self-distance matrices for small workloads and builds sparse kernels
#' directly from coordinates for large workloads or when self-distances have not
#' been materialized.
#'
#' @param object A `CoPro` object with multiple cell types
#' @param sigmaValues A vector of sigma values used for kernel calculation
#' @param lowerLimit The lower limit for the kernel function, default is 1e-7
#' @param upperQuantile The quantile used for clipping the kernel values, default is 0.85
#' @param normalizeKernel Whether to normalize the kernel matrix? Default = FALSE
#' @param minAveCellNeighor Minimum average number of neighbors. Default = 2
#' @param rowNormalizeKernel Whether to row-normalize kernel matrices. Default = FALSE
#' @param colNormalizeKernel Whether to column-normalize kernel matrices. Default = FALSE
#' @param verbose Whether to output progress information
#' @param overwrite Whether to overwrite existing kernel matrices. If FALSE,
#'  will add self-kernel matrices to existing cross-type kernels. Default = FALSE
#' @param normalizeMethod How the reference distance is estimated when
#'  `normalizeDistance = TRUE`: `"global"` (median nearest-neighbor distance
#'  over all cells, ignoring type labels), `"spacing"` (median nearest-partner
#'  distance per block, combined by median), or `"percentile"` (pre-1.2.0
#'  behavior). See [computeDistance()].
#' @param method One of `"auto"`, `"dense"`, or `"sparse"`. The sparse path
#'  constructs exact thresholded self-kernels directly from coordinates and
#'  does not require [computeSelfDistance()].
#' @param autoThreshold Cell-count threshold used by `method = "auto"`.
#'  Sparse construction is selected when a self-kernel dimension reaches this
#'  value, aggregate dense self-kernel entries reach its square, or required
#'  self-distance matrices are absent. Default 5000.
#' @param distType,xDistScale,yDistScale,zDistScale,normalizeDistance,normalizeTarget,truncateLowDist
#'  Distance options for the sparse path, matching [computeKernelMatrix()]:
#'  `NULL` inherits the geometry recorded by [computeDistance()] /
#'  [computeSelfDistance()], and a value that contradicts that record is an
#'  error.
#'
#' @return `CoPro` object with self-kernel matrices added to the kernelMatrices slot
#' @export
#' @rdname computeSelfKernel
#' @aliases computeSelfKernel,CoProSingle-method
#' @aliases computeSelfKernel,CoProMulti-method
#'
#' @examples
#' \dontrun{
#' # Assume you have a CoPro object with multiple cell types
#' # First compute cross-type distances and kernels
#' object <- computeDistance(object)
#' object <- computeKernelMatrix(object, sigmaValues = c(0.01, 0.05, 0.1))
#' 
#' # Then add self-distances and self-kernels
#' object <- computeSelfDistance(object)
#' object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))
#' 
#' # Now you have both cross-type and self-type kernel matrices
#' }
setGeneric(
  "computeSelfKernel",
  function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
           normalizeKernel = FALSE, minAveCellNeighor = 2,
           rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
           verbose = TRUE, overwrite = FALSE,
           method = c("auto", "dense", "sparse"), autoThreshold = 5000L,
           distType = NULL, xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
           normalizeDistance = NULL, normalizeMethod = NULL, normalizeTarget = NULL,
           truncateLowDist = NULL) standardGeneric("computeSelfKernel")
)

#' @rdname computeSelfKernel
#' @export
setMethod("computeSelfKernel", "CoProSingle", 
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   verbose = TRUE, overwrite = FALSE,
                   method = c("auto", "dense", "sparse"), autoThreshold = 5000L,
                   distType = NULL, xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
                   normalizeDistance = NULL, normalizeMethod = NULL, normalizeTarget = NULL,
                   truncateLowDist = NULL) {
            .computeSelfKernelDispatch(
              object, sigmaValues, lowerLimit, upperQuantile,
              normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
              colNormalizeKernel, method = match.arg(method),
              autoThreshold = autoThreshold, distType = distType,
              xDistScale = xDistScale, yDistScale = yDistScale,
              zDistScale = zDistScale, normalizeDistance = normalizeDistance,
              normalizeMethod = normalizeMethod,
              normalizeTarget = normalizeTarget, truncateLowDist = truncateLowDist,
              verbose = verbose, overwrite = overwrite, is_multi = FALSE
            )
          })

#' @rdname computeSelfKernel
#' @export
setMethod("computeSelfKernel", "CoProMulti", 
          function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
                   normalizeKernel = FALSE, minAveCellNeighor = 2,
                   rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
                   verbose = TRUE, overwrite = FALSE,
                   method = c("auto", "dense", "sparse"), autoThreshold = 5000L,
                   distType = NULL, xDistScale = NULL, yDistScale = NULL, zDistScale = NULL,
                   normalizeDistance = NULL, normalizeMethod = NULL, normalizeTarget = NULL,
                   truncateLowDist = NULL) {
            .computeSelfKernelDispatch(
              object, sigmaValues, lowerLimit, upperQuantile,
              normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
              colNormalizeKernel, method = match.arg(method),
              autoThreshold = autoThreshold, distType = distType,
              xDistScale = xDistScale, yDistScale = yDistScale,
              zDistScale = zDistScale, normalizeDistance = normalizeDistance,
              normalizeMethod = normalizeMethod,
              normalizeTarget = normalizeTarget, truncateLowDist = truncateLowDist,
              verbose = verbose, overwrite = overwrite, is_multi = TRUE
            )
          })

#' Core function for computing self-kernels (single slide)
#' @noRd
.computeSelfKernelCore <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                  normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                  colNormalizeKernel, verbose, overwrite) {
  
  # Get cell types
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) {
    stop("No cell types of interest specified")
  }
  
  if (length(cts) == 1) {
    warning("Only one cell type detected. Use computeKernelMatrix() instead for single cell type.")
    return(object)
  }
  
  # Check for conflicting normalization options
  if (rowNormalizeKernel && colNormalizeKernel) {
    stop("Cannot do both row-wise and column-wise normalization.")
  }
  
  # Validate sigma values
  if (!is.numeric(sigmaValues) || any(sigmaValues <= 0)) {
    stop("sigmaValues must be positive numeric values")
  }
  
  # Check if self-distance matrices exist
  missing_self_dists <- character(0)
  for (ct in cts) {
    flat_name <- .createDistMatrixName(ct, ct, slide = NULL)
    if (!flat_name %in% names(object@distances)) {
      missing_self_dists <- c(missing_self_dists, ct)
    }
  }
  
  if (length(missing_self_dists) > 0) {
    stop(paste("Self-distance matrices missing for cell types:", 
               paste(missing_self_dists, collapse = ", "), 
               ". Run computeSelfDistance() first."))
  }
  
  # Initialize or preserve existing kernels
  if (overwrite || length(object@kernelMatrices) == 0) {
    kernel_matrices <- list()
  } else {
    kernel_matrices <- object@kernelMatrices
  }
  
  if (verbose) {
    cat("Computing self-kernel matrices for", length(cts), "cell types\n")
    if (normalizeKernel) {
      cat("normalizeKernel is set to TRUE. Self-kernel matrices will be",
          "normalized so that median row sums will be 1\n")
    }
  }
  
  # Track sigma values to remove
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  sigmaValuesToRemove <- vector(mode = "logical", length = length(sigmaValues))
  names(sigmaValuesToRemove) <- sigma_names
  
  # Compute self-kernels for each sigma and cell type
  for (tt in seq_along(sigmaValues)) {
    sigma_val <- sigmaValues[tt]
    sigma_name <- sigma_names[tt]
    
    if (verbose) {
      cat("Processing sigma value:", sigma_val, "\n")
    }
    
    sigma_invalid <- FALSE
    
    for (ct in cts) {
      # Get self-distance matrix
      dist_flat_name <- .createDistMatrixName(ct, ct, slide = NULL)
      dist_mat <- object@distances[[dist_flat_name]]
      
      if (is.null(dist_mat)) {
        warning(paste("Self-distance matrix not found for cell type", ct, ". Skipping."))
        sigma_invalid <- TRUE
        next
      }
      
      # Compute kernel from distance
      kernel_current <- kernel_from_distance(
        sigma = sigma_val,
        dist_mat = dist_mat,
        lower_limit = lowerLimit
      )
      
      # Check if sigma should be removed
      should_remove <- .CheckSigmaValuesToRemove(
        kernel_current = kernel_current, lowerLimit = lowerLimit,
        minAveCellNeighor = minAveCellNeighor, sigma_choose = sigma_val,
        sigmaValues = sigmaValues, i = ct, j = ct
      )
      
      if (should_remove) {
        sigma_invalid <- TRUE
        # Still store the kernel even if marked for removal
        flat_name <- .createKernelMatrixName(sigma_val, ct, ct, slide = NULL)
        kernel_matrices[[flat_name]] <- kernel_current
        next
      }
      
      # Process kernel matrix (clipping and normalization)
      kernel_current <- .processKernelMatrix(kernel_current, lowerLimit, upperQuantile,
                                           normalizeKernel, rowNormalizeKernel, colNormalizeKernel)
      
      # Store in flat structure
      flat_name <- .createKernelMatrixName(sigma_val, ct, ct, slide = NULL)
      kernel_matrices[[flat_name]] <- kernel_current
      
      if (verbose) {
        cat("  Cell type:", ct, "- Kernel matrix:", nrow(kernel_current), "x", ncol(kernel_current), "\n")
      }
    }
    
    if (sigma_invalid) {
      sigmaValuesToRemove[sigma_name] <- TRUE
    }
  }
  
  # Clean up removed sigma values if needed
  if (any(sigmaValuesToRemove)) {
    if (verbose) {
      cat("Removing", sum(sigmaValuesToRemove), "sigma values due to invalid kernels\n")
    }
    
    sigmas_to_remove <- as.numeric(gsub("sigma_", "", names(sigmaValuesToRemove)[sigmaValuesToRemove]))
    
    # Remove flat entries that contain these sigma values (only for self-kernels)
    to_remove <- sapply(names(kernel_matrices), function(name) {
      parsed <- .parseKernelMatrixName(name)
      # Only remove self-kernels with invalid sigmas
      parsed$sigma %in% sigmas_to_remove && parsed$cellType1 == parsed$cellType2
    })
    
    if (any(to_remove)) {
      kernel_matrices <- kernel_matrices[!to_remove]
    }
  }
  
  object@kernelMatrices <- kernel_matrices
  return(object)
}

#' Core function for computing self-kernels (multi-slide)
#' @noRd
.computeSelfKernelCoreMulti <- function(object, sigmaValues, lowerLimit, upperQuantile,
                                       normalizeKernel, minAveCellNeighor, rowNormalizeKernel,
                                       colNormalizeKernel, verbose, overwrite) {
  
  # Get cell types and slides
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) {
    stop("No cell types of interest specified")
  }
  
  if (length(cts) == 1) {
    warning("Only one cell type detected. Use computeKernelMatrix() instead for single cell type.")
    return(object)
  }
  
  slides <- getSlideList(object)
  
  # Check for conflicting normalization options
  if (rowNormalizeKernel && colNormalizeKernel) {
    stop("Cannot do both row-wise and column-wise normalization.")
  }
  
  # Validate sigma values
  if (!is.numeric(sigmaValues) || any(sigmaValues <= 0)) {
    stop("sigmaValues must be positive numeric values")
  }
  
  # Check if self-distance matrices exist
  missing_self_dists <- character(0)
  for (ct in cts) {
    for (sID in slides) {
      flat_name <- .createDistMatrixName(ct, ct, slide = sID)
      if (!flat_name %in% names(object@distances)) {
        missing_self_dists <- c(missing_self_dists, paste(ct, sID, sep = "@"))
      }
    }
  }
  
  if (length(missing_self_dists) > 0) {
    stop(paste("Self-distance matrices missing for:", 
               paste(missing_self_dists, collapse = ", "), 
               ". Run computeSelfDistance() first."))
  }
  
  # Initialize or preserve existing kernels
  if (overwrite || length(object@kernelMatrices) == 0) {
    kernel_matrices_all <- list()
  } else {
    kernel_matrices_all <- object@kernelMatrices
  }
  
  if (verbose) {
    cat("Computing self-kernel matrices for", length(cts), "cell types across", length(slides), "slides\n")
    if (normalizeKernel) {
      cat("normalizeKernel is set to TRUE. Self-kernel matrices will be",
          "normalized so that median row sums will be 1\n")
    }
  }
  
  # Track sigma values to remove
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  sigmaValuesToRemove <- vector(mode = "logical", length = length(sigmaValues))
  names(sigmaValuesToRemove) <- sigma_names
  
  # Compute self-kernels for each sigma, slide, and cell type
  for (tt in seq_along(sigmaValues)) {
    sigma_val <- sigmaValues[tt]
    sigma_name <- sigma_names[tt]
    
    if (verbose) {
      cat("Processing sigma value:", sigma_val, "\n")
    }
    
    sigma_valid_across_slides <- TRUE
    
    for (sID in slides) {
      for (ct in cts) {
        # Get self-distance matrix for this slide and cell type
        dist_flat_name <- .createDistMatrixName(ct, ct, slide = sID)
        dist_mat <- object@distances[[dist_flat_name]]
        
        if (is.null(dist_mat)) {
          if (verbose) message(paste("Self-distance matrix not found for", ct, "in slide", sID, ". Skipping."))
          next
        }
        
        # Compute kernel from distance
        kernel_current <- kernel_from_distance(
          sigma = sigma_val,
          dist_mat = dist_mat,
          lower_limit = lowerLimit
        )
        
        # Check kernel validity
        if (!.checkKernelValidityMulti(kernel_current, lowerLimit, minAveCellNeighor,
                                      sigma_val, ct, ct, sID)) {
          sigma_valid_across_slides <- FALSE
        }
        
        # Process kernel matrix (clipping and normalization)
        kernel_current <- .processKernelMatrix(kernel_current, lowerLimit, upperQuantile,
                                             normalizeKernel, rowNormalizeKernel, colNormalizeKernel)
        
        # Store in flat structure
        flat_name <- .createKernelMatrixName(sigma_val, ct, ct, slide = sID)
        kernel_matrices_all[[flat_name]] <- kernel_current
        
        if (verbose) {
          cat("  Slide:", sID, ", Cell type:", ct, "- Kernel matrix:", nrow(kernel_current), "x", ncol(kernel_current), "\n")
        }
      }
    }
    
    # Mark sigma for removal if invalid across slides
    if (!sigma_valid_across_slides) {
      sigmaValuesToRemove[sigma_name] <- TRUE
      warning(paste("Removing sigma value", sigma_val, "as it was invalid for one or more self-kernels."))
    }
  }
  
  # Clean up removed sigma values if needed
  if (any(sigmaValuesToRemove)) {
    if (verbose) {
      cat("Removing", sum(sigmaValuesToRemove), "sigma values due to invalid self-kernels\n")
    }
    
    sigmas_to_remove <- as.numeric(gsub("sigma_", "", names(sigmaValuesToRemove)[sigmaValuesToRemove]))
    
    # Remove flat entries that contain these sigma values (only for self-kernels)
    to_remove <- sapply(names(kernel_matrices_all), function(name) {
      parsed <- .parseKernelMatrixName(name)
      # Only remove self-kernels with invalid sigmas
      parsed$sigma %in% sigmas_to_remove && parsed$cellType1 == parsed$cellType2
    })
    
    if (any(to_remove)) {
      kernel_matrices_all <- kernel_matrices_all[!to_remove]
    }
  }
  
  object@kernelMatrices <- kernel_matrices_all
  return(object)
}

#' Get Self-Distance Matrix
#'
#' Convenience function to retrieve self-distance matrices computed by computeSelfDistance.
#'
#' @param object A CoPro object
#' @param cellType Cell type name
#' @param slide Slide ID (for CoProMulti objects)
#' @param verbose Whether to print error messages
#'
#' @return Self-distance matrix for the specified cell type
#' @export
getSelfDistMat <- function(object, cellType, slide = NULL, verbose = TRUE) {
  if (inherits(object, "CoProMulti") && is.null(slide)) {
    stop("slide parameter is required for CoProMulti objects")
  }
  
  flat_name <- .createDistMatrixName(cellType, cellType, slide = slide)
  
  if (!flat_name %in% names(object@distances)) {
    msg <- if (is.null(slide)) {
      paste("Self-distance matrix not found for cell type:", cellType,
            ". Run computeSelfDistance() first.")
    } else {
      paste("Self-distance matrix not found for cell type:", cellType,
            "in slide:", slide, ". Run computeSelfDistance() first.")
    }
    if (verbose) stop(msg)
    warning(msg)
    return(NULL)
  }

  return(object@distances[[flat_name]])
}

#' Get Self-Kernel Matrix
#'
#' Convenience function to retrieve self-kernel matrices computed by computeSelfKernel.
#'
#' @param object A CoPro object
#' @param sigma Sigma value
#' @param cellType Cell type name
#' @param slide Slide ID (for CoProMulti objects)
#' @param verbose Whether to print error messages
#' @param materialize For encoded float32 kernels, whether to return a
#'   temporary standard double-precision sparse matrix.
#'
#' @return Self-kernel matrix for the specified parameters
#' @export
getSelfKernelMatrix <- function(object, sigma, cellType, slide = NULL,
                                verbose = TRUE, materialize = TRUE) {
  if (inherits(object, "CoProMulti") && is.null(slide)) {
    stop("slide parameter is required for CoProMulti objects")
  }
  
  flat_name <- .createKernelMatrixName(sigma, cellType, cellType, slide = slide)
  
  if (!flat_name %in% names(object@kernelMatrices)) {
    msg <- if (is.null(slide)) {
      paste("Self-kernel matrix not found for cell type:", cellType,
            "with sigma:", sigma, ". Run computeSelfKernel() first.")
    } else {
      paste("Self-kernel matrix not found for cell type:", cellType,
            "with sigma:", sigma, "in slide:", slide, ". Run computeSelfKernel() first.")
    }
    if (verbose) stop(msg)
    warning(msg)
    return(NULL)
  }

  kernel <- object@kernelMatrices[[flat_name]]
  if (materialize && .isFloat32SparseKernel(kernel)) {
    kernel <- asDoubleSparseMatrix(kernel)
  }
  kernel
}
