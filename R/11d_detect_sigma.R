# =============================================================================
# Data-driven sigma selection
# -----------------------------------------------------------------------------
# `sigma` is the Gaussian bandwidth in whatever units the coordinates arrived
# in -- microns, pixels, arbitrary Visium units -- so no single recommended
# value can travel between datasets. What DOES travel is how many neighbors a
# cell ends up talking to. For a cell `a` of type i, the kernel row sum
#
#     m_a(sigma) = sum_{b in type j} exp(-0.5 (d_ab / sigma)^2)
#
# is the effective number of type-j neighbors that cell a is coupled to. It is
# dimensionless, biologically interpretable, and strictly increasing in sigma,
# so "I want 5-20 effective neighbors" inverts to a sigma range on any dataset.
#
# Cost control: m_a is dominated by the nearest neighbors, so it is evaluated
# from a sample of anchor cells against their `nNeighbor` nearest partners
# rather than from all pairs. Truncating at K neighbors underestimates m by
# roughly exp(-K / m) (uniform-density approximation), i.e. < 0.2% at the
# K = 128 default for targets up to m = 20 -- far below the resolution that
# matters here, since downstream bandwidth selection refines the choice anyway.
#
# Anchors are individual cells spread by a deterministic stride, not contiguous
# regions. Both avoid the O(N^2) full computation, but spread anchors sample the
# tissue's density distribution without bias, whereas a handful of patches
# inherits whatever density those patches happened to have. The stride also
# keeps the estimate reproducible, which matters because the same estimator
# sets the distance scale on both the dense and the sparse kernel path.
# =============================================================================

#' Draw anchor cell indices for a neighbor sample
#'
#' Deterministic on purpose. The same estimator sets the distance scale on both
#' the dense and the sparse path, and those two must produce identical kernels;
#' a random subsample would make them differ by however much the two draws
#' differed. A systematic stride also makes every kernel reproducible without
#' the caller having to manage a seed.
#' @noRd
.anchorSample <- function(n, nAnchor) {
  if (n <= nAnchor) return(seq_len(n))
  unique(as.integer(round(seq(1, n, length.out = nAnchor))))
}

#' Bounding-box density of a coordinate matrix (points per unit area/volume)
#' @noRd
.boxDensity <- function(coords) {
  d <- ncol(coords)
  span <- apply(coords, 2L, function(z) diff(range(z)))
  span[!is.finite(span) | span <= 0] <- .Machine$double.eps
  nrow(coords) / prod(span)
}

#' Radius expected to enclose `k` neighbors at density `rho`
#' @noRd
.radiusForCount <- function(k, rho, d) {
  unit_ball <- if (d == 2L) pi else 4 * pi / 3
  (k / (rho * unit_ball))^(1 / d)
}

#' Distances from sampled anchors to their nearest partners
#'
#' Grows a fixed-radius search until the median anchor has `nNeighbor`
#' partners, then keeps the `nNeighbor` smallest distances per anchor.
#'
#' @param A,B coordinate matrices. `B` is the partner set; pass `B = A` with
#'   `within = TRUE` for the same-type case, where a cell is not its own
#'   neighbor.
#' @param anchors integer row indices into `A`.
#' @return numeric matrix `length(anchors)` x `nNeighbor` of ascending
#'   distances, padded with `Inf` where an anchor has fewer partners. `Inf`
#'   contributes exactly zero kernel weight, so padding needs no special
#'   handling downstream.
#' @noRd
.anchorNeighborDistances <- function(A, B, anchors, within, nNeighbor,
                                     maxExpand = 24L) {
  nA <- length(anchors)
  d <- ncol(A)
  out <- matrix(Inf, nrow = nA, ncol = nNeighbor)
  if (nA == 0L || nrow(B) == 0L) return(out)

  query <- A[anchors, , drop = FALSE]
  # The self pair is removed by original index below, so `want` already means
  # the requested number of non-self partners.
  want <- nNeighbor
  rho <- .boxDensity(B)
  r <- .radiusForCount(max(2 * want, 8), rho, d)
  span <- apply(rbind(query, B), 2L, function(z) diff(range(z)))
  r_max <- sqrt(sum(span^2))
  if (!is.finite(r_max) || r_max <= 0) return(out)
  r <- min(max(r, .Machine$double.eps), r_max)

  tri <- NULL
  for (step in seq_len(maxExpand)) {
    tri <- .frnnGrid(query, B, r)
    if (within) {
      # `.frnnGrid` was given B = A explicitly, so j indexes A: a pair is the
      # self pair exactly when the partner row is the anchor's own row.
      self <- anchors[tri$i] == tri$j
      if (any(self)) {
        tri$i <- tri$i[!self]; tri$j <- tri$j[!self]; tri$d <- tri$d[!self]
      }
    }
    counts <- tabulate(tri$i, nbins = nA)
    if (stats::median(counts) >= want || r >= r_max) break
    r <- min(2 * r, r_max)
  }
  if (length(tri$i) == 0L) return(out)

  ord <- order(tri$i, tri$d)
  i <- tri$i[ord]
  dist <- tri$d[ord]
  rank_in_anchor <- sequence(rle(i)$lengths)
  keep <- rank_in_anchor <= nNeighbor
  out[cbind(i[keep], rank_in_anchor[keep])] <- dist[keep]
  out
}

#' Median effective neighbor count implied by a bandwidth
#'
#' The effective count is the Gaussian kernel row sum, i.e. how much total
#' kernel weight a cell receives, not how many cells fall inside some radius.
#' @noRd
.medianEffectiveNeighbors <- function(distances, sigma) {
  if (!is.finite(sigma) || sigma <= 0) return(0)
  stats::median(rowSums(exp(-0.5 * (distances / sigma)^2)))
}

#' Invert the effective-neighbor curve for a bandwidth
#'
#' `.medianEffectiveNeighbors()` is non-decreasing in sigma, so bisection on
#' log(sigma) is unconditionally safe. Returns `NA_real_` when the target lies
#' outside what the truncated neighbor sample can represent.
#' @noRd
.solveSigmaForNeighbors <- function(distances, target, steps = 60L) {
  finite <- distances[is.finite(distances) & distances > 0]
  if (length(finite) == 0L || !is.finite(target) || target <= 0) {
    return(NA_real_)
  }
  lo <- min(finite) / 100
  hi <- max(finite)
  # Expand the bracket rather than assume it: a very clustered or very sparse
  # block can put the solution outside the raw distance range.
  for (i in seq_len(40L)) {
    if (.medianEffectiveNeighbors(distances, lo) <= target) break
    lo <- lo / 4
  }
  for (i in seq_len(40L)) {
    if (.medianEffectiveNeighbors(distances, hi) >= target) break
    hi <- hi * 4
  }
  if (.medianEffectiveNeighbors(distances, hi) < target) return(NA_real_)
  if (.medianEffectiveNeighbors(distances, lo) > target) return(NA_real_)

  log_lo <- log(lo); log_hi <- log(hi)
  for (i in seq_len(steps)) {
    mid <- (log_lo + log_hi) / 2
    if (.medianEffectiveNeighbors(distances, exp(mid)) < target) {
      log_lo <- mid
    } else {
      log_hi <- mid
    }
  }
  exp((log_lo + log_hi) / 2)
}

#' Median nearest-partner distance from an anchor neighbor sample
#'
#' The typical spacing between a cell and its closest partner of the other
#' type. Unlike a quantile of all pairwise distances, this does not depend on
#' the tissue's overall extent or shape.
#' @noRd
.medianNearestSpacing <- function(distances) {
  nearest <- distances[, 1L]
  nearest <- nearest[is.finite(nearest) & nearest > 0]
  if (length(nearest) == 0L) return(NA_real_)
  stats::median(nearest)
}

#' Typical nearest-partner spacing for one block
#'
#' The reference length used by `normalizeMethod = "spacing"`. Estimated from
#' an anchor sample rather than from every pair, and summarized by a median, so
#' it costs O(sample) and is unmoved by a handful of coincident or outlying
#' cells.
#'
#' Contrast with `.lowPercentileBlock()`, which answers a different question --
#' how small the smallest distances get, used to floor them. Sharing one number
#' for both jobs is what made the old scaling sensitive to the densest block in
#' the object.
#' @noRd
.blockNearestSpacing <- function(A, B, within, nAnchor = 500L) {
  partners <- if (within) A else B
  anchors <- .anchorSample(nrow(A), nAnchor)
  distances <- .anchorNeighborDistances(
    A, partners, anchors, within = within, nNeighbor = 8L
  )
  .medianNearestSpacing(distances)
}

#' Block spacing addressed the way the dense distance paths address a block
#'
#' Wraps [.blockNearestSpacing()] with the coordinate lookup, so the dense
#' single- and multi-slide paths reach the same estimator the sparse paths use
#' without each one repeating the `.getCoordinateMatrix()` calls.
#' @noRd
.blockSpacingForBlock <- function(object, distType, xDistScale, yDistScale,
                                  zDistScale, slideID, cellType1, cellType2) {
  A <- .getCoordinateMatrix(object, cellType1, distType,
                            xDistScale, yDistScale, zDistScale, slideID = slideID)
  within <- identical(cellType1, cellType2)
  B <- if (within) {
    A
  } else {
    .getCoordinateMatrix(object, cellType2, distType,
                         xDistScale, yDistScale, zDistScale, slideID = slideID)
  }
  .blockNearestSpacing(A, B, within = within)
}

#' Combine per-block reference distances into one scaling reference
#'
#' `"spacing"` takes the median across blocks: every cell-type pair gets a
#' say, and one unusually dense pair cannot compress the whole object.
#' `"percentile"` takes the minimum, reproducing the pre-1.2.0 behavior.
#' `"global"` has only one value per slide to combine, and combines it the same
#' way `"spacing"` does.
#' @noRd
.combineDistanceReference <- function(values, method) {
  values <- values[is.finite(values) & values > 0]
  if (length(values) == 0L) return(NA_real_)
  if (identical(method, "percentile")) min(values) else stats::median(values)
}

#' Typical spacing between cells, ignoring cell type: the `"global"` reference
#'
#' Every other reference is measured on the blocks a step happens to build, and
#' that is what made the unit depend on the step. A cross-type block measures
#' the gap between two compartments, so it moves with colocalization; a
#' within-type block measures the packing of one type, so it moves with that
#' type's abundance -- in 2-D, a type with a tenth of the cells sits about
#' `sqrt(10)` further from its own nearest neighbor. Either way, two steps
#' looking at different blocks derive different units for the same tissue.
#'
#' Pooling all cells of interest removes both dependencies: the reference is a
#' property of the tissue, so cross-type and within-type distances agree by
#' construction, in any order, with no pin needed. It answers "how far apart
#' are cells here", which is what a length unit should mean.
#'
#' What it deliberately does not do is equalize effective neighbor counts
#' across types -- a rare type still couples to fewer partners at a given
#' bandwidth. That is a question about sigma, not about units, and
#' [detectSigmaRange()] answers it per block.
#'
#' @param slideID Restrict to one slide. `NULL` on a multi-slide object
#'   measures each slide and takes the median, so one unusually dense section
#'   cannot set the unit for the rest.
#' @noRd
.globalSpacingReference <- function(object, distType, xDistScale = 1,
                                    yDistScale = 1, zDistScale = 1,
                                    slideID = NULL) {
  slides <- if (!is.null(slideID)) {
    as.list(slideID)
  } else if (isMultiSlide(object)) {
    as.list(getSlideList(object))
  } else {
    list(NULL)
  }
  values <- vapply(slides, function(s) {
    coords <- .getCoordinateMatrix(object, cellType = NULL, distType,
                                   xDistScale, yDistScale, zDistScale,
                                   slideID = s)
    if (nrow(coords) < 2L) return(NA_real_)
    .blockNearestSpacing(coords, coords, within = TRUE)
  }, numeric(1))
  .combineDistanceReference(values, "global")
}

#' Turn a combined reference distance into a multiplicative scale factor
#' @noRd
.distanceScaleFactor <- function(reference, normalizeTarget, method) {
  if (!is.finite(reference) || reference <= 0) {
    stop(sprintf(
      paste0("Cannot compute distance normalization: no valid %s reference ",
             "across cell-type blocks. Check for degenerate or coincident ",
             "coordinates, or set normalizeDistance = FALSE."),
      method
    ))
  }
  normalizeTarget / reference
}

#' Predicted kernel density for a block at a given support radius
#'
#' Estimates the local partner density from the sampled K-th neighbor distance
#' (`rho ~ K / unit_ball * d_K^d`) and converts it to the expected fraction of
#' the block that a radius-`R` search retains. Near saturation this
#' over-predicts, which is the safe direction: it biases the dense/sparse
#' decision toward dense, and dense is what is correct there.
#' @noRd
.predictedBlockDensity <- function(distances, radius, nPartner, d) {
  if (!is.finite(radius) || radius <= 0 || nPartner <= 0) return(NA_real_)
  finite <- is.finite(distances)
  per_anchor <- rowSums(finite)
  usable <- per_anchor > 0L
  if (!any(usable)) return(NA_real_)
  # Largest sampled distance per anchor, with the count that produced it.
  d_k <- apply(distances[usable, , drop = FALSE], 1L, function(z) {
    z <- z[is.finite(z)]
    z[length(z)]
  })
  k <- per_anchor[usable]
  expected <- k * (radius / d_k)^d
  min(1, stats::median(expected) / nPartner)
}

#' Predict what a kernel will cost before building it
#'
#' Answers the question `method = "auto"` needs: at the largest requested
#' sigma, what fraction of each block does the fixed-radius search actually
#' retain? A "sparse" kernel that retains most of its block is not sparse --
#' a `dgCMatrix` costs 12 bytes per stored entry against 8 for a dense double,
#' so past roughly two-thirds density sparse storage is strictly worse.
#'
#' Uses a smaller anchor sample than [detectSigmaRange()] because the answer
#' only has to be right to within a factor that changes a storage decision.
#'
#' @return list with `density` (per block), `nnz` (predicted represented
#'   nonzeros across all blocks), `denseEntries`, and `scaleFactor`.
#' @noRd
.kernelStorageProbe <- function(object, cts, is_multi, geometry, sigmaValues,
                                lowerLimit, nAnchor = 128L, nNeighbor = 64L) {
  blocks <- .float32KernelBlocks(object, cts, is_multi)
  samples <- .blockNeighborSamples(object, cts, is_multi, geometry,
                                   nAnchor, nNeighbor, verbose = FALSE)

  spacings <- vapply(samples, function(s) .medianNearestSpacing(s$distances),
                     numeric(1))
  # When distances are rescaled, the support radius in ORIGINAL units shrinks
  # by the same factor. A pinned factor is exact and is used as is; otherwise
  # the probe estimates it, exactly for "global" and approximately for the
  # block-based methods (where the typical spacing stands in for the low
  # percentile). Approximate is fine: this only picks a storage format.
  scale_factor <- if (!isTRUE(geometry$normalizeDistance)) {
    1
  } else {
    pin <- .pinnedScaleFactor(object)
    if (!is.null(pin)) {
      pin$factor
    } else {
      reference <- if (identical(geometry$normalizeMethod, "global")) {
        .globalSpacingReference(object, geometry$distType, geometry$xDistScale,
                                geometry$yDistScale, geometry$zDistScale)
      } else {
        .combineDistanceReference(spacings, geometry$normalizeMethod)
      }
      if (is.finite(reference) && reference > 0) {
        geometry$normalizeTarget / reference
      } else {
        1
      }
    }
  }
  radius <- .kernelSupportMultiplier(lowerLimit) * max(sigmaValues) / scale_factor

  density <- vapply(samples, function(s) {
    partners <- if (identical(s$cellType1, s$cellType2)) s$n1 else s$n2
    .predictedBlockDensity(s$distances, radius, partners, s$dim)
  }, numeric(1))

  sizes <- vapply(blocks, function(b) as.numeric(b$n1) * as.numeric(b$n2),
                  numeric(1))
  usable <- is.finite(density)
  list(
    density = density,
    maxDensity = if (any(usable)) max(density[usable]) else NA_real_,
    nnz = sum(density[usable] * sizes[usable]),
    denseEntries = sum(sizes),
    scaleFactor = scale_factor,
    dim = samples[[1]]$dim
  )
}

#' Neighbor sample for every kernel block of an object
#' @noRd
.blockNeighborSamples <- function(object, cts, is_multi, geometry,
                                  nAnchor, nNeighbor, verbose) {
  blocks <- .float32KernelBlocks(object, cts, is_multi)
  coord_cache <- new.env(parent = emptyenv())
  coords_for <- function(slide, cellType) {
    key <- paste(if (is.null(slide)) "" else slide, cellType, sep = "|")
    if (is.null(coord_cache[[key]])) {
      coord_cache[[key]] <- .getCoordinateMatrix(
        object, cellType, geometry$distType,
        geometry$xDistScale, geometry$yDistScale, geometry$zDistScale,
        slideID = slide
      )
    }
    coord_cache[[key]]
  }

  lapply(blocks, function(block) {
    A <- coords_for(block$slide, block$cellType1)
    B <- if (block$symmetric) A else coords_for(block$slide, block$cellType2)
    anchors <- .anchorSample(nrow(A), nAnchor)
    distances <- .anchorNeighborDistances(
      A, B, anchors, within = block$symmetric, nNeighbor = nNeighbor
    )
    if (verbose) {
      message(sprintf(
        "  %s -> %s%s: %d x %d cells, %d anchors",
        block$cellType1, block$cellType2,
        if (is.null(block$slide)) "" else paste0(" on ", block$slide),
        block$n1, block$n2, length(anchors)
      ))
    }
    c(block, list(distances = distances, nAnchorUsed = length(anchors),
                  dim = ncol(A)))
  })
}

#' Detect a usable sigma range from the data
#'
#' Chooses Gaussian bandwidths by the only quantity that is comparable across
#' datasets: how many neighbors a cell is effectively coupled to. For a cell
#' `a`, the effective neighbor count at bandwidth sigma is the kernel row sum
#' `m_a(sigma) = sum_b exp(-0.5 (d_ab / sigma)^2)`. This function finds the
#' sigma values at which the median cell reaches `minNeighbors` and
#' `maxNeighbors`, for every cell-type pair (and slide), and reports a shared
#' range plus a recommended grid.
#'
#' The returned bandwidths are in the coordinate units of `locationData`, so
#' they can be passed straight to [computeKernelMatrix()] with
#' `normalizeDistance = FALSE`. There is no need to rescale distances first:
#' selecting sigma from the data is what makes the analysis unit-independent.
#'
#' Speed comes from sampling. Only `nAnchor` cells per block are examined, each
#' against its `nNeighbor` nearest partners found by fixed-radius search, so
#' cost scales with the sample rather than with the number of pairs. The
#' estimate is deliberately approximate; the downstream bandwidth selection in
#' [computeNormalizedCorrelation()] refines the choice within the range.
#'
#' @param object A `CoProSingle` or `CoProMulti` object, after [subsetData()].
#' @param minNeighbors,maxNeighbors Target effective neighbor counts bracketing
#'   the useful range. Defaults of 5 and 20 keep the kernel above the point
#'   where cells are effectively isolated and below the point where every pair
#'   is coupled.
#' @param nSigma Number of bandwidths in the recommended grid, spaced
#'   logarithmically across the detected range.
#' @param nAnchor Anchor cells sampled per block.
#' @param nNeighbor Nearest partners retained per anchor. Truncation error in
#'   the effective count is about `exp(-nNeighbor / m)`, so the default of 128
#'   comfortably covers targets up to a few dozen neighbors.
#' @param distType,xDistScale,yDistScale,zDistScale Coordinate geometry. When
#'   omitted, any geometry already recorded on the object is used.
#' Anchors are chosen by a deterministic stride, not at random, so repeated
#' calls on the same object return the same bandwidths and the dense and sparse
#' kernel paths agree exactly on the distance scale.
#'
#' @param verbose Whether to report per-block progress.
#'
#' @return An object of class `CoProSigmaRange`: a list with `sigmaValues` (the
#'   recommended grid), `sigmaRange` (lower and upper bound), `feasible`
#'   (whether one range satisfies every block), and `blocks` (a per-block
#'   diagnostic `data.frame` with the bandwidths that block would need, its
#'   median nearest-partner spacing, and its effective neighbor count at the
#'   recommended bounds).
#' @family spatial-pipeline
#' @seealso [computeKernelMatrix()], [computeSparseKernelFloat32()]
#' @examples
#' toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
#' obj <- newCoProSingle(
#'   normalizedData = toy$normalizedData,
#'   locationData   = toy$locationData,
#'   metaData       = toy$metaData,
#'   cellTypes      = toy$cellTypes
#' )
#' obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
#' rng <- detectSigmaRange(obj, verbose = FALSE)
#' rng$sigmaValues
#' @export
#' @rdname detectSigmaRange
setGeneric(
  "detectSigmaRange",
  function(object, minNeighbors = 5, maxNeighbors = 20, nSigma = 5L,
           nAnchor = 500L, nNeighbor = 128L,
           distType = NULL, xDistScale = NULL, yDistScale = NULL,
           zDistScale = NULL, verbose = TRUE) {
    standardGeneric("detectSigmaRange")
  }
)

#' @noRd
.detectSigmaRangeCore <- function(object, minNeighbors, maxNeighbors, nSigma,
                                  nAnchor, nNeighbor, distType,
                                  xDistScale, yDistScale, zDistScale,
                                  verbose, is_multi) {
  if (!is.numeric(minNeighbors) || length(minNeighbors) != 1L ||
      !is.finite(minNeighbors) || minNeighbors <= 0) {
    stop("minNeighbors must be a positive finite scalar.")
  }
  if (!is.numeric(maxNeighbors) || length(maxNeighbors) != 1L ||
      !is.finite(maxNeighbors) || maxNeighbors <= minNeighbors) {
    stop("maxNeighbors must be a finite scalar greater than minNeighbors.")
  }
  nSigma <- as.integer(nSigma)
  nAnchor <- as.integer(nAnchor)
  nNeighbor <- as.integer(nNeighbor)
  if (is.na(nSigma) || nSigma < 1L) stop("nSigma must be a positive integer.")
  if (is.na(nAnchor) || nAnchor < 1L) stop("nAnchor must be a positive integer.")
  if (is.na(nNeighbor) || nNeighbor < 2L) {
    stop("nNeighbor must be an integer of at least 2.")
  }
  if (maxNeighbors > nNeighbor / 4) {
    warning(sprintf(
      paste0("maxNeighbors (%g) is large relative to nNeighbor (%d); the ",
             "truncated neighbor sample may underestimate the upper ",
             "bandwidth. Raise nNeighbor to at least %d."),
      maxNeighbors, nNeighbor, as.integer(4 * maxNeighbors)
    ))
  }

  geometry <- .resolveDistanceGeometry(
    object,
    requested = list(distType = distType, xDistScale = xDistScale,
                     yDistScale = yDistScale, zDistScale = zDistScale,
                     normalizeDistance = NULL, normalizeMethod = NULL, normalizeTarget = NULL,
                     truncateLowDist = NULL),
    what = "detectSigmaRange", verbose = verbose
  )
  if (identical(geometry$distType, "Morphology-Aware")) {
    stop(paste("detectSigmaRange() needs Euclidean coordinates.",
               "Pass distType = 'Euclidean2D' or 'Euclidean3D'."))
  }
  if (isTRUE(geometry$normalizeDistance)) {
    # The two numbers are separate on purpose: the scale factor fixes the unit
    # distances are reported in, sigma fixes how far a cell reaches. Reporting
    # the factor lets a caller convert between them instead of guessing.
    pin <- .pinnedScaleFactor(object)
    warning(paste(
      "This object records normalizeDistance = TRUE, but detectSigmaRange()",
      "reports bandwidths in raw coordinate units. Build kernels with",
      "normalizeDistance = FALSE to use them directly,",
      if (is.null(pin)) {
        "or multiply them by @distanceScaleFactor."
      } else {
        sprintf("or multiply them by the recorded scale factor (%g).",
                pin$factor)
      }
    ))
  }

  cts <- if (length(object@cellTypesOfInterest) != 0L) {
    object@cellTypesOfInterest
  } else {
    unique(object@cellTypesSub)
  }
  if (length(cts) == 0L) {
    stop("No cell types available. Run subsetData() first.")
  }

  if (verbose) {
    message(sprintf("Sampling neighbors for sigma detection (%s):",
                    geometry$distType))
  }
  samples <- .blockNeighborSamples(object, cts, is_multi, geometry,
                                   nAnchor, nNeighbor, verbose)

  rows <- lapply(samples, function(s) {
    sigma_lo <- .solveSigmaForNeighbors(s$distances, minNeighbors)
    sigma_hi <- .solveSigmaForNeighbors(s$distances, maxNeighbors)
    data.frame(
      slide = if (is.null(s$slide)) NA_character_ else s$slide,
      cellType1 = s$cellType1, cellType2 = s$cellType2,
      n1 = s$n1, n2 = s$n2, nAnchor = s$nAnchorUsed,
      medianSpacing = .medianNearestSpacing(s$distances),
      sigmaLo = sigma_lo, sigmaHi = sigma_hi,
      stringsAsFactors = FALSE
    )
  })
  block_table <- do.call(rbind, rows)

  resolved <- stats::complete.cases(block_table[, c("sigmaLo", "sigmaHi")])
  if (!any(resolved)) {
    stop(paste("Could not resolve a bandwidth for any cell-type block.",
               "Check for degenerate or duplicated coordinates."))
  }
  if (!all(resolved)) {
    warning(sprintf(
      "%d of %d blocks had too few sampled neighbors to resolve a bandwidth.",
      sum(!resolved), nrow(block_table)
    ))
  }

  lower <- max(block_table$sigmaLo[resolved])
  upper <- min(block_table$sigmaHi[resolved])
  feasible <- lower <= upper
  if (!feasible) {
    # No single bandwidth puts every block inside the target band, which is
    # what happens when a rare, scattered type sits alongside a dense one.
    # Falling back to the median block keeps a usable answer and names the
    # blocks that will sit outside it.
    lower <- stats::median(block_table$sigmaLo[resolved])
    upper <- stats::median(block_table$sigmaHi[resolved])
    if (lower > upper) {
      mid <- sqrt(lower * upper)
      lower <- mid; upper <- mid
    }
    starved <- block_table$sigmaLo[resolved] > upper
    saturated <- block_table$sigmaHi[resolved] < lower
    warning(sprintf(
      paste0("No single sigma gives every block %g-%g effective neighbors. ",
             "Reporting the median block's range instead; %d block(s) will ",
             "sit below the target and %d above it. See the `blocks` table."),
      minNeighbors, maxNeighbors, sum(starved), sum(saturated)
    ))
  }

  sigma_values <- if (nSigma == 1L) {
    sqrt(lower * upper)
  } else if (isTRUE(all.equal(lower, upper))) {
    rep(lower, nSigma)
  } else {
    exp(seq(log(lower), log(upper), length.out = nSigma))
  }
  sigma_values <- signif(unique(sigma_values), 3L)

  block_table$neighborsAtLower <- vapply(samples, function(s) {
    .medianEffectiveNeighbors(s$distances, lower)
  }, numeric(1))
  block_table$neighborsAtUpper <- vapply(samples, function(s) {
    .medianEffectiveNeighbors(s$distances, upper)
  }, numeric(1))

  structure(
    list(
      sigmaValues = sigma_values,
      sigmaRange = c(lower = lower, upper = upper),
      feasible = feasible,
      targetNeighbors = c(min = minNeighbors, max = maxNeighbors),
      blocks = block_table,
      geometry = geometry,
      samples = lapply(samples, function(s) {
        list(slide = s$slide, cellType1 = s$cellType1,
             cellType2 = s$cellType2, n1 = s$n1, n2 = s$n2,
             dim = s$dim, distances = s$distances)
      })
    ),
    class = "CoProSigmaRange"
  )
}

#' @rdname detectSigmaRange
#' @export
setMethod(
  "detectSigmaRange", "CoPro",
  function(object, minNeighbors = 5, maxNeighbors = 20, nSigma = 5L,
           nAnchor = 500L, nNeighbor = 128L,
           distType = NULL, xDistScale = NULL, yDistScale = NULL,
           zDistScale = NULL, verbose = TRUE) {
    .detectSigmaRangeCore(object, minNeighbors, maxNeighbors, nSigma,
                          nAnchor, nNeighbor, distType, xDistScale,
                          yDistScale, zDistScale, verbose,
                          is_multi = is(object, "CoProMulti"))
  }
)

#' Print a detected sigma range
#' @param x A `CoProSigmaRange` from [detectSigmaRange()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @method print CoProSigmaRange
#' @keywords internal
#' @export
print.CoProSigmaRange <- function(x, ...) {
  cat("CoPro sigma range (coordinate units)\n")
  cat(sprintf("  target effective neighbors: %g - %g\n",
              x$targetNeighbors[["min"]], x$targetNeighbors[["max"]]))
  cat(sprintf("  detected range: %.4g - %.4g%s\n",
              x$sigmaRange[["lower"]], x$sigmaRange[["upper"]],
              if (x$feasible) "" else "  [no range satisfies every block]"))
  cat(sprintf("  recommended sigmaValues: %s\n",
              paste(format(x$sigmaValues), collapse = ", ")))
  cat(sprintf("  blocks: %d\n", nrow(x$blocks)))
  show_cols <- c("slide", "cellType1", "cellType2", "n1", "n2",
                 "medianSpacing", "sigmaLo", "sigmaHi",
                 "neighborsAtLower", "neighborsAtUpper")
  print(utils::head(x$blocks[, show_cols, drop = FALSE], 10L))
  if (nrow(x$blocks) > 10L) {
    cat(sprintf("  ... %d more block(s)\n", nrow(x$blocks) - 10L))
  }
  invisible(x)
}
