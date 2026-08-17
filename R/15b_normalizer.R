## Whitening operators for the normalized-correlation denominator.
##
## computeNormalizedCorrelation() divides the bilinear statistic T = a' K b by
## the whitened Frobenius norm ||R_x^{1/2} K_c R_y^{1/2}||_F, which is exactly
## the null standard deviation of T when the score vectors have within-type
## covariance proportional to R_x and R_y. Choosing R_x and R_y therefore fixes
## how the criterion behaves across bandwidths, and the choice is not innocuous:
##
##   R = I                 under-counts the noise at large sigma, so the
##                         criterion drifts upward with sigma and the selected
##                         bandwidth rails to the top of the grid.
##   R = K(sigma)          (the matched-bandwidth self-kernel) over-counts it,
##                         and the selected bandwidth rails to the floor.
##   R = exp(-D^2/2l^2)    with l the within-type autocorrelation range of the
##                         scores is the operator the theory actually calls for,
##                         and flattens the null path.
##
## The range l is a property of the score field, not of the kernel bandwidth, so
## it is estimated once per cell type and reused across the sigma grid. It must
## be estimated from all PC score columns rather than from a fitted canonical
## score: the canonical direction is chosen to maximize the very association
## under test, so using it would leak signal into the denominator.
##
## Known limitation, measured rather than assumed. Fitting one range to the
## feature-averaged autocorrelation presumes the PCs share a correlation
## length. In simulation, where they do, this flattens the null path (drift
## 1.21x against 1.19x for the true operator and 3.30x un-whitened). On a real
## targeted panel they do not: per-PC ranges spanned 0.3-4.1 cell-spacings on
## colon MERFISH, the flat average over 40 PCs returned a sub-spacing range,
## and the selected bandwidth moved an order of magnitude across a plausible
## range of l. A variance-weighted or mixture aggregation over per-PC ranges
## would be the principled fix; until then treat the variogram normalizer as a
## diagnostic and use a permutation null when the selected bandwidth matters.

.NORMALIZER_MODES <- c("legacy", "unwhitened", "kernel", "variogram")

#' Default tuning for whitening-operator construction
#'
#' `distType` and the axis scales are not stored on the object, so they are
#' repeated here; they must match what was passed to [computeDistance()] or the
#' whitening operator will live in different units from the kernel.
#' @noRd
.normalizerControlDefaults <- function() {
  list(
    distType = "Euclidean2D",
    xDistScale = 1,
    yDistScale = 1,
    zDistScale = 1,
    ## variogram estimation
    maxCells = 2000L,       # subsample size for the autocorrelation estimate
    nBins = 25L,            # distance bins
    maxLagQuantile = 0.05,  # initial upper lag, as a quantile of pair distances
    minCorrelation = 0.05,  # bins below this are past the usable decay
    minBins = 5L,           # refine the lag window until this many bins survive
    ## whitening operator assembly
    lowerLimit = 1e-4,      # truncation of exp(-d^2/2l^2)
    range = NULL            # named numeric: supply ranges instead of estimating
  )
}

#' @noRd
.resolveNormalizerControl <- function(control) {
  defaults <- .normalizerControlDefaults()
  if (is.null(control) || length(control) == 0L) return(defaults)
  if (!is.list(control)) stop("`normalizerControl` must be a list.")
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    stop("Unknown `normalizerControl` entries: ",
         paste(unknown, collapse = ", "), ". Valid entries: ",
         paste(names(defaults), collapse = ", "), ".")
  }
  utils::modifyList(defaults, control)
}

#' Feature-averaged spatial autocorrelation range of a score matrix
#'
#' Bins the cell-by-cell score correlation by spatial lag and fits the Gaussian
#' form `r(d) = exp(-d^2 / 2 l^2)` by regressing `log r` on `d^2` through the
#' origin. Averaging over all score columns is what keeps the estimate free of
#' the canonical direction, and hence free of the association being tested.
#'
#' The lag window is refined downward until at least `minBins` bins sit above
#' `minCorrelation`, so a range much shorter than the domain is still resolved.
#'
#' @param scores Cell-by-PC score matrix, rows aligned with `coords`.
#' @param coords Cell-by-dimension coordinate matrix, in the same units as the
#'   distances the kernel was built from.
#' @return The fitted range, or `NA_real_` if the decay could not be fitted.
#' @importFrom stats dist lm quantile sd
#' @noRd
.estimateSpatialRange <- function(scores, coords, control) {
  n <- nrow(scores)
  if (n < 30L || ncol(scores) < 2L) return(NA_real_)
  if (nrow(coords) != n) {
    stop("Score matrix and coordinates disagree on cell count (",
         n, " vs ", nrow(coords), ").")
  }

  idx <- if (n > control$maxCells) {
    sort(sample.int(n, control$maxCells))
  } else {
    seq_len(n)
  }
  X <- scores[idx, , drop = FALSE]
  ## standardize each feature so the cross-product is a correlation across
  ## features; a constant column would otherwise divide by zero
  keep_col <- apply(X, 2, stats::sd) > 0
  if (sum(keep_col) < 2L) return(NA_real_)
  X <- scale(X[, keep_col, drop = FALSE])

  ## `dist()` enumerates pairs in the same order as the column-major lower
  ## triangle, so the two vectors below are aligned pair for pair.
  C <- tcrossprod(X) / ncol(X)
  scale0 <- mean(diag(C))
  if (!is.finite(scale0) || scale0 <= 0) return(NA_real_)
  lt <- lower.tri(C)
  cc <- C[lt] / scale0
  dd <- as.numeric(stats::dist(coords[idx, , drop = FALSE]))

  positive <- dd > 0
  if (!any(positive)) return(NA_real_)
  cc <- cc[positive]
  dd <- dd[positive]

  maxd <- stats::quantile(dd, control$maxLagQuantile, names = FALSE)
  for (attempt in seq_len(8L)) {
    if (!is.finite(maxd) || maxd <= 0) return(NA_real_)
    breaks <- seq(0, maxd, length.out = control$nBins + 1L)
    inwin <- dd <= maxd
    if (sum(inwin) < 50L) return(NA_real_)
    bin <- .bincode(dd[inwin], breaks, include.lowest = TRUE)
    ac <- tapply(cc[inwin], bin, mean)
    mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
    binid <- as.integer(names(ac))
    ok <- is.finite(ac) & ac >= control$minCorrelation
    if (sum(ok) >= control$minBins) {
      fit <- stats::lm(log(ac[ok]) ~ 0 + I(mid[binid[ok]]^2))
      slope <- unname(stats::coef(fit)[1])
      if (!is.finite(slope) || slope >= 0) return(NA_real_)
      range <- sqrt(-1 / (2 * slope))
      ## a range at or beyond the fitting window is not identified by it
      if (!is.finite(range) || range <= 0 || range > maxd) return(NA_real_)
      return(range)
    }
    ## too few informative bins: the decay happens inside the first bin or two,
    ## so tighten the window and try again
    maxd <- maxd / 2
  }
  NA_real_
}

#' Sparse Gaussian correlation operator at a fitted range
#'
#' Built through the same fixed-radius machinery as the kernels themselves, so
#' it stays sparse at realistic cell counts. `.frnnGrid()` omits self-pairs, so
#' the unit diagonal is added back explicitly -- `R` is a correlation operator
#' and must have `R_ii = 1`.
#'
#' @param coords Raw (unscaled) coordinates.
#' @param range Autocorrelation range, in normalized distance units.
#' @param scaling_factor Raw-to-normalized distance ratio.
#' @noRd
.buildWhiteningOperator <- function(coords, range, scaling_factor, lowerLimit) {
  n <- nrow(coords)
  if (n == 0L) return(NULL)
  triplets <- .buildBlockTriplets(
    A = coords, B = NULL, percentile = NA_real_,
    scaling_factor = scaling_factor, max_sigma = range,
    lowerLimit = lowerLimit, truncateLowDist = FALSE
  )
  if (is.null(triplets)) return(Matrix::Diagonal(n))
  R <- .sparseKernelFromTriplets(
    triplets, sigma = range, lowerLimit = lowerLimit, symmetric = TRUE
  )
  R + Matrix::Diagonal(n)
}

#' Per-cell-type whitening operators for the variogram normalizer
#'
#' @param scoreMats Named list of cell-by-PC score matrices for one slide.
#' @return Named list of operators, with attribute `"ranges"`.
#' @noRd
.variogramOperators <- function(object, scoreMats, cts, control, slide = NULL) {
  scaling_factor <- if (length(object@distanceScaleFactor) == 1L) {
    object@distanceScaleFactor
  } else {
    warning("`distanceScaleFactor` is not set on this object (it predates the ",
            "slot, or computeDistance() was never run); assuming 1. If the ",
            "distances were normalized, the whitening operator will be in the ",
            "wrong units -- rerun computeDistance().")
    1
  }

  operators <- setNames(vector("list", length(cts)), cts)
  ranges <- setNames(rep(NA_real_, length(cts)), cts)

  for (ct in cts) {
    scores <- scoreMats[[ct]]
    if (is.null(scores) || nrow(scores) == 0L) next
    coords <- .getCoordinateMatrix(
      object, ct, control$distType, control$xDistScale,
      control$yDistScale, control$zDistScale, slideID = slide
    )
    if (nrow(coords) != nrow(scores)) {
      stop("Cell count mismatch for cell type '", ct,
           if (is.null(slide)) "" else paste0("' on slide '", slide),
           "': ", nrow(scores), " score rows vs ", nrow(coords),
           " coordinates. The normalizer cannot be built.")
    }

    supplied <- if (is.null(control$range)) NULL else control$range[[ct]]
    range <- if (!is.null(supplied)) {
      if (!is.numeric(supplied) || length(supplied) != 1L || supplied <= 0) {
        stop("`normalizerControl$range[['", ct, "']]` must be a positive scalar.")
      }
      supplied
    } else {
      ## coordinates are raw; the range is reported in normalized units
      .estimateSpatialRange(scores, coords * scaling_factor, control)
    }

    if (!is.finite(range)) {
      warning("Could not estimate a spatial autocorrelation range for cell ",
              "type '", ct, "'",
              if (is.null(slide)) "" else paste0(" on slide '", slide, "'"),
              "; falling back to no whitening (R = I) for this cell type. ",
              "Supply one through normalizerControl$range if the scores are ",
              "known to be spatially smooth.")
      next
    }
    ranges[[ct]] <- range
    operators[[ct]] <- .buildWhiteningOperator(
      coords, range, scaling_factor, control$lowerLimit
    )
  }

  attr(operators, "ranges") <- ranges
  operators
}

#' Build the resolver that supplies R_x / R_y to .whitenedFrobNorm()
#'
#' Returns a list with `mode`, a `get(sigma, cellType, slide)` accessor, and
#' `ranges` (populated for the variogram normalizer). The accessor returns
#' `NULL` to mean "no whitening", which `.whitenedFrobNorm()` reads as R = I.
#'
#' @param scoreMats For a single-slide object, a named list of score matrices.
#'   For a multi-slide object, a list of such lists keyed by slide.
#' @noRd
.makeWhiteningResolver <- function(object, normalizer, control, scoreMats, cts,
                                   slides = NULL) {
  normalizer <- match.arg(normalizer, .NORMALIZER_MODES)
  control <- .resolveNormalizerControl(control)

  selfKernel <- function(sigma, cellType, slide = NULL) {
    tryCatch(
      getKernelMatrix(object, sigma = sigma, cellType1 = cellType,
                      cellType2 = cellType, slide = slide,
                      verbose = FALSE, materialize = FALSE),
      error = function(e) NULL
    )
  }

  restoreUnitDiagonal <- function(R) {
    if (is.null(R)) return(NULL)
    if (.isFloat32SparseKernel(R)) R <- asDoubleSparseMatrix(R)
    if (inherits(R, "sparseMatrix")) {
      d <- Matrix::diag(R)
      if (any(d != 1)) R <- R + Matrix::Diagonal(nrow(R), 1 - d)
    } else if (inherits(R, "Matrix")) {
      R <- as.matrix(R)
      diag(R) <- 1
    } else if (is.matrix(R)) {
      diag(R) <- 1
    }
    R
  }

  if (normalizer == "unwhitened") {
    return(list(mode = normalizer, ranges = NULL,
                get = function(sigma, cellType, slide = NULL) NULL))
  }

  if (normalizer == "legacy") {
    ## Preserve the historical selection rule (matched-sigma self-kernel when
    ## present, otherwise R = I), but repair the self-kernel's diagonal before
    ## treating it as a correlation operator. The stored analysis kernel keeps
    ## its zero diagonal; only the whitening copy receives R_ii = 1.
    getter <- function(sigma, cellType, slide = NULL) {
      restoreUnitDiagonal(selfKernel(sigma, cellType, slide))
    }
    return(list(mode = normalizer, ranges = NULL, get = getter))
  }

  if (normalizer == "kernel") {
    getter <- function(sigma, cellType, slide = NULL) {
      R <- selfKernel(sigma, cellType, slide)
      if (is.null(R)) {
        stop("normalizer = \"kernel\" needs the matched-sigma within-type ",
             "kernel for cell type '", cellType, "' at sigma = ", sigma,
             if (is.null(slide)) "" else paste0(" on slide '", slide, "'"),
             ", which this object does not contain. Run computeSelfKernel() ",
             "first, or choose normalizer = \"variogram\" (recommended) or ",
             "\"unwhitened\".")
      }
      ## .frnnGrid() and the dense path both leave a zero diagonal on a
      ## self-kernel, but a correlation operator has R_ii = 1. Restore it, so
      ## that the null variance this normalizer claims to compute is the one it
      ## actually computes.
      restoreUnitDiagonal(R)
    }
    return(list(mode = normalizer, ranges = NULL, get = getter))
  }

  ## variogram: operators are sigma-independent, so build them once
  if (is.null(slides)) {
    operators <- .variogramOperators(object, scoreMats, cts, control, slide = NULL)
    ranges <- attr(operators, "ranges")
    getter <- function(sigma, cellType, slide = NULL) operators[[cellType]]
  } else {
    per_slide <- setNames(vector("list", length(slides)), slides)
    ranges <- list()
    for (sID in slides) {
      mats <- scoreMats[[sID]]
      if (is.null(mats)) next
      present <- cts[vapply(cts, function(ct) !is.null(mats[[ct]]), logical(1))]
      if (length(present) == 0L) next
      per_slide[[sID]] <- .variogramOperators(object, mats, present, control,
                                              slide = sID)
      ranges[[sID]] <- attr(per_slide[[sID]], "ranges")
    }
    getter <- function(sigma, cellType, slide = NULL) {
      if (is.null(slide)) return(NULL)
      ops <- per_slide[[slide]]
      if (is.null(ops)) NULL else ops[[cellType]]
    }
  }

  list(mode = normalizer, ranges = ranges, get = getter)
}

#' Warn when a cross-type denominator was built from self-kernels
#'
#' `computeSelfDistance()` derives its own normalization factor from the
#' self-distance percentiles instead of reusing `@distanceScaleFactor`, so a
#' self-kernel indexed by "sigma = s" is generally at a different physical
#' bandwidth from the cross-kernel indexed by the same s. Anything that whitens
#' a *cross*-type kernel with those inherits the mismatch.
#'
#' A within-type analysis is exempt: there the whitening operator and the
#' analysis kernel are the same matrix, built by `computeKernelMatrix()` under
#' one scaling factor, so there is nothing to mismatch. Hence the count passed
#' in is of whitened cross-type pairs, not of whitened pairs.
#' @noRd
.warnSelfKernelUnits <- function(resolver, whitened_cross_pairs) {
  if (whitened_cross_pairs == 0L) return(invisible(NULL))
  if (!resolver$mode %in% c("legacy", "kernel")) return(invisible(NULL))
  warning(
    "The denominator was whitened with matched-sigma self-kernels. ",
    "computeSelfDistance() normalizes self-distances with its own scaling ",
    "factor rather than the one computeDistance() recorded, so those kernels ",
    "are generally at a different physical bandwidth from the cross-kernel at ",
    "the same nominal sigma. Prefer normalizer = \"variogram\".",
    call. = FALSE
  )
  invisible(NULL)
}

#' Human-readable summary of what a resolver actually applied
#' @noRd
.describeNormalizer <- function(resolver, whitened_pairs, total_pairs) {
  switch(
    resolver$mode,
    unwhitened = "unwhitened Frobenius norm ||K_c||_F (R = I)",
    kernel = "whitened Frobenius norm with R = matched-sigma self-kernel",
    variogram = paste0(
      "whitened Frobenius norm with R = exp(-d^2/2l^2), l estimated from the ",
      "PC scores"
    ),
    legacy = if (whitened_pairs == 0L) {
      paste0("unwhitened Frobenius norm ||K_c||_F (R = I) -- no within-type ",
             "kernels found on this object")
    } else if (whitened_pairs == total_pairs) {
      "whitened Frobenius norm with R = matched-sigma self-kernel"
    } else {
      sprintf(paste0("MIXED: %d of %d cell-type pairs whitened by a self-kernel",
                     ", the rest unwhitened. These are different statistics; ",
                     "set `normalizer` explicitly."),
              whitened_pairs, total_pairs)
    }
  )
}

#' Report the normalizer actually used
#'
#' The denominator of the normalized correlation determines how the criterion
#' behaves across the sigma grid, and it used to be decided implicitly by
#' whether the object happened to carry within-type kernels. This records the
#' resolved choice on the object so a stored result can be interpreted later.
#'
#' @param object A `CoPro` or `CoProMulti` object that has been through
#'   [computeNormalizedCorrelation()].
#' @return A list with the normalizer `mode`, a human-readable `description`,
#'   and the fitted autocorrelation `ranges` when `normalizer = "variogram"`
#'   was used. `NULL` if the object predates this record.
#' @family scores-and-correlation
#' @seealso [computeNormalizedCorrelation()]
#' @examples
#' # obj <- computeNormalizedCorrelation(obj, normalizer = "variogram")
#' # getNormalizerInfo(obj)
#' @export
getNormalizerInfo <- function(object) {
  if (!methods::is(object, "CoPro")) {
    stop("getNormalizerInfo() expects a CoPro-derived object.")
  }
  attr(object@normalizedCorrelation, "normalizer", exact = TRUE)
}
