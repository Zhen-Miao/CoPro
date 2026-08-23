## get PC matrices
.getAllPCMats <- function(allPCs, scalePCs) {

  if (length(allPCs) == 0) {
    stop("PCA results do not exist, run computePCA() first.")
  }

  PCmats <- setNames(
    vector("list", length = length(allPCs)),
    names(allPCs)
  )

  ## optionally, scale the PCs before running CCA
  if (scalePCs) {
    for (i in names(allPCs)) {
      pca_A_sd <- allPCs[[i]]$sdev
      PCmats[[i]] <- scale(allPCs[[i]]$x,
                           center = FALSE,
                           scale = pca_A_sd
      )
    }
  } else {
    for (i in names(allPCs)) {
      PCmats[[i]] <- allPCs[[i]]$x
    }
  }
  return(PCmats)
}


#' Fraction of nonzero entries per column
#'
#' `colSums(x != 0)` builds a full logical copy of the matrix -- 100 MB for a
#' 50k x 500 dense matrix, and an entire second sparse matrix for a
#' `dgCMatrix`. A compressed-column matrix already stores the count: the gaps
#' in its column pointer. `drop0()` first, because an explicitly stored zero
#' would otherwise be counted as a nonzero.
#'
#' A symmetric compressed-column matrix stores only one triangle, so its column
#' pointer counts roughly half the nonzeros. No caller passes one today; the
#' guard is here so that a future one falls back to the slow-but-correct path
#' instead of silently undercounting.
#'
#' `IterableMatrix` supports no comparison operator at all -- `x != 0` errors
#' with "comparison (!=) is possible only for atomic and list types" -- and the
#' two obvious BPCells substitutes are both wrong. `binarize(x)` is `x > 0`, so
#' it drops negative entries; `matrix_stats(col_stats = "nonzero")` counts
#' *stored* entries, so an explicitly stored zero counts as a nonzero (the
#' `drop0()` problem again). `binarize(x^2)` is `x^2 > 0`, which is `x != 0`
#' for either sign and every storage type (`^` evaluates in double, so integer
#' matrices cannot wrap), streams the matrix once, and never materializes a
#' logical copy. Only magnitudes below ~1.5e-162 square to zero, far outside
#' anything normalized expression produces.
#' @noRd
.columnNonzeroFraction <- function(x) {
  if (inherits(x, "CsparseMatrix") && !inherits(x, "symmetricMatrix")) {
    return(diff(Matrix::drop0(x)@p) / nrow(x))
  }
  if (.is_bpcells(x)) {
    return(colSums(BPCells::binarize(x^2)) / nrow(x))
  }
  colSums(x != 0) / nrow(x)
}

#' Column variances of a BPCells matrix without `BPCells::colVars()`
#'
#' `BPCells::colVars()` is an S3 generic whose `IterableMatrix` method is
#' exported but not registered, so a call from another namespace dispatches
#' only when BPCells is attached or when the Bioconductor package
#' MatrixGenerics is installed (BPCells then re-registers the method on that
#' S4 generic at load time). A library with neither -- a CI runner holding
#' hard dependencies only, or a user who never ran `library(BPCells)` --
#' fails with "no applicable method for 'colVars'". `matrix_stats()` is the
#' plain exported function the method itself calls, so this is the identical
#' streamed computation (sample variance with an `n - 1` denominator, `NaN`
#' for a single row) without the dispatch hazard. Names are kept so callers
#' see exactly what `colVars()` returned.
#' @noRd
.bpcellsColumnVariances <- function(x) {
  BPCells::matrix_stats(x, col_stats = "variance")$col_stats["variance", ]
}

#' Which columns are numerically unsafe to divide by
#'
#' The one degeneracy rule shared by every route into PCA and by the frozen
#' score reference: a column's scale cannot serve as a divisor when it is not
#' finite, is below `zero_sd_threshold`, or belongs to a gene detected in fewer
#' than `nz_proportion_threshold` of the cells. Callers pin the scale of every
#' flagged column at 1, which carries a degenerate gene through unscaled
#' rather than amplifying its noise or turning the column into `NaN`.
#'
#' The rule is written as the conjunction of the *safe* conditions and then
#' negated so that a missing value on either input comes out `TRUE`, never
#' `NA`: `NA | FALSE` is `NA`, and `x[NA] <- 1` silently skips that element,
#' which is how a per-site `which()` or `|` chain could let a column through.
#' @noRd
.unsafeScaleColumns <- function(scale_values, nonzero_fraction,
                                zero_sd_threshold = 1e-3,
                                nz_proportion_threshold = 0.01) {
  safe <- is.finite(scale_values) &
    scale_values >= zero_sd_threshold &
    is.finite(nonzero_fraction) &
    nonzero_fraction >= nz_proportion_threshold
  !safe
}

#' Per-column standard deviations
#'
#' `apply(x, 2, sd)` pays an R-level closure call per column;
#' `matrixStats::colSds()` is ~2.5x faster on a dense numeric matrix. It only
#' accepts a base matrix, so `dgCMatrix` and anything else keeps the `apply()`
#' path.
#'
#' The two use different variance algorithms and can disagree by 1 ulp
#' (observed: 1.11e-16 on 1 of 60 genes in the test fixture). That is enough to
#' flip the sign of a principal component coming out of `prcomp_irlba()`, which
#' sounds worse than it is: the CCA weight's coordinate on that PC flips with
#' it and the two cancel. Cell scores, gene scores, regression gene weights,
#' normalized correlations and the selected sigma are all invariant to it. A PC
#' or CCA axis sign is knife-edge under any implementation -- a different BLAS
#' or R build moves it too -- so the guarantee worth holding is the invariance,
#' not bit-identity. `test-pca-workflow.R` asserts both that invariance and
#' that this function is what actually ships.
#'
#' Names are set from `colnames()` rather than taken from `colSds()`, whose
#' `useNames` default has changed across matrixStats releases; `apply()` always
#' names its result, and the two paths must not differ in that.
#' @importFrom matrixStats colSds
#' @noRd
.columnSds <- function(x) {
  if (is.matrix(x) && is.numeric(x)) {
    sds <- matrixStats::colSds(x)
    names(sds) <- colnames(x)
    return(sds)
  }
  apply(x, 2, sd)
}

#' Guarded root-mean-square column scales for uncentered scaling
#'
#' `scale(x, center = FALSE, scale = TRUE)` divides every column by
#' `sqrt(colSums(x^2) / (n - 1))` and does nothing about a zero divisor: an
#' all-zero gene yields a column of `NaN`s, and a near-constant one gets its
#' noise amplified. This is the uncentered twin of the guard in
#' [center_scale_matrix_opt()] and `.sparse_pca_parameters()` -- same two
#' thresholds, same response of pinning the divisor at 1 -- so a degenerate
#' gene is carried through unscaled on every route into PCA rather than only
#' on the centered and the sparse ones.
#'
#' The divisor is the root-mean-square, not the standard deviation, because
#' that is what `scale()` and the uncentered sparse branch already use; only
#' guarded columns change value.
#' @noRd
.uncenteredColumnScales <- function(x,
                                    zero_sd_threshold = 1e-3,
                                    nz_proportion_threshold = 0.01) {
  n <- nrow(x)
  col_nz <- .columnNonzeroFraction(x)
  col_scales <- sqrt(pmax(as.numeric(colSums(x^2)), 0) / max(1, n - 1L))
  names(col_scales) <- colnames(x)
  unsafe <- .unsafeScaleColumns(col_scales, col_nz,
                                zero_sd_threshold, nz_proportion_threshold)
  col_scales[unsafe] <- 1.0
  col_scales
}

#' centering and scaling the matrix
#' @importFrom stats sd
#' @param matrix Input matrix to be column-centered
#'
#' @return centered and scaled matrix
#' @noRd
center_scale_matrix_opt <- function(input_matrix,
                                    zero_sd_threshold = 1e-3,
                                    nz_proportion_threshold = 0.01) {

  if (!.is_bpcells(input_matrix)) {
    # Original behavior for base matrix / Matrix::dgCMatrix
    col_means <- colMeans(input_matrix)
    col_sds <- .columnSds(input_matrix)
    col_nz <- .columnNonzeroFraction(input_matrix)

    unsafe <- .unsafeScaleColumns(col_sds, col_nz,
                                  zero_sd_threshold, nz_proportion_threshold)
    col_sds_safe <- col_sds
    col_sds_safe[unsafe] <- 1.0

    return(scale(input_matrix, center = col_means, scale = col_sds_safe))
  }

  # ---- BPCells path ----

  col_means <- colMeans(input_matrix)
  col_sds <- sqrt(.bpcellsColumnVariances(input_matrix))
  col_nz <- .columnNonzeroFraction(input_matrix)
  unsafe <- .unsafeScaleColumns(col_sds, col_nz,
                                zero_sd_threshold, nz_proportion_threshold)
  col_sds_safe <- col_sds
  col_sds_safe[unsafe] <- 1.0

  # Center then scale using BPCells broadcasting (no base::scale())
  centered <- BPCells::add_cols(input_matrix, -col_means)
  scaled   <- BPCells::multiply_cols(centered, 1 / col_sds_safe)

  scaled
}



#' Normalize a vector to unit length
#' @param v Input vector or matrix (if matrix, treated as column vector)
#' @return Normalized vector as column matrix
#' @noRd
normalize_vec <- function(v) {
  v_norm <- sqrt(sum(v^2))

  if (v_norm < 1e-12) {
    warning("Vector has very small norm, may cause numerical issues")
    return(matrix(0, nrow = length(v), ncol = 1))
  }

  normalized <- v / v_norm
  if (is.matrix(v)) {
    return(normalized)
  } else {
    return(matrix(normalized, ncol = 1))
  }
}

#' Normalize w under the CCA constraint w'(X'X)w = 1
#'
#' When \code{sdev2} is \code{NULL} (whitened PCs), this reduces to
#' \code{||w|| = 1}. Otherwise the weighted norm
#' \code{sqrt(sum(w^2 * sdev2))} is used so that the constraint becomes
#' \code{w' diag(sdev^2) w = 1}, which is equivalent to \code{||Xw|| = 1}
#' in the original (unwhitened) PC space.
#'
#' @param v Numeric vector or single-column matrix
#' @param sdev2 Squared sdev vector (same length as v), or NULL for standard norm
#' @return Normalized vector as column matrix
#' @noRd
normalize_vec_weighted <- function(v, sdev2 = NULL) {
  if (is.null(sdev2)) {
    return(normalize_vec(v))
  }
  .validateSdev2(sdev2, length(v))
  v_norm <- sqrt(sum(as.numeric(v)^2 * sdev2))
  if (v_norm < 1e-12) {
    warning("Vector has very small weighted norm, may cause numerical issues")
    return(matrix(0, nrow = length(v), ncol = 1))
  }
  normalized <- v / v_norm
  if (is.matrix(v)) {
    return(normalized)
  } else {
    return(matrix(normalized, ncol = 1))
  }
}

#' Apply D-inverse preconditioner to gradient, then normalize under D-norm
#'
#' For the generalized eigenvalue problem \code{Yw = λDw}, the correct power
#' iteration is: (1) compute gradient \code{v = Yw}, (2) apply
#' \code{D^{-1}}: \code{u = v / sdev2}, (3) normalize so that
#' \code{u'Du = 1}. This function performs steps 2--3 in one call.
#' When \code{sdev2} is \code{NULL}, falls back to standard normalization.
#'
#' @param v Gradient vector (single-column matrix)
#' @param sdev2 Squared sdev vector, or NULL for standard norm
#' @return Preconditioned and normalized vector as column matrix
#' @noRd
normalize_gradient_weighted <- function(v, sdev2 = NULL) {
  if (is.null(sdev2)) {
    return(normalize_vec(v))
  }
  .validateSdev2(sdev2, length(v))
  # D^{-1} v
  u <- as.numeric(v) / sdev2
  # Normalize under D-norm: sqrt(u' D u) = sqrt(sum(u^2 * sdev2))
  normalize_vec_weighted(matrix(u, ncol = 1), sdev2)
}

#' Validate a diagonal CCA metric
#' @noRd
.validateSdev2 <- function(sdev2, expected_length, label = "sdev2") {
  if (!is.numeric(sdev2) || length(sdev2) != expected_length ||
      any(!is.finite(sdev2)) || any(sdev2 <= 0)) {
    stop(label, " must contain one finite positive value per feature")
  }
  invisible(sdev2)
}

#' Validate parameters shared by the iterative optimizers
#'
#' Exported low-level optimizers can be called without going through
#' `runSkrCCA()`, so validation belongs at every public entry point rather than
#' only in the high-level pipeline.
#' @noRd
.validateOptimizerParams <- function(max_iter, tol, step_size,
                                     max_iter_name = "max_iter") {
  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1 ||
      max_iter != floor(max_iter)) {
    stop(max_iter_name, " must be a positive integer")
  }
  if (!is.numeric(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite number")
  }
  if (!is.numeric(step_size) || length(step_size) != 1L ||
      !is.finite(step_size) || step_size <= 0 || step_size > 1) {
    stop("step_size must be a single numeric value in (0, 1]")
  }
  invisible(NULL)
}

#' Update the object's usable sigma grid after building kernel blocks
#'
#' A kernel call may request only a subset of an existing bandwidth grid. Such
#' a call must remove bandwidths it proved invalid, but it must not make the
#' unrequested bandwidths unreachable. When no grid has been recorded yet, the
#' valid bandwidths from this call establish it.
#' @param object A `CoPro` object.
#' @param surviving Valid bandwidths requested by the current call.
#' @param invalid Bandwidths requested by the current call that produced no
#'   usable kernel.
#' @return The object with `@sigmaValues` updated.
#' @noRd
.pruneSigmaValues <- function(object, surviving, invalid = numeric()) {
  surviving <- as.numeric(surviving)
  invalid <- as.numeric(invalid)
  if (length(object@sigmaValues) == 0L) {
    object@sigmaValues <- surviving
  } else if (length(invalid) > 0L) {
    object@sigmaValues <- setdiff(object@sigmaValues, invalid)
  }
  object
}

#' Capture and restore the caller's global RNG state
#'
#' The returned value distinguishes an absent `.Random.seed` from an existing
#' seed. Seeded entry points use this with `on.exit()` so reproducibility does
#' not come at the cost of changing the caller's random stream.
#' @noRd
.captureRNGState <- function() {
  exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  list(
    exists = exists,
    seed = if (exists) get(".Random.seed", envir = .GlobalEnv,
                          inherits = FALSE) else NULL
  )
}

#' @noRd
.restoreRNGState <- function(state) {
  if (isTRUE(state$exists)) {
    assign(".Random.seed", state$seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

# Declare globals used across the package to quiet R CMD check NOTES
utils::globalVariables(c(
  ".computeSpatialCrossCorrelation",
  "getKernelMatrix",
  "getSlideList",
  "x",
  "y"
))

#' Check normalizedData for NA / NaN / Inf values
#'
#' Fails fast with a targeted error naming the offending column(s) rather
#' than letting the data flow into PCA or kernel computation where the
#' failure mode is cryptic. BPCells `IterableMatrix` inputs are checked with a
#' streaming column-mean pass. This reads the on-disk values once without
#' materializing the matrix, so the same finite-data contract applies to every
#' storage backend.
#'
#' @param x A matrix, sparse matrix, or BPCells IterableMatrix.
#' @return TRUE invisibly if valid, stops with an informative error otherwise.
#' @noRd
.validateNormalizedData <- function(x) {
  if (inherits(x, "IterableMatrix")) {
    if (nrow(x) == 0L || ncol(x) == 0L) return(invisible(TRUE))
    bad_cols <- which(!is.finite(as.numeric(colMeans(x))))
    if (length(bad_cols) > 0L) {
      preview <- utils::head(bad_cols, 5)
      stop(
        "normalizedData contains NA, NaN, or Inf values in ",
        length(bad_cols), " column(s) (e.g. ",
        paste(preview, collapse = ", "),
        "). Please remove or impute these before constructing a CoPro object.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }

  if (inherits(x, "sparseMatrix")) {
    vals <- x@x
    if (length(vals) == 0) return(invisible(TRUE))
    if (anyNA(vals) || any(!is.finite(vals))) {
      stop(
        "normalizedData contains NA, NaN, or Inf values. ",
        "Please remove or impute these before constructing a CoPro object.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }

  if (is.matrix(x)) {
    bad_cols <- which(colSums(!is.finite(x)) > 0)
    if (length(bad_cols) > 0) {
      preview <- utils::head(bad_cols, 5)
      stop(
        "normalizedData contains NA, NaN, or Inf values in ",
        length(bad_cols), " column(s) (e.g. ",
        paste(preview, collapse = ", "),
        "). Please remove or impute these before constructing a CoPro object.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }

  invisible(TRUE)
}

#' Validate that cell types and slide IDs don't contain pipe characters
#' @param cellTypes Character vector of cell type names
#' @param slideIDs Character vector of slide IDs (optional)
#' @return TRUE if valid, stops execution with error if invalid
#' @noRd
.validateSeparatorSafety <- function(cellTypes = NULL, slideIDs = NULL) {

  # Check cell types for pipe characters
  if (!is.null(cellTypes)) {
    cellTypes <- as.character(cellTypes)
    pipe_in_cellTypes <- grepl("\\|", cellTypes)
    if (any(pipe_in_cellTypes)) {
      problematic_types <- cellTypes[pipe_in_cellTypes]
      stop(paste("Cell type names cannot contain pipe characters (|).",
                 "Problematic cell types:", paste(problematic_types, collapse = ", "),
                 "\nPlease rename these cell types to avoid conflicts with internal naming."))
    }
  }

  # Check slide IDs for pipe characters
  if (!is.null(slideIDs)) {
    slideIDs <- as.character(slideIDs)
    pipe_in_slideIDs <- grepl("\\|", slideIDs)
    if (any(pipe_in_slideIDs)) {
      problematic_slides <- slideIDs[pipe_in_slideIDs]
      stop(paste("Slide IDs cannot contain pipe characters (|).",
                 "Problematic slide IDs:", paste(problematic_slides, collapse = ", "),
                 "\nPlease rename these slide IDs to avoid conflicts with internal naming."))
    }
  }

  return(TRUE)
}



#' Assign distance matrix manually
#'
#' @param object A `CoPro` object
#' @param distanceList A list object that contains all pairwise distances
#' between any two pairs of cells.
#'
#' @return A `CoPro` object with specified
#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
#'
setGeneric("assignDistanceManually",
           function(object,
                    distanceList) standardGeneric("assignDistanceManually")
)


#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
setMethod(
  "assignDistanceManually", "CoPro",
  function(object, distanceList) {
    if (!is.list(distanceList)) {
      stop(paste(
        "distanceList must be a nested list object with names",
        "specified by cell types"
      ))
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

    if (!identical(names(distanceList), cts)) {
      stop(paste(
        "The names of distanceList do not match cell types",
        "of interest"
      ))
    }

    for (i in cts) {
      if (!identical(names(distanceList[[i]]), cts)) {
        stop(paste("The names of distanceList[[", i,
          "]] do not match cell types ",
          "of interest",
          sep = ""
        ))
      }
    }

    object@distances <- distanceList
    return(object)
  }
)

#' Get slide IDs from CoPro object
#'
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @return For CoProMulti: the slideID vector; For CoProSingle: character(0)
#' @export
#' @rdname getSlideID
setGeneric("getSlideID", function(object) standardGeneric("getSlideID"))

#' @rdname getSlideID
setMethod("getSlideID", "CoProSingle", function(object) {
  character(0)
})

#' @rdname getSlideID
setMethod("getSlideID", "CoProMulti", function(object) {
  if ("slideID" %in% colnames(object@metaDataSub)) {
    return(object@metaDataSub$slideID)
  } else {
    return(character(0))
  }
})

#' Get slide list from CoPro object
#'
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @return For CoProMulti: the slideList vector; For CoProSingle: character(0)
#' @export
#' @rdname getSlideList
setGeneric("getSlideList", function(object) standardGeneric("getSlideList"))

#' @rdname getSlideList
setMethod("getSlideList", "CoProSingle", function(object) {
  character(0)
})

#' @rdname getSlideList
setMethod("getSlideList", "CoProMulti", function(object) {
  object@slideList
})

#' Check if object is multi-slide
#'
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @return Logical indicating if object contains multiple slides
#' @export
#' @rdname isMultiSlide
setGeneric("isMultiSlide", function(object) standardGeneric("isMultiSlide"))

#' @rdname isMultiSlide
setMethod("isMultiSlide", "CoProSingle", function(object) {
  FALSE
})

#' @rdname isMultiSlide
setMethod("isMultiSlide", "CoProMulti", function(object) {
  TRUE
})

# ============================================================================
# Per-slide adequacy
# ============================================================================

# Per-slide minimum number of cells per cell type. Below this the per-slide
# second-moment estimate from that slide is too noisy to divide by: the
# covariance has rank <= n-1, so for many features the estimate is dominated
# by sampling noise and the per-slide sigma floor starts to bite.
#
# Shared by both CCA spaces. Gene space applies it in .prepareGeneSpaceData();
# PCA space applies it in .dropDegenerateSlides() on the sumcor route. The
# threshold is a numerical-stability floor, not a statistical one: 10 keeps the
# slide-level second moment well defined while still admitting modest cell
# populations.
.min_cells_per_slide <- 10

#' Drop slides that cannot support a per-slide normalization
#'
#' Only the `sumcor` objectives divide by a per-slide quantity, so only they
#' need this. Under `sumcov` a thin slide simply contributes little to the
#' summed operator and is harmless, which is why the caller there reports
#' rather than drops -- dropping would silently change results computed before
#' this function existed.
#'
#' A slide is kept only when *every* requested cell type clears `minCells` on
#' it. A slide missing one type of a pair contributes nothing to that pair's
#' correlation but would still occupy a slot in the average, so partial slides
#' are dropped rather than partially used.
#'
#' @param X_list_all `X_list_all[[slide]][[cellType]]` cell-by-feature matrices.
#' @param cts Cell types that must all be adequately represented.
#' @param minCells Minimum cells per (slide, cell type).
#' @param nFeatures Optional feature count; when supplied, surviving slides
#'   with `n <= nFeatures` get a rank-deficiency warning. The per-slide Gram
#'   matrix is then singular, which is tolerable -- `sigma` is only zero if the
#'   weight lands in its null space -- but worth surfacing.
#' @param what Label used in messages to name the calling routine.
#' @param verbose Emit the per-slide report.
#' @return A list with `X_list_all` (filtered), `slides` (survivors), and
#'   `dropped`.
#' @noRd
.dropDegenerateSlides <- function(X_list_all, cts, minCells = .min_cells_per_slide,
                                  nFeatures = NULL, what = "sumcor",
                                  verbose = TRUE) {
  slides <- names(X_list_all)
  if (is.null(slides)) {
    stop("X_list_all must be a named list of slides.")
  }

  counts <- lapply(slides, function(s) {
    vapply(cts, function(ct) {
      X <- X_list_all[[s]][[ct]]
      if (is.null(X)) 0L else nrow(X)
    }, integer(1))
  })
  names(counts) <- slides

  keep <- vapply(slides, function(s) all(counts[[s]] >= minCells), logical(1))
  dropped <- slides[!keep]

  for (s in dropped) {
    thin <- cts[counts[[s]] < minCells]
    warning(sprintf(
      "%s: slide '%s' dropped -- cell type(s) %s have fewer than %d cells (%s).",
      what, s, paste(thin, collapse = ", "), minCells,
      paste(sprintf("%s=%d", thin, counts[[s]][thin]), collapse = ", ")
    ), call. = FALSE)
  }

  survivors <- slides[keep]
  if (length(survivors) == 0) {
    stop(sprintf(
      "%s: no slide has at least %d cells for every requested cell type. %s",
      what, minCells,
      "Use objective = \"sumcov\", or analyse the slides separately."
    ))
  }

  if (!is.null(nFeatures)) {
    for (s in survivors) {
      thin <- cts[counts[[s]] <= nFeatures]
      if (length(thin) > 0) {
        warning(sprintf(
          paste0("%s: slide '%s' has n <= nFeatures (%d) for cell type(s) %s, ",
                 "so its per-slide Gram matrix is rank deficient. The ",
                 "per-slide scale is still usable but is estimated from few ",
                 "cells."),
          what, s, nFeatures, paste(thin, collapse = ", ")
        ), call. = FALSE)
      }
    }
  }

  if (verbose && length(dropped) > 0) {
    message(sprintf("  %s: using %d of %d slides.",
                    what, length(survivors), length(slides)))
  }

  list(
    X_list_all = X_list_all[survivors],
    slides = survivors,
    dropped = dropped
  )
}
# Standard deviations below this point are numerically unsafe as divisors in
# score transfer and PCA back-projection. Replacing them by one drops the
# unstable rescaling while retaining the corresponding unscaled coefficient.
.COPRO_SD_GUARD <- 1e-5

#' Replace unsafe standard deviations before division
#' @noRd
.safeStandardDeviations <- function(x, threshold = .COPRO_SD_GUARD) {
  unsafe <- !is.finite(x) | x < threshold
  x[unsafe] <- 1
  x
}
