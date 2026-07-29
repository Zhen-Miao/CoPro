# =============================================================================
# Distance geometry: one recorded description of the coordinates that the
# distances -- and therefore the kernels -- were actually built on.
# -----------------------------------------------------------------------------
# Three code paths can build the kernel's coordinate basis: computeDistance()
# (dense), computeSparseKernel() / computeKernelMatrix(method = "sparse"), and
# the streaming gene-space CCA. Each used to re-derive `distType` and the
# per-axis scales from its own arguments, so they could disagree silently and
# produce kernels on coordinates the caller never asked for.
#
# The fix is a single record written by whichever call actually built the
# coordinates, stored in `@distanceGeometry`. Later calls read it and defer to
# it instead of re-defaulting, and downstream helpers that need the raw ->
# analysis coordinate map (.sigmaAwareBins, .recoverDistanceScaleFactor) read
# it instead of assuming raw, unscaled, 2-D `x,y`.
# =============================================================================

#' Fields that define the coordinate geometry, in report order
#' @noRd
.DISTANCE_GEOMETRY_FIELDS <- c(
  "distType", "xDistScale", "yDistScale", "zDistScale",
  "normalizeDistance", "normalizeMethod", "normalizeTarget", "truncateLowDist"
)

#' Hard defaults used when neither the caller nor the object supplies a value.
#'
#' `distType` is absent here because its default depends on the coordinates
#' (see `.defaultDistType()`).
#'
#' `normalizeDistance` defaults to `FALSE` as of CoPro 1.2.0. Rescaling
#' distances was a way to make one recommended sigma travel between datasets
#' recorded in different units; [detectSigmaRange()] now derives sigma from the
#' data itself, which achieves the same thing without moving the coordinates
#' the results are reported on. See `.noteNormalizeDistanceDefault()`.
#' @noRd
.DISTANCE_GEOMETRY_DEFAULTS <- list(
  xDistScale = 1, yDistScale = 1, zDistScale = 1,
  normalizeDistance = FALSE, normalizeMethod = "spacing",
  normalizeTarget = 0.01, truncateLowDist = TRUE
)

#' Supported distance-normalization strategies
#' @noRd
.NORMALIZE_METHODS <- c("spacing", "percentile")

#' Is this `normalizeDistance` the "reuse the recorded factor" instruction?
#'
#' The self-distance and self-kernel paths accept `normalizeDistance =
#' "inherit"`, which says to reuse the scaling factor [computeDistance()]
#' recorded rather than derive one from the within-type blocks. It describes
#' where the factor comes from, not a different coordinate basis, so it agrees
#' with any record by construction and must never read as a contradiction.
#' @noRd
.isInheritNormalize <- function(value) identical(value, "inherit")

#' The `normalizeDistance` to write into a record when the caller said "inherit"
#'
#' A record is a description of the coordinates the object's matrices live on.
#' Storing the instruction instead of its outcome would leave later steps
#' inheriting `"inherit"` from an object that has nothing left to inherit from.
#' @noRd
.recordedNormalizeDistance <- function(resolved, recorded) {
  if (!.isInheritNormalize(resolved)) return(resolved)
  if (!is.null(recorded[["normalizeDistance"]])) {
    recorded[["normalizeDistance"]]
  } else {
    # Unreachable in practice: inheriting with nothing recorded already errored
    # in .resolveSelfDistanceScaling(). TRUE is the honest description anyway --
    # the blocks did get scaled.
    TRUE
  }
}

#' Validate a normalizeMethod argument
#' @noRd
.checkNormalizeMethod <- function(method) {
  if (is.null(method)) return(NULL)
  if (!is.character(method) || length(method) != 1L ||
      !(method %in% .NORMALIZE_METHODS)) {
    stop(sprintf("normalizeMethod must be one of: %s.",
                 paste(.NORMALIZE_METHODS, collapse = ", ")))
  }
  method
}

#' Announce the changed `normalizeDistance` default once per session
#'
#' The default flipped from `TRUE` to `FALSE` in CoPro 1.2.0, which changes
#' results for anyone who relied on it implicitly. Silence would let an
#' existing script keep running while quietly analyzing different kernels, so
#' the change is announced the first time a step falls back to the default.
#' @noRd
.noteNormalizeDistanceDefault <- local({
  announced <- FALSE
  function(what) {
    if (announced) return(invisible(FALSE))
    announced <<- TRUE
    message(sprintf(
      paste0(
        "%s: normalizeDistance now defaults to FALSE (changed in CoPro ",
        "1.2.0); sigma is interpreted in raw coordinate units.\n",
        "  Use detectSigmaRange() to choose sigma for this dataset, or pass ",
        "normalizeDistance = TRUE to reproduce earlier results."
      ), what
    ))
    invisible(TRUE)
  }
})

#' Distance type implied by the coordinate columns when none is supplied
#' @noRd
.defaultDistType <- function(object) {
  if ("z" %in% tolower(colnames(object@locationDataSub))) {
    "Euclidean3D"
  } else {
    "Euclidean2D"
  }
}

#' Assemble a geometry record
#'
#' @param source Name of the function that built the coordinates. Recorded so
#'   that a mismatch report can name the step the user needs to re-run.
#' @noRd
.makeDistanceGeometry <- function(distType, xDistScale = 1, yDistScale = 1,
                                  zDistScale = 1, normalizeDistance = FALSE,
                                  normalizeMethod = "spacing",
                                  normalizeTarget = 0.01, truncateLowDist = TRUE,
                                  source = NA_character_) {
  list(
    distType = distType,
    xDistScale = xDistScale,
    yDistScale = yDistScale,
    zDistScale = zDistScale,
    normalizeDistance = normalizeDistance,
    normalizeMethod = normalizeMethod,
    normalizeTarget = normalizeTarget,
    truncateLowDist = truncateLowDist,
    source = source
  )
}

#' Read the recorded geometry, tolerating objects built before the slot existed
#'
#' Objects serialized by an earlier CoPro version have no `distanceGeometry`
#' slot, and `@` on a missing slot is an error, so probe `slotNames()` first.
#'
#' @return The recorded geometry list, or `list()` when nothing is recorded.
#' @noRd
.getDistanceGeometry <- function(object) {
  if (!("distanceGeometry" %in% methods::slotNames(object))) {
    return(list())
  }
  geom <- object@distanceGeometry
  if (!is.list(geom) || length(geom) == 0) {
    return(list())
  }
  geom
}

#' Per-axis coordinate scales implied by a geometry record
#'
#' @return Named numeric vector `c(x, y, z)`; `z` is `NA_real_` for a 2-D
#'   geometry, so callers can tell "no z axis" from "z scaled by 1".
#' @noRd
.geometryAxisScales <- function(geom) {
  pick <- function(field) {
    value <- geom[[field]]
    if (is.null(value) || !is.numeric(value) || length(value) != 1 ||
        !is.finite(value) || value <= 0) {
      .DISTANCE_GEOMETRY_DEFAULTS[[field]]
    } else {
      value
    }
  }
  is_3d <- identical(geom$distType, "Euclidean3D")
  c(x = pick("xDistScale"), y = pick("yDistScale"),
    z = if (is_3d) pick("zDistScale") else NA_real_)
}

#' Format one geometry field for a user-facing message
#' @noRd
.formatGeometryValue <- function(value) {
  if (is.null(value)) return("NULL")
  if (is.character(value)) return(sprintf("'%s'", value))
  if (is.logical(value)) return(as.character(value))
  format(value)
}

#' Fields where an explicitly requested value contradicts the recorded one
#'
#' @param recorded Recorded geometry (possibly `list()`).
#' @param requested Named list of caller-supplied values; `NULL` entries mean
#'   "not supplied" and never conflict.
#' @return Character vector of human-readable descriptions, one per conflict.
#' @noRd
.geometryConflicts <- function(recorded, requested) {
  if (length(recorded) == 0) return(character(0))
  conflicts <- character(0)
  for (field in .DISTANCE_GEOMETRY_FIELDS) {
    want <- requested[[field]]
    have <- recorded[[field]]
    if (is.null(want) || is.null(have)) next
    if (identical(field, "normalizeDistance") && .isInheritNormalize(want)) next
    if (!isTRUE(all.equal(want, have))) {
      conflicts <- c(conflicts, sprintf(
        "  %s: requested %s, but distances were built with %s",
        field, .formatGeometryValue(want), .formatGeometryValue(have)
      ))
    }
  }
  conflicts
}

#' Resolve the geometry a kernel step should use
#'
#' Precedence is caller-supplied -> recorded on the object -> package default.
#' A caller-supplied value that contradicts the record is an error rather than
#' a silent override: the recorded geometry is what the distances (and any
#' kernels already built from them) live on, so honoring the argument would
#' produce a mixed-basis object.
#'
#' @param object The CoPro object (supplies the coordinate-dependent
#'   `distType` default).
#' @param requested Named list of caller-supplied values, `NULL` where the
#'   caller did not supply one.
#' @param what Name of the calling step, used in messages.
#' @return A geometry record with every field filled in.
#' @noRd
.resolveDistanceGeometry <- function(object, requested, what = "computeKernelMatrix",
                                     verbose = TRUE) {
  recorded <- .getDistanceGeometry(object)

  conflicts <- .geometryConflicts(recorded, requested)
  if (length(conflicts) > 0) {
    stop(sprintf(
      paste0("%s was given distance arguments that contradict the geometry ",
             "the existing distances were built on by %s():\n%s\n",
             "Kernels must be built on the same coordinates as the distances. ",
             "Either drop the conflicting argument(s) so the recorded geometry ",
             "is inherited, or re-run %s() with the geometry you want."),
      what, if (is.null(recorded$source)) "an earlier step" else recorded$source,
      paste(conflicts, collapse = "\n"),
      if (is.null(recorded$source)) "computeDistance" else recorded$source
    ), call. = FALSE)
  }

  .checkNormalizeMethod(requested[["normalizeMethod"]])

  resolved <- list()
  defaulted <- character(0)
  for (field in .DISTANCE_GEOMETRY_FIELDS) {
    resolved[[field]] <- if (!is.null(requested[[field]])) {
      requested[[field]]
    } else if (!is.null(recorded[[field]])) {
      recorded[[field]]
    } else if (identical(field, "distType")) {
      .defaultDistType(object)
    } else {
      defaulted <- c(defaulted, field)
      .DISTANCE_GEOMETRY_DEFAULTS[[field]]
    }
  }
  if ("normalizeDistance" %in% defaulted) {
    .noteNormalizeDistanceDefault(what)
  }

  inherited <- length(recorded) > 0 &&
    any(vapply(.DISTANCE_GEOMETRY_FIELDS,
               function(field) is.null(requested[[field]]) &&
                 !is.null(recorded[[field]]),
               logical(1)))
  if (verbose && inherited) {
    message(sprintf(
      "%s: inheriting distance geometry from %s (distType = '%s', scales = %g / %g / %g).",
      what, if (is.null(recorded$source)) "the object" else paste0(recorded$source, "()"),
      resolved$distType, resolved$xDistScale, resolved$yDistScale,
      resolved$zDistScale
    ))
  }

  resolved$source <- if (length(recorded) > 0 && !is.null(recorded$source)) {
    recorded$source
  } else {
    what
  }
  resolved
}

#' Warn when a step builds coordinates that disagree with the recorded ones
#'
#' Used by entry points such as [computeSparseKernel()] whose arguments have
#' concrete defaults, so "not supplied" and "supplied the default value" cannot
#' be told apart. Erroring would punish the documented usage (calling it
#' *instead of* [computeDistance()], where nothing is recorded and nothing
#' fires); a warning still removes the silence.
#' @noRd
.warnDistanceGeometryMismatch <- function(object, requested, what) {
  recorded <- .getDistanceGeometry(object)
  conflicts <- .geometryConflicts(recorded, requested)
  if (length(conflicts) == 0) return(invisible(FALSE))
  warning(sprintf(
    paste0("%s is using different coordinates than the geometry recorded by ",
           "%s():\n%s\nThe result will use the arguments given to %s, so this ",
           "object will hold matrices built on two different coordinate bases."),
    what, if (is.null(recorded$source)) "an earlier step" else recorded$source,
    paste(conflicts, collapse = "\n"), what
  ), call. = FALSE)
  invisible(TRUE)
}

#' Get the coordinate geometry a CoPro object's distances and kernels use
#'
#' Returns the record written by whichever step built the object's coordinate
#' basis -- [computeDistance()], [computeSelfDistance()],
#' [computeSparseKernel()], or [computeKernelMatrix()] on its sparse path. Use
#' it whenever an analysis has to reproduce the geometry the kernels were fit
#' on, for instance to rebuild a spatial neighbor graph for a permutation null.
#'
#' The record persists after `computeKernelMatrix(dropDistances = TRUE)` (the
#' default) has cleared `@distances`, so a stored object still knows which
#' coordinates its kernels live on.
#'
#' @param object A `CoPro` object.
#' @return A named list with `distType`, `xDistScale`, `yDistScale`,
#'   `zDistScale`, `normalizeDistance`, `normalizeTarget`, `truncateLowDist`,
#'   and `source` (the function that wrote the record). `NULL` when no
#'   distance- or kernel-building step has run, or when `object` predates the
#'   `distanceGeometry` slot.
#' @family accessors
#' @seealso [computeDistance()], [computeKernelMatrix()]
#' @examples
#' toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
#' obj <- newCoProSingle(
#'   normalizedData = toy$normalizedData,
#'   locationData   = toy$locationData,
#'   metaData       = toy$metaData,
#'   cellTypes      = toy$cellTypes
#' )
#' obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
#' obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
#' getDistanceGeometry(obj)
#' @export
#' @rdname getDistanceGeometry
#' @aliases getDistanceGeometry,CoPro-method
setGeneric("getDistanceGeometry",
           function(object) standardGeneric("getDistanceGeometry"))

#' @rdname getDistanceGeometry
#' @export
setMethod("getDistanceGeometry", "CoPro", function(object) {
  geom <- .getDistanceGeometry(object)
  if (length(geom) == 0) NULL else geom
})
