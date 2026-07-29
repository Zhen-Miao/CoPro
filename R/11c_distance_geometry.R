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
  normalizeDistance = FALSE, normalizeMethod = "global",
  normalizeTarget = 0.01, truncateLowDist = TRUE
)

#' Supported distance-normalization strategies
#'
#' `"global"` is the default because it is the only one whose reference does
#' not depend on which blocks a step happens to build. See
#' `.globalSpacingReference()` for why that matters.
#' @noRd
.NORMALIZE_METHODS <- c("global", "spacing", "percentile")

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
                                  normalizeMethod = "global",
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

#' The distance scale factor already pinned on an object, if any
#'
#' Cross-type and within-type blocks reach normalization through different
#' entry points -- [computeDistance()] / [computeKernelMatrix()] for the cross
#' blocks, [computeSelfDistance()] / [computeSelfKernel()] for the self blocks
#' -- and each used to derive its own reference length from its own blocks. The
#' two references differ whenever cell types differ in abundance or in how
#' tightly they colocalize, so an object could end up holding cross distances
#' and self distances on two different units, with `@distanceScaleFactor`
#' describing only one of them.
#'
#' @return `NULL` when no normalization pass has been recorded, otherwise
#'   `list(factor, source)`.
#' @noRd
.pinnedScaleFactor <- function(object) {
  geom <- .getDistanceGeometry(object)
  if (!isTRUE(geom$normalizeDistance)) return(NULL)
  factor <- tryCatch(object@distanceScaleFactor, error = function(e) numeric(0))
  if (length(factor) != 1L || !is.finite(factor) || factor <= 0) return(NULL)
  list(
    factor = factor,
    source = if (is.null(geom$source)) "an earlier step" else geom$source
  )
}

#' Reuse the pinned scale factor so all blocks stay on one unit
#'
#' First normalization wins. A step that adds blocks to an already-normalized
#' object adopts that object's factor rather than deriving a second one, so
#' cross-type and within-type distances remain comparable and
#' `@distanceScaleFactor` describes every block in the object.
#'
#' There is no universally right shared reference -- see the `normalizeMethod`
#' discussion in [computeDistance()] -- so this pins consistency rather than
#' claiming optimality. [computeDistance()] rebuilds `@distances` from scratch
#' and therefore always re-derives; so does any additive step called with
#' `overwrite = TRUE`, which clears the pin via [.clearDistanceRecord()].
#' Those are the two documented ways to re-pin.
#'
#' @param computed The factor this step derived from its own blocks.
#' @return The factor to use.
#' @noRd
.adoptScaleFactor <- function(object, computed, what, verbose = TRUE) {
  pin <- .pinnedScaleFactor(object)
  if (is.null(pin)) return(computed)
  if (verbose && !isTRUE(all.equal(pin$factor, computed))) {
    message(sprintf(
      paste0("%s: using the distance scale factor already set by %s (%g) ",
             "instead of the %g its own blocks imply, so every block in this ",
             "object stays on one unit.\n",
             "  Pass overwrite = TRUE, or re-run computeDistance(), to ",
             "re-derive the scale."),
      what, if (identical(pin$source, "an earlier step")) pin$source
            else paste0(pin$source, "()"),
      pin$factor, computed
    ))
  }
  pin$factor
}

#' The scale factor a normalization pass should use
#'
#' One entry point for every step that rescales distances, so the choice of
#' reference, the conversion to a factor, and the reuse of an already-pinned
#' factor cannot drift apart between the dense, sparse, float32, self, and
#' streaming paths.
#'
#' Under `normalizeMethod = "global"` the reference comes from the cells, not
#' from the blocks this particular step builds, so every step computes the same
#' number and `adopt` never has anything to do. The pin still runs: it is what
#' keeps an object consistent when the record was written under one of the
#' block-based methods.
#'
#' @param blockValues Per-block reference distances gathered by the caller
#'   (spacings or low percentiles). Ignored by `"global"`.
#' @param adopt Whether an already-pinned factor wins. `FALSE` for
#'   [computeDistance()], which rebuilds `@distances` wholesale and is
#'   therefore the step that re-pins.
#' @param slideID Restrict a `"global"` reference to one slide; `NULL` combines
#'   across slides.
#' @return A positive scale factor.
#' @noRd
.normalizationScaleFactor <- function(object, blockValues, normalizeMethod,
                                      normalizeTarget, distType,
                                      xDistScale = 1, yDistScale = 1,
                                      zDistScale = 1, what = "computeDistance",
                                      verbose = TRUE, slideID = NULL,
                                      adopt = TRUE) {
  reference <- if (identical(normalizeMethod, "global")) {
    .globalSpacingReference(object, distType, xDistScale, yDistScale,
                            zDistScale, slideID = slideID)
  } else {
    .combineDistanceReference(blockValues, normalizeMethod)
  }
  computed <- .distanceScaleFactor(reference, normalizeTarget, normalizeMethod)
  if (!adopt) return(computed)
  .adoptScaleFactor(object, computed, what = what, verbose = verbose)
}

#' Forget the recorded geometry and pinned scale factor
#'
#' Called by the additive steps when `overwrite = TRUE`: the blocks the record
#' described are about to be discarded, so keeping it would make the step
#' inherit -- and be checked against -- a basis nothing in the object still
#' uses. After this the step re-resolves from its own arguments and re-pins.
#' @noRd
.clearDistanceRecord <- function(object) {
  if ("distanceGeometry" %in% methods::slotNames(object)) {
    object@distanceGeometry <- list()
  }
  object@distanceScaleFactor <- numeric(0)
  object
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
