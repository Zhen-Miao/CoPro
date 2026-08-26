#' Dependency graph for derived CoPro state
#'
#' Each edge points from a pipeline stage to state that becomes stale when that
#' stage changes. Keeping the graph here makes invalidation rules explicit and
#' prevents individual mutators from growing their own incomplete slot lists.
#' @noRd
.coProStateDependencies <- list(
  subset = c("integrated", "pca", "distance"),
  integrated = "pca",
  pca = "cca",
  distance = "kernel",
  self_distance = "self_kernel",
  kernel = c("kernel_cache", "cca"),
  self_kernel = c("kernel_cache", "correlation", "permutation"),
  kernel_cache = character(0),
  cca = c("correlation", "scores", "permutation"),
  additive_cca = c(
    "correlation", "regression_scores", "score_tests", "permutation"
  ),
  correlation = character(0),
  scores = c("regression_scores", "score_tests"),
  regression_scores = character(0),
  score_tests = character(0),
  permutation = character(0)
)

#' Object slots owned by each derived-state stage
#' @noRd
.coProStateSlots <- list(
  integrated = "integratedData",
  pca = c("pcaResults", "pcaGlobal", "nPCA", "scalePCs"),
  distance = c("distances", "distanceScaleFactor", "distanceGeometry"),
  kernel = c("kernelMatrices", "sigmaValues"),
  cca = c("skrCCAOut", "nCC"),
  correlation = c(
    "normalizedCorrelation", "bidirCorrelation", "sigmaValueChoice"
  ),
  scores = c("cellScores", "geneScores"),
  regression_scores = "geneScoresRegression",
  score_tests = "geneScoresTest",
  permutation = c(
    "nPermu", "skrCCAPermuOut", "cellPermu",
    "normalizedCorrelationPermu", "bidirCorrelationPermu", "conditionalPermu"
  )
)

.coProStateAttributes <- list(
  kernel_cache = "kernelNormalizerCache",
  permutation = c("permutationProvenance", "fairSigmaPermu")
)

.coProEmptySlotValues <- list(
  integratedData = list(),
  pcaResults = list(),
  pcaGlobal = list(),
  nPCA = numeric(0),
  scalePCs = logical(0),
  distances = list(),
  distanceScaleFactor = numeric(0),
  distanceGeometry = list(),
  kernelMatrices = list(),
  sigmaValues = numeric(0),
  skrCCAOut = list(),
  nCC = numeric(0),
  normalizedCorrelation = list(),
  bidirCorrelation = list(),
  sigmaValueChoice = numeric(0),
  cellScores = list(),
  geneScores = list(),
  geneScoresRegression = list(),
  geneScoresTest = list(),
  nPermu = numeric(0),
  skrCCAPermuOut = list(),
  cellPermu = list(),
  normalizedCorrelationPermu = list(),
  bidirCorrelationPermu = list(),
  conditionalPermu = list()
)

#' Find every transitive dependent of one or more changed stages
#' @noRd
.coProDependentStages <- function(changed) {
  unknown <- setdiff(changed, names(.coProStateDependencies))
  if (length(unknown) > 0L) {
    stop("Unknown CoPro state stage(s): ", paste(unknown, collapse = ", "))
  }

  pending <- unname(unlist(.coProStateDependencies[changed], use.names = FALSE))
  dependents <- character(0)
  while (length(pending) > 0L) {
    stage <- pending[[1L]]
    pending <- pending[-1L]
    if (stage %in% dependents) next
    dependents <- c(dependents, stage)
    pending <- c(
      pending,
      unname(.coProStateDependencies[[stage]])
    )
  }
  dependents
}

#' Remove only within-cell-type kernels
#'
#' Adding or replacing self-distances cannot make an existing cross-type kernel
#' stale. It does invalidate any self-kernel built from the old distances, so
#' preserve the independent cross-type blocks and prune the usable sigma grid to
#' the blocks that remain.
#' @noRd
.clearCoProSelfKernels <- function(object) {
  kernels <- object@kernelMatrices
  if (length(kernels) == 0L) return(object)

  parsed <- lapply(names(kernels), function(name) {
    tryCatch(.parseKernelMatrixName(name), error = function(e) NULL)
  })
  is_self <- vapply(parsed, function(info) {
    !is.null(info) && identical(info$cellType1, info$cellType2)
  }, logical(1))
  object@kernelMatrices <- kernels[!is_self]

  remaining <- parsed[!is_self]
  remaining_sigmas <- unique(vapply(remaining, function(info) {
    if (is.null(info)) NA_real_ else info$sigma
  }, numeric(1)))
  remaining_sigmas <- remaining_sigmas[is.finite(remaining_sigmas)]
  if (length(object@kernelMatrices) == 0L) {
    object@sigmaValues <- numeric(0)
  } else if (length(remaining_sigmas) > 0L) {
    object@sigmaValues <- object@sigmaValues[
      object@sigmaValues %in% remaining_sigmas
    ]
  }
  object
}

#' Invalidate state downstream of a changed pipeline stage
#'
#' The changed stage itself is left intact because the caller is responsible
#' for rebuilding it. Only its transitive dependents are cleared. This lets a
#' PCA rebuild retain spatial kernels, while a distance rebuild correctly drops
#' those kernels and every result derived from them.
#'
#' @param object A `CoPro` object.
#' @param changed One or more names from `.coProStateDependencies`.
#' @return `object` with stale downstream state removed.
#' @noRd
.invalidateCoProState <- function(object, changed) {
  if (!methods::is(object, "CoPro")) {
    stop("State invalidation requires a CoPro object.")
  }

  stages <- .coProDependentStages(changed)
  if ("self_kernel" %in% stages && !"kernel" %in% stages) {
    object <- .clearCoProSelfKernels(object)
  }

  regular_stages <- setdiff(stages, "self_kernel")
  slots <- unique(unlist(.coProStateSlots[regular_stages], use.names = FALSE))
  for (slot_name in slots) {
    methods::slot(object, slot_name) <- .coProEmptySlotValues[[slot_name]]
  }

  attributes <- unique(unlist(
    .coProStateAttributes[regular_stages], use.names = FALSE
  ))
  for (attribute_name in attributes) {
    attr(object, attribute_name) <- NULL
  }

  # computeGeneAndCellScores() mirrors scores into metadata for convenience.
  # Those columns are derived state too, not user metadata.
  if ("scores" %in% stages && ncol(object@metaDataSub) > 0L) {
    score_columns <- grepl("^cellScore_", colnames(object@metaDataSub))
    if (any(score_columns)) {
      object@metaDataSub <- object@metaDataSub[, !score_columns, drop = FALSE]
    }
  }

  object
}
