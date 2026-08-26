.seed_cca_dependents <- function(object) {
  object@skrCCAOut <- list(old_cca = list())
  object@nCC <- 2
  object@normalizedCorrelation <- list(old_correlation = 1)
  object@bidirCorrelation <- list(old_bidir = 1)
  object@sigmaValueChoice <- 99
  object@cellScores <- list(old_cell_scores = matrix(1))
  object@geneScores <- list(old_gene_scores = matrix(1))
  object@geneScoresRegression <- list(old_regression_scores = matrix(1))
  object@geneScoresTest <- list(old_score_tests = TRUE)
  object@nPermu <- 10
  object@skrCCAPermuOut <- list(old_permutation = list())
  object@cellPermu <- list(old_permutation = seq_len(nrow(object@metaDataSub)))
  object@normalizedCorrelationPermu <- list(old_permutation = 1)
  object@bidirCorrelationPermu <- list(old_permutation = 1)
  object@conditionalPermu <- list(old_conditional = TRUE)
  attr(object, "permutationProvenance") <- list(old = TRUE)
  attr(object, "fairSigmaPermu") <- list(old = TRUE)
  object@metaDataSub$cellScore_sigma_99_cc_index_1 <- 1
  object
}

.expect_cca_dependents_empty <- function(object) {
  slots <- c(
    "normalizedCorrelation", "bidirCorrelation", "sigmaValueChoice",
    "cellScores", "geneScores", "geneScoresRegression", "geneScoresTest",
    "nPermu", "skrCCAPermuOut", "cellPermu",
    "normalizedCorrelationPermu", "bidirCorrelationPermu", "conditionalPermu"
  )
  for (slot_name in slots) {
    expect_length(methods::slot(object, slot_name), 0)
  }
  expect_null(attr(object, "permutationProvenance", exact = TRUE))
  expect_null(attr(object, "fairSigmaPermu", exact = TRUE))
  expect_false(any(grepl("^cellScore_", colnames(object@metaDataSub))))
}

.state_test_object <- function(seed = 810) {
  object <- create_test_copro_single(
    n_cells = 120, n_genes = 20, n_cell_types = 3, seed = seed
  )
  subsetData(
    object,
    cellTypesOfInterest = c("CellTypeA", "CellTypeB", "CellTypeC"),
    saveOriginal = TRUE
  )
}

test_that("subsetData clears every state derived from the previous subset", {
  object <- .state_test_object()
  object@integratedData <- list(old_integration = matrix(1))
  object@pcaResults <- list(old_pca = list())
  object@pcaGlobal <- list(old_pca = list())
  object@nPCA <- 5
  object@scalePCs <- TRUE
  object@distances <- list(old_distance = matrix(1))
  object@distanceScaleFactor <- 2
  object@distanceGeometry <- list(source = "old")
  object@kernelMatrices <- list(old_kernel = matrix(1))
  object@sigmaValues <- 99
  attr(object, "kernelNormalizerCache") <- list(old_cache = TRUE)
  object <- .seed_cca_dependents(object)

  object <- subsetData(
    object,
    cellTypesOfInterest = c("CellTypeA", "CellTypeB"),
    saveOriginal = TRUE
  )

  upstream_slots <- c(
    "integratedData", "pcaResults", "pcaGlobal", "nPCA", "scalePCs",
    "distances", "distanceScaleFactor", "distanceGeometry",
    "kernelMatrices", "sigmaValues", "skrCCAOut", "nCC"
  )
  for (slot_name in upstream_slots) {
    expect_length(methods::slot(object, slot_name), 0)
  }
  expect_null(attr(object, "kernelNormalizerCache", exact = TRUE))
  .expect_cca_dependents_empty(object)
})

test_that("computePCA clears CCA state but preserves spatial state", {
  object <- .state_test_object(seed = 811)
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    verbose = FALSE
  )
  distances_before <- object@distances
  kernels_before <- object@kernelMatrices
  sigmas_before <- object@sigmaValues
  attr(object, "kernelNormalizerCache") <- list(valid_cache = TRUE)
  object <- .seed_cca_dependents(object)

  object <- suppressMessages(computePCA(object, nPCA = 5))

  expect_identical(object@distances, distances_before)
  expect_identical(object@kernelMatrices, kernels_before)
  expect_identical(object@sigmaValues, sigmas_before)
  expect_identical(
    attr(object, "kernelNormalizerCache", exact = TRUE),
    list(valid_cache = TRUE)
  )
  expect_length(object@pcaGlobal, 3)
  expect_length(object@skrCCAOut, 0)
  expect_length(object@nCC, 0)
  .expect_cca_dependents_empty(object)
})

test_that("computeDistance clears kernels and their downstream state only", {
  object <- .state_test_object(seed = 812)
  object <- suppressMessages(computePCA(object, nPCA = 5))
  pca_before <- object@pcaGlobal
  object@kernelMatrices <- list(old_kernel = matrix(1))
  object@sigmaValues <- 99
  attr(object, "kernelNormalizerCache") <- list(old_cache = TRUE)
  object <- .seed_cca_dependents(object)

  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )

  expect_identical(object@pcaGlobal, pca_before)
  expect_length(object@distances, 3)
  expect_length(object@kernelMatrices, 0)
  expect_length(object@sigmaValues, 0)
  expect_null(attr(object, "kernelNormalizerCache", exact = TRUE))
  expect_length(object@skrCCAOut, 0)
  expect_length(object@nCC, 0)
  .expect_cca_dependents_empty(object)
})

test_that("computeKernelMatrix clears CCA outputs while preserving its inputs", {
  object <- .state_test_object(seed = 813)
  object <- suppressMessages(computePCA(object, nPCA = 5))
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  pca_before <- object@pcaGlobal
  distances_before <- object@distances
  attr(object, "kernelNormalizerCache") <- list(old_cache = TRUE)
  object <- .seed_cca_dependents(object)

  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    verbose = FALSE
  )

  expect_identical(object@pcaGlobal, pca_before)
  expect_identical(object@distances, distances_before)
  expect_gt(length(object@kernelMatrices), 0)
  expect_identical(object@sigmaValues, 3)
  expect_null(attr(object, "kernelNormalizerCache", exact = TRUE))
  expect_length(object@skrCCAOut, 0)
  expect_length(object@nCC, 0)
  .expect_cca_dependents_empty(object)
})

test_that("runSkrCCA replaces CCA and clears all results derived from it", {
  object <- .state_test_object(seed = 814)
  object <- suppressMessages(computePCA(object, nPCA = 5))
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    verbose = FALSE
  )
  object <- .seed_cca_dependents(object)

  object <- suppressMessages(runSkrCCA(
    object, scalePCs = TRUE, nCC = 1, maxIter = 30
  ))

  expect_gt(length(object@skrCCAOut), 0)
  expect_identical(object@nCC, 1)
  .expect_cca_dependents_empty(object)
})

test_that("computeSelfDistance preserves cross-kernels and drops stale self-kernels", {
  object <- .state_test_object(seed = 815)
  object <- suppressMessages(computePCA(object, nPCA = 5))
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    verbose = FALSE
  )
  cross_names <- names(object@kernelMatrices)
  self_name <- CoPro:::.createKernelMatrixName(
    3, "CellTypeA", "CellTypeA"
  )
  n_a <- sum(object@cellTypesSub == "CellTypeA")
  object@kernelMatrices[[self_name]] <- diag(n_a)
  object <- .seed_cca_dependents(object)
  cca_before <- object@skrCCAOut
  scores_before <- lapply(
    c("cellScores", "geneScores", "geneScoresRegression", "geneScoresTest"),
    function(slot_name) methods::slot(object, slot_name)
  )

  object <- computeSelfDistance(object, verbose = FALSE, overwrite = FALSE)

  expect_true(all(cross_names %in% names(object@kernelMatrices)))
  expect_false(self_name %in% names(object@kernelMatrices))
  expect_identical(object@sigmaValues, 3)
  expect_identical(object@skrCCAOut, cca_before)
  expect_identical(object@nCC, 2)
  score_slots <- c(
    "cellScores", "geneScores", "geneScoresRegression", "geneScoresTest"
  )
  for (index in seq_along(score_slots)) {
    expect_identical(
      methods::slot(object, score_slots[[index]]), scores_before[[index]]
    )
  }
  expect_true(any(grepl("^cellScore_", colnames(object@metaDataSub))))

  invalidated_slots <- c(
    "normalizedCorrelation", "bidirCorrelation", "sigmaValueChoice",
    "nPermu", "skrCCAPermuOut", "cellPermu",
    "normalizedCorrelationPermu", "bidirCorrelationPermu", "conditionalPermu"
  )
  for (slot_name in invalidated_slots) {
    expect_length(methods::slot(object, slot_name), 0)
  }
  expect_null(attr(object, "kernelNormalizerCache", exact = TRUE))
  expect_null(attr(object, "permutationProvenance", exact = TRUE))
  expect_null(attr(object, "fairSigmaPermu", exact = TRUE))
})

test_that("additive self-kernels preserve CCA weights and scores", {
  object <- .state_test_object(seed = 816)
  object <- suppressMessages(computePCA(object, nPCA = 5))
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    verbose = FALSE
  )
  object <- .seed_cca_dependents(object)
  cca_before <- object@skrCCAOut
  scores_before <- object@cellScores

  object <- suppressWarnings(computeSelfKernel(
    object, sigmaValues = 3, method = "sparse", overwrite = FALSE,
    minAveCellNeighor = 1, verbose = FALSE
  ))

  expect_identical(object@skrCCAOut, cca_before)
  expect_identical(object@nCC, 2)
  expect_identical(object@cellScores, scores_before)
  expect_gt(sum(vapply(names(object@kernelMatrices), function(name) {
    parsed <- CoPro:::.parseKernelMatrixName(name)
    identical(parsed$cellType1, parsed$cellType2)
  }, logical(1))), 0)
  expect_length(object@normalizedCorrelation, 0)
  expect_length(object@normalizedCorrelationPermu, 0)
})

test_that("additive float32 kernels preserve CCA state at untouched sigmas", {
  object <- .state_test_object(seed = 817)
  object <- suppressMessages(computePCA(object, nPCA = 5))
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    verbose = FALSE
  )
  existing_kernels <- names(object@kernelMatrices)
  object <- .seed_cca_dependents(object)
  cca_before <- object@skrCCAOut
  preserved_slots <- c(
    "cellScores", "geneScores", "geneScoresRegression", "geneScoresTest"
  )
  scores_before <- lapply(
    preserved_slots, function(slot_name) methods::slot(object, slot_name)
  )

  object <- suppressWarnings(computeSparseKernelFloat32(
    object, sigmaValues = 4, overwrite = FALSE,
    minAveCellNeighor = 1, verbose = FALSE
  ))

  expect_true(all(existing_kernels %in% names(object@kernelMatrices)))
  expect_setequal(object@sigmaValues, c(3, 4))
  expect_identical(object@skrCCAOut, cca_before)
  expect_identical(object@nCC, 2)
  for (index in seq_along(preserved_slots)) {
    expect_identical(
      methods::slot(object, preserved_slots[[index]]), scores_before[[index]]
    )
  }
  expect_true(any(grepl("^cellScore_", colnames(object@metaDataSub))))

  invalidated_slots <- c(
    "normalizedCorrelation", "bidirCorrelation", "sigmaValueChoice",
    "nPermu", "skrCCAPermuOut", "cellPermu",
    "normalizedCorrelationPermu", "bidirCorrelationPermu", "conditionalPermu"
  )
  for (slot_name in invalidated_slots) {
    expect_length(methods::slot(object, slot_name), 0)
  }
  expect_null(attr(object, "kernelNormalizerCache", exact = TRUE))
  expect_null(attr(object, "permutationProvenance", exact = TRUE))
  expect_null(attr(object, "fairSigmaPermu", exact = TRUE))

  # overwrite = TRUE remains a full kernel rebuild that clears CCA state.
  object <- .seed_cca_dependents(object)
  object <- suppressWarnings(computeSparseKernelFloat32(
    object, sigmaValues = 4, overwrite = TRUE,
    minAveCellNeighor = 1, verbose = FALSE
  ))
  expect_length(object@skrCCAOut, 0)
  expect_identical(object@sigmaValues, 4)
  .expect_cca_dependents_empty(object)
})

test_that("float32 replacement invalidates CCA state at touched sigmas", {
  object <- .state_test_object(seed = 818)
  object <- suppressMessages(computePCA(object, nPCA = 5))
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    verbose = FALSE
  )
  kernel_name <- names(object@kernelMatrices)[[1L]]
  kernel_before <- as.matrix(object@kernelMatrices[[kernel_name]])
  object <- .seed_cca_dependents(object)

  object <- suppressWarnings(computeSparseKernelFloat32(
    object, sigmaValues = 3, rowNormalizeKernel = TRUE,
    overwrite = FALSE, minAveCellNeighor = 1, verbose = FALSE
  ))

  kernel_after <- as.matrix(asDoubleSparseMatrix(
    object@kernelMatrices[[kernel_name]]
  ))
  expect_false(isTRUE(all.equal(kernel_after, kernel_before)))
  expect_identical(object@sigmaValues, 3)
  expect_length(object@skrCCAOut, 0)
  expect_length(object@nCC, 0)
  .expect_cca_dependents_empty(object)
})

test_that("single-type self-kernel no-ops preserve all derived state", {
  object <- create_test_copro_single(
    n_cells = 80, n_genes = 20, n_cell_types = 1, seed = 819
  )
  object <- subsetData(
    object, cellTypesOfInterest = "CellTypeA", saveOriginal = TRUE
  )
  object <- suppressMessages(computePCA(object, nPCA = 5))
  object <- computeDistance(
    object, distType = "Euclidean2D", normalizeDistance = FALSE,
    verbose = FALSE
  )
  object <- computeKernelMatrix(
    object, sigmaValues = 3, method = "dense", dropDistances = FALSE,
    minAveCellNeighor = 1, verbose = FALSE
  )
  object <- .seed_cca_dependents(object)
  attr(object, "kernelNormalizerCache") <- list(valid_cache = TRUE)

  for (overwrite in c(FALSE, TRUE)) {
    expect_warning(
      result <- computeSelfKernel(
        object, sigmaValues = 3, method = "dense", overwrite = overwrite,
        minAveCellNeighor = 1, verbose = FALSE
      ),
      "Only one cell type detected"
    )
    expect_identical(result, object)
  }
})
