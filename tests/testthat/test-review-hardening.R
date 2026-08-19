test_that("unsafe standard deviations are replaced before inverse scaling", {
  observed <- .safeStandardDeviations(c(0, 1e-6, 1e-5, 2, NA, Inf))
  expect_identical(observed, c(1, 1, 1e-5, 2, 1, 1))

  pca <- list(rotation = diag(2), sdev = c(1e-6, 2))
  weights <- matrix(c(3, 4), ncol = 1)
  scores <- .computeGeneScores(
    weights, pca, scalePCs = TRUE, nCC = 1, gene_names = c("g1", "g2")
  )
  expect_equal(unname(scores[, 1]), c(3, 2))
})

test_that("degenerate bidirectional correlations are NA rather than perfect", {
  expect_warning(
    cross <- .computeSpatialCrossCorrelation(
      1, 2, matrix(1, 1, 1), normalize_K = "none", filter_kernel = FALSE
    ),
    "at least two"
  )
  expect_identical(cross, NA_real_)

  expect_warning(
    vectorized <- .computeAllCCCorrelations(
      matrix(1, 1, 1), matrix(2, 1, 1), matrix(1, 1, 1), "none"
    ),
    "at least two"
  )
  expect_true(is.na(vectorized[[1]]))

  expect_warning(
    self <- .computeSpatialSelfCorrelation(
      1, matrix(1, 1, 1), normalize_K = "none", filter_kernel = FALSE
    ),
    "at least two"
  )
  expect_identical(self, NA_real_)
})

test_that("self correlation matches the symmetric direction for unnormalized kernels", {
  scores <- matrix(c(-2, -1, 1, 3), ncol = 1)
  kernel <- matrix(c(
    0, 1, 0.2, 0,
    1, 0, 0.5, 0.1,
    0.2, 0.5, 0, 1,
    0, 0.1, 1, 0
  ), 4, 4, byrow = TRUE)
  expected <- cor(as.vector(crossprod(kernel, scores)), as.vector(scores))
  observed <- .computeSpatialSelfCorrelation(
    scores, kernel, normalize_K = "none", filter_kernel = FALSE
  )
  expect_equal(observed, expected)
})

test_that("Sinkhorn scaling warns when its iteration budget is exhausted", {
  expect_warning(
    sinkhorn_knopp(matrix(c(1, 2, 3, 7), 2, 2), tol = 0, max_iter = 1),
    "did not converge"
  )
})

test_that("self correlation averages both directions after asymmetric Sinkhorn scaling", {
  scores <- matrix(c(-2, -1, 1, 3), ncol = 1)
  scaled_kernel <- matrix(c(
    0, 0.8, 0.1, 0,
    1.2, 0, 0.4, 0.2,
    0.3, 0.7, 0, 0.9,
    0, 0.1, 1.1, 0
  ), 4, 4, byrow = TRUE)
  testthat::local_mocked_bindings(
    sinkhorn_knopp = function(A, ...) scaled_kernel,
    .package = "CoPro"
  )

  expected <- mean(c(
    cor(as.vector(crossprod(scaled_kernel, scores)), as.vector(scores)),
    cor(as.vector(scaled_kernel %*% scores), as.vector(scores))
  ))
  observed <- .computeSpatialSelfCorrelation(
    scores, diag(4), normalize_K = "sinkhorn_knopp", filter_kernel = FALSE
  )
  expect_equal(observed, expected)
})

test_that("spatial partitioning handles constant coordinate axes", {
  one_axis <- data.frame(x = rep(4, 30), y = seq_len(30))
  labels <- .partitionByLocation(one_axis, n = 3, maxCell = 7)
  expect_length(labels, 30)
  expect_lte(max(table(labels)), 7)

  coincident <- data.frame(x = rep(4, 23), y = rep(9, 23))
  labels <- .partitionByLocation(coincident, n = 2, maxCell = 6)
  expect_length(labels, 23)
  expect_lte(max(table(labels)), 6)
})

test_that("quantile normalization rejects non-finite inputs", {
  expect_error(
    quantile_normalize(
      matrix(c(1, 2, 3, NA), 2, 2), matrix(1:4, 2, 2), verbose = FALSE
    ),
    "finite"
  )
})

test_that("transfer gene selection rejects an empty overlap", {
  expect_warning(
    expect_error(.trans_gene_sel(c("a", "b"), c("c", "d")), "no overlapping"),
    "mismatch"
  )
})

test_that("all-zero regression p-values remain finite", {
  observed <- get.onesided.p.value(c(1, -1), c(0, 0))
  expect_true(all(is.finite(observed)))
  expect_gt(observed[[1]], 0)
  expect_lte(observed[[1]], 0.5)
  expect_gte(observed[[2]], 0.5)
})

test_that("sigma keys are stable under the package's numeric round trip", {
  sigmas <- c(0.1, pi, 1.3466659723781049, 12345.6789012345)
  keys <- .sigmaName(sigmas)
  parsed <- as.numeric(sub("^sigma_", "", keys))
  expect_identical(.sigmaName(parsed), keys)
  expect_identical(.sigmaName(0.1), "sigma_0.1")
  original <- vapply(
    sigmas, .kernelNormalizerKey, character(1),
    cellType1 = "A", cellType2 = "B"
  )
  round_tripped <- vapply(
    parsed, .kernelNormalizerKey, character(1),
    cellType1 = "A", cellType2 = "B"
  )
  expect_identical(round_tripped, original)
})

test_that("aggregate correlation supports sum and RSS denominators", {
  slide_denominators <- c(2, 3, 5)
  expect_equal(
    .aggregateNormCorr(slide_denominators, slide_denominators, "sum"),
    1
  )
  expect_equal(
    .aggregateNormCorr(slide_denominators, slide_denominators, "rss"),
    sum(slide_denominators) / sqrt(sum(slide_denominators^2))
  )
  per_slide <- c(0.2, 0.7, 0.9)
  observed <- .aggregateNormCorr(
    per_slide * slide_denominators, slide_denominators, "sum"
  )
  expect_gte(observed, min(per_slide))
  expect_lte(observed, max(per_slide))
})

test_that("sparse zero-distance handling matches an empty dense support", {
  A <- matrix(c(0, 0), nrow = 1)
  B <- rbind(c(0, 0), c(100, 100))
  expect_warning(
    triplets <- .buildBlockTriplets(
      A, B, percentile = NA_real_, scaling_factor = 1,
      max_sigma = 0.01, lowerLimit = 1e-7, truncateLowDist = FALSE
    ),
    "Zero distances"
  )
  expect_null(triplets)
})

test_that("failed sigma weights leave initialized scores as NA", {
  q <- function(x) suppressMessages(suppressWarnings(x))
  obj <- create_test_copro_single(
    n_cells = 100, n_genes = 20, n_cell_types = 2, seed = 812
  )
  obj <- q(subsetData(obj, c("CellTypeA", "CellTypeB")))
  obj <- q(computePCA(obj, nPCA = 5))
  obj <- q(computeKernelMatrix(
    obj, sigmaValues = c(0.05, 0.1), method = "sparse",
    normalizeDistance = TRUE, verbose = FALSE
  ))
  obj <- q(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  failed_name <- .sigmaName(obj@sigmaValues[[1]])
  obj@skrCCAOut[failed_name] <- list(NULL)
  expect_warning(
    scored <- computeGeneAndCellScores(obj),
    "weights are unavailable"
  )
  failed_sigma <- obj@sigmaValues[[1]]
  failed_cell_key <- .createCellScoresName(failed_sigma, "CellTypeA")
  failed_gene_key <- .createGeneScoresName(failed_sigma, "CellTypeA")
  expect_true(all(is.na(scored@cellScores[[failed_cell_key]])))
  expect_true(all(is.na(scored@geneScores[[failed_gene_key]])))
})

test_that("skipped normalized-correlation pairs remain NA", {
  q <- function(x) suppressMessages(suppressWarnings(x))
  obj <- create_test_copro_single(
    n_cells = 100, n_genes = 20, n_cell_types = 2, seed = 813
  )
  obj <- q(subsetData(obj, c("CellTypeA", "CellTypeB")))
  obj <- q(computePCA(obj, nPCA = 5))
  obj <- q(computeKernelMatrix(
    obj, sigmaValues = 0.1, method = "sparse",
    normalizeDistance = TRUE, verbose = FALSE
  ))
  obj <- q(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))
  obj@skrCCAOut[[.sigmaName(0.1)]][["CellTypeA"]] <- NULL

  seen <- character()
  correlated <- withCallingHandlers(
    suppressMessages(computeNormalizedCorrelation(obj)),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("results missing", seen)))
  expect_true(any(grepl("All correlations were NA", seen)))
  expect_true(all(is.na(correlated@normalizedCorrelation[[1]]$normalizedCorrelation)))
  expect_length(correlated@sigmaValueChoice, 0)
})

test_that("regression gene scores reject mismatched cell names", {
  q <- function(x) suppressMessages(suppressWarnings(x))
  obj <- create_test_copro_single(
    n_cells = 100, n_genes = 20, n_cell_types = 2, seed = 814
  )
  obj <- q(subsetData(obj, c("CellTypeA", "CellTypeB")))
  obj <- q(computePCA(obj, nPCA = 5))
  obj <- q(computeKernelMatrix(
    obj, sigmaValues = 0.1, method = "sparse",
    normalizeDistance = TRUE, verbose = FALSE
  ))
  obj <- q(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))
  obj <- q(computeGeneAndCellScores(obj))
  key <- .createCellScoresName(0.1, "CellTypeA")
  rownames(obj@cellScores[[key]])[[1]] <- "not-an-expression-cell"
  expect_error(
    computeRegressionGeneScores(obj, verbose = FALSE),
    "Cell-score names do not match"
  )
})

test_that("mixed-effect gene testing reports a missing optional dependency", {
  expect_error(
    .requireLmerTest("lmerTest sentinel", available = FALSE),
    "lmerTest sentinel"
  )
  testthat::local_mocked_bindings(
    .requireLmerTest = function(...) stop(
      "testGeneMixedEffect() requires the optional package 'lmerTest'."
    ),
    .package = "CoPro"
  )
  expect_error(
    testGeneMixedEffect(NULL, NULL, NULL, NULL),
    "requires the optional package 'lmerTest'"
  )
})

test_that("gene GLM returns NA for an aliased gene coefficient", {
  obj <- create_test_copro_single(
    n_cells = 100, n_genes = 20, n_cell_types = 2, seed = 815
  )
  obj <- suppressMessages(subsetData(obj, c("CellTypeA", "CellTypeB")))
  in_type <- obj@cellTypesSub == "CellTypeA"
  cells <- rownames(obj@normalizedDataSub)[in_type]
  obj@sigmaValues <- 0.1
  obj@nCC <- 1L
  scores <- matrix(
    seq_along(cells), ncol = 1,
    dimnames = list(cells, "CC_1")
  )
  obj@cellScores <- setNames(
    list(scores), .createCellScoresName(0.1, "CellTypeA")
  )
  obj@normalizedDataSub[in_type, 1] <- 1

  observed <- testGeneGLM(
    obj, sigmaName = "sigma_0.1", cellTypeChoice = "CellTypeA",
    covariates = character(), CCChoice = "CC_1", frm = "y ~ x"
  )
  expect_true(all(is.na(observed[1, c("Estimate", "Pr(>|t|)")])))
})
