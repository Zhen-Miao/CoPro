# Tests for the conditional (sequential step-down) permutation test and the
# sigma-aware bin sizing used by the spatial permutation null.
#
# The conditional test controls two multiplicities at once: sigma selection
# (fair-sigma max-statistic) and canonical-axis ordering (closed step-down).
# The key correctness property is that deflating the FIXED observed
# CC1..CC(k-1) directions reproduces the stored observed CC_k exactly, and that
# the k = 1 conditional test is identical to runSkrCCAPermu_FairSigma().

# ---------------------------------------------------------------------------
# Helper: full single-slide pipeline with several sigma values
# ---------------------------------------------------------------------------
.cond_pipeline <- function(nCC = 3, seed = 42) {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = seed)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1, 0.2),
                             dropDistances = FALSE, verbose = FALSE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE,
                    scalePCs = TRUE)
  suppressWarnings(
    obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = nCC, maxIter = 200)
  )
  obj <- computeNormalizedCorrelation(obj)
  obj
}

test_that("distance scale factor is recovered and sigma-aware bins are sane", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 2)

  sf <- .recoverDistanceScaleFactor(obj)
  expect_true(is.finite(sf) && sf > 0)

  # sf must be the constant mapping raw distance -> stored normalized distance
  D <- obj@distances[[1]]
  loc <- as.data.frame(obj@locationDataSub)
  i <- 1L; j <- 2L
  a <- rownames(D)[i]; b <- colnames(D)[j]
  d_raw <- sqrt((loc[a, "x"] - loc[b, "x"])^2 + (loc[a, "y"] - loc[b, "y"])^2)
  expect_equal(unname(D[i, j]) / d_raw, sf, tolerance = 1e-6)

  bins <- suppressWarnings(
    .sigmaAwareBins(obj, sigma = obj@sigmaValueChoice, verbose = FALSE))
  expect_gte(bins$num_bins_x, 2L)
  expect_gte(bins$num_bins_y, 2L)
  expect_true(is.finite(bins$scale_factor))
})

test_that("distance scale factor survives @distances being dropped", {
  skip_on_cran()
  # computeKernelMatrix(dropDistances = TRUE) is the DEFAULT and clears
  # @distances after building the kernel. The scale factor persisted in
  # @distanceScaleFactor at computeDistance time must still be recoverable.
  obj <- .cond_pipeline(nCC = 2)
  sf_before <- .recoverDistanceScaleFactor(obj)
  expect_true(is.finite(sf_before) && sf_before > 0)

  obj@distances <- list()                     # simulate dropDistances = TRUE
  sf_after <- .recoverDistanceScaleFactor(obj)
  expect_true(is.finite(sf_after) && sf_after > 0)
  expect_equal(sf_after, sf_before, tolerance = 1e-12)

  # sigma-aware bins must NOT fall back to the hard-coded 10 x 10 grid
  bins <- .sigmaAwareBins(obj, sigma = obj@sigmaValueChoice, verbose = FALSE)
  expect_true(is.finite(bins$scale_factor))
  expect_equal(bins$scale_factor, sf_before, tolerance = 1e-12)
})

test_that(".fitConditionalAxis reproduces stored observed CC_k exactly", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 3)
  cts <- obj@cellTypesOfInterest
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  sname <- paste0("sigma_", obj@sigmaValueChoice)
  Wobs <- obj@skrCCAOut[[sname]]
  df <- obj@normalizedCorrelation[[sname]]

  for (k in seq_len(obj@nCC)) {
    fit <- .fitConditionalAxis(PCmats, obj@kernelMatrices,
                               obj@sigmaValueChoice, cts,
                               W_lower = Wobs, k_minus_1 = k - 1)
    obs_k <- df$normalizedCorrelation[df$CC_index == k]
    expect_equal(fit$ncorr, obs_k, tolerance = 1e-6,
                 info = paste("CC", k))
  }
})

test_that("precomputed permutation scoring inputs preserve results", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 2)
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  Wobs <- obj@skrCCAOut[[paste0("sigma_", sigma)]]
  Y0 <- compute_Y_resi(PCmats, obj@kernelMatrices, sigma, cts)
  kernel_info <- .get_ncorr_kernel_info(obj@kernelMatrices, sigma, cts)
  w1 <- lapply(Wobs, function(w) w[, 1, drop = FALSE])

  uncached <- .compute_ncorr_quick(
    PCmats, w1, obj@kernelMatrices, sigma, cts
  )
  cached <- .compute_ncorr_quick(
    PCmats, w1, obj@kernelMatrices, sigma, cts,
    kernel_info = kernel_info, Y_resi = Y0
  )
  expect_equal(cached, uncached, tolerance = 1e-12)

  fit_uncached <- .fitConditionalAxis(
    PCmats, obj@kernelMatrices, sigma, cts, k_minus_1 = 0
  )
  fit_cached <- .fitConditionalAxis(
    PCmats, obj@kernelMatrices, sigma, cts, k_minus_1 = 0,
    Y_resi = Y0, kernel_info = kernel_info
  )
  for (ct in cts) {
    expect_equal(abs(fit_cached$w[[ct]]), abs(fit_uncached$w[[ct]]),
                 tolerance = 1e-5)
  }
  expect_equal(fit_cached$ncorr, fit_uncached$ncorr, tolerance = 1e-10)
})

test_that("conditional k = 1 p-value equals fair-sigma p-value", {
  skip_on_cran()
  obj1 <- .cond_pipeline(nCC = 1)
  seed <- 202

  set.seed(seed)
  fair <- suppressMessages(suppressWarnings(
    runSkrCCAPermu_FairSigma(obj1, nPermu = 25, permu_method = "bin",
                             verbose = FALSE)))
  pf_perm <- vapply(fair@normalizedCorrelationPermu,
                    function(d) d$normalizedCorrelation[1], numeric(1))
  pf <- calculate_pvalue(fair, cc_index = 1)$p_value

  set.seed(seed)
  cond <- suppressMessages(suppressWarnings(
    runSkrCCAPermu_Conditional(obj1, nPermu = 25, permu_method = "bin",
                               verbose = FALSE)))
  pc_perm <- cond@conditionalPermu$perm_stats[, 1]
  pc <- calculate_pvalue_stepdown(cond)$p_raw[1]

  # same permutations + identical RNG consumption -> bitwise-identical null
  expect_equal(unname(pf_perm), unname(pc_perm), tolerance = 1e-10)
  expect_equal(pf, pc, tolerance = 1e-12)
})

test_that("calculate_pvalue uses the same max-over-pairs statistic in every row", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 1)
  observed <- max(getNormCorr(obj)$normalizedCorrelation)
  make_null <- function(values) {
    data.frame(
      sigmaValues = obj@sigmaValueChoice,
      cellType1 = c("A", "A"),
      cellType2 = c("B", "C"),
      CC_index = c(1L, 1L),
      normalizedCorrelation = values,
      stringsAsFactors = FALSE
    )
  }
  obj@normalizedCorrelationPermu <- list(
    permu_1 = make_null(c(observed - 2, observed + 2)),
    permu_2 = make_null(c(observed - 2, observed - 1))
  )

  result <- calculate_pvalue(obj, cc_index = 1, alternative = "greater")
  expect_equal(unname(result$permu_values), c(observed + 2, observed - 1))
  expect_equal(result$p_value, 2 / 3)
  expect_identical(result$pair_aggregation, "max")
})

test_that("quick permutation scoring does not silently use the first pair", {
  expect_error(
    .compute_ncorr_quick(list(), list(), list(), 0.1, c("A", "B", "C")),
    "one predeclared cell-type pair"
  )
})

test_that("step-down p-values are monotone with stop-at-first-NS", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 3)
  set.seed(11)
  cond <- suppressMessages(suppressWarnings(
    runSkrCCAPermu_Conditional(obj, nPermu = 25, permu_method = "bin",
                               verbose = FALSE)))
  pa <- calculate_pvalue_stepdown(cond)

  # closed step-down: non-decreasing and dominates the raw p-values
  expect_true(all(diff(pa$p_stepdown) >= -1e-12))
  expect_true(all(pa$p_stepdown >= pa$p_raw - 1e-12))

  # Phipson-Smyth: p-values are strictly positive and respect the MC floor
  expect_true(all(pa$p_raw > 0))
  expect_true(all(pa$p_raw >= pa$mc_floor - 1e-12))
  expect_equal(unique(pa$mc_floor), 1 / (cond@conditionalPermu$nPermu + 1))

  # significance forms a leading prefix (once false, stays false)
  sigv <- pa$significant
  ff <- which(!sigv)
  if (length(ff) > 0) expect_true(all(!sigv[ff[1]:length(sigv)]))
  expect_equal(attr(pa, "n_significant_axes"), sum(sigv))
})

test_that("explicit num_bins override bypasses sigma-aware sizing", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 2)
  set.seed(5)
  cond <- suppressMessages(suppressWarnings(
    runSkrCCAPermu_Conditional(obj, nPermu = 5, permu_method = "bin",
                               num_bins_x = 8L, num_bins_y = 8L,
                               verbose = FALSE)))
  nb <- cond@conditionalPermu$num_bins
  expect_equal(unname(nb[["x"]]), 8L)
  expect_equal(unname(nb[["y"]]), 8L)
})

test_that("calculate_pvalue_stepdown errors before the test is run", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 2)
  expect_error(calculate_pvalue_stepdown(obj), "runSkrCCAPermu_Conditional")
})

test_that("Y-deflation equals Freedman-Lane data residualization (whitened PCs)", {
  # Deflating Y = X_i^T K X_j by the observed CC1..CC(k-1) weight directions is
  # algebraically identical to residualizing the data by the observed lower
  # canonical variates and recomputing Y, when X^T X = c I (scalePCs = TRUE).
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 3)
  cts <- obj@cellTypesOfInterest
  PCm <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  A <- PCm[[cts[1]]]; B <- PCm[[cts[2]]]
  K <- get_kernel_matrix_flat(obj@kernelMatrices, obj@sigmaValueChoice,
                              cts[1], cts[2], slide = NULL)
  W <- obj@skrCCAOut[[paste0("sigma_", obj@sigmaValueChoice)]]
  leading_sv <- function(Y) irlba::irlba(Y, nv = 1, tol = 1e-9)$d[1]
  Y0 <- crossprod(A, K %*% B)

  for (k in 2:obj@nCC) {
    Ydef <- Y0
    Ar <- A; Br <- B
    for (qq in seq_len(k - 1)) {
      w1 <- W[[cts[1]]][, qq, drop = FALSE]
      w2 <- W[[cts[2]]][, qq, drop = FALSE]
      Ydef <- Ydef - as.numeric(t(w1) %*% Ydef %*% w2) * (w1 %*% t(w2))
      a <- A %*% w1; b <- B %*% w2
      Ar <- Ar - a %*% (crossprod(a, Ar) / sum(a^2))
      Br <- Br - b %*% (crossprod(b, Br) / sum(b^2))
    }
    Yres <- crossprod(Ar, K %*% Br)
    expect_equal(leading_sv(Ydef), leading_sv(Yres), tolerance = 1e-6,
                 info = paste("CC", k))
  }
})

test_that("conditional null projects fixed axes from a permuted operator", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 2)
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCm <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  # A deterministic non-identity permutation makes the observed directions no
  # longer singular vectors of the current operator.
  PCm[[cts[2]]] <- PCm[[cts[2]]][rev(seq_len(nrow(PCm[[cts[2]]]))), , drop = FALSE]
  Y0 <- compute_Y_resi(PCm, obj@kernelMatrices, sigma, cts)
  W <- obj@skrCCAOut[[paste0("sigma_", sigma)]]

  projected <- apply_deflation(
    Y0, W, 1L, cts, deflation = "projection"
  )[[cts[1]]][[cts[2]]]
  rank1 <- apply_deflation(
    Y0, W, 1L, cts, deflation = "rank1"
  )[[cts[1]]][[cts[2]]]
  w1 <- W[[cts[1]]][, 1L, drop = FALSE]
  w2 <- W[[cts[2]]][, 1L, drop = FALSE]

  expect_lt(sqrt(sum(crossprod(w1, projected)^2)), 1e-8)
  expect_lt(sqrt(sum((projected %*% w2)^2)), 1e-8)
  expect_gt(sqrt(sum(crossprod(w1, rank1)^2)) +
              sqrt(sum((rank1 %*% w2)^2)), 1e-4)
})

test_that("fixed-sigma p-values never maximize the observed statistic over other sigmas", {
  skip_on_cran()
  obj <- .cond_pipeline(nCC = 1)
  fixed_sigma <- min(obj@sigmaValues)
  observed_fixed <- max(
    getNormCorr(obj)$normalizedCorrelation[
      getNormCorr(obj)$sigmaValues == fixed_sigma
    ]
  )
  attr(obj, "permutationProvenance") <- list(
    method = "fixed_sigma", sigma_values = fixed_sigma,
    sigma_aggregation = "fixed", pair_aggregation = "max",
    sigma_predeclared = TRUE, selection_adjusted = TRUE
  )
  obj@normalizedCorrelationPermu <- lapply(c(-1, 1), function(offset) {
    data.frame(
      sigmaValues = fixed_sigma,
      cellType1 = "CellTypeA", cellType2 = "CellTypeB", CC_index = 1L,
      normalizedCorrelation = observed_fixed + offset
    )
  })

  result <- calculate_pvalue(obj)
  expect_equal(result$observed, observed_fixed)
  expect_identical(result$sigma_aggregation, "fixed")
  expect_true(result$selection_adjusted)
})

test_that("permutation fits preserve unscaled-PC CCA metrics", {
  skip_on_cran()
  obj <- create_test_copro_single(n_cells = 90, n_cell_types = 2, seed = 77)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE)
  suppressWarnings(
    obj <- computePCA(obj, nPCA = 7, scalePCs = FALSE)
  )
  obj <- runSkrCCA(obj, scalePCs = FALSE, nCC = 2)
  obj <- computeNormalizedCorrelation(obj)
  set.seed(8)
  obj <- suppressWarnings(runSkrCCAPermu(
    obj, nPermu = 10, permu_method = "global", sigma = 0.1,
    verbose = FALSE
  ))

  metrics <- lapply(obj@pcaGlobal, function(x) x$sdev^2)
  first <- obj@skrCCAPermuOut[[1L]]
  for (ct in obj@cellTypesOfInterest) {
    expect_equal(
      crossprod(first[[ct]], metrics[[ct]] * first[[ct]]),
      diag(2), tolerance = 1e-9
    )
  }
})

test_that("calculate_pvalue() provenance stays bound to its null after a later conditional run", {
  skip_on_cran()
  # Regression: permutationProvenance was a single object-level attribute.
  # runSkrCCAPermu_Conditional() overwrote it without touching the base-path
  # null in @normalizedCorrelationPermu, so a subsequent calculate_pvalue()
  # paired that null with conditional provenance -- silently changing the
  # reported p-value and suppressing the sigma-selection warning. Provenance is
  # now bound to the null itself.
  obj <- .cond_pipeline(nCC = 2)

  set.seed(1)
  obj <- suppressMessages(suppressWarnings(
    runSkrCCAPermu(obj, nPermu = 25, permu_method = "global", verbose = FALSE)))
  obj <- suppressMessages(computeNormalizedCorrelationPermu(obj))
  pv_before <- suppressWarnings(calculate_pvalue(obj, cc_index = 1))

  null_before <- lapply(obj@normalizedCorrelationPermu,
                        function(x) x$normalizedCorrelation)

  cond <- suppressMessages(suppressWarnings(
    runSkrCCAPermu_Conditional(obj, nPermu = 25, permu_method = "global",
                               verbose = FALSE)))

  # The conditional run overwrites the shared object-level provenance...
  expect_identical(attr(cond, "permutationProvenance")$method,
                   "conditional_stepdown")
  # ...but leaves the base-path null it does not own untouched...
  expect_equal(lapply(cond@normalizedCorrelationPermu,
                      function(x) x$normalizedCorrelation),
               null_before)
  # ...so calculate_pvalue() returns exactly the same result as before it ran.
  pv_after <- suppressWarnings(calculate_pvalue(cond, cc_index = 1))
  expect_equal(pv_after$p_value, pv_before$p_value)
  expect_identical(pv_after$selection_adjusted, pv_before$selection_adjusted)
  expect_equal(pv_after$sigma_values, pv_before$sigma_values)
})

test_that("n_cores > 1 falls back to sequential under load_all instead of crashing", {
  skip_on_cran()
  # Regression: the PSOCK bootstrap used dirname(find.package("CoPro")) +
  # library(CoPro) in fresh workers, which errors when CoPro is only
  # devtools::load_all()'ed (this project's documented dev workflow, and how
  # this test itself runs). It must now fall back to sequential with a warning,
  # not abort inside clusterCall().
  if (!is.null(.installedCoProLibrary())) {
    skip("CoPro is installed; the load_all() PSOCK fallback path is not exercised.")
  }
  obj <- .cond_pipeline(nCC = 2)

  set.seed(7)
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(
      runSkrCCAPermu_Conditional(obj, nPermu = 12, permu_method = "global",
                                 n_cores = 2, verbose = FALSE)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("sequential", warns)))
  expect_s4_class(res, "CoPro")
  p_raw <- res@conditionalPermu$per_axis$p_raw
  expect_length(p_raw, 2L)
  expect_true(all(is.finite(p_raw) & p_raw > 0 & p_raw <= 1))
})
