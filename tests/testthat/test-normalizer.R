## Tests for the `normalizer` argument of computeNormalizedCorrelation().
##
## The denominator ||R_x^{1/2} K_c R_y^{1/2}||_F is the null SD of the bilinear
## statistic, so the whitening operators decide how the criterion drifts across
## the sigma grid. These tests pin the back-compatible path, the explicit
## failure of the self-kernel path, the range estimator, and the property the
## variogram normalizer exists for: a flatter null path.

## Two cell types whose expression is spatially smooth at a known range and,
## unless `shared = TRUE`, independent across the two types (i.e. the null).
make_smooth_object <- function(n_per_type = 200, n_genes = 60, range = 0.05,
                               shared = FALSE, seed = 1,
                               subset_to = c("TypeA", "TypeB"), n_pca = 8) {
  set.seed(seed)
  n <- 2 * n_per_type
  loc <- data.frame(x = runif(n), y = runif(n))
  ids <- paste0("cell_", seq_len(n))
  rownames(loc) <- ids

  smoother <- exp(-0.5 * (as.matrix(dist(cbind(loc$x, loc$y))) / range)^2)
  smoother <- smoother / sqrt(rowSums(smoother^2))

  is_a <- seq_len(n) <= n_per_type
  field <- function() smoother %*% matrix(rnorm(n * n_genes), n, n_genes)
  expr <- matrix(0, n, n_genes)
  if (shared) {
    common <- field()
    expr[is_a, ] <- common[is_a, ]
    expr[!is_a, ] <- common[!is_a, ]
  } else {
    expr[is_a, ] <- field()[is_a, ]
    expr[!is_a, ] <- field()[!is_a, ]
  }
  expr <- expr + matrix(rnorm(n * n_genes, sd = 0.2), n, n_genes)
  dimnames(expr) <- list(ids, paste0("Gene", seq_len(n_genes)))

  meta <- data.frame(cell_type = ifelse(is_a, "TypeA", "TypeB"), row.names = ids)
  obj <- newCoProSingle(expr, loc, meta, meta$cell_type)
  obj <- subsetData(obj, cellTypesOfInterest = subset_to)
  obj <- computePCA(obj, nPCA = n_pca)
  obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  obj
}

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

cc1 <- function(object, sigma_name) {
  df <- object@normalizedCorrelation[[sigma_name]]
  df$normalizedCorrelation[df$CC_index == 1]
}

test_that("legacy equals unwhitened when the object has no within-type kernels", {
  obj <- quiet(make_smooth_object(n_per_type = 120, n_genes = 40))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = c(0.02, 0.05), verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  ## computeKernelMatrix() builds only cross-type kernels, so this is the
  ## situation every standard two-type workflow is actually in.
  expect_false(any(grepl("TypeA\\|TypeA", names(obj@kernelMatrices))))

  legacy <- quiet(computeNormalizedCorrelation(obj))
  unwh <- quiet(computeNormalizedCorrelation(obj, normalizer = "unwhitened"))

  expect_equal(cc1(legacy, "sigma_0.05"), cc1(unwh, "sigma_0.05"))
  expect_equal(legacy@sigmaValueChoice, unwh@sigmaValueChoice)
  expect_match(getNormalizerInfo(legacy)$description, "no within-type kernels")
})

test_that("normalizer = 'kernel' fails loudly instead of silently unwhitening", {
  obj <- quiet(make_smooth_object(n_per_type = 120, n_genes = 40))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  expect_error(
    quiet(computeNormalizedCorrelation(obj, normalizer = "kernel")),
    "within-type kernel"
  )
})

test_that("legacy silently whitens once self-kernels are present", {
  obj <- quiet(make_smooth_object(n_per_type = 120, n_genes = 40))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(computeSelfDistance(obj, distType = "Euclidean2D", verbose = FALSE))
  obj <- quiet(computeSelfKernel(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  skip_if_not(any(grepl("TypeA\\|TypeA", names(obj@kernelMatrices))),
              "self-kernels were not stored under the expected keys")

  legacy <- quiet(computeNormalizedCorrelation(obj))
  unwh <- quiet(computeNormalizedCorrelation(obj, normalizer = "unwhitened"))

  ## This is the reproducibility hazard the argument exists to surface: the
  ## same call returns a different statistic depending on the object's history.
  expect_false(isTRUE(all.equal(cc1(legacy, "sigma_0.05"),
                                cc1(unwh, "sigma_0.05"))))
  expect_match(getNormalizerInfo(legacy)$description, "self-kernel")
  expect_match(getNormalizerInfo(unwh)$description, "R = I")

  ## A cross-type self-kernel still carries the physical-units caveat, but its
  ## diagonal is repaired in the private whitening copy before use.
  seen <- character(0)
  withCallingHandlers(
    suppressMessages(computeNormalizedCorrelation(obj)),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("different physical bandwidth", seen)))
  expect_false(any(grepl("all-zero diagonal", seen)))
})

test_that("a zero-diagonal stored kernel gets a unit-diagonal whitening copy", {
  ## Self-kernels are stored with a zero diagonal on purpose -- correct in the
  ## numerator, where a cell should not be its own neighbour. But R_x is a
  ## correlation operator, so R_ii = 1; zeroing it forces trace(R) = 0, which
  ## for a nonzero symmetric matrix means negative eigenvalues exist.
  obj <- quiet(make_smooth_object(n_per_type = 150, n_genes = 40,
                                  subset_to = "TypeA", n_pca = 6))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  K <- as.matrix(getKernelMatrix(obj, sigma = 0.05, cellType1 = "TypeA",
                                 cellType2 = "TypeA", verbose = FALSE))
  expect_equal(sum(diag(K)), 0)
  eigenvalues <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(eigenvalues), 0)

  ## Both the legacy and explicit self-kernel modes now build a valid
  ## correlation operator without mutating the stored analysis kernel.
  expect_no_warning(
    suppressMessages(computeNormalizedCorrelation(obj)))
  expect_no_warning(
    suppressMessages(computeNormalizedCorrelation(obj, normalizer = "kernel")))
  expect_no_warning(
    suppressMessages(computeNormalizedCorrelation(obj, normalizer = "unwhitened")))
})

test_that("legacy and kernel modes use the same repaired self-kernel", {
  obj <- quiet(make_smooth_object(n_per_type = 150, n_genes = 40,
                                  subset_to = "TypeA", n_pca = 6))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  legacy <- quiet(computeNormalizedCorrelation(obj))
  fixed <- quiet(computeNormalizedCorrelation(obj, normalizer = "kernel"))
  expect_equal(cc1(legacy, "sigma_0.05"), cc1(fixed, "sigma_0.05"))
})

test_that("a within-type analysis is exempt from the self-kernel units warning", {
  ## With one cell type the whitening operator IS the analysis kernel, built by
  ## computeKernelMatrix() under a single scaling factor, so there is no
  ## mismatch to warn about.
  obj <- quiet(make_smooth_object(n_per_type = 120, n_genes = 40,
                                  subset_to = "TypeA", n_pca = 6))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  ## The units warning is about mixing a self-kernel built by
  ## computeSelfDistance() with a cross-kernel built by computeDistance().
  ## There is no cross-kernel in a within-type analysis.
  seen <- character(0)
  withCallingHandlers(
    suppressMessages(computeNormalizedCorrelation(obj)),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("different physical bandwidth", seen)))
})

test_that("the variogram normalizer recovers a planted autocorrelation range", {
  ## Smoothing white noise with bandwidth b gives a field of range b * sqrt(2).
  smoothing <- 0.04
  expected <- smoothing * sqrt(2)

  obj <- quiet(make_smooth_object(n_per_type = 300, n_genes = 60,
                                  range = smoothing, seed = 7))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))
  out <- quiet(computeNormalizedCorrelation(obj, normalizer = "variogram"))

  info <- getNormalizerInfo(out)
  expect_identical(info$mode, "variogram")
  expect_true(all(is.finite(info$ranges)))

  ## ranges are reported in normalized distance units
  expected_norm <- expected * obj@distanceScaleFactor
  expect_equal(unname(info$ranges[["TypeA"]]), expected_norm, tolerance = 0.25)
  expect_equal(unname(info$ranges[["TypeB"]]), expected_norm, tolerance = 0.25)
})

test_that("a supplied range bypasses estimation and changes the denominator", {
  obj <- quiet(make_smooth_object(n_per_type = 150, n_genes = 40))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  small <- quiet(computeNormalizedCorrelation(
    obj, normalizer = "variogram",
    normalizerControl = list(range = c(TypeA = 0.01, TypeB = 0.01))))
  large <- quiet(computeNormalizedCorrelation(
    obj, normalizer = "variogram",
    normalizerControl = list(range = c(TypeA = 0.08, TypeB = 0.08))))

  expect_equal(unname(getNormalizerInfo(small)$ranges[["TypeA"]]), 0.01)
  expect_equal(unname(getNormalizerInfo(large)$ranges[["TypeA"]]), 0.08)
  ## a wider whitening operator means a larger null SD, hence a smaller ratio
  expect_lt(abs(cc1(large, "sigma_0.05")), abs(cc1(small, "sigma_0.05")))
})

test_that("normalizerControl rejects unknown entries", {
  obj <- quiet(make_smooth_object(n_per_type = 120, n_genes = 40))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  expect_error(
    quiet(computeNormalizedCorrelation(
      obj, normalizer = "variogram",
      normalizerControl = list(rangeGuess = 0.02))),
    "Unknown .normalizerControl. entries"
  )
})

test_that("the variogram normalizer flattens the null path across sigma", {
  ## The defect being fixed: with R = I the criterion drifts upward with sigma
  ## under a pure null, so bandwidth selection rails to the top of the grid.
  obj <- quiet(make_smooth_object(n_per_type = 250, n_genes = 60,
                                  range = 0.04, shared = FALSE, seed = 11))
  sigmas <- c(0.015, 0.03, 0.06, 0.12)
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = sigmas, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  path <- function(o) {
    vapply(paste0("sigma_", sigmas), function(s) abs(cc1(o, s)), numeric(1))
  }
  unwh_path <- path(quiet(computeNormalizedCorrelation(obj, normalizer = "unwhitened")))
  vario_out <- quiet(computeNormalizedCorrelation(obj, normalizer = "variogram"))
  skip_if_not(all(is.finite(getNormalizerInfo(vario_out)$ranges)),
              "range estimation did not converge on this draw")
  vario_path <- path(vario_out)

  drift <- function(p) max(p) / min(p)
  expect_gt(drift(unwh_path), drift(vario_path))
  ## and the unwhitened path is the specific failure mode: monotone increasing
  expect_equal(which.max(unwh_path), length(sigmas), ignore_attr = TRUE)
})

test_that("the variogram normalizer works on a multi-slide object", {
  obj <- create_test_copro_multi(n_cells_per_slide = 90, n_slides = 2,
                                 n_genes = 40, n_cell_types = 2, seed = 5)
  cts <- sort(unique(as.character(obj@cellTypes)))[1:2]
  obj <- quiet(subsetData(obj, cellTypesOfInterest = cts))
  obj <- quiet(computePCA(obj, nPCA = 5))
  ## the helper lays cells out over a 10 x 10 field, so sigma = 0.05 only means
  ## a populated kernel once the distances are rescaled
  obj <- quiet(computeDistance(obj, distType = "Euclidean2D",
                               normalizeDistance = TRUE, verbose = FALSE))
  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))

  out <- quiet(computeNormalizedCorrelation(
    obj, normalizer = "variogram",
    ## white expression on this helper, so the range is supplied rather than fit
    normalizerControl = list(range = setNames(rep(0.05, length(cts)), cts))))

  info <- getNormalizerInfo(out)
  expect_identical(info$mode, "variogram")
  expect_named(info$ranges, getSlideList(obj))
  expect_true(nrow(out@normalizedCorrelation[["sigma_0.05"]]) > 0)
})
