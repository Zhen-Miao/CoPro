# Tests for detectSigmaRange(), the spacing-based distance normalization, and
# the storage-format decision that method = "auto" now makes.

.sigma_obj_single <- function(n_per_type = 400, extent = 1000, unit = 1,
                              types = c("CellTypeA", "CellTypeB"), seed = 11) {
  set.seed(seed)
  n <- n_per_type * length(types)
  ids <- paste0("cell_", seq_len(n))
  expr <- matrix(rnorm(n * 20, 2, 1), n, 20,
                 dimnames = list(ids, paste0("Gene", seq_len(20))))
  expr[expr < 0] <- 0
  loc <- data.frame(x = runif(n, 0, extent) * unit,
                    y = runif(n, 0, extent) * unit, row.names = ids)
  ct <- rep(types, each = n_per_type)
  meta <- data.frame(cell_type = ct, row.names = ids)
  obj <- newCoProSingle(normalizedData = expr, locationData = loc,
                        metaData = meta, cellTypes = ct)
  subsetData(obj, cellTypesOfInterest = types)
}

.sigma_obj_multi <- function(extents = c(S1 = 1000, S2 = 1000), n_per_type = 300,
                             types = c("CellTypeA", "CellTypeB"), seed = 12) {
  set.seed(seed)
  parts <- lapply(names(extents), function(s) {
    n <- n_per_type * length(types)
    ids <- paste0(s, "::cell_", seq_len(n))
    list(
      expr = matrix(pmax(rnorm(n * 20, 2, 1), 0), n, 20,
                    dimnames = list(ids, paste0("Gene", seq_len(20)))),
      loc = data.frame(x = runif(n, 0, extents[[s]]),
                       y = runif(n, 0, extents[[s]]), row.names = ids),
      ct = rep(types, each = n_per_type), slide = rep(s, n)
    )
  })
  expr <- do.call(rbind, lapply(parts, `[[`, "expr"))
  loc <- do.call(rbind, lapply(parts, `[[`, "loc"))
  ct <- unlist(lapply(parts, `[[`, "ct"))
  slide <- unlist(lapply(parts, `[[`, "slide"))
  meta <- data.frame(cell_type = ct, row.names = rownames(expr))
  obj <- newCoProMulti(normalizedData = expr, locationData = loc,
                       metaData = meta, cellTypes = ct, slideID = slide)
  subsetData(obj, cellTypesOfInterest = types)
}

# True median effective neighbor count, computed the slow exact way.
.true_effective_neighbors <- function(obj, ct_i, ct_j, sigma) {
  loc <- obj@locationDataSub
  A <- as.matrix(loc[obj@cellTypesSub == ct_i, c("x", "y")])
  B <- as.matrix(loc[obj@cellTypesSub == ct_j, c("x", "y")])
  D <- fields::rdist(A, B)
  stats::median(rowSums(exp(-0.5 * (D / sigma)^2)))
}

test_that("detected bandwidths hit the requested effective neighbor counts", {
  obj <- .sigma_obj_single()
  rng <- detectSigmaRange(obj, minNeighbors = 5, maxNeighbors = 20,
                          verbose = FALSE)

  expect_s3_class(rng, "CoProSigmaRange")
  expect_true(rng$feasible)
  expect_true(rng$sigmaRange[["lower"]] < rng$sigmaRange[["upper"]])

  low <- .true_effective_neighbors(obj, "CellTypeA", "CellTypeB",
                                   rng$sigmaRange[["lower"]])
  high <- .true_effective_neighbors(obj, "CellTypeA", "CellTypeB",
                                    rng$sigmaRange[["upper"]])
  # The estimate is sampled and neighbor-truncated on purpose, so allow 15%.
  expect_equal(low, 5, tolerance = 0.15)
  expect_equal(high, 20, tolerance = 0.15)
})

test_that("detected bandwidths scale with the coordinate units", {
  microns <- detectSigmaRange(.sigma_obj_single(unit = 1), verbose = FALSE)
  pixels <- detectSigmaRange(.sigma_obj_single(unit = 10), verbose = FALSE)

  # Same tissue, coordinates x10 -> every bandwidth x10. This is the property
  # that makes a recommended sigma travel between datasets.
  expect_equal(pixels$sigmaRange[["lower"]] / microns$sigmaRange[["lower"]],
               10, tolerance = 0.05)
  expect_equal(pixels$sigmaRange[["upper"]] / microns$sigmaRange[["upper"]],
               10, tolerance = 0.05)
})

test_that("detection is deterministic and does not disturb the RNG", {
  obj <- .sigma_obj_single()
  set.seed(99)
  before <- .Random.seed
  a <- detectSigmaRange(obj, verbose = FALSE)
  b <- detectSigmaRange(obj, verbose = FALSE)
  expect_identical(a$sigmaValues, b$sigmaValues)
  expect_identical(a$sigmaRange, b$sigmaRange)
  expect_identical(before, .Random.seed)
})

test_that("the recommended grid spans the detected range and is usable", {
  obj <- .sigma_obj_single()
  rng <- detectSigmaRange(obj, nSigma = 4L, verbose = FALSE)

  expect_length(rng$sigmaValues, 4L)
  expect_true(all(diff(rng$sigmaValues) > 0))
  expect_gte(min(rng$sigmaValues), rng$sigmaRange[["lower"]] * 0.99)
  expect_lte(max(rng$sigmaValues), rng$sigmaRange[["upper"]] * 1.01)

  # The whole point: they can be handed straight to computeKernelMatrix.
  built <- computeKernelMatrix(obj, sigmaValues = rng$sigmaValues,
                               verbose = FALSE)
  expect_length(built@sigmaValues, 4L)
  expect_gt(length(built@kernelMatrices), 0L)
})

test_that("per-block diagnostics cover every cell-type pair and slide", {
  obj <- .sigma_obj_multi()
  rng <- detectSigmaRange(obj, verbose = FALSE)

  expect_equal(nrow(rng$blocks), 2L)             # one pair x two slides
  expect_setequal(rng$blocks$slide, c("S1", "S2"))
  expect_true(all(rng$blocks$medianSpacing > 0))
  expect_true(all(is.finite(rng$blocks$sigmaLo)))
  expect_true(all(rng$blocks$neighborsAtLower <= rng$blocks$neighborsAtUpper))
})

test_that("a single cell type is handled through the within-type block", {
  obj <- .sigma_obj_single(types = "CellTypeA", n_per_type = 500)
  rng <- detectSigmaRange(obj, verbose = FALSE)
  expect_equal(nrow(rng$blocks), 1L)
  expect_equal(rng$blocks$cellType1, rng$blocks$cellType2)
  expect_true(is.finite(rng$sigmaRange[["lower"]]))
})

test_that("detection reports rather than hides an infeasible shared range", {
  # One slide 100x denser than the other: no single bandwidth can give both
  # blocks 5-20 effective neighbors.
  obj <- .sigma_obj_multi(extents = c(S1 = 2000, S2 = 60))
  expect_warning(rng <- detectSigmaRange(obj, verbose = FALSE),
                 "No single sigma")
  expect_false(rng$feasible)
  expect_true(all(is.finite(rng$sigmaValues)))
})

test_that("invalid neighbor targets are rejected", {
  obj <- .sigma_obj_single(n_per_type = 100)
  expect_error(detectSigmaRange(obj, minNeighbors = 0, verbose = FALSE),
               "minNeighbors")
  expect_error(detectSigmaRange(obj, minNeighbors = 10, maxNeighbors = 5,
                                verbose = FALSE), "maxNeighbors")
})

# ---------------------------------------------------------------------------
# normalizeMethod
# ---------------------------------------------------------------------------

test_that("normalizeDistance now defaults to FALSE and is recorded", {
  obj <- .sigma_obj_single(n_per_type = 150)
  out <- computeKernelMatrix(obj, sigmaValues = 40, verbose = FALSE)
  geom <- getDistanceGeometry(out)
  expect_false(geom$normalizeDistance)
  expect_equal(geom$normalizeMethod, "spacing")
})

test_that("spacing normalization resists a single much denser block", {
  scale_of <- function(obj, method) {
    computeKernelMatrix(obj, sigmaValues = 0.5, normalizeDistance = TRUE,
                        normalizeMethod = method, method = "float32",
                        verbose = FALSE)@distanceScaleFactor
  }
  balanced <- .sigma_obj_multi(extents = c(S1 = 1000, S2 = 1000))
  skewed <- .sigma_obj_multi(extents = c(S1 = 1000, S2 = 100))

  spacing_shift <- scale_of(skewed, "spacing") / scale_of(balanced, "spacing")
  pct_shift <- scale_of(skewed, "percentile") / scale_of(balanced, "percentile")

  # "percentile" takes the minimum across blocks, so the densest block sets the
  # unit for the whole object; "spacing" takes the median and moves far less.
  expect_lt(abs(log(spacing_shift)), abs(log(pct_shift)))
})

test_that("dense and sparse paths agree under both normalization methods", {
  obj <- .sigma_obj_single(n_per_type = 150, extent = 10)
  for (method in c("spacing", "percentile")) {
    d <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, normalizeMethod = method,
                         verbose = FALSE)
    dense <- computeKernelMatrix(d, sigmaValues = 0.1, method = "dense",
                                 dropDistances = FALSE, verbose = FALSE)
    sparse <- computeKernelMatrix(obj, sigmaValues = 0.1, method = "sparse",
                                  distType = "Euclidean2D",
                                  normalizeDistance = TRUE,
                                  normalizeMethod = method, verbose = FALSE)
    expect_equal(dense@distanceScaleFactor, sparse@distanceScaleFactor,
                 tolerance = 1e-10, info = method)
    expect_equal(
      as.matrix(getKernelMatrix(dense, 0.1, "CellTypeA", "CellTypeB",
                                verbose = FALSE)),
      as.matrix(getKernelMatrix(sparse, 0.1, "CellTypeA", "CellTypeB",
                                verbose = FALSE)),
      tolerance = 1e-10, ignore_attr = TRUE, info = method)
  }
})

test_that("an unknown normalizeMethod is rejected", {
  obj <- .sigma_obj_single(n_per_type = 100)
  expect_error(
    computeKernelMatrix(obj, sigmaValues = 40, normalizeDistance = TRUE,
                        normalizeMethod = "median", verbose = FALSE),
    "normalizeMethod")
})

# ---------------------------------------------------------------------------
# storage-format decision
# ---------------------------------------------------------------------------

test_that("auto picks float32 above the threshold and dense below", {
  obj <- .sigma_obj_single(n_per_type = 400)
  rng <- detectSigmaRange(obj, verbose = FALSE)
  sig <- rng$sigmaValues[1]

  big <- computeKernelMatrix(obj, sigmaValues = sig, method = "auto",
                             autoThreshold = 10L, verbose = FALSE)
  expect_s3_class(getKernelMatrix(big, sig, "CellTypeA", "CellTypeB",
                                  verbose = FALSE, materialize = FALSE),
                  "CoProFloat32SparseMatrix")

  small <- computeKernelMatrix(obj, sigmaValues = sig, method = "auto",
                               autoThreshold = 100000L, verbose = FALSE)
  expect_true(isS4(getKernelMatrix(small, sig, "CellTypeA", "CellTypeB",
                                   verbose = FALSE)))
})

test_that("method = 'sparse' still selects the float64 representation", {
  obj <- .sigma_obj_single(n_per_type = 300)
  out <- computeKernelMatrix(obj, sigmaValues = 40, method = "sparse",
                             verbose = FALSE)
  K <- getKernelMatrix(out, 40, "CellTypeA", "CellTypeB", verbose = FALSE)
  expect_true(inherits(K, "sparseMatrix"))
  expect_false(inherits(K, "CoProFloat32SparseMatrix"))
})

test_that("auto warns when the kernel would not actually be sparse", {
  obj <- .sigma_obj_single(n_per_type = 400, extent = 1000)
  # A bandwidth this large couples essentially every pair, so a fixed-radius
  # kernel retains the whole block and saves nothing.
  expect_warning(
    computeKernelMatrix(obj, sigmaValues = 600, method = "auto",
                        autoThreshold = 10L, verbose = FALSE),
    "predicted to be")
})

test_that("auto stays quiet at a bandwidth detection recommends", {
  obj <- .sigma_obj_single(n_per_type = 400, extent = 1000)
  rng <- detectSigmaRange(obj, verbose = FALSE)
  expect_no_warning(
    computeKernelMatrix(obj, sigmaValues = rng$sigmaValues[1],
                        method = "auto", autoThreshold = 10L, verbose = FALSE))
})

test_that("the density probe tracks the density actually built", {
  obj <- .sigma_obj_single(n_per_type = 400, extent = 1000)
  geom <- CoPro:::.resolveDistanceGeometry(
    obj, requested = list(distType = "Euclidean2D"),
    what = "test", verbose = FALSE)

  for (sig in c(20, 40, 80)) {
    probe <- CoPro:::.kernelStorageProbe(
      obj, c("CellTypeA", "CellTypeB"), is_multi = FALSE,
      geometry = geom, sigmaValues = sig, lowerLimit = 1e-7)
    built <- computeKernelMatrix(obj, sigmaValues = sig, method = "sparse",
                                 verbose = FALSE)
    K <- getKernelMatrix(built, sig, "CellTypeA", "CellTypeB", verbose = FALSE)
    actual <- length(K@x) / (nrow(K) * ncol(K))
    # A storage decision only needs the right ballpark, and the estimator is
    # deliberately biased toward over-predicting (which favors dense).
    expect_gt(probe$maxDensity, actual * 0.5)
    expect_lt(probe$maxDensity, min(1, actual * 2.5) + 1e-8)
  }
})
