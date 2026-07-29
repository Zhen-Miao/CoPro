## Self-distances used to be normalized on their own low-distance percentile,
## independently of the factor computeDistance() derived from the cross-type
## distances. The two differ, so "sigma = s" meant a different physical
## bandwidth for a self-kernel than for a cross-kernel. These tests pin the
## new "inherit" mode, the warning on the historical mode, and the fact that
## building self-kernels no longer overwrites the recorded factor.
##
## `normalizeDistance` defaults to FALSE as of 1.2.0, and with no scaling there
## is no factor to disagree about, so every test that is about the scaling asks
## for it explicitly.

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

make_two_type_object <- function(n = 300, n_genes = 40, seed = 2) {
  set.seed(seed)
  ids <- paste0("cell_", seq_len(n))
  loc <- data.frame(x = runif(n), y = runif(n), row.names = ids)
  expr <- matrix(rnorm(n * n_genes), n, n_genes,
                 dimnames = list(ids, paste0("Gene", seq_len(n_genes))))
  meta <- data.frame(cell_type = rep(c("TypeA", "TypeB"), each = n / 2),
                     row.names = ids)
  obj <- newCoProSingle(expr, loc, meta, meta$cell_type)
  obj <- quiet(subsetData(obj, cellTypesOfInterest = c("TypeA", "TypeB")))
  obj <- quiet(computePCA(obj, nPCA = 5))
  quiet(computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                        verbose = FALSE))
}

self_key <- function(obj) {
  grep("TypeA[|]TypeA", names(obj@distances), value = TRUE)[1]
}

test_that('normalizeDistance = "inherit" reuses the recorded cross-type factor', {
  obj <- make_two_type_object()
  recorded <- obj@distanceScaleFactor
  expect_true(is.finite(recorded) && recorded > 0)

  own <- quiet(computeSelfDistance(obj, normalizeDistance = TRUE, verbose = FALSE))
  inherited <- quiet(computeSelfDistance(obj, normalizeDistance = "inherit",
                                         verbose = FALSE))

  key <- self_key(inherited)
  expect_false(is.na(key))
  d_own <- own@distances[[key]]
  d_inh <- inherited@distances[[key]]
  finite <- is.finite(d_own) & is.finite(d_inh)

  ## the two modes differ by exactly the ratio of the two scaling factors
  ratio <- median(d_own[finite] / d_inh[finite])
  expect_gt(abs(ratio - 1), 0.05)

  ## and "inherit" puts self-distances on the same scale as the cross-distances
  raw <- quiet(computeSelfDistance(obj, normalizeDistance = FALSE, verbose = FALSE))
  d_raw <- raw@distances[[key]]
  expect_equal(median(d_inh[finite] / d_raw[finite]), recorded, tolerance = 1e-8)
})

test_that("the historical mode warns when it disagrees with the recorded factor", {
  obj <- make_two_type_object()
  expect_warning(
    suppressMessages(computeSelfDistance(obj, normalizeDistance = TRUE,
                                         verbose = FALSE)),
    "differs from the factor computeDistance"
  )
})

test_that("normalizeDistance = FALSE leaves the distances unscaled", {
  obj <- make_two_type_object()
  raw <- quiet(computeSelfDistance(obj, normalizeDistance = FALSE, verbose = FALSE))
  key <- self_key(raw)

  coords <- cbind(obj@locationDataSub$x[obj@cellTypesSub == "TypeA"],
                  obj@locationDataSub$y[obj@cellTypesSub == "TypeA"])
  expected <- as.matrix(dist(coords))
  got <- raw@distances[[key]]
  off <- row(got) != col(got)
  ## truncateLowDist floors the very smallest distances, so compare the bulk
  expect_equal(median(got[off] / expected[off]), 1, tolerance = 1e-6)
})

test_that('"inherit" errors when no factor has been recorded', {
  obj <- make_two_type_object()
  obj@distanceScaleFactor <- numeric(0)
  expect_error(
    quiet(computeSelfDistance(obj, normalizeDistance = "inherit", verbose = FALSE)),
    "reuses the scaling factor recorded"
  )
})

test_that("normalizeDistance rejects values that are neither logical nor 'inherit'", {
  obj <- make_two_type_object()
  expect_error(
    quiet(computeSelfDistance(obj, normalizeDistance = "yes", verbose = FALSE)),
    'must be TRUE, FALSE, or "inherit"'
  )
})

test_that('"inherit" is not read as a contradiction of the recorded geometry', {
  ## `normalizeDistance = "inherit"` says where the scaling factor comes from,
  ## not which coordinates were used, so the geometry record it is checked
  ## against must not treat it as a different basis -- and must not end up
  ## storing the instruction, which a later step would have nothing to inherit
  ## from.
  obj <- make_two_type_object()
  expect_identical(getDistanceGeometry(obj)$normalizeDistance, TRUE)

  expect_no_warning(
    suppressMessages(computeSelfDistance(obj, normalizeDistance = "inherit",
                                         verbose = FALSE))
  )

  with_self <- quiet(computeSelfKernel(obj, sigmaValues = 0.05, method = "sparse",
                                       normalizeDistance = "inherit",
                                       verbose = FALSE))
  expect_identical(getDistanceGeometry(with_self)$normalizeDistance, TRUE)
  expect_true(length(with_self@kernelMatrices) > 0)
})

test_that("building self-kernels does not overwrite the recorded scale factor", {
  ## The sparse self-kernel path derives its own factor. Recording it would
  ## leave the object claiming a scale its cross-kernels were not built at,
  ## and the variogram normalizer reads that slot to size its operators.
  obj <- make_two_type_object()
  recorded <- obj@distanceScaleFactor

  obj <- quiet(computeKernelMatrix(obj, sigmaValues = 0.05, verbose = FALSE))
  expect_equal(obj@distanceScaleFactor, recorded)

  with_self <- quiet(computeSelfKernel(obj, sigmaValues = 0.05, method = "sparse",
                                       verbose = FALSE))
  expect_equal(with_self@distanceScaleFactor, recorded)

  with_self_auto <- quiet(computeSelfKernel(obj, sigmaValues = 0.05, verbose = FALSE))
  expect_equal(with_self_auto@distanceScaleFactor, recorded)
})
