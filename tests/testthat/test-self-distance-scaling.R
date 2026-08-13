## Self-distances used to be normalized on their own low-distance percentile,
## independently of the factor computeDistance() derived from the cross-type
## distances. The two differ, so "sigma = s" meant a different physical
## bandwidth for a self-kernel than for a cross-kernel. These tests pin the
## "inherit" mode and the fact that building self-kernels never overwrites the
## recorded factor.
##
## `normalizeDistance = TRUE` no longer derives a second unit on top of a
## recorded one: it adopts the pinned factor, so it cannot disagree with the
## cross-type scale and there is no disagreement warning to assert. That path is
## covered in test-distance-scale-consistency.R. What remains distinctive about
## `"inherit"` is that it refuses to guess when nothing has been pinned, where
## `TRUE` would derive a unit from the within-type blocks.
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

  ## TRUE adopts the pinned factor rather than deriving a second unit, so the
  ## two modes now land on the same scale instead of differing by the ratio of
  ## two independently derived factors.
  expect_equal(median(d_own[finite] / d_inh[finite]), 1, tolerance = 1e-8)

  ## and "inherit" puts self-distances on the same scale as the cross-distances.
  ## Getting unscaled blocks out of an object whose record says TRUE needs
  ## overwrite = TRUE: contradicting the record is otherwise an error, since the
  ## blocks already in the object are on the recorded unit.
  raw <- quiet(computeSelfDistance(obj, normalizeDistance = FALSE,
                                   overwrite = TRUE, verbose = FALSE))
  d_raw <- raw@distances[[key]]
  expect_equal(median(d_inh[finite] / d_raw[finite]), recorded, tolerance = 1e-8)
})

test_that('"inherit" still differs from TRUE when nothing has been pinned', {
  ## The one case where the two modes part company: with no recorded factor
  ## there is nothing to adopt, so TRUE measures the within-type blocks while
  ## "inherit" refuses. Uses "spacing", since "global" reads the unit off the
  ## cells and would agree with the cross-type scale regardless.
  obj <- make_two_type_object()
  obj@distanceScaleFactor <- numeric(0)
  obj@distanceGeometry <- list()

  derived <- quiet(computeSelfDistance(obj, normalizeDistance = TRUE,
                                       normalizeMethod = "spacing",
                                       verbose = FALSE))
  expect_true(is.finite(derived@distanceScaleFactor))
  expect_gt(derived@distanceScaleFactor, 0)

  expect_error(
    quiet(computeSelfDistance(obj, normalizeDistance = "inherit", verbose = FALSE)),
    "reuses the scaling factor recorded"
  )
})

test_that("normalizeDistance = FALSE leaves the distances unscaled", {
  obj <- make_two_type_object()
  ## overwrite = TRUE: see above -- the fixture's record says TRUE, and only a
  ## step that discards the blocks that record describes may contradict it.
  raw <- quiet(computeSelfDistance(obj, normalizeDistance = FALSE,
                                   overwrite = TRUE, verbose = FALSE))
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
  is_self <- vapply(names(with_self_auto@kernelMatrices), function(key) {
    parsed <- .parseKernelMatrixName(key)
    identical(parsed$cellType1, parsed$cellType2)
  }, logical(1))
  expect_true(all(vapply(
    with_self_auto@kernelMatrices[is_self],
    inherits, logical(1), "CoProFloat32SparseMatrix"
  )))
})

test_that("invalid self-kernel bandwidths are pruned from sigmaValues", {
  obj <- make_two_type_object(n = 120)
  dense_obj <- quiet(computeSelfDistance(
    obj, normalizeDistance = "inherit", verbose = FALSE
  ))

  for (method in c("dense", "sparse", "float32")) {
    input <- if (method == "dense") dense_obj else obj
    input@sigmaValues <- c(1e-10, 1)
    out <- quiet(computeSelfKernel(
      input, sigmaValues = c(1e-10, 1), method = method,
      minAveCellNeighor = 1, normalizeDistance = "inherit",
      verbose = FALSE
    ))
    expect_equal(out@sigmaValues, 1, info = method)
    expect_false(any(vapply(names(out@kernelMatrices), function(key) {
      identical(.parseKernelMatrixName(key)$sigma, 1e-10)
    }, logical(1))), info = method)
  }
})

test_that("a narrow self-kernel call preserves a broader cross-kernel grid", {
  broad <- c(0.05, 0.1, 0.2)
  obj <- make_two_type_object(n = 180)
  cross <- quiet(computeKernelMatrix(
    obj, sigmaValues = broad, method = "dense", dropDistances = FALSE,
    minAveCellNeighor = 1, verbose = FALSE
  ))
  dense_input <- quiet(computeSelfDistance(
    cross, normalizeDistance = "inherit", verbose = FALSE
  ))

  routes <- list(
    dense = quiet(computeSelfKernel(
      dense_input, sigmaValues = 0.1, method = "dense",
      minAveCellNeighor = 1, normalizeDistance = "inherit", verbose = FALSE
    )),
    sparse = quiet(computeSelfKernel(
      cross, sigmaValues = 0.1, method = "sparse",
      minAveCellNeighor = 1, normalizeDistance = "inherit", verbose = FALSE
    )),
    float32 = quiet(computeSelfKernel(
      cross, sigmaValues = 0.1, method = "float32",
      minAveCellNeighor = 1, normalizeDistance = "inherit", verbose = FALSE
    ))
  )

  for (method in names(routes)) {
    out <- routes[[method]]
    expect_equal(out@sigmaValues, broad, info = method)
    for (sigma in broad) {
      kernel <- getKernelMatrix(
        out, sigma, "TypeA", "TypeB", verbose = FALSE
      )
      expect_false(is.null(kernel), info = paste(method, sigma))
    }
  }
})

test_that("dense multi-slide self kernels preserve unrequested sigmas", {
  broad <- c(0.05, 0.1, 0.2)
  obj <- create_test_copro_multi(
    n_cells_per_slide = 90, n_slides = 2, n_genes = 15,
    n_cell_types = 2, seed = 813
  )
  obj <- quiet(subsetData(
    obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB")
  ))
  obj <- quiet(computeDistance(
    obj, normalizeDistance = TRUE, verbose = FALSE
  ))
  obj <- quiet(computeKernelMatrix(
    obj, sigmaValues = broad, method = "dense", dropDistances = FALSE,
    minAveCellNeighor = 1, verbose = FALSE
  ))
  obj <- quiet(computeSelfDistance(
    obj, normalizeDistance = "inherit", verbose = FALSE
  ))
  out <- quiet(computeSelfKernel(
    obj, sigmaValues = 0.1, method = "dense",
    minAveCellNeighor = 1, normalizeDistance = "inherit", verbose = FALSE
  ))

  expect_equal(out@sigmaValues, broad)
})
