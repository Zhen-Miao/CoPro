# Cross-type and within-type blocks reach normalization through different entry
# points, and each used to derive its own reference length from its own blocks.
# Those references differ whenever cell types differ in abundance or in how
# tightly they colocalize, so an object could hold cross distances and self
# distances on two different units with @distanceScaleFactor describing only
# one of them. First normalization now wins: later steps adopt the pinned
# factor. See the follow-up to issues #39 / PR #41.

# Two types with deliberately different within-type spacing, so a self-derived
# reference and a cross-derived reference cannot coincide by accident.
.uneven_density_obj <- function(n_dense = 200, n_sparse = 200, seed = 11) {
  set.seed(seed)
  ids <- c(paste0("a", seq_len(n_dense)), paste0("b", seq_len(n_sparse)))
  loc <- data.frame(
    x = c(runif(n_dense, 0, 3), runif(n_sparse, 0, 30)),
    y = c(runif(n_dense, 0, 3), runif(n_sparse, 0, 30)),
    row.names = ids
  )
  labels <- c(rep("A", n_dense), rep("B", n_sparse))
  expr <- matrix(rnorm(length(ids) * 40), length(ids), 40,
                 dimnames = list(ids, paste0("g", seq_len(40))))
  obj <- newCoProSingle(
    normalizedData = expr, locationData = loc,
    metaData = data.frame(cell_type = labels, row.names = ids),
    cellTypes = labels
  )
  computePCA(subsetData(obj, cellTypesOfInterest = c("A", "B")), nPCA = 6)
}

# The factor a stored block was actually multiplied by, read back from the
# block itself. Uses a high quantile: truncateLowDist floors the smallest
# distances, which moves low quantiles for reasons unrelated to scaling.
.recovered_factor <- function(object, type1, type2, coords, labels) {
  q90 <- function(d) {
    v <- as.numeric(d)
    stats::quantile(v[is.finite(v) & v > 0], 0.9, names = FALSE)
  }
  stored <- object@distances[[.createDistMatrixName(type1, type2, slide = NULL)]]
  raw <- fields::rdist(coords[labels == type1, , drop = FALSE],
                       coords[labels == type2, , drop = FALSE])
  q90(stored) / q90(raw)
}

.raw_coords <- function(object) as.matrix(object@locationDataSub[, c("x", "y")])


test_that("computeSelfDistance adopts the factor computeDistance pinned", {
  obj <- .uneven_density_obj()
  coords <- .raw_coords(obj)
  labels <- obj@cellTypesSub

  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, verbose = FALSE))

  cross <- .recovered_factor(s, "A", "B", coords, labels)
  self <- .recovered_factor(s, "A", "A", coords, labels)
  expect_equal(cross, self, tolerance = 1e-8)
  expect_equal(s@distanceScaleFactor, d@distanceScaleFactor)
  # ... and the cross blocks were not silently rescaled a second time.
  expect_equal(cross, d@distanceScaleFactor, tolerance = 1e-8)
})

test_that("the self-derived factor really would have differed", {
  # Guards the test above: if both references happened to agree, adopting
  # would be a no-op and the test would pass without exercising anything.
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)
  own <- suppressMessages(computeSelfDistance(d, normalizeDistance = TRUE,
                                              verbose = FALSE, overwrite = TRUE))

  expect_gt(abs(own@distanceScaleFactor / d@distanceScaleFactor - 1), 0.1)
})

test_that("normalizeTarget and normalizeMethod are inherited, not re-defaulted", {
  obj <- .uneven_density_obj()
  coords <- .raw_coords(obj)
  labels <- obj@cellTypesSub

  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       normalizeTarget = 0.05, normalizeMethod = "percentile",
                       verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, verbose = FALSE))

  geom <- getDistanceGeometry(s)
  expect_equal(geom$normalizeTarget, 0.05)
  expect_equal(geom$normalizeMethod, "percentile")
  expect_equal(geom$source, "computeDistance")
  expect_equal(.recovered_factor(s, "A", "A", coords, labels),
               d@distanceScaleFactor, tolerance = 1e-8)
})

test_that("computeSelfDistance inherits normalizeDistance = FALSE", {
  obj <- .uneven_density_obj()
  coords <- .raw_coords(obj)
  labels <- obj@cellTypesSub

  d <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, verbose = FALSE))

  expect_false(getDistanceGeometry(s)$normalizeDistance)
  expect_equal(s@distanceScaleFactor, 1)
  expect_equal(.recovered_factor(s, "A", "B", coords, labels), 1, tolerance = 1e-8)
  expect_equal(.recovered_factor(s, "A", "A", coords, labels), 1, tolerance = 1e-8)
})

test_that("asking computeSelfDistance to contradict the record is an error", {
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)

  expect_error(computeSelfDistance(d, normalizeDistance = FALSE, verbose = FALSE),
               "contradict the geometry")
  expect_error(computeSelfDistance(d, yDistScale = 2.5, verbose = FALSE),
               "contradict the geometry")
})

test_that("overwrite = TRUE re-derives the factor from the self blocks", {
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, normalizeDistance = TRUE,
                                            verbose = FALSE, overwrite = TRUE))

  # The cross blocks the old pin described are gone, so nothing is left on the
  # old unit and the record now names the step that wrote it.
  expect_false(.createDistMatrixName("A", "B", slide = NULL) %in% names(s@distances))
  expect_equal(getDistanceGeometry(s)$source, "computeSelfDistance")
  expect_false(isTRUE(all.equal(s@distanceScaleFactor, d@distanceScaleFactor)))
  # A cleared record also means the geometry may be changed freely again.
  expect_no_error(suppressMessages(
    computeSelfDistance(d, yDistScale = 2.5, normalizeDistance = TRUE,
                        verbose = FALSE, overwrite = TRUE)
  ))
})

test_that("computeSelfKernel does not overwrite the pinned factor", {
  obj <- .uneven_density_obj()

  k <- suppressMessages(computeKernelMatrix(obj, sigmaValues = 0.05, method = "sparse",
                                            normalizeDistance = TRUE, verbose = FALSE))
  ks <- suppressMessages(computeSelfKernel(k, sigmaValues = 0.05, method = "sparse",
                                           verbose = FALSE))

  expect_equal(ks@distanceScaleFactor, k@distanceScaleFactor)
  # Cross-type kernels built under the old factor are untouched, so the slot
  # keeps describing them.
  expect_equal(
    as.matrix(getKernelMatrix(ks, 0.05, "A", "B", verbose = FALSE)),
    as.matrix(getKernelMatrix(k, 0.05, "A", "B", verbose = FALSE))
  )
})

test_that("the pin also flows the other way, self-kernel first", {
  obj <- .uneven_density_obj()

  ks <- suppressMessages(computeSelfKernel(obj, sigmaValues = 0.05, method = "sparse",
                                           normalizeDistance = TRUE, verbose = FALSE))
  k <- suppressMessages(computeKernelMatrix(ks, sigmaValues = 0.05, method = "sparse",
                                            verbose = FALSE))

  expect_equal(k@distanceScaleFactor, ks@distanceScaleFactor)
})

test_that("dense and sparse still agree when nothing is pinned", {
  # Adoption must not paper over a genuine dense/sparse disagreement: both
  # start from a bare object here, so both derive their own factor.
  obj <- .uneven_density_obj()

  dense <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                           verbose = FALSE)
  sparse <- suppressMessages(computeKernelMatrix(obj, sigmaValues = 0.05,
                                                 method = "sparse",
                                                 normalizeDistance = TRUE,
                                                 verbose = FALSE))
  expect_equal(sparse@distanceScaleFactor, dense@distanceScaleFactor,
               tolerance = 1e-10)
})

test_that("multi-slide self-distances adopt the pinned factor", {
  obj <- create_test_copro_multi(n_cells_per_slide = 120, n_slides = 2,
                                 n_genes = 30, n_cell_types = 2, seed = 3)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, verbose = FALSE))

  expect_equal(s@distanceScaleFactor, d@distanceScaleFactor)
  expect_true(getDistanceGeometry(s)$normalizeDistance)
  expect_equal(getDistanceGeometry(s)$source, "computeDistance")
})

test_that("a self-distance-only workflow still normalizes on its own blocks", {
  obj <- .uneven_density_obj()
  s <- suppressMessages(computeSelfDistance(obj, normalizeDistance = TRUE,
                                            verbose = FALSE))

  expect_equal(getDistanceGeometry(s)$source, "computeSelfDistance")
  expect_true(is.finite(s@distanceScaleFactor) && s@distanceScaleFactor > 0)
  expect_equal(getDistanceGeometry(s)$normalizeTarget, 0.01)
})

test_that("computeSelfDistance honors an explicit normalizeTarget when unpinned", {
  obj <- .uneven_density_obj()
  a <- suppressMessages(computeSelfDistance(obj, normalizeDistance = TRUE,
                                            verbose = FALSE))
  b <- suppressMessages(computeSelfDistance(obj, normalizeDistance = TRUE,
                                            normalizeTarget = 0.05, verbose = FALSE))

  expect_equal(getDistanceGeometry(b)$normalizeTarget, 0.05)
  expect_equal(b@distanceScaleFactor, a@distanceScaleFactor * 5, tolerance = 1e-8)
})
