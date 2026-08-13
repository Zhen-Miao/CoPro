# The coordinate geometry a CoPro object's distances and kernels were built on
# is recorded in @distanceGeometry, and every path that builds coordinates
# defers to it. Before this record existed, computeKernelMatrix(method =
# "sparse", normalizeDistance = TRUE) re-derived distType and the per-axis scales from its own defaults
# and silently built kernels on coordinates computeDistance() was never told
# about. See issue #39.

.two_type_obj <- function(n_cells = 300, seed = 1) {
  obj <- create_test_copro_single(n_cells = n_cells, n_genes = 40,
                                  n_cell_types = 2, seed = seed)
  subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
}

# Same object but with a z coordinate, so distType can differ from the value
# the sparse path would infer from the presence of that column.
.two_type_obj_3d <- function(n_cells = 300, seed = 1) {
  d <- generate_test_data_single(n_cells = n_cells, n_genes = 40,
                                 n_cell_types = 2, seed = seed)
  d$locationData$z <- runif(n_cells, 0, 10)
  obj <- newCoProSingle(normalizedData = d$normalizedData,
                        locationData = d$locationData,
                        metaData = d$metaData,
                        cellTypes = d$cellTypes)
  subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
}


# --- the record itself --------------------------------------------------------

test_that("computeDistance records the geometry it used", {
  obj <- .two_type_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", xDistScale = 1.5,
                       yDistScale = 2.5, normalizeTarget = 0.02,
                       truncateLowDist = FALSE, verbose = FALSE, normalizeDistance = TRUE)

  geom <- getDistanceGeometry(d)
  expect_equal(geom$distType, "Euclidean2D")
  expect_equal(geom$xDistScale, 1.5)
  expect_equal(geom$yDistScale, 2.5)
  expect_true(geom$normalizeDistance)
  expect_equal(geom$normalizeMethod, "global")
  expect_equal(geom$normalizeTarget, 0.02)
  expect_false(geom$truncateLowDist)
  expect_equal(geom$source, "computeDistance")
})

test_that("the geometry record survives dropDistances = TRUE", {
  obj <- .two_type_obj()
  d <- computeDistance(obj, yDistScale = 2.5, verbose = FALSE, normalizeDistance = TRUE)
  k <- computeKernelMatrix(d, sigmaValues = 0.1, method = "dense",
                           dropDistances = TRUE, verbose = FALSE, normalizeDistance = TRUE)

  expect_length(k@distances, 0L)
  expect_equal(getDistanceGeometry(k)$yDistScale, 2.5)
})

test_that("getDistanceGeometry returns NULL before any distance step", {
  expect_null(getDistanceGeometry(.two_type_obj()))
})

test_that("multi-slide computeDistance records geometry and scale factor", {
  obj <- create_test_copro_multi(n_cells_per_slide = 80, n_slides = 2,
                                 n_genes = 30, n_cell_types = 2, seed = 4)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  d <- computeDistance(obj, distType = "Euclidean2D", yDistScale = 2,
                       verbose = FALSE, normalizeDistance = TRUE)

  expect_equal(getDistanceGeometry(d)$yDistScale, 2)
  expect_equal(getDistanceGeometry(d)$source, "computeDistance")
  expect_length(d@distanceScaleFactor, 1L)
  expect_true(is.finite(d@distanceScaleFactor) && d@distanceScaleFactor > 0)
  expect_equal(.recoverDistanceScaleFactor(d), d@distanceScaleFactor)
})


# --- the sparse path now inherits instead of re-defaulting ---------------------

test_that("sparse kernels honor per-axis scales set by computeDistance", {
  obj <- .two_type_obj()
  sigma <- 0.1
  d <- computeDistance(obj, distType = "Euclidean2D", yDistScale = 2.5,
                       verbose = FALSE, normalizeDistance = TRUE)

  dense <- computeKernelMatrix(d, sigmaValues = sigma, method = "dense",
                               dropDistances = FALSE, verbose = FALSE, normalizeDistance = TRUE)
  sparse <- computeKernelMatrix(d, sigmaValues = sigma, method = "sparse",
                                verbose = FALSE, normalizeDistance = TRUE)

  Kd <- getKernelMatrix(dense, sigma, "CellTypeA", "CellTypeB", verbose = FALSE)
  Ks <- getKernelMatrix(sparse, sigma, "CellTypeA", "CellTypeB", verbose = FALSE)
  expect_equal(as.matrix(Ks), as.matrix(Kd), tolerance = 1e-8,
               ignore_attr = TRUE)

  # Guard against the test passing because yDistScale had no effect at all.
  flat <- computeKernelMatrix(
    computeDistance(obj, distType = "Euclidean2D", verbose = FALSE, normalizeDistance = TRUE),
    sigmaValues = sigma, method = "dense", dropDistances = FALSE,
    verbose = FALSE, normalizeDistance = TRUE)
  Kf <- getKernelMatrix(flat, sigma, "CellTypeA", "CellTypeB", verbose = FALSE)
  expect_false(isTRUE(all.equal(as.matrix(Kd), as.matrix(Kf),
                                tolerance = 1e-8, check.attributes = FALSE)))
})

test_that('cross-type sparse kernels implement normalizeDistance = "inherit"', {
  obj <- .two_type_obj(n_cells = 120)
  pinned <- computeDistance(
    obj, normalizeDistance = TRUE, normalizeMethod = "spacing",
    verbose = FALSE
  )
  factor_before <- pinned@distanceScaleFactor

  inherited <- suppressWarnings(computeSparseKernel(
    pinned, sigmaValues = 0.1, minAveCellNeighor = 1,
    normalizeDistance = "inherit", normalizeMethod = "spacing",
    verbose = FALSE
  ))
  expect_identical(inherited@distanceScaleFactor, factor_before)
  expect_true(getDistanceGeometry(inherited)$normalizeDistance)

  unpinned <- .two_type_obj(n_cells = 120)
  expect_error(
    computeSparseKernel(
      unpinned, sigmaValues = 1, minAveCellNeighor = 1,
      normalizeDistance = "inherit", verbose = FALSE
    ),
    "has none"
  )
})

test_that("sparse kernels stay 2-D when computeDistance asked for 2-D", {
  obj <- .two_type_obj_3d()
  sigma <- 0.1
  d <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE, normalizeDistance = TRUE)

  dense <- computeKernelMatrix(d, sigmaValues = sigma, method = "dense",
                               dropDistances = FALSE, verbose = FALSE, normalizeDistance = TRUE)
  sparse <- computeKernelMatrix(d, sigmaValues = sigma, method = "sparse",
                                verbose = FALSE, normalizeDistance = TRUE)

  expect_equal(getDistanceGeometry(sparse)$distType, "Euclidean2D")
  expect_equal(
    as.matrix(getKernelMatrix(sparse, sigma, "CellTypeA", "CellTypeB",
                              verbose = FALSE)),
    as.matrix(getKernelMatrix(dense, sigma, "CellTypeA", "CellTypeB",
                              verbose = FALSE)),
    tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("sparse kernels inherit a 3-D geometry too", {
  obj <- .two_type_obj_3d()
  sigma <- 0.1
  d <- computeDistance(obj, distType = "Euclidean3D", zDistScale = 0.5,
                       verbose = FALSE, normalizeDistance = TRUE)

  dense <- computeKernelMatrix(d, sigmaValues = sigma, method = "dense",
                               dropDistances = FALSE, verbose = FALSE, normalizeDistance = TRUE)
  sparse <- computeKernelMatrix(d, sigmaValues = sigma, method = "sparse",
                                verbose = FALSE, normalizeDistance = TRUE)

  expect_equal(getDistanceGeometry(sparse)$distType, "Euclidean3D")
  expect_equal(getDistanceGeometry(sparse)$zDistScale, 0.5)
  expect_equal(
    as.matrix(getKernelMatrix(sparse, sigma, "CellTypeA", "CellTypeB",
                              verbose = FALSE)),
    as.matrix(getKernelMatrix(dense, sigma, "CellTypeA", "CellTypeB",
                              verbose = FALSE)),
    tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("method = 'auto' above the threshold inherits the same geometry", {
  obj <- .two_type_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", yDistScale = 2.5,
                       verbose = FALSE, normalizeDistance = TRUE)
  # autoThreshold below the block size forces the sparse route.
  auto <- computeKernelMatrix(d, sigmaValues = 0.1, method = "auto",
                              autoThreshold = 10L, verbose = FALSE, normalizeDistance = TRUE)
  dense <- computeKernelMatrix(d, sigmaValues = 0.1, method = "dense",
                               dropDistances = FALSE, verbose = FALSE, normalizeDistance = TRUE)

  expect_equal(getDistanceGeometry(auto)$yDistScale, 2.5)
  # method = "auto" resolves to the float32 sparse route above the threshold,
  # so agreement with dense is bounded by float32 precision, not double.
  expect_equal(
    as.matrix(getKernelMatrix(auto, 0.1, "CellTypeA", "CellTypeB",
                              verbose = FALSE)),
    as.matrix(getKernelMatrix(dense, 0.1, "CellTypeA", "CellTypeB",
                              verbose = FALSE)),
    tolerance = 2e-6, ignore_attr = TRUE)
})

test_that("sparse kernels without a prior computeDistance keep their defaults", {
  obj <- .two_type_obj()
  sparse <- computeKernelMatrix(obj, sigmaValues = 0.1, method = "sparse",
                                yDistScale = 2.5, verbose = FALSE, normalizeDistance = TRUE)
  geom <- getDistanceGeometry(sparse)
  expect_equal(geom$yDistScale, 2.5)
  expect_equal(geom$source, "computeKernelMatrix")
})


# --- contradictions are reported, not silently resolved -----------------------

test_that("kernel arguments that contradict the record are an error", {
  obj <- .two_type_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", yDistScale = 2.5,
                       verbose = FALSE, normalizeDistance = TRUE)

  expect_error(
    computeKernelMatrix(d, sigmaValues = 0.1, method = "sparse",
                        yDistScale = 1, verbose = FALSE, normalizeDistance = TRUE),
    "yDistScale")
  expect_error(
    computeKernelMatrix(d, sigmaValues = 0.1, method = "dense",
                        distType = "Euclidean3D", verbose = FALSE, normalizeDistance = TRUE),
    "distType")
  # Restating the recorded value is not a contradiction.
  expect_no_error(
    computeKernelMatrix(d, sigmaValues = 0.1, method = "sparse",
                        yDistScale = 2.5, verbose = FALSE, normalizeDistance = TRUE))
})

test_that("computeSparseKernel warns rather than errors on a mismatch", {
  obj <- .two_type_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", yDistScale = 2.5,
                       verbose = FALSE, normalizeDistance = TRUE)

  expect_warning(
    computeSparseKernel(d, sigmaValues = 0.1, distType = "Euclidean2D",
                        verbose = FALSE, normalizeDistance = TRUE),
    "different coordinates")
  # Called as intended -- instead of computeDistance() -- nothing fires.
  expect_no_warning(
    computeSparseKernel(obj, sigmaValues = 0.1, distType = "Euclidean2D",
                        verbose = FALSE, normalizeDistance = TRUE))
})


# --- downstream helpers read the record ---------------------------------------

test_that(".recoverDistanceScaleFactor is exact under anisotropic scaling", {
  obj <- .two_type_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", yDistScale = 2.5,
                       normalizeDistance = FALSE, verbose = FALSE)
  # normalizeDistance = FALSE leaves @distanceScaleFactor at 1, exercising the
  # coordinate-probing fallback rather than the stored scalar.
  d@distanceScaleFactor <- numeric(0)

  expect_equal(.recoverDistanceScaleFactor(d), 1, tolerance = 1e-10)
})

test_that(".sigmaAwareBins sizes each axis by its own coordinate scale", {
  obj <- .two_type_obj()
  sigma <- 0.5
  # normalizeDistance = FALSE pins the global scale factor at 1, so the bin
  # counts follow directly from the raw extents and the per-axis scales.
  iso <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = FALSE, verbose = FALSE)
  aniso <- computeDistance(obj, distType = "Euclidean2D", yDistScale = 2.5,
                           normalizeDistance = FALSE, verbose = FALSE)

  loc <- as.data.frame(obj@locationDataSub)
  ext_x <- diff(range(loc$x))
  ext_y <- diff(range(loc$y))

  bins_iso <- .sigmaAwareBins(iso, sigma, verbose = FALSE)
  bins_aniso <- .sigmaAwareBins(aniso, sigma, verbose = FALSE)

  expect_equal(bins_iso$num_bins_x, as.integer(floor(ext_x / (2 * sigma))))
  expect_equal(bins_iso$num_bins_y, as.integer(floor(ext_y / (2 * sigma))))

  # Stretching y by 2.5 before distances are taken makes a 2*sigma patch 2.5x
  # narrower in raw y units, so the y axis gains bins and x is untouched.
  expect_equal(bins_aniso$num_bins_x, bins_iso$num_bins_x)
  expect_equal(bins_aniso$num_bins_y,
               as.integer(floor(2.5 * ext_y / (2 * sigma))))

  # Sizing y as if the coordinates were unscaled -- the old behavior -- would
  # give a patch 2.5x too wide on the analysis scale, past the 4*sigma
  # guardrail. The geometry-aware grid stays inside it on both axes.
  expect_lt((ext_y / bins_aniso$num_bins_y) * 2.5, 4 * sigma)
  expect_gt((ext_y / bins_iso$num_bins_y) * 2.5, 4 * sigma)
})


# --- objects predating the slot ------------------------------------------------

test_that("an object with no distanceGeometry slot degrades gracefully", {
  # Emulate an object serialized before the slot existed: a plain list has no
  # such slot, and .getDistanceGeometry() must not error on the missing name.
  legacy <- methods::setClass("LegacyCoProProbe",
                              representation(locationDataSub = "data.frame"))
  probe <- legacy(locationDataSub = data.frame(x = 1:3, y = 1:3))
  expect_equal(.getDistanceGeometry(probe), list())
  expect_equal(unname(.geometryAxisScales(list())[c("x", "y")]), c(1, 1))
  expect_true(is.na(.geometryAxisScales(list())[["z"]]))
})
