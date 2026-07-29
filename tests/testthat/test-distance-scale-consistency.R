# Cross-type and within-type blocks reach normalization through different entry
# points, and each used to derive its own reference length from its own blocks.
# Those references differ whenever cell types differ in abundance or in how
# tightly they colocalize, so an object could hold cross distances and self
# distances on two different units with @distanceScaleFactor describing only
# one of them. See the follow-up to issues #39 / PR #41.
#
# Two mechanisms answer this, and both are covered here:
#   * normalizeMethod = "global" (the default) reads the reference off the
#     cells rather than off the blocks, so the two families agree with no
#     coordination at all -- see the last section.
#   * the pin, for the block-based methods, where the first normalization wins
#     and later steps adopt its factor. Those tests pass normalizeMethod
#     explicitly, since under the default there would be nothing to adopt.

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
                       normalizeMethod = "spacing", verbose = FALSE)
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
                       normalizeMethod = "spacing", verbose = FALSE)
  own <- suppressMessages(computeSelfDistance(d, normalizeDistance = TRUE,
                                              normalizeMethod = "spacing",
                                              distType = "Euclidean2D",
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
                       normalizeMethod = "spacing", verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, normalizeDistance = TRUE,
                                            normalizeMethod = "spacing",
                                            distType = "Euclidean2D",
                                            verbose = FALSE, overwrite = TRUE))

  # The cross blocks the old pin described are gone, so nothing is left on the
  # old unit and the record now names the step that wrote it.
  expect_false(.createDistMatrixName("A", "B", slide = NULL) %in% names(s@distances))
  expect_equal(getDistanceGeometry(s)$source, "computeSelfDistance")
  expect_false(isTRUE(all.equal(s@distanceScaleFactor, d@distanceScaleFactor)))
  # A cleared record also means the geometry may be changed freely again.
  expect_no_error(suppressMessages(
    computeSelfDistance(d, yDistScale = 2.5, normalizeDistance = TRUE,
                        normalizeMethod = "spacing", distType = "Euclidean2D",
                        verbose = FALSE, overwrite = TRUE)
  ))
})

test_that("overwrite = TRUE re-derives without being told to keep normalizing", {
  # The record is cleared so the pin is dropped, but it is read for defaults
  # first: otherwise normalizeDistance falls back to the 1.2.0 default of FALSE
  # and overwrite silently turns normalization off instead of re-deriving it.
  # The test above passes normalizeDistance = TRUE explicitly and cannot see it.
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       normalizeMethod = "spacing", verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, verbose = FALSE, overwrite = TRUE))

  expect_true(getDistanceGeometry(s)$normalizeDistance)
  expect_equal(getDistanceGeometry(s)$normalizeMethod, "spacing")
  expect_true(is.finite(s@distanceScaleFactor) && s@distanceScaleFactor > 0)
  # Genuinely re-derived from the self blocks, not inherited and not 1.
  expect_false(isTRUE(all.equal(s@distanceScaleFactor, d@distanceScaleFactor)))
  expect_false(isTRUE(all.equal(s@distanceScaleFactor, 1)))
})

test_that("a normalization that could not happen is not recorded as one", {
  # Every (slide, cellType) block sits at the <= 5 cell skip threshold, so no
  # block is built and no reference can be measured. Recording the request
  # anyway left normalizeDistance = TRUE beside the untouched scaling_factor of
  # 1, which .pinnedScaleFactor() then served to every later step as a real pin.
  types <- paste0("T", seq_len(4))
  n_slides <- 4
  per_block <- 5
  set.seed(7)
  slides <- rep(paste0("S", seq_len(n_slides)), each = length(types) * per_block)
  labels <- rep(rep(types, each = per_block), times = n_slides)
  n <- length(labels)
  ids <- paste0("c", seq_len(n))
  obj <- newCoProMulti(
    normalizedData = matrix(rnorm(n * 30), n, 30,
                            dimnames = list(ids, paste0("g", seq_len(30)))),
    locationData = data.frame(x = runif(n, 0, 10), y = runif(n, 0, 10),
                              row.names = ids),
    metaData = data.frame(cell_type = labels, slide = slides, row.names = ids),
    cellTypes = labels, slideID = slides
  )
  obj <- subsetData(obj, cellTypesOfInterest = types)

  expect_warning(
    s <- suppressMessages(computeSelfDistance(obj, normalizeDistance = TRUE,
                                              verbose = FALSE)),
    "no valid reference"
  )
  expect_length(s@distances, 0)
  expect_false(getDistanceGeometry(s)$normalizeDistance)
  expect_null(.pinnedScaleFactor(s))
})

test_that("computeSparseKernelFloat32 leaves a pinned factor alone", {
  # It wrote the slot unconditionally, so a non-normalizing run overwrote the
  # factor computeDistance() had pinned with its own inert 1 -- and
  # .recoverDistanceScaleFactor() reads that slot to map analysis coordinates
  # back to raw ones for the permutation null.
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)
  expect_gt(d@distanceScaleFactor, 0)

  f32 <- suppressMessages(suppressWarnings(
    computeSparseKernelFloat32(d, sigmaValues = 0.5, verbose = FALSE)
  ))
  expect_equal(f32@distanceScaleFactor, d@distanceScaleFactor)
  expect_equal(.recoverDistanceScaleFactor(f32), d@distanceScaleFactor)
})

test_that("normalizeTarget is validated on every path that accepts it", {
  obj <- .uneven_density_obj()
  for (bad in list(-0.01, 0, NA_real_, Inf, c(0.01, 0.02), "0.01")) {
    expect_error(
      suppressMessages(computeSelfDistance(obj, normalizeDistance = TRUE,
                                           normalizeTarget = bad, verbose = FALSE)),
      "normalizeTarget must be a positive finite scalar"
    )
    expect_error(
      computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                      normalizeTarget = bad, verbose = FALSE),
      "normalizeTarget must be a positive finite scalar"
    )
  }
})

test_that("the remediation message names a remedy that path actually has", {
  # computeSelfDistance() can re-derive, because overwrite clears the record.
  # The kernel steps cannot -- re-deriving would leave the kernels they build
  # on a different unit from the ones already in the object -- so their message
  # must not offer overwrite = TRUE, which there only replaces kernel matrices.
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       normalizeMethod = "spacing", verbose = FALSE)

  self_msg <- paste(capture.output(
    invisible(computeSelfDistance(d, verbose = TRUE)), type = "message"
  ), collapse = "\n")
  expect_match(self_msg, "Pass overwrite = TRUE to re-derive the scale")

  k <- suppressMessages(computeKernelMatrix(obj, sigmaValues = 0.05,
                                            method = "sparse",
                                            normalizeDistance = TRUE,
                                            normalizeMethod = "spacing",
                                            verbose = FALSE))
  kern_msg <- paste(capture.output(
    invisible(computeSelfKernel(k, sigmaValues = 0.05, method = "sparse",
                                verbose = TRUE)), type = "message"
  ), collapse = "\n")
  expect_match(kern_msg, "replaces\\s+kernel matrices only")
  expect_false(grepl("Pass overwrite = TRUE", kern_msg))

  # ... and the claim the message makes is true: overwrite really is a no-op
  # for the factor on the kernel path.
  ks <- suppressMessages(computeSelfKernel(k, sigmaValues = 0.05, method = "sparse",
                                           verbose = FALSE, overwrite = TRUE))
  expect_equal(ks@distanceScaleFactor, k@distanceScaleFactor)
})

test_that("an object predating @distanceGeometry keeps its pinned factor", {
  # .pinnedScaleFactor() gated on the record alone, so an object serialized
  # before 1.2.0 -- valid factor, no record -- re-derived, which is exactly the
  # two-unit bug this machinery exists to prevent. The slot probe has to ask
  # the instance too: slotNames() reads the class definition and answers TRUE
  # for such an object, then @ raises on the missing slot.
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       normalizeMethod = "spacing", verbose = FALSE)
  legacy <- d
  attr(legacy, "distanceGeometry") <- NULL

  expect_true("distanceGeometry" %in% methods::slotNames(legacy))
  expect_false(methods::.hasSlot(legacy, "distanceGeometry"))
  expect_no_error(.getDistanceGeometry(legacy))
  expect_length(.getDistanceGeometry(legacy), 0)

  pin <- .pinnedScaleFactor(legacy)
  expect_equal(pin$factor, d@distanceScaleFactor)

  s <- suppressMessages(suppressWarnings(
    computeSelfDistance(legacy, normalizeDistance = TRUE,
                        normalizeMethod = "spacing", verbose = FALSE)
  ))
  expect_equal(s@distanceScaleFactor, d@distanceScaleFactor)
})

test_that("computeSelfKernel does not overwrite the pinned factor", {
  obj <- .uneven_density_obj()

  k <- suppressMessages(computeKernelMatrix(obj, sigmaValues = 0.05, method = "sparse",
                                            normalizeDistance = TRUE,
                                            normalizeMethod = "spacing",
                                            verbose = FALSE))
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
                                           normalizeDistance = TRUE,
                                           normalizeMethod = "spacing",
                                           verbose = FALSE))
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
                       normalizeMethod = "spacing", verbose = FALSE)
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


# ---------------------------------------------------------------------------
# normalizeMethod = "global": the reference is a property of the tissue, so
# consistency does not depend on the pin, on ordering, or on which blocks a
# step happens to build.
# ---------------------------------------------------------------------------

test_that("global is the default reference when normalizeDistance = TRUE", {
  obj <- .uneven_density_obj()
  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)
  expect_equal(getDistanceGeometry(d)$normalizeMethod, "global")
})

test_that("global gives cross-only and self-only objects the same unit", {
  # Neither object ever sees the other's factor: these are two independent
  # runs on the same coordinates, so agreement here is the estimator's doing,
  # not the pin's.
  obj <- .uneven_density_obj()
  cross <- computeDistance(obj, distType = "Euclidean2D",
                           normalizeDistance = TRUE, verbose = FALSE)
  self <- suppressMessages(computeSelfDistance(obj, distType = "Euclidean2D",
                                               normalizeDistance = TRUE,
                                               verbose = FALSE))

  expect_equal(self@distanceScaleFactor, cross@distanceScaleFactor,
               tolerance = 1e-12)
  # The block-based alternatives really do disagree on this object, which is
  # what makes the equality above worth asserting.
  spacing_cross <- computeDistance(obj, distType = "Euclidean2D",
                                   normalizeDistance = TRUE,
                                   normalizeMethod = "spacing",
                                   verbose = FALSE)@distanceScaleFactor
  spacing_self <- suppressMessages(computeSelfDistance(
    obj, distType = "Euclidean2D", normalizeDistance = TRUE,
    normalizeMethod = "spacing", verbose = FALSE))@distanceScaleFactor
  expect_gt(abs(spacing_self / spacing_cross - 1), 0.1)
})

test_that("global survives overwrite = TRUE, which re-derives from scratch", {
  obj <- .uneven_density_obj()
  coords <- .raw_coords(obj)
  labels <- obj@cellTypesSub

  d <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                       verbose = FALSE)
  s <- suppressMessages(computeSelfDistance(d, normalizeDistance = TRUE,
                                            verbose = FALSE, overwrite = TRUE))

  expect_equal(s@distanceScaleFactor, d@distanceScaleFactor, tolerance = 1e-12)
  expect_equal(.recovered_factor(s, "A", "A", coords, labels),
               d@distanceScaleFactor, tolerance = 1e-8)
})

test_that("the global reference ignores cell-type proportions", {
  # Same coordinates, relabeled: 50/50 versus 90/10. A per-block reference
  # moves with the split because each type's own spacing does; a reference
  # measured on the pooled cells cannot see the labels at all.
  factor_for <- function(frac_a) {
    set.seed(4)
    n <- 400
    ids <- paste0("c", seq_len(n))
    loc <- data.frame(x = runif(n, 0, 20), y = runif(n, 0, 20), row.names = ids)
    labels <- rep(c("A", "B"), c(round(n * frac_a), n - round(n * frac_a)))
    expr <- matrix(rnorm(n * 30), n, 30,
                   dimnames = list(ids, paste0("g", seq_len(30))))
    obj <- newCoProSingle(expr, loc,
                          data.frame(cell_type = labels, row.names = ids),
                          labels)
    obj <- computePCA(subsetData(obj, cellTypesOfInterest = c("A", "B")),
                      nPCA = 5)
    vapply(c("global", "spacing"), function(m) {
      computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE,
                      normalizeMethod = m,
                      verbose = FALSE)@distanceScaleFactor
    }, numeric(1))
  }
  balanced <- factor_for(0.5)
  skewed <- factor_for(0.9)

  expect_equal(skewed[["global"]], balanced[["global"]], tolerance = 1e-12)
  expect_gt(abs(log(skewed[["spacing"]] / balanced[["spacing"]])),
            abs(log(skewed[["global"]] / balanced[["global"]])) + 0.05)
})

test_that("global falls back to percentile for Morphology-Aware distances", {
  skip_if_not_installed("igraph")
  obj <- .uneven_density_obj(n_dense = 60, n_sparse = 60)
  expect_warning(
    d <- computeDistance(obj, distType = "Morphology-Aware",
                         normalizeDistance = TRUE, verbose = FALSE),
    "'global' is not defined"
  )
  expect_equal(getDistanceGeometry(d)$normalizeMethod, "percentile")
})
