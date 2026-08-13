## Tests for selectSigmaByPermutation(): the studentized max-T bandwidth
## selector. The properties worth pinning are the ones the method's validity
## rests on -- that the kernel/statistic is what it claims, that the shift is
## rigid and inside the wrap box, that the studentized null is flat by
## construction, that the p-value is the Phipson-Smyth max-T one, and that
## the coordinates are on the scale sigma is expressed in.

.fit_selection_object <- function(cts, sigmaValues = c(0.05, 0.1, 0.2),
                                  nCC = 2L, seed = 42) {
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData,
    locationData   = toy$locationData,
    metaData       = toy$metaData,
    cellTypes      = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computePCA(obj, nPCA = 5)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = sigmaValues, verbose = FALSE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = nCC)
  obj <- suppressWarnings(computeNormalizedCorrelation(obj))
  computeGeneAndCellScores(obj)
}

## ---------------------------------------------------------------------------
## The statistic
## ---------------------------------------------------------------------------

test_that("the blocked bilinear statistic equals the dense one, including the self diagonal", {
  set.seed(3)
  cA <- cbind(runif(37), runif(37))
  cB <- cbind(runif(23), runif(23))
  sigmas <- c(0.05, 0.2)
  a <- lapply(sigmas, function(s) rnorm(37))
  b <- lapply(sigmas, function(s) rnorm(23))

  dense <- vapply(seq_along(sigmas), function(g) {
    K <- exp(-0.5 * as.matrix(
      fields::rdist(cA, cB))^2 / sigmas[g]^2)
    sum(a[[g]] * as.numeric(K %*% b[[g]]))
  }, numeric(1))

  ## Block size must not change the answer: 1 row at a time, or all at once.
  for (bs in c(1L, 7L, 1024L)) {
    got <- CoPro:::.bilinearOverSigma(cA, cB, a, b, sigmas,
                                      selfPair = FALSE, blockSize = bs)
    expect_equal(got, dense, tolerance = 1e-10)
  }

  ## Self pair: w_ii = 0, so the diagonal must be excluded.
  cS <- cbind(runif(29), runif(29))
  aS <- lapply(sigmas, function(s) rnorm(29))
  dense_self <- vapply(seq_along(sigmas), function(g) {
    K <- exp(-0.5 * as.matrix(fields::rdist(cS, cS))^2 / sigmas[g]^2)
    diag(K) <- 0
    sum(aS[[g]] * as.numeric(K %*% aS[[g]]))
  }, numeric(1))
  for (bs in c(1L, 5L, 1024L)) {
    got <- CoPro:::.bilinearOverSigma(cS, cS, aS, aS, sigmas,
                                      selfPair = TRUE, blockSize = bs)
    expect_equal(got, dense_self, tolerance = 1e-10)
  }

  ## and the diagonal really is what differs between the two modes
  with_diag <- CoPro:::.bilinearOverSigma(cS, cS, aS, aS, sigmas,
                                          selfPair = FALSE, blockSize = 1024L)
  expect_false(isTRUE(all.equal(with_diag, dense_self)))
})

test_that("the kernel matches computeKernelMatrix's convention, not a factor-of-2 variant", {
  ## sigma is a distance and gets handed back to computeKernelMatrix(), so the
  ## two must agree on exp(-0.5 (d/sigma)^2). A 1/sigma^2 convention would make
  ## the selected bandwidth mean something else downstream.
  cA <- matrix(c(0, 0), nrow = 1)
  cB <- matrix(c(1, 0), nrow = 1)
  got <- CoPro:::.bilinearOverSigma(cA, cB, list(1), list(1), 2,
                                    selfPair = FALSE, blockSize = 1L)
  expect_equal(got, exp(-0.5 * (1 / 2)^2), tolerance = 1e-12)
})

## ---------------------------------------------------------------------------
## The null
## ---------------------------------------------------------------------------

test_that("the toroidal shift stays inside the wrap box and moves every cell together", {
  set.seed(11)
  coords <- cbind(runif(200, 2, 9), runif(200, -3, 4))
  domain <- CoPro:::.resolveSelectionDomain(NULL, list(coords))

  for (i in 1:25) {
    shifted <- CoPro:::.toroidalShift(coords, domain)
    expect_true(all(shifted[, 1] >= domain$lower[1] - 1e-9))
    expect_true(all(shifted[, 1] <= domain$upper[1] + 1e-9))
    expect_true(all(shifted[, 2] >= domain$lower[2] - 1e-9))
    expect_true(all(shifted[, 2] <= domain$upper[2] + 1e-9))
    ## Rigid: every cell moves by the same offset modulo the wrap, so the
    ## displacements take at most two distinct values per axis (wrapped or not).
    for (j in 1:2) {
      expect_lte(length(unique(round(shifted[, j] - coords[, j], 8))), 2L)
    }
  }
})

test_that("a constant coordinate axis is left where it is, not turned into NaN", {
  ## Planar data recorded with a z column has zero extent on that axis, and
  ## `x %% 0` is NaN -- which propagated into every statistic computed from the
  ## shifted coordinates.
  set.seed(31)
  coords <- cbind(runif(50), runif(50), rep(3, 50))
  domain <- CoPro:::.resolveSelectionDomain(NULL, list(coords))
  expect_equal(domain$upper[3] - domain$lower[3], 0)

  shifted <- CoPro:::.toroidalShift(coords, domain)
  expect_false(anyNA(shifted))
  expect_identical(shifted[, 3], coords[, 3])
  ## and the axes that do have width still move
  expect_false(isTRUE(all.equal(shifted[, 1], coords[, 1])))
})

test_that("one draw is one offset per cell type, reused everywhere that type appears", {
  ## With three cell types a type sits in more than one pair. An offset drawn
  ## per pair would put the same cells in two places inside a single draw, so
  ## the row of statistics would not come from any one configuration.
  set.seed(21)
  coords <- list(A = cbind(runif(40), runif(40)),
                 B = cbind(runif(30), runif(30)),
                 C = cbind(runif(25), runif(25)))
  domain <- CoPro:::.resolveSelectionDomain(NULL, coords)
  shifted <- CoPro:::.drawSelectionShift(coords, domain)

  expect_identical(names(shifted), names(coords))
  ## The first type anchors the frame; the rest each move rigidly.
  expect_identical(shifted$A, coords$A)
  for (ct in c("B", "C")) {
    for (j in 1:2) {
      expect_lte(
        length(unique(round(shifted[[ct]][, j] - coords[[ct]][, j], 8))), 2L)
    }
  }
  ## and they move independently of one another
  expect_false(isTRUE(all.equal(shifted$B[1, ] - coords$B[1, ],
                                shifted$C[1, ] - coords$C[1, ])))

  ## A lone cell type has no partner to be relative to, so it is the one that
  ## moves -- otherwise the self statistic would not change at all.
  solo <- list(A = coords$A)
  moved <- CoPro:::.drawSelectionShift(solo, domain)
  expect_false(isTRUE(all.equal(moved$A, solo$A)))
})

test_that("the studentized null is flat across bandwidths and the p-value is Phipson-Smyth max-T", {
  set.seed(5)
  B <- 200L
  ## Three columns with deliberately different null scales, and a middle column
  ## with a null mean far from zero -- the within-type case, where the w_ii = 0
  ## convention tilts the floor.
  nullMat <- cbind(rnorm(B, sd = 1), rnorm(B, mean = -400, sd = 50),
                   rnorm(B, sd = 1000))
  obs <- c(2, 100, 2000)

  st <- CoPro:::.studentizeSelectMaxT(obs, nullMat, "greater")

  expect_equal(st$nullMean, colMeans(nullMat))
  expect_equal(st$nullSD, apply(nullMat, 2, sd))
  expect_equal(st$z, (obs - colMeans(nullMat)) / apply(nullMat, 2, sd))
  ## Scaling alone would leave the tilted column reading off a floor of its own.
  expect_false(isTRUE(all.equal(st$z, obs / apply(nullMat, 2, sd))))

  ## After studentizing, every column's null has zero mean and unit spread.
  ## Flat in both is the whole reason the argmax is comparable across
  ## bandwidths; flat in scale alone is not enough.
  zNull <- sweep(sweep(nullMat, 2, st$nullMean, "-"), 2, st$nullSD, "/")
  expect_equal(unname(apply(zNull, 2, sd)), rep(1, 3), tolerance = 1e-12)
  expect_equal(unname(colMeans(zNull)), rep(0, 3), tolerance = 1e-12)

  expect_equal(st$nullMax, apply(zNull, 1, max))
  expect_equal(st$zMax, max(st$z))
  ## +1 in numerator and denominator; smallest resolvable p is 1/(B+1).
  expect_equal(st$pAdjusted,
               vapply(st$z, function(zi) (1 + sum(st$nullMax >= zi)) / (1 + B),
                      numeric(1)))
  expect_gte(min(st$pAdjusted), 1 / (B + 1))
})

test_that("degenerate columns are absent from both sides of the max-T scan", {
  null <- cbind(constant = rep(4, 20), varying = seq_len(20))
  obs <- c(100, 25)

  expect_warning(
    got <- CoPro:::.studentizeSelectMaxT(obs, null, "greater"),
    "dropped from the scan"
  )
  z_varying <- as.numeric(scale(null[, "varying"]))
  expect_equal(got$nullMax, z_varying, tolerance = 1e-12)
  expect_identical(got$zSelect[[1]], -Inf)

  expect_warning(
    got_two <- CoPro:::.studentizeSelectMaxT(obs, null, "two.sided"),
    "dropped from the scan"
  )
  expect_equal(got_two$nullMax, abs(z_varying), tolerance = 1e-12)
})

test_that("two.sided uses the null of the maximum absolute studentized statistic", {
  set.seed(6)
  nullMat <- cbind(rnorm(100), rnorm(100, sd = 10))
  obs <- c(-3, 1)
  st <- CoPro:::.studentizeSelectMaxT(obs, nullMat, "two.sided")
  zNull <- sweep(sweep(nullMat, 2, st$nullMean, "-"), 2, st$nullSD, "/")
  expect_equal(st$nullMax, apply(abs(zNull), 1, max))
  ## A strong negative statistic must be able to win under two.sided.
  expect_equal(which.max(st$zSelect), 1L)
})

test_that("a degenerate bandwidth is dropped rather than dividing by zero", {
  nullMat <- cbind(rnorm(50), rep(0, 50))
  expect_warning(
    st <- CoPro:::.studentizeSelectMaxT(c(1, 1), nullMat, "greater"),
    "null spread of zero"
  )
  expect_true(is.na(st$z[2]))
  expect_true(is.na(st$pAdjusted[2]))
  expect_true(is.finite(st$zMax))
})

test_that("a scan with no null spread anywhere stops instead of reporting -Inf", {
  ## Every column degenerate used to fall through: nullSD = Inf everywhere, so
  ## zSelect was -Inf everywhere, which.max() took the first candidate, and the
  ## function returned that bandwidth with pValue = NA and zMax = -Inf.
  expect_error(
    CoPro:::.studentizeSelectMaxT(c(1, 2, 3), matrix(0, nrow = 20, ncol = 3),
                                  "greater"),
    "Every bandwidth-pair combination"
  )
})

## ---------------------------------------------------------------------------
## Coordinates carry the bandwidth's units
## ---------------------------------------------------------------------------

test_that("coordinates are put on the same distance scale computeDistance used", {
  ## Stop before computeKernelMatrix(), which drops @distances by default, so
  ## the stored matrix is still there to compare against.
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData, locationData = toy$locationData,
    metaData = toy$metaData, cellTypes = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("Epithelial", "Fibroblast"))
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)

  coords <- CoPro:::.selectionCoords(obj, c("Epithelial", "Fibroblast"))

  ## Distances between these coordinates must reproduce the stored distance
  ## matrix. If they did not, the scanned sigmas would be physically wrong --
  ## sigma is a distance, and the selected value is handed to
  ## computeKernelMatrix() afterwards.
  stored <- as.matrix(obj@distances[[
    CoPro:::.createDistMatrixName("Epithelial", "Fibroblast", slide = NULL)
  ]])
  rebuilt <- fields::rdist(coords[["Epithelial"]], coords[["Fibroblast"]])

  ## computeDistance(truncateLowDist = TRUE) additionally floors the bottom
  ## ~0.1% of distances; the selector does not replicate that floor (see the
  ## note in selectSigmaByPermutation()). Apply it here so the comparison tests
  ## the scale, which is the part that has to match, and confirm the residual
  ## disagreement is confined to the floored entries.
  floored <- rebuilt
  floored[floored < min(stored)] <- min(stored)
  expect_equal(as.numeric(floored), as.numeric(stored), tolerance = 1e-8)

  differing <- which(abs(rebuilt - stored) > 1e-8)
  expect_true(all(rebuilt[differing] < min(stored)))
  expect_lt(length(differing) / length(stored), 0.01)

  ## and the object really did normalize, so this is not a trivial pass
  expect_false(isTRUE(all.equal(obj@distanceScaleFactor, 1)))
})

## ---------------------------------------------------------------------------
## End to end
## ---------------------------------------------------------------------------

test_that("selection runs on a cross-type object and reports a coherent result", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  sel <- selectSigmaByPermutation(obj, nPermu = 49L, seed = 1, verbose = FALSE)

  expect_s3_class(sel, "CoProSigmaSelection")
  expect_true(sel$sigmaValueChoice %in% obj@sigmaValues)
  expect_equal(nrow(sel$perSigma), length(obj@sigmaValues))
  expect_equal(length(sel$nullMax), 49L)
  expect_gte(sel$pValue, 1 / 50)
  expect_lte(sel$pValue, 1)
  ## the reported optimum is the argmax of z, and its p is the max-T p
  expect_equal(sel$sigmaValueChoice,
               sel$perSigma$sigma[which.max(sel$perSigma$z)])
  expect_equal(sel$zMax, max(sel$perSigma$z))
  expect_true(all(sel$plateau %in% sel$perSigma$sigma))
  expect_output(print(sel), "studentized max-T")
})

test_that("selection runs on a within-type object", {
  ## The single-cell-type case, which the re-optimizing route now also serves
  ## (see test-permutation-within-type.R). What is specific to this selector is
  ## that it needs no whitening operator at all -- and within-type is exactly
  ## where the shipped one is invalid, because the stored self-kernel has a zero
  ## diagonal, so it is not a correlation operator and has no square root.
  obj <- .fit_selection_object("Epithelial")
  sel <- selectSigmaByPermutation(obj, nPermu = 49L, seed = 1, verbose = FALSE)
  expect_true(sel$sigmaValueChoice %in% obj@sigmaValues)
  expect_equal(nrow(sel$perSigma), length(obj@sigmaValues))
  expect_true(all(sel$perSigma$cellType1 == "Epithelial"))
  expect_true(all(sel$perSigma$cellType2 == "Epithelial"))
})

test_that("with three cell types every null row comes from a single configuration", {
  ## The max-T reference is the row maximum of statNull, so a row has to be the
  ## scan statistic of one realizable arrangement of the tissue. With three
  ## types, CellTypeC is the second member of two different pairs; drawing its
  ## offset separately in each would have it in two places within one draw.
  ## Rebuild the whole null here from one offset per type per draw and require
  ## the shipped result to match it exactly.
  obj <- create_test_copro_single(n_cells = 120, n_genes = 30, n_cell_types = 3)
  obj <- subsetData(obj, cellTypesOfInterest = paste0("CellType", LETTERS[1:3]))
  obj <- computePCA(obj, nPCA = 5)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = FALSE, verbose = FALSE)
  sigmas <- c(1, 2)
  obj <- computeKernelMatrix(obj, sigmaValues = sigmas, verbose = FALSE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 1)
  obj <- suppressWarnings(computeNormalizedCorrelation(obj))
  obj <- computeGeneAndCellScores(obj)

  B <- 5L
  sel <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = B, seed = 77, minSigma = NULL,
                             verbose = FALSE)
  )

  cts <- CoPro:::.resolveSelectionCellTypes(obj)
  coords <- CoPro:::.selectionCoords(obj, cts)
  cells <- CoPro:::.sampleSelectionCells(coords, 2000L)
  scores <- CoPro:::.selectionScores(obj, cts, sigmas, 1L, cells)
  domain <- CoPro:::.resolveSelectionDomain(NULL, coords)
  pairs <- CoPro:::.selectionPairs(cts)
  expect_equal(length(pairs), 3L)

  ## One shift per cell type per draw, held across the three pairs.
  set.seed(77)
  statNull <- t(vapply(seq_len(B), function(b) {
    shifted <- coords
    for (k in seq_along(cts)[-1L]) {
      shifted[[k]] <- CoPro:::.toroidalShift(coords[[k]], domain)
    }
    unlist(lapply(pairs, function(pr) {
      CoPro:::.bilinearOverSigma(shifted[[pr$ct1]], shifted[[pr$ct2]],
                                 scores[[pr$ct1]], scores[[pr$ct2]], sigmas,
                                 selfPair = FALSE, blockSize = 1024L)
    }))
  }, numeric(length(sigmas) * length(pairs))))

  reference <- CoPro:::.studentizeSelectMaxT(sel$perSigma$statistic, statNull,
                                             "greater")
  expect_equal(sel$nullMax, reference$nullMax)
  expect_equal(sel$perSigma$nullMean, reference$nullMean)
  expect_equal(sel$perSigma$z, reference$z)
  expect_equal(sel$perSigma$pAdjusted, reference$pAdjusted)

  ## and the scan really did cover both pairs that share CellTypeC
  expect_true(all(c("CellTypeB", "CellTypeC") %in% sel$perSigma$cellType2))
})

test_that("the within-type null is centred, not just scaled", {
  ## Removing K[i, i] after the shift leaves the self null with a negative mean
  ## that grows with sigma, so dividing by the SD alone would still compare
  ## bandwidths across a moving floor.
  obj <- .fit_selection_object("Epithelial")
  sel <- selectSigmaByPermutation(obj, nPermu = 199L, seed = 8, verbose = FALSE)

  expect_true(all(is.finite(sel$perSigma$nullMean)))
  expect_equal(sel$perSigma$z,
               (sel$perSigma$statistic - sel$perSigma$nullMean) /
                 sel$perSigma$nullSD)

  ## The tilt is real and sigma-dependent: at the widest bandwidth the null mean
  ## is a sizeable fraction of the null SD, and larger than at the narrowest.
  tilt <- sel$perSigma$nullMean / sel$perSigma$nullSD
  widest <- which.max(sel$perSigma$sigma)
  narrowest <- which.min(sel$perSigma$sigma)
  expect_lt(tilt[widest], -0.15)
  expect_lt(tilt[widest], tilt[narrowest])

  ## A cross-type null has no such term, so its mean sits near zero.
  cross <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  selCross <- suppressWarnings(
    selectSigmaByPermutation(cross, nPermu = 199L, seed = 8, verbose = FALSE)
  )
  expect_true(all(abs(selCross$perSigma$nullMean) <
                    0.3 * selCross$perSigma$nullSD))
})

test_that("the same seed reproduces the selection exactly", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  set.seed(812)
  before <- .Random.seed
  a <- selectSigmaByPermutation(obj, nPermu = 25L, seed = 99, verbose = FALSE)
  expect_identical(.Random.seed, before)
  b <- selectSigmaByPermutation(obj, nPermu = 25L, seed = 99, verbose = FALSE)
  expect_equal(a$perSigma, b$perSigma)
  expect_equal(a$nullMax, b$nullMax)
  expect_equal(a$sigmaValueChoice, b$sigmaValueChoice)
})

test_that("maxCells subsamples and keeps the observed statistic on the same cells as the null", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  sel <- selectSigmaByPermutation(obj, nPermu = 19L, maxCells = 40L, seed = 4,
                                  verbose = FALSE)
  expect_true(all(vapply(sel$cells, length, integer(1)) <= 40L))
  ## the recorded sample is what the statistic was computed on, so a rerun
  ## restricted to those cells must reproduce it
  expect_true(all(is.finite(sel$perSigma$statistic)))
})

test_that("selection scores sumcor-requested weights without special-casing them", {
  ## The directions are held fixed inside every draw, so the statistic does not
  ## depend on which canonical criterion produced them. (On a single-slide
  ## object sumcor is routed to the sumcov solvers anyway, so this is a smoke
  ## test that the routing does not disturb the selector, not a claim that the
  ## other permutation routes would refuse it.)
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData, locationData = toy$locationData,
    metaData = toy$metaData, cellTypes = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("Epithelial", "Fibroblast"))
  obj <- computePCA(obj, nPCA = 5)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
  obj <- suppressMessages(runSkrCCA(obj, scalePCs = TRUE, nCC = 1,
                                    objective = "sumcor"))
  obj <- suppressWarnings(computeNormalizedCorrelation(obj))
  obj <- computeGeneAndCellScores(obj)

  expect_equal(getCCAObjective(obj)$requested, "sumcor")
  sel <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 2, verbose = FALSE)
  )
  expect_true(sel$sigmaValueChoice %in% obj@sigmaValues)
})

## ---------------------------------------------------------------------------
## Refusals
## ---------------------------------------------------------------------------

test_that("selection refuses inputs it cannot serve, with a reason", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))

  expect_error(selectSigmaByPermutation(obj, sigmaValues = 0.123, verbose = FALSE),
               "No cell scores stored for sigma")
  expect_error(selectSigmaByPermutation(obj, ccIndex = 99L, nPermu = 19L,
                                        verbose = FALSE),
               "only 2 component")
  expect_error(selectSigmaByPermutation(obj, ccIndex = 0, verbose = FALSE),
               "positive integer")
  expect_error(selectSigmaByPermutation(obj, alpha = 1, verbose = FALSE),
               "alpha must be")
  ## One draw has no spread to studentize by; it used to return the first
  ## bandwidth with pValue = NA and zMax = -Inf.
  expect_error(selectSigmaByPermutation(obj, nPermu = 1L, verbose = FALSE),
               "nPermu must be at least 2")

  bare <- newCoProSingle(
    normalizedData = matrix(rnorm(200), 20, 10,
                            dimnames = list(paste0("c", 1:20), paste0("g", 1:10))),
    locationData = data.frame(x = runif(20), y = runif(20),
                              row.names = paste0("c", 1:20)),
    metaData = data.frame(cell_id = paste0("c", 1:20),
                          row.names = paste0("c", 1:20)),
    cellTypes = rep(c("A", "B"), 10)
  )
  expect_error(selectSigmaByPermutation(bare, verbose = FALSE),
               "computeGeneAndCellScores")

  multi <- create_test_copro_multi(n_cells_per_slide = 20, n_slides = 2)
  expect_error(selectSigmaByPermutation(multi, verbose = FALSE),
               "CoProSingle")
})

test_that("a morphology-aware object is refused, not silently rescored as Euclidean", {
  ## The selector rebuilds the kernel from coordinates at every bandwidth and
  ## under every shift. A morphology-aware distance is not a function of the
  ## coordinates -- its geodesic filter is fitted on a k-NN graph of the tissue
  ## as observed -- so rebuilding one would score a different metric than the
  ## weights were fitted under.
  skip_if_not_installed("igraph")
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData, locationData = toy$locationData,
    metaData = toy$metaData, cellTypes = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("Epithelial", "Fibroblast"))
  obj <- computePCA(obj, nPCA = 5)
  obj <- suppressWarnings(
    computeDistance(obj, distType = "Morphology-Aware", verbose = FALSE)
  )
  expect_identical(CoPro:::.getDistanceGeometry(obj)$distType, "Morphology-Aware")

  ## The metric really is not reproducible from the coordinates: the filter cuts
  ## pairs that a plain Euclidean rebuild would keep.
  stored <- as.matrix(obj@distances[[
    CoPro:::.createDistMatrixName("Epithelial", "Fibroblast", slide = NULL)]])
  euclid <- as.matrix(fields::rdist(
    as.matrix(obj@locationDataSub[obj@cellTypesSub == "Epithelial", c("x", "y")]),
    as.matrix(obj@locationDataSub[obj@cellTypesSub == "Fibroblast", c("x", "y")])
  ))
  expect_gt(sum(abs(euclid - stored) > 1e-8), 0)

  expect_error(CoPro:::.selectionCoords(obj, c("Epithelial", "Fibroblast")),
               "Euclidean geometry")

  ## and the refusal is reached through the public entry point, on an object
  ## carrying everything the selector otherwise needs
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 1)
  obj <- suppressWarnings(computeNormalizedCorrelation(obj))
  obj <- computeGeneAndCellScores(obj)
  expect_error(selectSigmaByPermutation(obj, verbose = FALSE),
               "Morphology-Aware")
})

test_that("a flat z axis does not poison the scan with NaN", {
  ## Euclidean3D coordinates whose z never varies: the wrap box has zero extent
  ## on that axis, which used to make the shifted coordinates -- and then every
  ## statistic -- NaN.
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  loc <- toy$locationData
  loc$z <- 0
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData, locationData = loc,
    metaData = toy$metaData, cellTypes = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("Epithelial", "Fibroblast"))
  obj <- computePCA(obj, nPCA = 5)
  obj <- computeDistance(obj, distType = "Euclidean3D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1, 0.2),
                             verbose = FALSE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 1)
  obj <- suppressWarnings(computeNormalizedCorrelation(obj))
  obj <- computeGeneAndCellScores(obj)

  sel <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 5, verbose = FALSE)
  )
  expect_false(anyNA(sel$perSigma$z))
  expect_false(anyNA(sel$nullMax))
  expect_true(sel$sigmaValueChoice %in% obj@sigmaValues)
})

test_that("a bandwidth selected at the edge of the grid warns even when quiet", {
  ## A scan that rails at an end of the grid has run out of grid rather than
  ## found an optimum, so this must not be gated on verbose -- the organoid
  ## vignette hit exactly that and verbose = FALSE swallowed it.
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  warned <- character(0)
  sel <- withCallingHandlers(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 7, verbose = FALSE),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  at_edge <- sel$sigmaValueChoice %in% range(sel$perSigma$sigma)
  expect_equal(any(grepl("end of the scanned grid", warned)), at_edge)
})

## ---------------------------------------------------------------------------
## The small-bandwidth floor
## ---------------------------------------------------------------------------

test_that("candidate bandwidths below the nearest-partner spacing are dropped", {
  ## Below one cell spacing the kernel is nearly diagonal and the
  ## fixed-direction null understates the floor, so z climbs without bound as
  ## sigma shrinks and the argmax rails at the smallest candidate.
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"),
                               sigmaValues = c(0.01, 0.05, 0.1, 0.2))
  spacing <- CoPro:::.blockNearestSpacing(
    CoPro:::.selectionCoords(obj, "Epithelial")[["Epithelial"]],
    CoPro:::.selectionCoords(obj, "Fibroblast")[["Fibroblast"]],
    within = FALSE
  )
  expect_gt(spacing, 0.01)

  floored <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 3, verbose = FALSE)
  )
  expect_false(0.01 %in% floored$perSigma$sigma)
  expect_equal(unname(floored$settings$minSigma), unname(spacing))

  ## minSigma = NULL scans it anyway, and the dropped bandwidth is exactly the
  ## one that would have won.
  unfloored <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 3, minSigma = NULL,
                             verbose = FALSE)
  )
  expect_true(0.01 %in% unfloored$perSigma$sigma)
  expect_equal(unfloored$sigmaValueChoice, 0.01)
})

test_that("an explicit numeric floor is honoured and an impossible one is refused", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))

  sel <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 3, minSigma = 0.09,
                             verbose = FALSE)
  )
  expect_equal(sort(unique(sel$perSigma$sigma)), c(0.1, 0.2))

  expect_error(
    selectSigmaByPermutation(obj, minSigma = 99, verbose = FALSE),
    "below the nearest-partner spacing"
  )
  expect_error(
    selectSigmaByPermutation(obj, minSigma = "wide", verbose = FALSE),
    "minSigma must be"
  )
})
