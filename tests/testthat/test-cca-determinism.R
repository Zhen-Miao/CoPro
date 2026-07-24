# The skrCCA optimizers must not depend on the RNG state.
#
# Regression: the block-relaxation path (three or more cell types, and every
# conditional higher axis) initialized its weight vectors with `irlba()`.
# `irlba()` draws a random starting vector -- and consumes the RNG stream even
# when handed one explicitly -- so the initial direction varied between runs,
# between sessions, and between a sequential run and a PSOCK worker. That
# carried through to the converged axis: values moved at the power-iteration
# tolerance and, because a power iteration converges to whichever sign its
# start points at, the SIGN of CC2+ weight vectors flipped at random. Gene
# weights are read directionally, so a random sign is not cosmetic.
#
# Every operator involved is nPC x nPC or a Gram matrix of that size, so exact
# LAPACK factorizations are used instead: deterministic, and cheaper here than
# a Krylov method.

.det_pipeline <- function(n_cell_types = 3, nCC = 2, scalePCs = TRUE,
                          seed = 5) {
  obj <- create_test_copro_single(n_cells = 150,
                                  n_cell_types = n_cell_types, seed = seed)
  cts <- paste0("CellType", LETTERS[seq_len(n_cell_types)])
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computeDistance(obj, normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1),
                             dropDistances = FALSE, verbose = FALSE)
  suppressWarnings(computePCA(obj, nPCA = 8, center = TRUE, scale. = TRUE,
                              scalePCs = scalePCs))
}

# Evaluate `expr` from a specific RNG state.
.from_seed <- function(seed, expr) {
  set.seed(seed)
  force(expr)
}

test_that("initialize_weights_svd is exact and RNG-independent", {
  set.seed(3)
  X <- matrix(rnorm(200 * 6), nrow = 200)
  X_list <- list(A = X)

  runs <- lapply(c(1, 2, 99), function(s)
    .from_seed(s, initialize_weights_svd(X_list, "A")[["A"]]))
  expect_equal(runs[[2]], runs[[1]])
  expect_equal(runs[[3]], runs[[1]])

  # ...and it really is the leading right singular vector
  expect_equal(abs(as.numeric(runs[[1]])),
               abs(svd(X)$v[, 1]), tolerance = 1e-10)
  expect_equal(sum(runs[[1]]^2), 1, tolerance = 1e-12)
})

test_that("initialize_next_component is exact and RNG-independent", {
  set.seed(3)
  Y <- matrix(rnorm(9 * 7), nrow = 9)
  Y_resi <- list(A = list(A = NULL, B = Y), B = list(A = t(Y), B = NULL))

  runs <- lapply(c(1, 2, 99), function(s)
    .from_seed(s, initialize_next_component(Y_resi, c("A", "B"))))
  expect_equal(runs[[2]], runs[[1]])
  expect_equal(runs[[3]], runs[[1]])

  # first type takes the leading LEFT singular vector, the second the RIGHT one
  decomp <- svd(Y)
  expect_equal(abs(as.numeric(runs[[1]][["A"]])), abs(decomp$u[, 1]),
               tolerance = 1e-10)
  expect_equal(abs(as.numeric(runs[[1]][["B"]])), abs(decomp$v[, 1]),
               tolerance = 1e-10)
})

test_that("near-tied singular values no longer make the start random", {
  # The worst case for a randomly started Krylov method: when the top two
  # singular values are nearly equal, different starts converge to entirely
  # different directions rather than differing in the last digits.
  set.seed(7)
  basis <- svd(matrix(rnorm(15 * 15), 15))
  Y <- basis$u %*%
    diag(c(1, 1 - 1e-9, seq(0.9, 0.1, length.out = 13))) %*% t(basis$v)
  Y_resi <- list(A = list(A = NULL, B = Y), B = list(A = t(Y), B = NULL))

  runs <- lapply(1:6, function(s)
    .from_seed(s, initialize_next_component(Y_resi, c("A", "B"))[["A"]]))
  for (i in 2:6) expect_identical(runs[[i]], runs[[1]])
})

test_that("three-type runSkrCCA gives identical weights from any RNG state", {
  skip_on_cran()
  for (scalePCs in c(TRUE, FALSE)) {
    obj <- .det_pipeline(scalePCs = scalePCs)
    fit <- function(s) .from_seed(s, suppressWarnings(
      runSkrCCA(obj, scalePCs = scalePCs, nCC = 2)@skrCCAOut
    ))
    reference <- fit(1)
    # Signs included: a randomly seeded start used to flip CC2+ outright.
    expect_identical(fit(2), reference,
                     info = paste("scalePCs =", scalePCs))
    expect_identical(fit(99), reference,
                     info = paste("scalePCs =", scalePCs))
  }
})

test_that("multi-slide weight initialization is RNG-independent", {
  set.seed(3)
  X_list_all <- list(
    slide1 = list(A = matrix(rnorm(60 * 5), 60)),
    slide2 = list(A = matrix(rnorm(40 * 5), 40))
  )
  runs <- lapply(c(1, 2, 99), function(s)
    .from_seed(s, initialize_weights_multi_slide(X_list_all, "A")[["A"]]))
  expect_identical(runs[[2]], runs[[1]])
  expect_identical(runs[[3]], runs[[1]])

  stacked <- rbind(X_list_all$slide1$A, X_list_all$slide2$A)
  expect_equal(abs(as.numeric(runs[[1]])), abs(svd(stacked)$v[, 1]),
               tolerance = 1e-10)
})

test_that("conditional higher axes are exact and RNG-independent", {
  skip_on_cran()
  obj <- .det_pipeline(n_cell_types = 2, nCC = 3)
  obj <- suppressWarnings(runSkrCCA(obj, scalePCs = TRUE, nCC = 3))
  obj <- computeNormalizedCorrelation(obj)

  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  Wobs <- obj@skrCCAOut[[paste0("sigma_", sigma)]]
  Y0 <- compute_Y_resi(PCmats, obj@kernelMatrices, sigma, cts)
  kinfo <- .get_ncorr_kernel_info(obj@kernelMatrices, sigma, cts)
  df <- obj@normalizedCorrelation[[paste0("sigma_", sigma)]]

  for (k in 2:3) {
    fits <- lapply(c(1, 2, 99), function(s) .from_seed(s, .fitConditionalAxis(
      PCmats, obj@kernelMatrices, sigma, cts, W_lower = Wobs,
      k_minus_1 = k - 1, Y_resi = Y0, kernel_info = kinfo
    )))
    expect_identical(fits[[2]], fits[[1]], info = paste("CC", k))
    expect_identical(fits[[3]], fits[[1]], info = paste("CC", k))

    # The two-type conditional axis is now solved by the same exact SVD that
    # produced the observed axis, so the null and observed statistics agree to
    # numerical precision rather than to the power-iteration tolerance.
    expect_equal(fits[[1]]$ncorr,
                 df$normalizedCorrelation[df$CC_index == k],
                 tolerance = 1e-10, info = paste("CC", k))
  }
})
