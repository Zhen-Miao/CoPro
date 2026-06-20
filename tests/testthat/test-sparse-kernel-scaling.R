# Scaling smoke test for the sparse kernel path. Confirms that, at a size where
# the dense n x n kernel would be infeasible, the sparse path (a) is selected by
# method = "auto", (b) returns a genuinely sparse dgCMatrix, and (c) completes
# quickly. Skipped on CI to keep CI fast.

test_that("sparse path scales to large data without forming dense matrices", {
  skip_on_ci()

  n <- 22000L  # > default autoThreshold (20000) -> method = "auto" picks sparse
  obj <- create_test_copro_single(n_cells = n, n_genes = 15,
                                  n_cell_types = 1, seed = 1)
  obj <- subsetData(obj, cellTypesOfInterest = "CellTypeA")

  t0 <- Sys.time()
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), method = "auto",
                             distType = "Euclidean2D",
                             normalizeDistance = FALSE, verbose = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  K <- getKernelMatrix(obj, sigma = 0.1, cellType1 = "CellTypeA",
                       cellType2 = "CellTypeA", verbose = FALSE)

  expect_s4_class(K, "dgCMatrix")
  expect_equal(dim(K), c(n, n))
  # The dense kernel would have n^2 = 4.84e8 entries (~3.9 GB of doubles);
  # the sparse kernel must hold only a small fraction of that.
  expect_lt(length(K@x), 0.05 * as.numeric(n) * as.numeric(n))
  # ... and it should finish quickly (well under a minute on CI hardware).
  expect_lt(elapsed, 120)
})

# A real-data end-to-end check (network + piggyback) is intentionally not part
# of the automated suite because the release assets and their structure are not
# guaranteed in headless/offline environments. Real datasets are exercised by
# the package vignettes via copro_download_data().
