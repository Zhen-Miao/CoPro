# Test self-distance and self-kernel functions
test_that("Self-distance and self-kernel functions work correctly", {
  skip_if_not_installed("CoPro")
  
  # This is a placeholder test - in practice you would need a proper CoPro object
  # with multiple cell types to test these functions
  
  # Test that functions exist and have correct signatures
  expect_true(exists("computeSelfDistance"))
  expect_true(exists("computeSelfKernel"))
  expect_true(exists("getSelfDistMat"))
  expect_true(exists("getSelfKernelMatrix"))
  
  # Test that generics are properly defined
  expect_true(isGeneric("computeSelfDistance"))
  expect_true(isGeneric("computeSelfKernel"))
  
  # Additional tests would require creating mock CoPro objects
  # or using real data, which is beyond the scope of this implementation
})

test_that("sparse multitype self-kernels match dense self-kernels", {
  obj <- create_test_copro_single(n_cells = 260, n_genes = 20,
                                  n_cell_types = 2, seed = 611)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  sigmas <- c(0.05, 0.1)

  dense <- computeSelfDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  dense <- computeSelfKernel(dense, sigmaValues = sigmas, method = "dense",
                             verbose = FALSE)
  sparse <- computeSelfKernel(obj, sigmaValues = sigmas, method = "sparse",
                              distType = "Euclidean2D", verbose = FALSE)

  for (sigma in sigmas) {
    for (ct in obj@cellTypesOfInterest) {
      Kd <- getSelfKernelMatrix(dense, sigma, ct, verbose = FALSE)
      Ks <- getSelfKernelMatrix(sparse, sigma, ct, verbose = FALSE)
      expect_s4_class(Ks, "dsCMatrix")
      expect_true(all(Matrix::diag(Ks) == 0))
      expect_equal(as.matrix(Ks), Kd, tolerance = 1e-8,
                   ignore_attr = TRUE)
    }
  }
})

test_that("self-kernel auto dispatches to sparse without dense distances", {
  obj <- create_test_copro_single(n_cells = 180, n_genes = 15,
                                  n_cell_types = 2, seed = 612)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  # Even below a deliberately high size threshold, missing dense self-distance
  # inputs make the fused sparse route the only one-pass path.
  out <- computeSelfKernel(obj, sigmaValues = 0.1, autoThreshold = 1e9,
                           distType = "Euclidean2D", verbose = FALSE)
  for (ct in obj@cellTypesOfInterest) {
    expect_s4_class(
      getSelfKernelMatrix(out, 0.1, ct, verbose = FALSE),
      "dsCMatrix"
    )
  }
})

test_that("sparse multislide self-kernels match dense self-kernels", {
  obj <- create_test_copro_multi(n_cells_per_slide = 180, n_slides = 2,
                                 n_genes = 15, n_cell_types = 2, seed = 613)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  dense <- computeSelfDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  dense <- computeSelfKernel(dense, sigmaValues = 0.1, method = "dense",
                             verbose = FALSE)
  sparse <- computeSelfKernel(obj, sigmaValues = 0.1, method = "sparse",
                              distType = "Euclidean2D", verbose = FALSE)

  for (slide in getSlideList(obj)) {
    for (ct in obj@cellTypesOfInterest) {
      Kd <- getSelfKernelMatrix(dense, 0.1, ct, slide, verbose = FALSE)
      Ks <- getSelfKernelMatrix(sparse, 0.1, ct, slide, verbose = FALSE)
      expect_s4_class(Ks, "dsCMatrix")
      expect_equal(as.matrix(Ks), Kd, tolerance = 1e-8,
                   ignore_attr = TRUE)
    }
  }
})

test_that("sparse self-kernels persist the distance scaling factor", {
  # The fused sparse self-kernel path builds directly from coordinates, so no
  # dense @distances survive to reconstruct the normalizer from. It must record
  # @distanceScaleFactor itself, or .recoverDistanceScaleFactor() (used by the
  # resampling utilities) silently degrades to NA.
  obj <- create_test_copro_single(n_cells = 260, n_genes = 20,
                                  n_cell_types = 2, seed = 614)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  sparse <- computeSelfKernel(obj, sigmaValues = 0.1, method = "sparse",
                              distType = "Euclidean2D", normalizeDistance = TRUE,
                              verbose = FALSE)
  expect_length(sparse@distanceScaleFactor, 1L)
  expect_true(is.finite(sparse@distanceScaleFactor) &&
                sparse@distanceScaleFactor > 0)

  # The downstream consumer must recover the stored factor, not fall through
  # to NA now that @distances is absent (the fused sparse path keeps none).
  recovered <- CoPro:::.recoverDistanceScaleFactor(sparse)
  expect_false(is.na(recovered))
  expect_equal(recovered, sparse@distanceScaleFactor)

  # Multi-slide path persists it too.
  objm <- create_test_copro_multi(n_cells_per_slide = 180, n_slides = 2,
                                  n_genes = 15, n_cell_types = 2, seed = 615)
  objm <- subsetData(objm, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  sparsem <- computeSelfKernel(objm, sigmaValues = 0.1, method = "sparse",
                               distType = "Euclidean2D", verbose = FALSE)
  expect_length(sparsem@distanceScaleFactor, 1L)
  expect_true(is.finite(sparsem@distanceScaleFactor) &&
                sparsem@distanceScaleFactor > 0)
})
