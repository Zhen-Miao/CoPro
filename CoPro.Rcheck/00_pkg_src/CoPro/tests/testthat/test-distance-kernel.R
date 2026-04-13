# Tests for distance and kernel computation

test_that("computeDistance works for two cell types", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  # Compute distance
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  
  expect_true(length(obj@distances) > 0)
  
  # Check distance matrix structure (using flat names)
  dist_names <- names(obj@distances)
  expect_true(any(grepl("CellTypeA", dist_names)))
  expect_true(any(grepl("CellTypeB", dist_names)))
})

test_that("computeDistance works for single cell type", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  
  # Compute distance
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  
  expect_true(length(obj@distances) > 0)
})

test_that("computeDistance normalizes correctly", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  obj_norm <- computeDistance(obj, distType = "Euclidean2D", 
                              normalizeDistance = TRUE, verbose = FALSE)
  obj_no_norm <- computeDistance(obj, distType = "Euclidean2D", 
                                 normalizeDistance = FALSE, verbose = FALSE)
  
  # Get first distance matrix
  dist_name <- names(obj_norm@distances)[1]
  
  # Check that normalized distances are different from non-normalized
  expect_false(identical(obj_norm@distances[[dist_name]], 
                        obj_no_norm@distances[[dist_name]]))
})

test_that("computeDistance requires subsetData to be called first", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  # Don't call subsetData - no cell types of interest set and no locationDataSub
  
  # Without subsetData, either:
  # 1. no cell type of interest specified, or
  # 2. locationDataSub doesn't have x, y columns
  expect_error(
    computeDistance(obj, distType = "Euclidean2D")
  )
})

test_that("computeKernelMatrix works with multiple sigma values", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  
  # Compute kernel with multiple sigma values
  sigmas <- c(0.05, 0.1, 0.2)
  obj <- computeKernelMatrix(obj, sigmaValues = sigmas, verbose = FALSE)
  
  expect_equal(obj@sigmaValues, sigmas)
  expect_true(length(obj@kernelMatrices) > 0)
  
  # Check kernel matrices are created for each sigma
  for (s in sigmas) {
    kernel_name <- paste0("kernel|sigma", s, "|CellTypeA|CellTypeB")
    expect_true(kernel_name %in% names(obj@kernelMatrices) ||
                paste0("kernel|sigma", s, "|CellTypeB|CellTypeA") %in% names(obj@kernelMatrices))
  }
})

test_that("computeKernelMatrix works for single cell type", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  
  expect_true(length(obj@kernelMatrices) > 0)
  kernel_name <- "kernel|sigma0.1|CellTypeA|CellTypeA"
  expect_true(kernel_name %in% names(obj@kernelMatrices))
})

test_that("computeKernelMatrix requires distance computation first", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  expect_error(
    computeKernelMatrix(obj, sigmaValues = c(0.1)),
    "run computeDistance"
  )
})

test_that("getKernelMatrix accessor works", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  
  # Test accessor
  K <- getKernelMatrix(obj, sigma = 0.1, 
                       cellType1 = "CellTypeA", cellType2 = "CellTypeB",
                       verbose = FALSE)
  
  expect_true(is.matrix(K))
  expect_true(nrow(K) > 0)
  expect_true(ncol(K) > 0)
})

test_that("getKernelMatrix handles symmetric access", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  
  # Get both orderings
  K_ab <- getKernelMatrix(obj, sigma = 0.1, 
                          cellType1 = "CellTypeA", cellType2 = "CellTypeB",
                          verbose = FALSE)
  K_ba <- getKernelMatrix(obj, sigma = 0.1, 
                          cellType1 = "CellTypeB", cellType2 = "CellTypeA",
                          verbose = FALSE)
  
  # K_ba should be transpose of K_ab
  expect_equal(K_ba, t(K_ab), tolerance = 1e-10)
})

test_that("getKernelMatrix validates sigma value", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  
  # Invalid sigma
  expect_error(
    getKernelMatrix(obj, sigma = 0.5, 
                    cellType1 = "CellTypeA", cellType2 = "CellTypeB"),
    "not in object sigmaValues"
  )
})

test_that("computeDistance works for CoProMulti", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2, 
                                  n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  
  # Should have distances for each slide
  dist_names <- names(obj@distances)
  expect_true(any(grepl("Slide1", dist_names)))
  expect_true(any(grepl("Slide2", dist_names)))
})

test_that("computeKernelMatrix works for CoProMulti", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2, 
                                  n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  
  # Should have kernels for each slide
  kernel_names <- names(obj@kernelMatrices)
  expect_true(any(grepl("Slide1", kernel_names)))
  expect_true(any(grepl("Slide2", kernel_names)))
})

