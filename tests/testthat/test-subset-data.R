# Tests for data subsetting functionality

test_that("subsetData works for CoProSingle with two cell types", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 3, seed = 42)
  
  # Subset to two cell types
  obj_sub <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  expect_equal(obj_sub@cellTypesOfInterest, c("CellTypeA", "CellTypeB"))
  expect_true(all(obj_sub@cellTypesSub %in% c("CellTypeA", "CellTypeB")))
  expect_equal(nrow(obj_sub@normalizedDataSub), length(obj_sub@cellTypesSub))
  expect_equal(nrow(obj_sub@locationDataSub), length(obj_sub@cellTypesSub))
  expect_equal(nrow(obj_sub@metaDataSub), length(obj_sub@cellTypesSub))
})

test_that("subsetData works for CoProSingle with one cell type", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 3, seed = 42)
  
  # Subset to one cell type
  obj_sub <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  
  expect_equal(obj_sub@cellTypesOfInterest, c("CellTypeA"))
  expect_true(all(obj_sub@cellTypesSub == "CellTypeA"))
})

test_that("subsetData works for CoProMulti", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2, 
                                  n_cell_types = 3, seed = 42)
  
  # Subset to two cell types
  obj_sub <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  expect_equal(obj_sub@cellTypesOfInterest, c("CellTypeA", "CellTypeB"))
  expect_true(all(obj_sub@cellTypesSub %in% c("CellTypeA", "CellTypeB")))
  
  # Check that slideID is still properly set in metadata
  expect_true("slideID" %in% colnames(obj_sub@metaDataSub))
})

test_that("subsetData rejects invalid cell types", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  
  expect_error(
    subsetData(obj, cellTypesOfInterest = c("InvalidType")),
    "not in cellTypes"
  )
})

test_that("subsetData rejects empty cell type selection", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  
  expect_error(
    subsetData(obj, cellTypesOfInterest = character(0)),
    "at least one cell type"
  )
})

test_that("subsetData preserves original data when saveOriginal = TRUE", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  original_nrow <- nrow(obj@normalizedData)
  
  obj_sub <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"), 
                        saveOriginal = TRUE)
  
  expect_equal(nrow(obj_sub@normalizedData), original_nrow)
})

test_that("subsetData clears original data when saveOriginal = FALSE", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  
  obj_sub <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"), 
                        saveOriginal = FALSE)
  
  expect_equal(nrow(obj_sub@normalizedData), 1)  # matrix(0) is 1x1
})

test_that("subsetData rejects if too few cells remain", {
  obj <- create_test_copro_single(n_cells = 20, n_cell_types = 4, seed = 42)
  
  # With 20 cells split across 4 types, some types will have < 10 cells
  # This might pass or fail depending on random distribution
  # Use a type we know will have few cells
  expect_error(
    subsetData(obj, cellTypesOfInterest = c("CellTypeD")),
    "Fewer than"
  )
})

