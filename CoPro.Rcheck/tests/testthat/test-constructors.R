# Tests for CoPro object constructors

test_that("newCoProSingle creates valid object", {
  test_data <- generate_test_data_single(n_cells = 100, n_genes = 50, seed = 42)
  
  obj <- newCoProSingle(
    normalizedData = test_data$normalizedData,
    locationData = test_data$locationData,
    metaData = test_data$metaData,
    cellTypes = test_data$cellTypes
  )
  
  expect_s4_class(obj, "CoProSingle")
  expect_equal(nrow(obj@normalizedData), 100)
  expect_equal(ncol(obj@normalizedData), 50)
  expect_equal(length(obj@cellTypes), 100)
  expect_equal(length(obj@geneList), 50)
})

test_that("newCoProSingle validates dimension mismatch", {
  test_data <- generate_test_data_single(n_cells = 100, n_genes = 50, seed = 42)
  
  # Mismatched cell types length

  expect_error(
    newCoProSingle(
      normalizedData = test_data$normalizedData,
      locationData = test_data$locationData,
      metaData = test_data$metaData,
      cellTypes = test_data$cellTypes[1:50]  # Wrong length
    ),
    "input data do not match dimensionality"
  )
})

test_that("newCoProSingle requires x and y columns", {
  test_data <- generate_test_data_single(n_cells = 100, n_genes = 50, seed = 42)
  
  # Remove x column
  bad_location <- test_data$locationData
  colnames(bad_location) <- c("a", "b")
  
  expect_error(
    newCoProSingle(
      normalizedData = test_data$normalizedData,
      locationData = bad_location,
      metaData = test_data$metaData,
      cellTypes = test_data$cellTypes
    ),
    "x, y"
  )
})

test_that("newCoProSingle validates rownames match", {
  test_data <- generate_test_data_single(n_cells = 100, n_genes = 50, seed = 42)
  
  # Change rownames of metaData
  bad_meta <- test_data$metaData
  rownames(bad_meta) <- paste0("wrong_", seq_len(100))
  
  expect_error(
    newCoProSingle(
      normalizedData = test_data$normalizedData,
      locationData = test_data$locationData,
      metaData = bad_meta,
      cellTypes = test_data$cellTypes
    ),
    "cell barcodes match"
  )
})

test_that("newCoProSingle rejects cell types with pipe characters", {
  test_data <- generate_test_data_single(n_cells = 100, n_genes = 50, seed = 42)
  
  # Add pipe character to cell type
  bad_cell_types <- test_data$cellTypes
  bad_cell_types[1] <- "Type|A"
  
  expect_error(
    newCoProSingle(
      normalizedData = test_data$normalizedData,
      locationData = test_data$locationData,
      metaData = test_data$metaData,
      cellTypes = bad_cell_types
    ),
    "pipe characters"
  )
})

test_that("newCoProMulti creates valid object", {
  test_data <- generate_test_data_multi(
    n_cells_per_slide = 50, 
    n_slides = 2, 
    n_genes = 30, 
    seed = 42
  )
  
  obj <- newCoProMulti(
    normalizedData = test_data$normalizedData,
    locationData = test_data$locationData,
    metaData = test_data$metaData,
    cellTypes = test_data$cellTypes,
    slideID = test_data$slideID
  )
  
  expect_s4_class(obj, "CoProMulti")
  expect_equal(nrow(obj@normalizedData), 100)
  expect_equal(ncol(obj@normalizedData), 30)
  expect_equal(length(obj@slideList), 2)
  expect_true("slideID" %in% colnames(obj@metaData))
})

test_that("newCoProMulti warns with single slide", {
  test_data <- generate_test_data_multi(
    n_cells_per_slide = 50, 
    n_slides = 1,
    n_genes = 30, 
    seed = 42
  )
  
  expect_warning(
    newCoProMulti(
      normalizedData = test_data$normalizedData,
      locationData = test_data$locationData,
      metaData = test_data$metaData,
      cellTypes = test_data$cellTypes,
      slideID = test_data$slideID
    ),
    "only one unique slide ID"
  )
})

test_that("newCoProMulti checks for unique cell IDs", {
  # Test that the function properly validates cell IDs
  # Note: R itself prevents duplicate row names in data.frames,

  # but the function still checks for duplicates in normalizedData
  
  test_data <- generate_test_data_multi(
    n_cells_per_slide = 50, 
    n_slides = 2, 
    n_genes = 30, 
    seed = 42
  )
  
  # Create a matrix with duplicate rownames (matrices allow this)
  old_names <- rownames(test_data$normalizedData)
  new_names <- old_names
  new_names[1] <- new_names[51]  # Create duplicate
  
  bad_matrix <- test_data$normalizedData
  rownames(bad_matrix) <- new_names  # This works for matrices
  
  # The error may come from R's data.frame or from newCoProMulti
  # Either way, duplicate IDs should cause an error
  expect_error(
    newCoProMulti(
      normalizedData = bad_matrix,
      locationData = test_data$locationData,
      metaData = test_data$metaData,
      cellTypes = test_data$cellTypes,
      slideID = test_data$slideID
    )
  )
})

test_that("isMultiSlide returns correct values", {
  single_obj <- create_test_copro_single()
  multi_obj <- create_test_copro_multi()
  
  expect_false(isMultiSlide(single_obj))
  expect_true(isMultiSlide(multi_obj))
})

test_that("getSlideList returns correct values", {
  single_obj <- create_test_copro_single()
  multi_obj <- create_test_copro_multi()
  
  expect_length(getSlideList(single_obj), 0)
  expect_length(getSlideList(multi_obj), 2)
  expect_equal(getSlideList(multi_obj), c("Slide1", "Slide2"))
})

