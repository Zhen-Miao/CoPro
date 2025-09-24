test_that("plotG12Functions works with valid input", {
  # Skip if required packages not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("spatstat.random")
  
  # Create synthetic test data
  set.seed(123)
  n_cells <- 60
  coords <- data.frame(
    x = runif(n_cells, 0, 100),
    y = runif(n_cells, 0, 100)
  )
  cell_types <- sample(c('TypeA', 'TypeB', 'TypeC'), n_cells, replace = TRUE)
  meta_data <- data.frame(cellID = paste0('Cell_', 1:n_cells))
  
  # Create expression data
  n_genes <- 20
  expr_data <- matrix(rnorm(n_cells * n_genes), nrow = n_cells, ncol = n_genes)
  rownames(expr_data) <- paste0('Cell_', 1:n_cells)
  colnames(expr_data) <- paste0('Gene_', 1:n_genes)
  
  # Create test object
  test_obj <- newCoProSingle(
    normalizedData = expr_data,
    locationData = coords,
    metaData = meta_data,
    cellTypes = cell_types
  )
  
  # Test basic functionality
  expect_no_error({
    result <- plotG12Functions(
      test_obj,
      r_um_range = c(5, 25),
      nsim = 10,  # Use fewer simulations for speed
      verbose = FALSE
    )
  })
  
  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("plot", "data", "summary"))
  
  # Check plot is a ggplot object
  expect_s3_class(result$plot, "ggplot")
  
  # Check data structure
  expect_s3_class(result$data, "data.frame")
  expect_true("r_um" %in% colnames(result$data))
  expect_true("g12_obs" %in% colnames(result$data))
  expect_true("pair_label" %in% colnames(result$data))
  
  # Check summary structure
  expect_s3_class(result$summary, "data.frame")
  expect_true("mean_g12" %in% colnames(result$summary))
  expect_true("max_g12" %in% colnames(result$summary))
})

test_that("plotG12Functions handles edge cases", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  # Create minimal test data
  set.seed(456)
  n_cells <- 10
  coords <- data.frame(
    x = runif(n_cells, 0, 10),
    y = runif(n_cells, 0, 10)
  )
  cell_types <- rep(c('TypeA', 'TypeB'), each = 5)
  meta_data <- data.frame(cellID = paste0('Cell_', 1:n_cells))
  
  expr_data <- matrix(rnorm(n_cells * 10), nrow = n_cells, ncol = 10)
  rownames(expr_data) <- paste0('Cell_', 1:n_cells)
  colnames(expr_data) <- paste0('Gene_', 1:10)
  
  test_obj <- newCoProSingle(
    normalizedData = expr_data,
    locationData = coords,
    metaData = meta_data,
    cellTypes = cell_types
  )
  
  # Test with insufficient points (should give warning but not error)
  expect_warning({
    result <- plotG12Functions(
      test_obj,
      min_points_per_type = 20,  # More than we have
      verbose = FALSE
    )
  })
})

test_that("plotG12Functions parameter validation works", {
  # Test invalid r_um_range
  expect_error(
    plotG12Functions(NULL, r_um_range = c(10, 5)),  # Invalid range
    "r_um_range must be a vector of length 2"
  )
  
  # Test invalid confidence_level
  expect_error(
    plotG12Functions(NULL, confidence_level = 1.5),  # Invalid level
    "confidence_level must be between 0 and 1"
  )
  
  # Test invalid object type
  expect_error(
    plotG12Functions("not_an_object"),
    "object must be a CoProSingle or CoProMulti object"
  )
})
