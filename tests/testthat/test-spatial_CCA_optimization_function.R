# Tests for spatial CCA optimization functions

test_that("optimize_bilinear converges for simple case", {
  # Create simple test data
  set.seed(42)
  n_cells <- 50
  n_features <- 10
  
  # Create two cell type matrices
  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features),
    TypeB = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )
  
  # Create a simple kernel matrix
  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)
  
  # Create flat kernel structure
  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = K
  )
  
  # Run optimization
  result <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1, 
                              max_iter = 100, tol = 1e-5)
  
  expect_true(is.list(result))
  expect_true("TypeA" %in% names(result))
  expect_true("TypeB" %in% names(result))
  expect_true(is.matrix(result$TypeA))
  expect_true(is.matrix(result$TypeB))
  expect_equal(ncol(result$TypeA), 1)
  expect_equal(ncol(result$TypeB), 1)
})

test_that("optimize_bilinear works for single cell type (within)", {
  set.seed(42)
  n_cells <- 50
  n_features <- 10
  
  # Create single cell type matrix
  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )
  
  # Create kernel matrix
  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)
  diag(K) <- 0  # Set diagonal to 0 for within-type
  
  # Create flat kernel structure
  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeA" = K
  )
  
  # Run optimization
  result <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1, 
                              max_iter = 100, tol = 1e-5)
  
  expect_true(is.list(result))
  expect_true("TypeA" %in% names(result))
  expect_equal(ncol(result$TypeA), 1)
})

test_that("optimize_bilinear_n computes multiple components", {
  set.seed(42)
  n_cells <- 50
  n_features <- 10
  
  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features),
    TypeB = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )
  
  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)
  
  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = K
  )
  
  # Get first component
  w_list_1 <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1, 
                                max_iter = 100, tol = 1e-5)
  
  # Get additional components
  w_list_n <- optimize_bilinear_n(
    X_list = X_list,
    flat_kernels = flat_kernels,
    sigma = 0.1,
    w_list = w_list_1,
    cellTypesOfInterest = c("TypeA", "TypeB"),
    nCC = 3,
    max_iter = 100,
    tol = 1e-5
  )
  
  expect_equal(ncol(w_list_n$TypeA), 3)
  expect_equal(ncol(w_list_n$TypeB), 3)
})

test_that("weight vectors are normalized", {
  set.seed(42)
  n_cells <- 50
  n_features <- 10
  
  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features),
    TypeB = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )
  
  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)
  
  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = K
  )
  
  result <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1, 
                              max_iter = 100, tol = 1e-5)
  
  # Check unit norm
  norm_A <- sqrt(sum(result$TypeA^2))
  norm_B <- sqrt(sum(result$TypeB^2))
  
  expect_equal(norm_A, 1, tolerance = 1e-6)
  expect_equal(norm_B, 1, tolerance = 1e-6)
})

test_that("optimize_bilinear converges with step_size < 1", {
  set.seed(42)
  n_cells <- 50
  n_features <- 10

  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features),
    TypeB = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )

  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)

  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = K
  )

  # Run with damped step size
  result <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1,
                              max_iter = 500, tol = 1e-5, step_size = 0.5)

  expect_true(is.list(result))
  expect_equal(ncol(result$TypeA), 1)
  expect_equal(ncol(result$TypeB), 1)

  # Weights should still be unit norm
  expect_equal(sqrt(sum(result$TypeA^2)), 1, tolerance = 1e-6)
  expect_equal(sqrt(sum(result$TypeB^2)), 1, tolerance = 1e-6)
})

test_that("step_size = 1 gives same result as default", {
  set.seed(42)
  n_cells <- 50
  n_features <- 10

  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features),
    TypeB = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )

  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)

  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = K
  )

  result_default <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1,
                                      max_iter = 100, tol = 1e-5)
  result_step1 <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1,
                                    max_iter = 100, tol = 1e-5, step_size = 1)

  expect_equal(result_default$TypeA, result_step1$TypeA)
  expect_equal(result_default$TypeB, result_step1$TypeB)
})

test_that("optimize_bilinear_n works with step_size < 1", {
  set.seed(42)
  n_cells <- 50
  n_features <- 10

  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features),
    TypeB = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )

  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)

  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = K
  )

  # Get first component with damped step size
  w_list_1 <- optimize_bilinear(X_list, flat_kernels, sigma = 0.1,
                                max_iter = 500, tol = 1e-5, step_size = 0.5)

  # Get additional components with damped step size
  w_list_n <- optimize_bilinear_n(
    X_list = X_list,
    flat_kernels = flat_kernels,
    sigma = 0.1,
    w_list = w_list_1,
    cellTypesOfInterest = c("TypeA", "TypeB"),
    nCC = 3,
    max_iter = 500,
    tol = 1e-5,
    step_size = 0.5
  )

  expect_equal(ncol(w_list_n$TypeA), 3)
  expect_equal(ncol(w_list_n$TypeB), 3)

  # All columns should be unit norm
  for (col in 1:3) {
    expect_equal(sqrt(sum(w_list_n$TypeA[, col]^2)), 1, tolerance = 1e-6)
    expect_equal(sqrt(sum(w_list_n$TypeB[, col]^2)), 1, tolerance = 1e-6)
  }
})

test_that("invalid step_size values are rejected", {
  set.seed(42)
  n_cells <- 50
  n_features <- 10

  X_list <- list(
    TypeA = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features),
    TypeB = matrix(rnorm(n_cells * n_features), nrow = n_cells, ncol = n_features)
  )

  dist_mat <- as.matrix(dist(matrix(runif(n_cells * 2), ncol = 2)))
  K <- exp(-0.5 * (dist_mat / 0.1)^2)

  flat_kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = K
  )

  expect_error(optimize_bilinear(X_list, flat_kernels, sigma = 0.1, step_size = 0),
               "step_size must be a single numeric value")
  expect_error(optimize_bilinear(X_list, flat_kernels, sigma = 0.1, step_size = -1),
               "step_size must be a single numeric value")
  expect_error(optimize_bilinear(X_list, flat_kernels, sigma = 0.1, step_size = 2),
               "step_size must be a single numeric value")
  expect_error(optimize_bilinear(X_list, flat_kernels, sigma = 0.1, step_size = "a"),
               "step_size must be a single numeric value")
})
