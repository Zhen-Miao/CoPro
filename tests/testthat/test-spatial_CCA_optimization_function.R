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

test_that("one-cell-type optimization selects largest algebraic eigenvalues", {
  # A four-cycle is a valid symmetric, non-negative, hollow proximity kernel.
  # On the centered PC subspace below, Y has eigenvalues 0 and -2. Power
  # iteration selects -2 because it has largest magnitude, but maximizing w'Yw
  # must select the largest algebraic eigenvalue, 0.
  K <- matrix(c(
    0, 1, 0, 1,
    1, 0, 1, 0,
    0, 1, 0, 1,
    1, 0, 1, 0
  ), 4, 4, byrow = TRUE)
  X <- cbind(
    c(1, -1, 1, -1) / 2,
    c(1, 0, -1, 0) / sqrt(2)
  )
  X_list <- list(TypeA = X)
  kernels <- list("kernel|sigma0.1|TypeA|TypeA" = K)

  result <- optimize_bilinear(X_list, kernels, sigma = 0.1)
  Y <- crossprod(X, K %*% X)

  expect_equal(as.numeric(crossprod(result$TypeA, Y %*% result$TypeA)),
               0, tolerance = 1e-12)
})

test_that("one-cell-type optimization returns all exact eigen axes", {
  set.seed(20260722)
  n <- 60
  p <- 7
  X <- matrix(rnorm(n * p), n, p)
  K0 <- matrix(rnorm(n * n), n, n)
  K <- (K0 + t(K0)) / 2
  X_list <- list(TypeA = X)
  kernels <- list("kernel|sigma0.1|TypeA|TypeA" = K)
  Y <- compute_symmetric_Y(X, K)
  expected <- eigen(Y, symmetric = TRUE)

  first <- optimize_bilinear(X_list, kernels, sigma = 0.1)
  result <- optimize_bilinear_n(
    X_list, kernels, sigma = 0.1, w_list = first,
    cellTypesOfInterest = "TypeA", nCC = 4
  )

  expect_equal(abs(crossprod(result$TypeA,
                             expected$vectors[, 1:4, drop = FALSE])),
               diag(4), tolerance = 1e-10)
  expect_equal(crossprod(result$TypeA), diag(4), tolerance = 1e-10)
})

test_that("one-cell-type eigen solver uses symmetric part and weighted metric", {
  set.seed(20260723)
  p <- 6
  Y <- matrix(rnorm(p * p), p, p)  # deliberately asymmetric
  metric <- runif(p, 0.5, 3)
  Y_resi <- list(TypeA = list(TypeA = Y))
  result <- solve_one_type_eigen(
    Y_resi, "TypeA", nCC = 3, sdev2_list = list(TypeA = metric)
  )

  inv_sqrt_d <- 1 / sqrt(metric)
  Ys <- (Y + t(Y)) / 2
  transformed <- sweep(sweep(Ys, 1, inv_sqrt_d, "*"),
                       2, inv_sqrt_d, "*")
  expected <- eigen(transformed, symmetric = TRUE)$vectors[, 1:3, drop = FALSE]
  mapped <- sweep(result$TypeA, 1, sqrt(metric), "*")

  expect_equal(abs(crossprod(mapped, expected)), diag(3), tolerance = 1e-10)
  expect_equal(crossprod(result$TypeA, metric * result$TypeA),
               diag(3), tolerance = 1e-10)
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

test_that("two-cell-type skrCCA axes equal one exact SVD", {
  set.seed(20260720)
  n_a <- 55
  n_b <- 43
  p <- 8
  X_list <- list(
    TypeA = matrix(rnorm(n_a * p), n_a, p),
    TypeB = matrix(rnorm(n_b * p), n_b, p)
  )
  K <- matrix(rnorm(n_a * n_b), n_a, n_b)
  kernels <- list("kernel|sigma0.1|TypeA|TypeB" = K)
  Y <- crossprod(X_list$TypeA, K %*% X_list$TypeB)
  expected <- svd(Y, nu = 4, nv = 4)

  first <- optimize_bilinear(
    X_list, kernels, sigma = 0.1, max_iter = 500, tol = 1e-10
  )
  result <- optimize_bilinear_n(
    X_list, kernels, sigma = 0.1, w_list = first,
    cellTypesOfInterest = names(X_list), nCC = 4,
    max_iter = 500, tol = 1e-10
  )

  expect_equal(abs(crossprod(result$TypeA, expected$u)), diag(4),
               tolerance = 1e-8)
  expect_equal(abs(crossprod(result$TypeB, expected$v)), diag(4),
               tolerance = 1e-8)
  expect_equal(crossprod(result$TypeA), diag(4), tolerance = 1e-10)
  expect_equal(crossprod(result$TypeB), diag(4), tolerance = 1e-10)
})

test_that("two-cell-type SVD respects weighted CCA constraints", {
  set.seed(20260721)
  n <- 50
  p <- 7
  X_list <- list(
    TypeA = matrix(rnorm(n * p), n, p),
    TypeB = matrix(rnorm(n * p), n, p)
  )
  K <- matrix(rnorm(n * n), n, n)
  kernels <- list("kernel|sigma0.1|TypeA|TypeB" = K)
  metrics <- list(TypeA = runif(p, 0.5, 3), TypeB = runif(p, 0.5, 3))

  first <- optimize_bilinear(
    X_list, kernels, sigma = 0.1, sdev2_list = metrics
  )
  result <- optimize_bilinear_n(
    X_list, kernels, sigma = 0.1, w_list = first,
    cellTypesOfInterest = names(X_list), nCC = 3,
    sdev2_list = metrics
  )

  expect_equal(crossprod(result$TypeA, metrics$TypeA * result$TypeA),
               diag(3), tolerance = 1e-10)
  expect_equal(crossprod(result$TypeB, metrics$TypeB * result$TypeB),
               diag(3), tolerance = 1e-10)
})

test_that("multi-cell-type higher axes use orthogonal projection", {
  set.seed(20260722)
  n <- 45
  p <- 8
  cell_types <- c("TypeA", "TypeB", "TypeC")
  X_list <- setNames(lapply(cell_types, function(x) {
    matrix(rnorm(n * p), n, p)
  }), cell_types)
  kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = matrix(rnorm(n * n), n, n),
    "kernel|sigma0.1|TypeA|TypeC" = matrix(rnorm(n * n), n, n),
    "kernel|sigma0.1|TypeB|TypeC" = matrix(rnorm(n * n), n, n)
  )

  first <- optimize_bilinear(
    X_list, kernels, sigma = 0.1, max_iter = 1000, tol = 1e-8
  )
  result <- optimize_bilinear_n(
    X_list, kernels, sigma = 0.1, w_list = first,
    cellTypesOfInterest = cell_types, nCC = 3,
    max_iter = 1000, tol = 1e-8
  )

  for (ct in cell_types) {
    expect_equal(crossprod(result[[ct]]), diag(3), tolerance = 1e-6)
  }
})

test_that("multiset optimizer supports different feature counts by cell type", {
  set.seed(20260727)
  n <- 36
  dims <- c(TypeA = 5L, TypeB = 7L, TypeC = 9L)
  X_list <- lapply(dims, function(p) matrix(rnorm(n * p), n, p))
  kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = matrix(rnorm(n * n), n, n),
    "kernel|sigma0.1|TypeA|TypeC" = matrix(rnorm(n * n), n, n),
    "kernel|sigma0.1|TypeB|TypeC" = matrix(rnorm(n * n), n, n)
  )

  first <- optimize_bilinear(
    X_list, kernels, sigma = 0.1, max_iter = 1000, tol = 1e-8
  )
  result <- optimize_bilinear_n(
    X_list, kernels, sigma = 0.1, w_list = first,
    cellTypesOfInterest = names(dims), nCC = 3,
    max_iter = 1000, tol = 1e-8
  )

  expect_equal(vapply(result, nrow, integer(1)), dims)
  expect_equal(unname(vapply(result, ncol, integer(1))), rep(3L, 3))
  for (ct in names(dims)) {
    expect_equal(crossprod(result[[ct]]), diag(3), tolerance = 1e-6)
  }
})

test_that("scalePCs is a pure reparametrization for >2 cell types", {
  # scalePCs = TRUE feeds whitened X~ = X diag(sdev)^-1 with sdev2_list = NULL;
  # scalePCs = FALSE feeds raw X with sdev2_list = D = sdev^2. These are the
  # same data in two coordinate systems (w~ = diag(sdev) w), so every canonical
  # axis must be identical up to the map w = diag(sdev)^-1 w~. This holds only
  # because apply_deflation() uses the weighted (oblique) projection for the
  # >2-type weighted case; a rank-one fallback there breaks it on CC2+.
  set.seed(20260724)
  n <- 40
  p <- 6
  cell_types <- c("TypeA", "TypeB", "TypeC")
  X_raw <- setNames(lapply(cell_types, function(x) matrix(rnorm(n * p), n, p)),
                    cell_types)
  kernels <- list(
    "kernel|sigma0.1|TypeA|TypeB" = matrix(rnorm(n * n), n, n),
    "kernel|sigma0.1|TypeA|TypeC" = matrix(rnorm(n * n), n, n),
    "kernel|sigma0.1|TypeB|TypeC" = matrix(rnorm(n * n), n, n)
  )
  D <- setNames(lapply(cell_types, function(x) runif(p, 0.5, 3)), cell_types)
  X_white <- setNames(
    lapply(cell_types, function(ct) sweep(X_raw[[ct]], 2, sqrt(D[[ct]]), "/")),
    cell_types
  )

  run <- function(Xin, sdev2) {
    w1 <- optimize_bilinear(Xin, kernels, sigma = 0.1,
                            max_iter = 20000, tol = 1e-11, sdev2_list = sdev2)
    optimize_bilinear_n(Xin, kernels, sigma = 0.1, w_list = w1,
                        cellTypesOfInterest = cell_types, nCC = 3,
                        max_iter = 20000, tol = 1e-11, sdev2_list = sdev2)
  }
  w_white <- run(X_white, NULL)
  w_raw <- run(X_raw, D)

  # Compare per axis, sign-aligned globally across all cell types (the only sign
  # symmetry of the multi-set objective is flipping every block together).
  for (cc in seq_len(3)) {
    mapped <- unlist(lapply(cell_types, function(ct)
      sweep(w_white[[ct]], 1, sqrt(D[[ct]]), "/")[, cc]))
    raw <- unlist(lapply(cell_types, function(ct) w_raw[[ct]][, cc]))
    s <- sign(sum(mapped * raw))
    expect_equal(mapped, s * raw, tolerance = 1e-6)
  }

  # Raw-coordinate axes satisfy the weighted CCA constraint w' D w = I.
  for (ct in cell_types) {
    expect_equal(crossprod(w_raw[[ct]], D[[ct]] * w_raw[[ct]]),
                 diag(3), tolerance = 1e-8)
  }
})

test_that("multi-slide two-type SVD uses the sum of slide operators", {
  set.seed(20260723)
  p <- 6
  slides <- c("slide1", "slide2")
  cell_types <- c("TypeA", "TypeB")
  X_all <- list(
    slide1 = list(
      TypeA = matrix(rnorm(35 * p), 35, p),
      TypeB = matrix(rnorm(28 * p), 28, p)
    ),
    slide2 = list(
      TypeA = matrix(rnorm(30 * p), 30, p),
      TypeB = matrix(rnorm(32 * p), 32, p)
    )
  )
  K1 <- matrix(rnorm(35 * 28), 35, 28)
  K2 <- matrix(rnorm(30 * 32), 30, 32)
  kernels <- list(
    "kernel|sigma0.1|slide1|TypeA|TypeB" = K1,
    "kernel|sigma0.1|slide2|TypeA|TypeB" = K2
  )
  Y_sum <- crossprod(X_all$slide1$TypeA, K1 %*% X_all$slide1$TypeB) +
    crossprod(X_all$slide2$TypeA, K2 %*% X_all$slide2$TypeB)
  expected <- svd(Y_sum, nu = 3, nv = 3)

  first <- optimize_bilinear_multi_slides(
    X_all, kernels, sigma = 0.1, slides = slides,
    max_iter = 500, tol = 1e-10
  )
  result <- optimize_bilinear_n_multi_slides(
    X_all, kernels, sigma = 0.1, slides = slides, w_list = first,
    cellTypesOfInterest = cell_types, nCC = 3,
    max_iter = 500, tol = 1e-10
  )

  expect_equal(abs(crossprod(result$TypeA, expected$u)), diag(3),
               tolerance = 1e-8)
  expect_equal(abs(crossprod(result$TypeB, expected$v)), diag(3),
               tolerance = 1e-8)
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
