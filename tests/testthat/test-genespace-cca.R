# Tests for gene-space average per-slide CCA optimization

# Helper: build synthetic per-slide covariance matrices from expression + kernel
make_slide_covmats <- function(Z_by_slide, K_by_slide, slides, cell_types) {
  C_self <- setNames(vector("list", length(slides)), slides)
  C_cross <- setNames(vector("list", length(slides)), slides)

  for (s in slides) {
    C_self[[s]] <- setNames(vector("list", length(cell_types)), cell_types)
    C_cross[[s]] <- list()

    for (ct in cell_types) {
      Z <- Z_by_slide[[s]][[ct]]
      C_self[[s]][[ct]] <- crossprod(Z) / nrow(Z)
    }

    pairs <- combn(cell_types, 2, simplify = FALSE)
    for (pair in pairs) {
      ct_i <- pair[1]
      ct_j <- pair[2]
      key <- paste0(ct_i, "-", ct_j)
      Z_i <- Z_by_slide[[s]][[ct_i]]
      Z_j <- Z_by_slide[[s]][[ct_j]]
      K_ij <- K_by_slide[[s]][[key]]
      C_cross[[s]][[key]] <- crossprod(Z_i, K_ij %*% Z_j) / sqrt(nrow(Z_i) * nrow(Z_j))
    }
  }

  list(C_self = C_self, C_cross = C_cross)
}

# Helper: generate synthetic data with a planted cross-cell-type spatial signal
make_synthetic_data <- function(n_slides = 3, n_genes = 20, n_cells = 50,
                                cell_types = c("TypeA", "TypeB"),
                                signal_strength = 2,
                                seed = 42) {
  set.seed(seed)
  slides <- paste0("slide", seq_len(n_slides))
  Z_by_slide <- setNames(vector("list", n_slides), slides)
  K_by_slide <- setNames(vector("list", n_slides), slides)

  # Known loading vectors: signal concentrated in first 5 genes
  loading_A <- rep(0, n_genes)
  loading_A[1:5] <- c(1, 0.8, 0.6, 0.4, 0.2)
  loading_A <- loading_A / sqrt(sum(loading_A^2))

  loading_B <- rep(0, n_genes)
  loading_B[1:5] <- c(0.2, 0.4, 0.6, 0.8, 1)
  loading_B <- loading_B / sqrt(sum(loading_B^2))

  # Distinct loading for TypeC: signal in genes 6-10
  loading_C <- rep(0, n_genes)
  if (n_genes >= 10) {
    loading_C[6:10] <- c(1, 0.7, 0.5, 0.3, 0.1)
  } else {
    loading_C[1:min(5, n_genes)] <- rev(loading_A[1:min(5, n_genes)])
  }
  loading_C <- loading_C / sqrt(sum(loading_C^2))

  for (s in slides) {
    Z_by_slide[[s]] <- list()
    K_by_slide[[s]] <- list()

    coords <- matrix(runif(n_cells * 2), ncol = 2)
    dist_mat <- as.matrix(dist(coords))
    K <- exp(-0.5 * (dist_mat / 0.3)^2)

    # Smooth spatial factor shared across cell types
    spatial_factor <- as.numeric(K %*% rnorm(n_cells))
    spatial_factor <- (spatial_factor - mean(spatial_factor)) / sd(spatial_factor)

    for (ct in cell_types) {
      if (ct == cell_types[1]) {
        loading <- loading_A
      } else if (length(cell_types) >= 3 && ct == cell_types[3]) {
        loading <- loading_C
      } else {
        loading <- loading_B
      }
      noise <- matrix(rnorm(n_cells * n_genes), nrow = n_cells, ncol = n_genes)
      Z <- noise + signal_strength * outer(spatial_factor, loading)
      Z <- scale(Z)
      Z[is.nan(Z)] <- 0
      Z_by_slide[[s]][[ct]] <- Z
    }

    pairs <- combn(cell_types, 2, simplify = FALSE)
    for (pair in pairs) {
      key <- paste0(pair[1], "-", pair[2])
      K_by_slide[[s]][[key]] <- K
    }
  }

  covmats <- make_slide_covmats(Z_by_slide, K_by_slide, slides, cell_types)
  list(
    Z_by_slide = Z_by_slide,
    K_by_slide = K_by_slide,
    C_self = covmats$C_self,
    C_cross = covmats$C_cross,
    slides = slides,
    cell_types = cell_types,
    n_genes = n_genes,
    loading_A = loading_A,
    loading_B = loading_B,
    loading_C = loading_C
  )
}

test_that("optimize_genespace_avg_corr converges and recovers planted signal", {
  dat <- make_synthetic_data(seed = 42)

  result <- expect_no_warning(optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 3000,
    tol = 1e-6,
    verbose = FALSE
  ))

  expect_true(is.list(result))
  expect_equal(sort(names(result)), sort(dat$cell_types))

  for (ct in dat$cell_types) {
    expect_true(is.matrix(result[[ct]]))
    expect_equal(nrow(result[[ct]]), dat$n_genes)
    expect_equal(ncol(result[[ct]]), 1)
    expect_equal(sqrt(sum(result[[ct]]^2)), 1, tolerance = 1e-8)
  }

  # Objective should be positive (non-trivial correlation found)
  obj <- .compute_p1b_objective(
    result, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  expect_gt(obj, 0)

  # Signal genes (1:5) should carry more weight than noise genes (6:20)
  for (ct in dat$cell_types) {
    signal_weight <- mean(abs(result[[ct]][1:5, 1]))
    noise_weight <- mean(abs(result[[ct]][6:dat$n_genes, 1]))
    expect_gt(signal_weight, noise_weight)
  }
})

test_that("optimize_genespace_avg_corr_n produces orthogonal components", {
  dat <- make_synthetic_data(seed = 42)

  w1 <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 3000,
    tol = 1e-6,
    verbose = FALSE
  )

  w_all <- optimize_genespace_avg_corr_n(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    w_list = w1,
    nCC = 2,
    max_iter = 3000,
    tol = 1e-6,
    verbose = FALSE
  )

  for (ct in dat$cell_types) {
    expect_equal(ncol(w_all[[ct]]), 2)
    # Both columns should be unit norm
    expect_equal(sqrt(sum(w_all[[ct]][, 1]^2)), 1, tolerance = 1e-8)
    expect_equal(sqrt(sum(w_all[[ct]][, 2]^2)), 1, tolerance = 1e-8)
    # Gram-Schmidt should produce near-exact orthogonality
    dot <- abs(sum(w_all[[ct]][, 1] * w_all[[ct]][, 2]))
    expect_lt(dot, 1e-4)
  }
})

test_that("P1b objective is invariant to per-slide covariance scaling", {
  dat <- make_synthetic_data(seed = 123)

  w1 <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 1000,
    tol = 1e-5,
    verbose = FALSE
  )

  # Scale one slide's covariance matrices by 10x — per-slide sigma
  # normalization should absorb this, leaving the objective unchanged
  C_self_scaled <- dat$C_self
  C_cross_scaled <- dat$C_cross
  scale_factor <- 10
  target_slide <- dat$slides[1]

  for (ct in dat$cell_types) {
    C_self_scaled[[target_slide]][[ct]] <-
      C_self_scaled[[target_slide]][[ct]] * scale_factor
  }
  for (key in names(C_cross_scaled[[target_slide]])) {
    C_cross_scaled[[target_slide]][[key]] <-
      C_cross_scaled[[target_slide]][[key]] * scale_factor
  }

  obj_original <- .compute_p1b_objective(
    w1, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  obj_scaled <- .compute_p1b_objective(
    w1, C_self_scaled, C_cross_scaled, dat$slides, dat$cell_types
  )

  # Uniform scaling of both C_self and C_cross cancels in the
  # per-slide normalized correlation: rho = w'C_cross w / (sigma * sigma)
  expect_equal(obj_original, obj_scaled, tolerance = 1e-10)

  # Stronger test: run optimizer independently on scaled data and check
  # that it recovers the same weights (batch-robustness property)
  w_scaled <- optimize_genespace_avg_corr(
    C_self_slide = C_self_scaled,
    C_cross_slide = C_cross_scaled,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 3000,
    tol = 1e-6,
    verbose = FALSE
  )
  for (ct in dat$cell_types) {
    cosine <- abs(sum(w1[[ct]] * w_scaled[[ct]])) /
      (sqrt(sum(w1[[ct]]^2)) * sqrt(sum(w_scaled[[ct]]^2)))
    expect_gt(cosine, 0.99)
  }
})

test_that("nCC validation works", {
  dat <- make_synthetic_data()

  w1 <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 500,
    tol = 1e-4,
    verbose = FALSE
  )

  expect_error(
    optimize_genespace_avg_corr_n(
      C_self_slide = dat$C_self,
      C_cross_slide = dat$C_cross,
      slides = dat$slides,
      cell_types = dat$cell_types,
      w_list = w1,
      nCC = 1,
      verbose = FALSE
    ),
    "must be greater"
  )
})

test_that("three cell types work correctly", {
  dat <- make_synthetic_data(cell_types = c("TypeA", "TypeB", "TypeC"))

  result <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 1000,
    tol = 1e-5,
    verbose = FALSE
  )

  expect_equal(length(result), 3)
  for (ct in dat$cell_types) {
    expect_equal(nrow(result[[ct]]), dat$n_genes)
    expect_equal(ncol(result[[ct]]), 1)
    expect_equal(sqrt(sum(result[[ct]]^2)), 1, tolerance = 1e-8)
  }

  # Signal recovery: TypeA signal in genes 1:5, TypeC signal in genes 6:10
  signal_A <- mean(abs(result[["TypeA"]][1:5, 1]))
  noise_A <- mean(abs(result[["TypeA"]][11:dat$n_genes, 1]))
  expect_gt(signal_A, noise_A)

  signal_C <- mean(abs(result[["TypeC"]][6:10, 1]))
  noise_C <- mean(abs(result[["TypeC"]][11:dat$n_genes, 1]))
  expect_gt(signal_C, noise_C)

  # Objective should be positive
  obj <- .compute_p1b_objective(
    result, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  expect_gt(obj, 0)
})

test_that("single cell type gives informative error", {
  C_self <- list(slide1 = list(TypeA = diag(5)))
  C_cross <- list(slide1 = list())

  expect_error(
    optimize_genespace_avg_corr(
      C_self_slide = C_self,
      C_cross_slide = C_cross,
      slides = "slide1",
      cell_types = "TypeA",
      verbose = FALSE
    ),
    "at least 2 cell types"
  )
})

test_that("runGeneSpaceCCA integration test with CoProMulti object", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D")
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE)

  obj <- runGeneSpaceCCA(obj, sigma = 0.1, nCC = 2,
                         max_iter = 500, tol = 1e-4, verbose = FALSE)

  # Gene scores populated with no NAs
  expect_gt(length(obj@geneScores), 0)
  gs_key <- names(obj@geneScores)[1]
  expect_equal(ncol(obj@geneScores[[gs_key]]), 2)
  expect_false(any(is.na(obj@geneScores[[gs_key]])))

  # Cell scores populated with no NAs and non-zero variance
  expect_gt(length(obj@cellScores), 0)
  cs_key <- names(obj@cellScores)[1]
  expect_equal(ncol(obj@cellScores[[cs_key]]), 2)
  expect_false(any(is.na(obj@cellScores[[cs_key]])))
  for (cc in seq_len(2)) {
    scores_cc <- obj@cellScores[[cs_key]][, cc]
    expect_gt(sd(scores_cc, na.rm = TRUE), 0)
  }

  # CCA output populated
  expect_true(paste0("sigma_", 0.1) %in% names(obj@skrCCAOut))
})

test_that("runGeneSpaceCCA validates sigma against available values", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D")
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE)

  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.5, verbose = FALSE),
    "not found in object@sigmaValues"
  )
})

test_that("runGeneSpaceCCA validates nCC as integer", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D")
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE)

  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 1.5, verbose = FALSE),
    "positive integer"
  )
})

test_that("runGeneSpaceCCA on CoProSingle gives informative error", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_single(n_cells = 60, n_genes = 30, seed = 42)

  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.1),
    "requires a CoProMulti object"
  )
})

test_that("reverse-key lookup in .get_C_cross works", {
  C_cross_s <- list("B-A" = matrix(1:4, 2, 2))

  result <- .get_C_cross(C_cross_s, "A", "B")
  expect_equal(result, t(C_cross_s[["B-A"]]))

  expect_error(
    .get_C_cross(C_cross_s, "A", "C"),
    "Cross-covariance not found"
  )
})
