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

# Helper: generate synthetic spatial data with planted signal
make_synthetic_data <- function(n_slides = 3, n_genes = 20, n_cells = 50,
                                cell_types = c("TypeA", "TypeB"),
                                seed = 42) {
  set.seed(seed)
  slides <- paste0("slide", seq_len(n_slides))
  Z_by_slide <- setNames(vector("list", n_slides), slides)
  K_by_slide <- setNames(vector("list", n_slides), slides)

  for (s in slides) {
    Z_by_slide[[s]] <- list()
    K_by_slide[[s]] <- list()

    # Shared spatial coordinates per slide
    coords <- matrix(runif(n_cells * 2), ncol = 2)
    dist_mat <- as.matrix(dist(coords))
    K <- exp(-0.5 * (dist_mat / 0.3)^2)

    for (ct in cell_types) {
      Z <- matrix(rnorm(n_cells * n_genes), nrow = n_cells, ncol = n_genes)
      Z <- scale(Z)
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
    n_genes = n_genes
  )
}

test_that("optimize_genespace_avg_corr converges", {
  dat <- make_synthetic_data()

  result <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 3000,
    tol = 1e-6,
    verbose = FALSE
  )

  expect_true(is.list(result))
  expect_equal(sort(names(result)), sort(dat$cell_types))

  for (ct in dat$cell_types) {
    expect_true(is.matrix(result[[ct]]))
    expect_equal(nrow(result[[ct]]), dat$n_genes)
    expect_equal(ncol(result[[ct]]), 1)
    # Unit norm
    expect_equal(sqrt(sum(result[[ct]]^2)), 1, tolerance = 1e-8)
  }
})

test_that("optimize_genespace_avg_corr_n produces orthogonal components", {
  dat <- make_synthetic_data()

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
    # CC1 and CC2 should be approximately orthogonal
    dot <- abs(sum(w_all[[ct]][, 1] * w_all[[ct]][, 2]))
    expect_lt(dot, 0.05)
  }
})

test_that("P1b objective is invariant to slide-level mean shifts", {
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

  obj_original <- .compute_p1b_objective(
    w1, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )

  # Add a large mean shift to one slide's data and recompute covmats
  Z_shifted <- dat$Z_by_slide
  shift <- rep(5, dat$n_genes)
  for (ct in dat$cell_types) {
    Z_shifted[["slide1"]][[ct]] <- sweep(Z_shifted[["slide1"]][[ct]], 2, shift, "+")
    # Re-standardize (which removes the mean shift)
    Z_shifted[["slide1"]][[ct]] <- scale(Z_shifted[["slide1"]][[ct]])
  }

  covmats_shifted <- make_slide_covmats(
    Z_shifted, dat$K_by_slide, dat$slides, dat$cell_types
  )

  w1_shifted <- optimize_genespace_avg_corr(
    C_self_slide = covmats_shifted$C_self,
    C_cross_slide = covmats_shifted$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 1000,
    tol = 1e-5,
    verbose = FALSE
  )

  obj_shifted <- .compute_p1b_objective(
    w1_shifted, covmats_shifted$C_self, covmats_shifted$C_cross,
    dat$slides, dat$cell_types
  )

  # Objectives should be very similar (mean shift absorbed by per-slide standardization)
  expect_equal(obj_original, obj_shifted, tolerance = 0.1)
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

  # nCC <= existing components should error

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
  }
})
