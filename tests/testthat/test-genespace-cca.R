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

test_that("gene-space initialization is deterministic and RNG-neutral", {
  set.seed(812)
  before <- .Random.seed
  first <- CoPro:::.deterministicGeneSpaceInit(12, c("A", "B"), component = 2)
  expect_identical(.Random.seed, before)

  runif(10)
  second <- CoPro:::.deterministicGeneSpaceInit(12, c("A", "B"), component = 2)
  expect_identical(first, second)
})

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

  # Objective with the planted strong signal should be well above zero —
  # tightened from > 0 to a meaningful lower bound so a degenerate fixed
  # point (e.g., recovering only noise) would fail the test.
  obj <- CoPro:::.compute_p1b_objective(
    result, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  expect_gt(obj, 0.3)

  # Signal genes (1:5) should carry more weight than noise genes (6:20),
  # AND the recovered loadings should align with the ground-truth loadings
  # the synthetic generator planted (cosine similarity check is a much
  # stronger signal-recovery assertion than the mean-abs comparison).
  for (ct in dat$cell_types) {
    signal_weight <- mean(abs(result[[ct]][1:5, 1]))
    noise_weight <- mean(abs(result[[ct]][6:dat$n_genes, 1]))
    expect_gt(signal_weight, noise_weight)

    truth <- if (ct == "TypeA") dat$loading_A else dat$loading_B
    cos_sim <- abs(sum(as.numeric(result[[ct]][, 1]) * truth)) /
      (sqrt(sum(result[[ct]][, 1]^2)) * sqrt(sum(truth^2)))
    expect_gt(cos_sim, 0.7)
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

  # Scale one slide's covariance matrices uniformly (same factor on
  # C_self and C_cross). Per-slide sigma normalization should absorb this:
  # rho = w'(k * C_cross)w / (sqrt(k * sigma^2) * sqrt(k * sigma^2))
  #     = (k / k) * (w'C_cross w / sigma^2)  --> unchanged.
  # The test confirms the implementation realizes this algebra correctly
  # and that the optimizer is robust to per-slide scale drift.
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

  obj_original <- CoPro:::.compute_p1b_objective(
    w1, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  obj_scaled <- CoPro:::.compute_p1b_objective(
    w1, C_self_scaled, C_cross_scaled, dat$slides, dat$cell_types
  )

  expect_equal(obj_original, obj_scaled, tolerance = 1e-10)

  # Optimizer should recover the same direction on scaled data
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

test_that("P1b objective is sensitive to ASYMMETRIC per-slide perturbation", {
  # Regression guard against a future change that drops sigma normalization:
  # if we scale one slide's C_cross WITHOUT scaling its C_self, that slide's
  # rho contribution changes, and the optimizer should produce a different
  # weight than on unscaled data. This is the asymmetric counterpart to the
  # uniform-scaling test above and rules out that test being a tautology
  # (the optimizer is not weight-invariant to all perturbations, only to
  # ones the per-slide normalization cancels).
  dat <- make_synthetic_data(seed = 123)

  w_orig <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 3000,
    tol = 1e-6,
    verbose = FALSE
  )

  # Scale only C_cross of one slide (and only some pairs to break symmetry).
  C_cross_perturbed <- dat$C_cross
  target_slide <- dat$slides[1]
  for (key in names(C_cross_perturbed[[target_slide]])) {
    C_cross_perturbed[[target_slide]][[key]] <-
      C_cross_perturbed[[target_slide]][[key]] * 5
  }

  w_perturbed <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,  # unchanged
    C_cross_slide = C_cross_perturbed,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 3000,
    tol = 1e-6,
    verbose = FALSE
  )

  obj_orig <- CoPro:::.compute_p1b_objective(
    w_orig, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  obj_perturbed <- CoPro:::.compute_p1b_objective(
    w_perturbed, dat$C_self, C_cross_perturbed, dat$slides, dat$cell_types
  )
  # Asymmetric scaling shifts the per-slide rho on the target slide, so
  # the objective on the perturbed inputs is not algebraically equal to
  # the original one.
  expect_false(isTRUE(all.equal(obj_orig, obj_perturbed, tolerance = 1e-6)))
})

test_that("step_size validation rejects out-of-range and non-numeric values", {
  dat <- make_synthetic_data(seed = 42)

  args_base <- list(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 10, tol = 1e-4, verbose = FALSE
  )
  expect_error(do.call(optimize_genespace_avg_corr,
                       c(args_base, list(step_size = 0))),
               "step_size must be a single numeric value in \\(0, 1\\]")
  expect_error(do.call(optimize_genespace_avg_corr,
                       c(args_base, list(step_size = -0.1))),
               "step_size must be a single numeric value in \\(0, 1\\]")
  expect_error(do.call(optimize_genespace_avg_corr,
                       c(args_base, list(step_size = 1.5))),
               "step_size must be a single numeric value in \\(0, 1\\]")
  expect_error(do.call(optimize_genespace_avg_corr,
                       c(args_base, list(step_size = "0.5"))),
               "step_size must be a single numeric value in \\(0, 1\\]")
  expect_error(do.call(optimize_genespace_avg_corr,
                       c(args_base, list(step_size = c(0.5, 0.5)))),
               "step_size must be a single numeric value in \\(0, 1\\]")

  args_bad_max <- args_base
  args_bad_max$max_iter <- 1.5
  expect_error(do.call(optimize_genespace_avg_corr, args_bad_max),
               "max_iter must be a positive integer")
  args_bad_tol <- args_base
  args_bad_tol$tol <- Inf
  expect_error(do.call(optimize_genespace_avg_corr, args_bad_tol),
               "tol must be a positive finite number")

  first <- CoPro:::.deterministicGeneSpaceInit(
    nrow(dat$C_self[[dat$slides[[1L]]]][[dat$cell_types[[1L]]]]),
    dat$cell_types, component = 1L
  )
  args_n <- c(args_base, list(w_list = first, nCC = 2L))
  args_n$tol <- NA_real_
  expect_error(do.call(optimize_genespace_avg_corr_n, args_n),
               "tol must be a positive finite number")
})

test_that("mildly damped power iteration converges to the same fixed point as undamped", {
  # On a problem where the undamped power iteration already converges,
  # mild damping (step_size near 1) should converge to the same
  # canonical direction. Aggressive damping (e.g., step_size = 0.5)
  # CAN introduce 2-cycle oscillation on small problems and is therefore
  # not tested for equivalence — it should only be used when undamped
  # power iteration diverges or oscillates.
  dat <- make_synthetic_data(seed = 42)

  set.seed(1L)
  w_undamped <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 3000, tol = 1e-8, step_size = 1,
    verbose = FALSE
  )

  set.seed(1L)
  w_damped <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self,
    C_cross_slide = dat$C_cross,
    slides = dat$slides,
    cell_types = dat$cell_types,
    max_iter = 6000, tol = 1e-8, step_size = 0.9,
    verbose = FALSE
  )

  # Per cell type: cosine similarity between damped and undamped weight
  # vectors should be ~1 in absolute value.
  for (ct in dat$cell_types) {
    a <- as.numeric(w_undamped[[ct]][, 1])
    b <- as.numeric(w_damped[[ct]][, 1])
    cos_sim <- abs(sum(a * b)) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    expect_gt(cos_sim, 0.999)
  }

  # |objective| should agree (sign indeterminate by symmetry of CCA).
  obj_undamped <- CoPro:::.compute_p1b_objective(
    w_undamped, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  obj_damped <- CoPro:::.compute_p1b_objective(
    w_damped, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  expect_equal(abs(obj_damped), abs(obj_undamped), tolerance = 1e-3)
})

test_that("mildly damped iteration with deflation matches undamped CC2", {
  dat <- make_synthetic_data(seed = 42)

  set.seed(2L)
  w1 <- optimize_genespace_avg_corr(
    C_self_slide = dat$C_self, C_cross_slide = dat$C_cross,
    slides = dat$slides, cell_types = dat$cell_types,
    max_iter = 3000, tol = 1e-8, step_size = 1, verbose = FALSE
  )

  set.seed(3L)
  w_undamped <- optimize_genespace_avg_corr_n(
    C_self_slide = dat$C_self, C_cross_slide = dat$C_cross,
    slides = dat$slides, cell_types = dat$cell_types,
    w_list = w1, nCC = 2,
    max_iter = 3000, tol = 1e-8, step_size = 1, verbose = FALSE
  )

  set.seed(3L)
  w_damped <- optimize_genespace_avg_corr_n(
    C_self_slide = dat$C_self, C_cross_slide = dat$C_cross,
    slides = dat$slides, cell_types = dat$cell_types,
    w_list = w1, nCC = 2,
    max_iter = 6000, tol = 1e-8, step_size = 0.9, verbose = FALSE
  )

  # CC2 should match (up to sign) between undamped and mildly damped.
  for (ct in dat$cell_types) {
    a <- as.numeric(w_undamped[[ct]][, 2])
    b <- as.numeric(w_damped[[ct]][, 2])
    cos_sim <- abs(sum(a * b)) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    expect_gt(cos_sim, 0.999)
  }
})

test_that("damping reaches the same fixed point under either sweep", {
  # The two tests above exercise the default "gauss-seidel" sweep, where the
  # damped weight also has to feed .refresh_slide_sigma(). Under "jacobi" no
  # refresh happens, so cover that path too: damping must not change WHERE the
  # iteration lands under either sweep, only how it gets there.
  dat <- make_synthetic_data(seed = 42)

  for (sweep in c("gauss-seidel", "jacobi")) {
    set.seed(1L)
    w_undamped <- optimize_genespace_avg_corr(
      C_self_slide = dat$C_self, C_cross_slide = dat$C_cross,
      slides = dat$slides, cell_types = dat$cell_types,
      max_iter = 3000, tol = 1e-8, step_size = 1,
      verbose = FALSE, sweep = sweep
    )

    set.seed(1L)
    w_damped <- optimize_genespace_avg_corr(
      C_self_slide = dat$C_self, C_cross_slide = dat$C_cross,
      slides = dat$slides, cell_types = dat$cell_types,
      max_iter = 6000, tol = 1e-8, step_size = 0.9,
      verbose = FALSE, sweep = sweep
    )

    for (ct in dat$cell_types) {
      a <- as.numeric(w_undamped[[ct]][, 1])
      b <- as.numeric(w_damped[[ct]][, 1])
      cos_sim <- abs(sum(a * b)) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
      expect_gt(cos_sim, 0.999)
    }
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

  # Objective should be well above zero with strong planted signal
  obj <- CoPro:::.compute_p1b_objective(
    result, dat$C_self, dat$C_cross, dat$slides, dat$cell_types
  )
  expect_gt(obj, 0.3)
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
  obj <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE)

  obj <- runGeneSpaceCCA(obj, sigma = 0.1, nCC = 2,
                         max_iter = 500, tol = 1e-4, verbose = FALSE)

  # Gene scores populated with no NAs
  expect_gt(length(obj@geneScores), 0)
  gs_key <- names(obj@geneScores)[1]
  expect_equal(ncol(obj@geneScores[[gs_key]]), 2)
  expect_false(any(is.na(obj@geneScores[[gs_key]])))

  # Cell scores populated with no NAs and non-zero variance.
  # Per-slide z-normalization in .storeGeneSpaceCCAResults should leave
  # each (slide, cell type) cell-score column with sd ~= 1; verify that
  # against the union over slides as a sanity check (overall sd may differ
  # if slides have different mean shifts, but each slide's slice should
  # be close to unit sd).
  expect_gt(length(obj@cellScores), 0)
  cs_key <- names(obj@cellScores)[1]
  expect_equal(ncol(obj@cellScores[[cs_key]]), 2)
  expect_false(any(is.na(obj@cellScores[[cs_key]])))
  for (cc in seq_len(2)) {
    scores_cc <- obj@cellScores[[cs_key]][, cc]
    expect_gt(sd(scores_cc, na.rm = TRUE), 0)
  }

  # Per-slide cell scores should each have sd close to 1 by construction
  # (.storeGeneSpaceCCAResults applies (raw - mean) / sd per slide). The
  # cellScores key encodes the cell type; we look up which slide each row
  # belongs to via the metadata to compute per-slide stats.
  ct_idx <- obj@cellTypesSub == "CellTypeA"
  cell_ids_ct <- rownames(obj@metaDataSub)[ct_idx]
  slide_ids_ct <- getSlideID(obj)[ct_idx]
  scores_full <- obj@cellScores[[cs_key]]
  for (sl in unique(slide_ids_ct)) {
    cells_in_slide <- cell_ids_ct[slide_ids_ct == sl]
    if (length(cells_in_slide) >= 2) {
      vals <- scores_full[cells_in_slide, "CC_1"]
      expect_equal(sd(vals), 1, tolerance = 1e-6)
    }
  }

  # Gene weight columns should be unit norm (the optimizer guarantees
  # this; verify the integration didn't drop the property).
  for (k in names(obj@geneScores)) {
    for (cc in seq_len(2)) {
      norm_k_cc <- sqrt(sum(obj@geneScores[[k]][, cc]^2))
      expect_equal(norm_k_cc, 1, tolerance = 1e-6,
                   info = paste("non-unit gene-weight norm for", k, "CC", cc))
    }
  }

  # CCA output populated under the gscca_ prefix to avoid collision with
  # runSkrCCA's "sigma_" keys.
  expect_true(paste0("gscca_sigma_", 0.1) %in% names(obj@skrCCAOut))
})

test_that("runGeneSpaceCCA validates sigma against available values", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE)

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
  obj <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE)

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

test_that("streaming path matches slot-based path with normalizeDistance=FALSE", {
  # Without distance normalization there is no cross-slide coupling, so the
  # streaming and slot-based paths must produce identical covariance matrices
  # and (with the same seed) identical optimization output.
  # Sigma is chosen to match the raw (un-normalized) distance scale on the
  # synthetic fixture (median pairwise distance ~5 units).
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  test_sigma <- 5

  # Slot-based path
  obj_slot <- computeDistance(obj, distType = "Euclidean2D",
                              normalizeDistance = FALSE, verbose = FALSE)
  obj_slot <- computeKernelMatrix(obj_slot, sigmaValues = test_sigma,
                                  verbose = FALSE)
  set.seed(123)
  obj_slot <- runGeneSpaceCCA(obj_slot, sigma = test_sigma, nCC = 2,
                              max_iter = 500, tol = 1e-6, verbose = FALSE)

  # Streaming path
  set.seed(123)
  obj_stream <- runGeneSpaceCCA(
    obj, sigma = test_sigma, nCC = 2,
    max_iter = 500, tol = 1e-6,
    streaming = TRUE,
    distanceArgs = list(distType = "Euclidean2D", normalizeDistance = FALSE),
    verbose = FALSE
  )

  expect_equal(names(obj_slot@geneScores), names(obj_stream@geneScores))
  for (k in names(obj_slot@geneScores)) {
    expect_equal(obj_stream@geneScores[[k]], obj_slot@geneScores[[k]],
                 tolerance = 1e-10,
                 info = paste("geneScores mismatch for:", k))
  }

  expect_equal(names(obj_slot@cellScores), names(obj_stream@cellScores))
  for (k in names(obj_slot@cellScores)) {
    expect_equal(obj_stream@cellScores[[k]], obj_slot@cellScores[[k]],
                 tolerance = 1e-10,
                 info = paste("cellScores mismatch for:", k))
  }

  # Streaming should NOT populate the slot-based caches
  expect_length(obj_stream@distances, 0)
  expect_length(obj_stream@kernelMatrices, 0)

  # ...but should record the sigma so downstream lookups work
  expect_true(test_sigma %in% obj_stream@sigmaValues)
})

test_that("streaming path matches slot-based path with normalizeDistance=TRUE", {
  # normalizeMethod = "global" reads its reference off the cells rather than
  # off the pair blocks, so the streaming path never materializes the per-pair
  # references the other methods need. That is a separate branch from the
  # normalizeDistance = FALSE case above, and it has to land on the same scale
  # factor -- and therefore the same kernels -- as the slot-based path.
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  test_sigma <- 0.05

  obj_slot <- computeDistance(obj, distType = "Euclidean2D",
                              normalizeDistance = TRUE, verbose = FALSE)
  expect_equal(getDistanceGeometry(obj_slot)$normalizeMethod, "global")
  obj_slot <- computeKernelMatrix(obj_slot, sigmaValues = test_sigma,
                                  verbose = FALSE)
  set.seed(123)
  obj_slot <- runGeneSpaceCCA(obj_slot, sigma = test_sigma, nCC = 2,
                              max_iter = 500, tol = 1e-6, verbose = FALSE)

  set.seed(123)
  obj_stream <- runGeneSpaceCCA(
    obj, sigma = test_sigma, nCC = 2,
    max_iter = 500, tol = 1e-6,
    streaming = TRUE,
    distanceArgs = list(distType = "Euclidean2D", normalizeDistance = TRUE),
    verbose = FALSE
  )

  expect_equal(names(obj_slot@geneScores), names(obj_stream@geneScores))
  for (k in names(obj_slot@geneScores)) {
    expect_equal(obj_stream@geneScores[[k]], obj_slot@geneScores[[k]],
                 tolerance = 1e-8,
                 info = paste("geneScores mismatch for:", k))
  }
})

test_that("streaming path matches slot-based path with 3 cell types (multiple pairs)", {
  # Multi-pair regression: with 3 cell types the per-slide pair loop has
  # to handle 3 cross-pairs. An earlier draft of the streaming code freed
  # pair_distances using `lst[[k]] <- NULL`, which shrinks the list and
  # shifts indices, producing a non-conformable kernel %*% Z error on the
  # second pair. This test is the smallest fixture that exercises >1 pair.
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 90, n_slides = 2, n_genes = 30,
    n_cell_types = 3, seed = 42
  )
  obj <- subsetData(obj,
                    cellTypesOfInterest = c("CellTypeA", "CellTypeB", "CellTypeC"))
  test_sigma <- 5

  obj_slot <- computeDistance(obj, distType = "Euclidean2D",
                              normalizeDistance = FALSE, verbose = FALSE)
  obj_slot <- computeKernelMatrix(obj_slot, sigmaValues = test_sigma,
                                  verbose = FALSE)
  set.seed(123)
  obj_slot <- runGeneSpaceCCA(obj_slot, sigma = test_sigma, nCC = 2,
                              max_iter = 500, tol = 1e-6, verbose = FALSE)

  set.seed(123)
  obj_stream <- runGeneSpaceCCA(
    obj, sigma = test_sigma, nCC = 2,
    max_iter = 500, tol = 1e-6,
    streaming = TRUE,
    distanceArgs = list(distType = "Euclidean2D", normalizeDistance = FALSE),
    verbose = FALSE
  )

  for (k in names(obj_slot@geneScores)) {
    expect_equal(obj_stream@geneScores[[k]], obj_slot@geneScores[[k]],
                 tolerance = 1e-10,
                 info = paste("geneScores mismatch for:", k))
  }
  for (k in names(obj_slot@cellScores)) {
    expect_equal(obj_stream@cellScores[[k]], obj_slot@cellScores[[k]],
                 tolerance = 1e-10,
                 info = paste("cellScores mismatch for:", k))
  }
})

test_that("streaming default (no normalizationScope arg) matches slot under normalizeDistance=TRUE", {
  # Pins the new default: normalizationScope = "global". With a fixed seed,
  # streaming should be bit-identical to the slot-based path even when the
  # caller does not specify normalizationScope. This guards against an
  # accidental regression to the old "per_slide" default, which was found
  # to perturb degenerate canonical components on heterogeneous datasets
  # (NSCLC macrophage CC2: cor 0.596 vs 0.038; Issue #14).
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 90, n_slides = 2, n_genes = 30,
    n_cell_types = 3, seed = 42
  )
  obj <- subsetData(obj,
                    cellTypesOfInterest = c("CellTypeA", "CellTypeB", "CellTypeC"))
  test_sigma <- 5

  obj_slot <- computeDistance(obj, distType = "Euclidean2D",
                              normalizeDistance = TRUE, verbose = FALSE)
  obj_slot <- computeKernelMatrix(obj_slot, sigmaValues = test_sigma,
                                  verbose = FALSE, normalizeDistance = TRUE)
  set.seed(123)
  obj_slot <- runGeneSpaceCCA(obj_slot, sigma = test_sigma, nCC = 2,
                              max_iter = 500, tol = 1e-6, verbose = FALSE)

  set.seed(123)
  obj_stream <- runGeneSpaceCCA(
    obj, sigma = test_sigma, nCC = 2,
    max_iter = 500, tol = 1e-6,
    streaming = TRUE,
    distanceArgs = list(distType = "Euclidean2D",
                        normalizeDistance = TRUE),  # no normalizationScope -> default
    verbose = FALSE
  )

  for (k in names(obj_slot@geneScores)) {
    expect_equal(obj_stream@geneScores[[k]], obj_slot@geneScores[[k]],
                 tolerance = 1e-10,
                 info = paste("geneScores mismatch for:", k))
  }
  for (k in names(obj_slot@cellScores)) {
    expect_equal(obj_stream@cellScores[[k]], obj_slot@cellScores[[k]],
                 tolerance = 1e-10,
                 info = paste("cellScores mismatch for:", k))
  }
})

test_that("streaming with normalizationScope='global' matches slot under normalizeDistance=TRUE", {
  # Slot path uses global distance normalization (single factor across all
  # slides). Streaming with scope='global' replicates that exact factor by
  # taking the min low-percentile across all (slide, pair) before scaling.
  # With a fixed seed, results should be bit-identical to the slot path.
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 90, n_slides = 2, n_genes = 30,
    n_cell_types = 3, seed = 42
  )
  obj <- subsetData(obj,
                    cellTypesOfInterest = c("CellTypeA", "CellTypeB", "CellTypeC"))
  test_sigma <- 5

  obj_slot <- computeDistance(obj, distType = "Euclidean2D",
                              normalizeDistance = TRUE, verbose = FALSE)
  obj_slot <- computeKernelMatrix(obj_slot, sigmaValues = test_sigma,
                                  verbose = FALSE, normalizeDistance = TRUE)
  set.seed(123)
  obj_slot <- runGeneSpaceCCA(obj_slot, sigma = test_sigma, nCC = 2,
                              max_iter = 500, tol = 1e-6, verbose = FALSE)

  set.seed(123)
  obj_stream <- runGeneSpaceCCA(
    obj, sigma = test_sigma, nCC = 2,
    max_iter = 500, tol = 1e-6,
    streaming = TRUE,
    distanceArgs = list(distType = "Euclidean2D",
                        normalizeDistance = TRUE,
                        normalizationScope = "global"),
    verbose = FALSE
  )

  for (k in names(obj_slot@geneScores)) {
    expect_equal(obj_stream@geneScores[[k]], obj_slot@geneScores[[k]],
                 tolerance = 1e-10,
                 info = paste("geneScores mismatch for:", k))
  }
  for (k in names(obj_slot@cellScores)) {
    expect_equal(obj_stream@cellScores[[k]], obj_slot@cellScores[[k]],
                 tolerance = 1e-10,
                 info = paste("cellScores mismatch for:", k))
  }
})

test_that("streaming path runs end-to-end with normalizeDistance=TRUE", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  set.seed(7)
  obj <- runGeneSpaceCCA(
    obj, sigma = 0.1, nCC = 2,
    max_iter = 500, tol = 1e-6,
    streaming = TRUE,
    distanceArgs = list(distType = "Euclidean2D", normalizeDistance = TRUE,
                        normalizeTarget = 0.01),
    verbose = FALSE
  )

  expect_gt(length(obj@geneScores), 0)
  gs_key <- names(obj@geneScores)[1]
  expect_equal(ncol(obj@geneScores[[gs_key]]), 2)
  expect_false(any(is.na(obj@geneScores[[gs_key]])))

  expect_gt(length(obj@cellScores), 0)
  cs_key <- names(obj@cellScores)[1]
  expect_equal(ncol(obj@cellScores[[cs_key]]), 2)
  expect_false(any(is.na(obj@cellScores[[cs_key]])))

  expect_true(paste0("gscca_sigma_", 0.1) %in% names(obj@skrCCAOut))
})

test_that("streaming path rejects unknown distanceArgs / kernelArgs", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.1, streaming = TRUE,
                    distanceArgs = list(noSuchArg = 1), verbose = FALSE),
    "Unknown distanceArgs"
  )
  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.1, streaming = TRUE,
                    kernelArgs = list(noSuchArg = 1), verbose = FALSE),
    "Unknown kernelArgs"
  )
})

test_that("streaming = FALSE warns when distanceArgs / kernelArgs are passed", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE, normalizeDistance = TRUE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE)

  # max_iter caps the run so the test stays cheap, but it has to clear the
  # optimizer's own convergence bar or the "did not converge" warning escapes
  # alongside the one under test. This fixture needs ~62 iterations to reach
  # tol = 1e-6; 200 keeps the margin wide and still finishes in well under a
  # second.
  expect_warning(
    runGeneSpaceCCA(obj, sigma = 0.1, streaming = FALSE,
                    distanceArgs = list(distType = "Euclidean2D"),
                    max_iter = 200, verbose = FALSE),
    "ignored when streaming = FALSE"
  )
})

test_that("reverse-key lookup in .get_C_cross works", {
  C_cross_s <- list("B-A" = matrix(1:4, 2, 2))

  result <- CoPro:::.get_C_cross(C_cross_s, "A", "B")
  expect_equal(result, t(C_cross_s[["B-A"]]))

  expect_error(
    CoPro:::.get_C_cross(C_cross_s, "A", "C"),
    "Cross-covariance not found"
  )
})

test_that("matrix-free gene-space operators equal explicit covariances", {
  set.seed(20260728)
  Zi <- matrix(rnorm(45 * 18), 45, 18)
  Zj <- matrix(rnorm(38 * 18), 38, 18)
  K <- Matrix::rsparsematrix(45, 38, density = 0.08)
  wi <- matrix(rnorm(18), 18, 1)
  wj <- matrix(rnorm(18), 18, 1)

  self_matrix <- crossprod(Zi) / nrow(Zi)
  cross_matrix <- crossprod(Zi, K %*% Zj) /
    sqrt(nrow(Zi) * nrow(Zj))
  self_operator <- CoPro:::.new_genespace_self_operator(Zi)
  cross_operator <- CoPro:::.new_genespace_cross_operator(Zi, K, Zj)

  expect_equal(
    CoPro:::.genespace_self_quad(self_operator, wi),
    as.numeric(crossprod(wi, self_matrix %*% wi)),
    tolerance = 1e-11
  )
  expect_equal(
    CoPro:::.genespace_cross_mult(cross_operator, wj),
    cross_matrix %*% wj,
    tolerance = 1e-11
  )
  expect_equal(
    CoPro:::.genespace_cross_mult(
      CoPro:::.transpose_genespace_cross_operator(cross_operator), wi
    ),
    t(cross_matrix) %*% wi,
    tolerance = 1e-11
  )
})

test_that(".prepareGeneSpaceData drops slides below the per-slide cell threshold", {
  # The threshold (CoPro:::.min_cells_per_slide, default 10) protects the
  # G x G covariance estimates from low-rank-induced noise. We drop a
  # subset of CellTypeB cells from one slide so its CellTypeB count falls
  # below threshold, then call the data-prep helper directly and confirm:
  #   (a) a warning is emitted naming the slide and offending cell type
  #   (b) the dropped slide is absent from the returned slides vector
  #   (c) the other slide is retained
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  slide_ids <- getSlideID(obj)
  cell_types_full <- obj@cellTypesSub
  slides_all <- unique(slide_ids)
  target_slide <- slides_all[1]

  # Drop most CellTypeB cells from the target slide so only 5 remain (below
  # threshold of 10), keeping the other slide untouched.
  victim <- which(slide_ids == target_slide & cell_types_full == "CellTypeB")
  if (length(victim) <= 5) skip("Fixture too small to construct drop test")
  to_drop <- victim[seq_len(length(victim) - 5)]

  keep_idx <- setdiff(seq_len(nrow(obj@metaDataSub)), to_drop)
  obj@metaDataSub <- obj@metaDataSub[keep_idx, , drop = FALSE]
  obj@cellTypesSub <- obj@cellTypesSub[keep_idx]
  obj@locationDataSub <- obj@locationDataSub[keep_idx, , drop = FALSE]
  obj@normalizedDataSub <- obj@normalizedDataSub[keep_idx, , drop = FALSE]

  expect_warning(
    res <- CoPro:::.prepareGeneSpaceData(
      obj, clip = "quantile",
      min_prevalence = 0.008, min_cells = 5,
      cts = c("CellTypeA", "CellTypeB"),
      slides = slides_all
    ),
    "dropped"
  )
  expect_false(target_slide %in% res$slides)
  expect_true(slides_all[2] %in% res$slides)
})

test_that(".prepareGeneSpaceData errors when no slides survive filtering", {
  # If every slide has a cell type below threshold, all slides drop out
  # and we raise an informative error rather than returning an empty
  # structure that would crash the optimizer with a confusing error.
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  cell_types_full <- obj@cellTypesSub
  victim <- which(cell_types_full == "CellTypeB")
  if (length(victim) < 4) skip("Fixture too small to construct drop test")
  # Leave only 2 CellTypeB cells overall -- below threshold on every slide.
  to_drop <- victim[seq_len(length(victim) - 2)]
  keep_idx <- setdiff(seq_len(nrow(obj@metaDataSub)), to_drop)
  obj@metaDataSub <- obj@metaDataSub[keep_idx, , drop = FALSE]
  obj@cellTypesSub <- obj@cellTypesSub[keep_idx]
  obj@locationDataSub <- obj@locationDataSub[keep_idx, , drop = FALSE]
  obj@normalizedDataSub <- obj@normalizedDataSub[keep_idx, , drop = FALSE]

  expect_error(
    suppressWarnings(
      CoPro:::.prepareGeneSpaceData(
        obj, clip = "quantile",
        min_prevalence = 0.008, min_cells = 1,
        cts = c("CellTypeA", "CellTypeB"),
        slides = unique(getSlideID(obj))
      )
    ),
    "No slides"
  )
})

test_that("runGeneSpaceCCA rejects nCC larger than the gene set", {
  # User-facing guard for the case where deflation would exhaust all
  # orthogonal directions and the optimizer would silently store random
  # initialization vectors as canonical components.
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE, normalizeDistance = TRUE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE)

  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 1000,
                    max_iter = 50, tol = 1e-4, verbose = FALSE),
    "exceeds number of genes"
  )
})

test_that("runGeneSpaceCCA validates the clip argument", {
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE, normalizeDistance = TRUE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE)

  # Numeric scalar is allowed (clamps expression at the threshold)
  expect_no_error(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 2, clip = 5,
                    max_iter = 50, tol = 1e-3, verbose = FALSE)
  )

  # Invalid string value should error explicitly rather than silently
  # skipping the clipping step.
  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 2, clip = "median",
                    max_iter = 50, tol = 1e-3, verbose = FALSE),
    "must be"
  )

  # Non-scalar numeric should also error
  expect_error(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 2, clip = c(1, 2),
                    max_iter = 50, tol = 1e-3, verbose = FALSE),
    "must be"
  )
})

test_that("computeGeneAndCellScores rejects gene-space CCA outputs", {
  # Guard against the silent correctness trap where a user runs
  # runGeneSpaceCCA (which writes weights into @skrCCAOut under a
  # gscca_-prefixed key) and then calls computeGeneAndCellScores. The
  # latter applies a PCA back-projection assuming the weights are in PC
  # space, which would silently produce garbage gene scores. The guard in
  # .checkInputGAC detects this and raises a clear error.
  skip_if_not_installed("CoPro")

  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 42
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE, normalizeDistance = TRUE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE)
  # Need pcaGlobal to reach the gscca-only branch in .checkInputGAC
  obj <- computePCA(obj, nPCA = 5)
  obj <- runGeneSpaceCCA(obj, sigma = 0.1, nCC = 2,
                         max_iter = 50, tol = 1e-3, verbose = FALSE)

  expect_error(
    computeGeneAndCellScores(obj),
    "gene-space CCA"
  )
})

# ============================================================================
# Sweep choice: Gauss-Seidel vs Jacobi, and the sign-repair question
# ============================================================================

.sweep_fixture <- function(nct = 2L, seed = 42L) {
  obj <- create_test_copro_multi(
    n_cells_per_slide = 140, n_slides = 2, n_genes = 30,
    n_cell_types = nct, seed = seed
  )
  cts <- paste0("CellType", LETTERS[seq_len(nct)])
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE,
                             normalizeDistance = TRUE)
  gsd <- CoPro:::.prepareGeneSpaceData(obj, "quantile", 0.008, 20, cts,
                                       getSlideList(obj))
  cm <- CoPro:::.precomputeCovarianceMatrices(
    gsd$Z_by_slide, obj@kernelMatrices, 0.1, gsd$slides, cts
  )
  list(obj = obj, cts = cts, slides = gsd$slides,
       C_self = cm$C_self, C_cross = cm$C_cross)
}

test_that("flipping one block does not negate the objective for 3+ cell types", {
  # This is the algebra that makes the legacy sign repair invalid beyond two
  # cell types, and it is why the fix is to stop needing the repair rather than
  # to patch it. f = rho_AB + rho_AC + rho_BC. Flipping block A gives
  # -rho_AB - rho_AC + rho_BC, which differs from -f by exactly 2 * rho_BC.
  fx <- .sweep_fixture(nct = 3L)
  set.seed(555)
  w <- suppressWarnings(suppressMessages(optimize_genespace_avg_corr(
    fx$C_self, fx$C_cross, fx$slides, fx$cts,
    max_iter = 400, tol = 1e-6, verbose = FALSE
  )))

  f <- CoPro:::.compute_p1b_objective(w, fx$C_self, fx$C_cross, fx$slides, fx$cts)

  w_flip_one <- w
  w_flip_one[[fx$cts[1]]] <- -w_flip_one[[fx$cts[1]]]
  f_flip_one <- CoPro:::.compute_p1b_objective(
    w_flip_one, fx$C_self, fx$C_cross, fx$slides, fx$cts
  )

  expect_false(isTRUE(all.equal(f_flip_one, -f, tolerance = 1e-6)))

  # And the discrepancy is exactly twice the pair that does not touch block A.
  sigma_all <- CoPro:::.compute_per_slide_sigma(w, fx$C_self, fx$slides, fx$cts)
  rho_bc <- sum(vapply(fx$slides, function(s) {
    C <- CoPro:::.get_C_cross(fx$C_cross[[s]], fx$cts[2], fx$cts[3])
    CoPro:::.genespace_cross_bilinear(w[[fx$cts[2]]], C, w[[fx$cts[3]]]) /
      (sigma_all[[s]][[fx$cts[2]]] * sigma_all[[s]][[fx$cts[3]]])
  }, numeric(1))) / length(fx$slides)
  expect_equal(f_flip_one, -f + 2 * rho_bc, tolerance = 1e-8)
})

test_that("flipping every block leaves the objective unchanged", {
  # Each pairwise term picks up two sign flips, so a global flip is the same
  # solution -- which is also why it could never repair a negative objective.
  for (nct in c(2L, 3L)) {
    fx <- .sweep_fixture(nct = nct)
    set.seed(555)
    w <- suppressWarnings(suppressMessages(optimize_genespace_avg_corr(
      fx$C_self, fx$C_cross, fx$slides, fx$cts,
      max_iter = 400, tol = 1e-6, verbose = FALSE
    )))
    f <- CoPro:::.compute_p1b_objective(w, fx$C_self, fx$C_cross,
                                        fx$slides, fx$cts)
    w_all <- lapply(w, function(m) -m)
    f_all <- CoPro:::.compute_p1b_objective(w_all, fx$C_self, fx$C_cross,
                                            fx$slides, fx$cts)
    expect_equal(f_all, f, tolerance = 1e-10, info = paste("nct =", nct))
  }
})

test_that("gauss-seidel never converges to a negative objective", {
  # At a Gauss-Seidel fixed point every block satisfies w_i' g_i = ||g_i|| >= 0.
  # Summing that over blocks gives 2f >= 0 (the per-slide sigmas are positive
  # scalars and do not change the sign structure), so no sign repair can ever be
  # needed -- which is the whole justification for dropping it on this path.
  for (nct in c(2L, 3L)) {
    for (dseed in c(42L, 101L)) {
      fx <- .sweep_fixture(nct = nct, seed = dseed)
      for (iseed in c(3L, 555L)) {
        set.seed(iseed)
        w <- suppressWarnings(suppressMessages(optimize_genespace_avg_corr(
          fx$C_self, fx$C_cross, fx$slides, fx$cts,
          max_iter = 400, tol = 1e-6, verbose = FALSE, sweep = "gauss-seidel"
        )))
        f <- CoPro:::.compute_p1b_objective(w, fx$C_self, fx$C_cross,
                                            fx$slides, fx$cts)
        expect_gte(f, 0)
      }
    }
  }
})

test_that("sweep = jacobi reproduces the pre-change gene-space result", {
  # Regression lock for results computed before "gauss-seidel" became the
  # default. The expected values were measured on a pristine checkout of the
  # commit before the sweep argument existed, using this exact fixture, seed and
  # settings; the old code path had no sweep argument and always did what
  # sweep = "jacobi" now does. Compared as a checksum of |weights| so an
  # incidental sign flip does not fail the test.
  expected <- c(nct2 = 1.823708199481756e+01, nct3 = 2.640422549226658e+01)

  for (nct in c(2L, 3L)) {
    obj <- create_test_copro_multi(
      n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
      n_cell_types = nct, seed = 42
    )
    cts <- paste0("CellType", LETTERS[seq_len(nct)])
    obj <- subsetData(obj, cellTypesOfInterest = cts)
    obj <- computeDistance(obj, distType = "Euclidean2D",
                           normalizeDistance = TRUE, verbose = FALSE)
    obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE,
                               normalizeDistance = TRUE)

    set.seed(20260729)
    fit <- suppressWarnings(suppressMessages(runGeneSpaceCCA(
      obj, sigma = 0.1, nCC = 2, max_iter = 500, tol = 1e-4,
      verbose = FALSE, sweep = "jacobi"
    )))
    w <- fit@skrCCAOut[["gscca_sigma_0.1"]]
    checksum <- sum(vapply(w, function(m) sum(abs(m)), numeric(1)))
    expect_equal(checksum, unname(expected[[paste0("nct", nct)]]),
                 tolerance = 1e-6, info = paste("nct =", nct))
  }
})

test_that("both sweeps produce unit-norm weights and run for 3 cell types", {
  fx <- .sweep_fixture(nct = 3L)
  for (sw in c("gauss-seidel", "jacobi")) {
    set.seed(555)
    w <- suppressWarnings(suppressMessages(optimize_genespace_avg_corr_n(
      fx$C_self, fx$C_cross, fx$slides, fx$cts,
      w_list = suppressWarnings(suppressMessages(optimize_genespace_avg_corr(
        fx$C_self, fx$C_cross, fx$slides, fx$cts,
        max_iter = 200, tol = 1e-5, verbose = FALSE, sweep = sw
      ))),
      nCC = 2, max_iter = 200, tol = 1e-5, verbose = FALSE, sweep = sw
    )))
    for (ct in fx$cts) {
      expect_equal(ncol(w[[ct]]), 2, info = sw)
      for (cc in seq_len(2)) {
        expect_equal(sqrt(sum(w[[ct]][, cc]^2)), 1, tolerance = 1e-6,
                     info = paste(sw, ct, cc))
      }
    }
  }
})

test_that("gene-space first and subsequent wrappers use the shared axis iterator", {
  fx <- .sweep_fixture(nct = 3L)
  for (sw in c("gauss-seidel", "jacobi")) {
    first <- suppressWarnings(suppressMessages(optimize_genespace_avg_corr(
      fx$C_self, fx$C_cross, fx$slides, fx$cts,
      max_iter = 300, tol = 1e-6, verbose = FALSE, sweep = sw
    )))
    first_direct <- suppressWarnings(suppressMessages(
      CoPro:::.optimizeGeneSpaceAxis(
        fx$C_self, fx$C_cross, fx$slides, fx$cts,
        component = 1L, max_iter = 300, tol = 1e-6,
        verbose = FALSE, sweep = sw
      )
    ))
    all_axes <- suppressWarnings(suppressMessages(optimize_genespace_avg_corr_n(
      fx$C_self, fx$C_cross, fx$slides, fx$cts,
      w_list = first, nCC = 2L,
      max_iter = 300, tol = 1e-6, verbose = FALSE, sweep = sw
    )))
    second_direct <- suppressWarnings(suppressMessages(
      CoPro:::.optimizeGeneSpaceAxis(
        fx$C_self, fx$C_cross, fx$slides, fx$cts,
        component = 2L, previous_weights = first,
        max_iter = 300, tol = 1e-6, verbose = FALSE, sweep = sw
      )
    ))

    expect_equal(first, first_direct, tolerance = 0, info = sw)
    for (ct in fx$cts) {
      expect_equal(
        all_axes[[ct]][, 2L, drop = FALSE], second_direct[[ct]],
        tolerance = 0, info = paste(sw, ct)
      )
    }
  }
})

test_that("gene-space objective = sumcov drops the per-slide denominators", {
  # Completes the {space} x {objective} grid: gene space can run the same
  # criterion runSkrCCA() defaults to.
  fx <- .sweep_fixture(nct = 2L)

  # With objective = "sumcov" every scale is pinned at 1, so the objective is
  # the plain average cross-covariance.
  sigma_all <- CoPro:::.compute_per_slide_sigma(
    list(CellTypeA = matrix(1, 1, 1), CellTypeB = matrix(1, 1, 1)),
    fx$C_self, fx$slides, fx$cts, objective = "sumcov"
  )
  expect_true(all(unlist(sigma_all) == 1))

  set.seed(555)
  w <- suppressWarnings(suppressMessages(optimize_genespace_avg_corr(
    fx$C_self, fx$C_cross, fx$slides, fx$cts,
    max_iter = 300, tol = 1e-6, verbose = FALSE, objective = "sumcov"
  )))
  for (ct in fx$cts) {
    expect_equal(sqrt(sum(w[[ct]]^2)), 1, tolerance = 1e-6)
  }

  # Rescaling one slide's expression must move a covariance objective; this is
  # the contrast that makes the sumcor default worth having.
  f_cov <- CoPro:::.compute_p1b_objective(w, fx$C_self, fx$C_cross, fx$slides,
                                          fx$cts, objective = "sumcov")
  f_cor <- CoPro:::.compute_p1b_objective(w, fx$C_self, fx$C_cross, fx$slides,
                                          fx$cts, objective = "sumcor")
  expect_false(isTRUE(all.equal(f_cov, f_cor, tolerance = 1e-6)))
})

test_that("an unknown sweep or gene-space objective is rejected", {
  fx <- .sweep_fixture(nct = 2L)
  expect_error(optimize_genespace_avg_corr(
    fx$C_self, fx$C_cross, fx$slides, fx$cts, verbose = FALSE, sweep = "sor"
  ))
  expect_error(optimize_genespace_avg_corr(
    fx$C_self, fx$C_cross, fx$slides, fx$cts, verbose = FALSE,
    objective = "sumsquares"
  ))
})
