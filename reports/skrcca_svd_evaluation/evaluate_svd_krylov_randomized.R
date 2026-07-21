#!/usr/bin/env Rscript

# Targeted numerical evaluation for the reviewer question about replacing
# CoPro/skrCCA's optimizer with an SVD, a Krylov method, or randomized SVD.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(Matrix)
})

set.seed(20260720)

out_dir <- "reports/skrcca_svd_evaluation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_pca <- 40L
n_cc <- 4L
sigma <- 0.01
tol <- 1e-8
max_iter <- 2000L

time_reps <- function(fun, reps = 20L) {
  calibration <- capture.output({
    calibration_time <- system.time(calibration_value <- fun())
  })
  invisible(calibration)
  calibration_elapsed <- unname(calibration_time[["elapsed"]])
  inner_reps <- if (calibration_elapsed < 0.002) {
    20L
  } else if (calibration_elapsed < 0.01) {
    20L
  } else {
    1L
  }
  elapsed <- numeric(reps)
  value <- calibration_value
  for (i in seq_len(reps)) {
    gc(FALSE)
    captured <- capture.output({
      timing <- system.time({
        for (j in seq_len(inner_reps)) value <- fun()
      })
    })
    invisible(captured)
    elapsed[i] <- unname(timing[["elapsed"]]) / inner_reps
  }
  list(
    value = value,
    median_sec = median(elapsed),
    mean_sec = mean(elapsed),
    min_sec = min(elapsed),
    times = elapsed,
    inner_reps = inner_reps
  )
}

copy_Y <- function(Y) {
  lapply(Y, function(row) lapply(row, function(x) {
    if (is.null(x)) NULL else x + 0
  }))
}

aggregate_Y <- function(Y_all, cell_types) {
  out <- setNames(vector("list", length(cell_types)), cell_types)
  for (i in cell_types) {
    out[[i]] <- setNames(vector("list", length(cell_types)), cell_types)
    for (j in cell_types) {
      if (i == j && length(cell_types) > 1L) next
      entries <- lapply(Y_all, function(Y) Y[[i]][[j]])
      if (!all(vapply(entries, is.null, logical(1)))) {
        out[[i]][[j]] <- Reduce("+", entries)
      }
    }
  }
  out
}

objective_one <- function(Y, w, cell_types = names(w)) {
  if (length(cell_types) == 1L) {
    ct <- cell_types[[1L]]
    return(drop(crossprod(w[[ct]], Y[[ct]][[ct]] %*% w[[ct]])))
  }
  pairs <- combn(cell_types, 2L)
  sum(vapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1L, k]
    j <- pairs[2L, k]
    drop(crossprod(w[[i]], Y[[i]][[j]] %*% w[[j]]))
  }, numeric(1)))
}

objective_axes <- function(Y, W, cell_types = names(W)) {
  vapply(seq_len(ncol(W[[1L]])), function(k) {
    wk <- lapply(W, function(x) x[, k, drop = FALSE])
    objective_one(Y, wk, cell_types)
  }, numeric(1))
}

joint_objective <- function(Y, W, cell_types = names(W)) {
  sum(objective_axes(Y, W, cell_types))
}

max_orthogonality_error <- function(W) {
  max(vapply(W, function(x) {
    max(abs(crossprod(x) - diag(ncol(x))))
  }, numeric(1)))
}

kkt_residual <- function(Y, w, cell_types = names(w)) {
  max(vapply(cell_types, function(i) {
    if (length(cell_types) == 1L) {
      gradient <- Y[[i]][[i]] %*% w[[i]]
    } else {
      gradient <- Reduce("+", lapply(setdiff(cell_types, i), function(j) {
        Y[[i]][[j]] %*% w[[j]]
      }))
    }
    lambda_i <- drop(crossprod(w[[i]], gradient))
    as.numeric(norm(gradient - lambda_i * w[[i]], "F") /
                 max(norm(gradient, "F"), .Machine$double.eps))
  }, numeric(1)))
}

align_columns <- function(estimate, truth) {
  out <- estimate
  for (k in seq_len(ncol(truth))) {
    if (sum(out[, k] * truth[, k]) < 0) out[, k] <- -out[, k]
  }
  out
}

two_type_svd <- function(Y, cell_types, n_cc, sdev2_list = NULL,
                         engine = c("base", "lapack", "irlba", "randomized"),
                         oversample = 10L, n_power = 2L) {
  engine <- match.arg(engine)
  stopifnot(length(cell_types) == 2L)
  i <- cell_types[[1L]]
  j <- cell_types[[2L]]
  A <- as.matrix(Y[[i]][[j]])

  if (!is.null(sdev2_list)) {
    inv_i <- 1 / sqrt(sdev2_list[[i]])
    inv_j <- 1 / sqrt(sdev2_list[[j]])
    A <- sweep(sweep(A, 1L, inv_i, "*"), 2L, inv_j, "*")
  }

  if (engine == "base") {
    dec <- svd(A, nu = n_cc, nv = n_cc)
    u <- dec$u[, seq_len(n_cc), drop = FALSE]
    v <- dec$v[, seq_len(n_cc), drop = FALSE]
    d <- dec$d[seq_len(n_cc)]
  } else if (engine == "lapack") {
    dec <- La.svd(A, nu = n_cc, nv = n_cc)
    u <- dec$u[, seq_len(n_cc), drop = FALSE]
    v <- t(dec$vt)[, seq_len(n_cc), drop = FALSE]
    d <- dec$d[seq_len(n_cc)]
  } else if (engine == "irlba") {
    dec <- irlba::irlba(A, nu = n_cc, nv = n_cc)
    u <- dec$u
    v <- dec$v
    d <- dec$d
  } else {
    target <- min(ncol(A), n_cc + oversample)
    omega <- matrix(rnorm(ncol(A) * target), nrow = ncol(A))
    Q <- qr.Q(qr(A %*% omega), complete = FALSE)
    if (n_power > 0L) {
      for (q in seq_len(n_power)) {
        Q_right <- qr.Q(qr(crossprod(A, Q)), complete = FALSE)
        Q <- qr.Q(qr(A %*% Q_right), complete = FALSE)
      }
    }
    small <- svd(crossprod(Q, A), nu = n_cc, nv = n_cc)
    u <- Q %*% small$u[, seq_len(n_cc), drop = FALSE]
    v <- small$v[, seq_len(n_cc), drop = FALSE]
    d <- small$d[seq_len(n_cc)]
  }

  if (!is.null(sdev2_list)) {
    u <- sweep(u, 1L, 1 / sqrt(sdev2_list[[i]]), "*")
    v <- sweep(v, 1L, 1 / sqrt(sdev2_list[[j]]), "*")
  }

  list(weights = setNames(list(u, v), cell_types), values = d)
}

current_axes_from_Y <- function(Y, X_list, cell_types, n_cc,
                                deflation = c("rank1", "projection"),
                                sdev2_list = NULL) {
  deflation <- match.arg(deflation)
  w <- CoPro:::initialize_weights_svd(X_list, cell_types)
  w <- CoPro:::bilinear_w_from_Y_resi(
    w_list_new = w,
    Y_resi = Y,
    n_features = nrow(Y[[cell_types[[1L]]]][[cell_types[[2L]]]]),
    max_iter = max_iter,
    tol = tol,
    step_size = 1,
    sdev2_list = sdev2_list
  )
  if (n_cc == 1L) return(w)

  residual <- copy_Y(Y)
  for (k in seq_len(n_cc - 1L)) {
    residual <- CoPro:::apply_deflation(
      residual, w, k, cell_types, sdev2_list = sdev2_list,
      deflation = deflation
    )
    init <- CoPro:::initialize_next_component(residual, cell_types)
    next_w <- CoPro:::bilinear_w_from_Y_resi(
      w_list_new = init,
      Y_resi = residual,
      n_features = nrow(Y[[cell_types[[1L]]]][[cell_types[[2L]]]]),
      max_iter = max_iter,
      tol = tol,
      step_size = 1,
      sdev2_list = sdev2_list
    )
    for (ct in cell_types) w[[ct]] <- cbind(w[[ct]], next_w[[ct]])
  }
  w
}

# Reproduce the pre-SVD two-type production route: form Y for CC1, iterate,
# form Y again for later axes, then use sequential rank-one deflation.
legacy_two_type_end_to_end <- function(X_list, flat_kernels, sigma, cell_types,
                                       n_cc) {
  Y_first <- CoPro:::compute_Y_resi(
    X_list, flat_kernels, sigma, cell_types
  )
  w <- CoPro:::initialize_weights_svd(X_list, cell_types)
  w <- CoPro:::bilinear_w_from_Y_resi(
    w, Y_first, ncol(X_list[[1L]]), max_iter = max_iter, tol = tol
  )
  if (n_cc == 1L) return(w)

  residual <- CoPro:::compute_Y_resi(
    X_list, flat_kernels, sigma, cell_types
  )
  for (k in seq_len(n_cc - 1L)) {
    residual <- CoPro:::apply_deflation(
      residual, w, k, cell_types, deflation = "rank1"
    )
    init <- CoPro:::initialize_next_component(residual, cell_types)
    next_w <- CoPro:::bilinear_w_from_Y_resi(
      init, residual, ncol(X_list[[1L]]), max_iter = max_iter, tol = tol
    )
    for (ct in cell_types) w[[ct]] <- cbind(w[[ct]], next_w[[ct]])
  }
  w
}

build_block_operator <- function(Y, cell_types) {
  p <- nrow(Y[[cell_types[[1L]]]][[cell_types[[2L]]]])
  m <- length(cell_types)
  B <- matrix(0, nrow = m * p, ncol = m * p)
  blocks <- split(seq_len(m * p), rep(seq_len(m), each = p))
  for (a in seq_len(m - 1L)) {
    for (b in (a + 1L):m) {
      i <- cell_types[[a]]
      j <- cell_types[[b]]
      B[blocks[[a]], blocks[[b]]] <- as.matrix(Y[[i]][[j]])
      B[blocks[[b]], blocks[[a]]] <- as.matrix(Y[[j]][[i]])
    }
  }
  list(matrix = B, blocks = blocks)
}

block_eigen_weights <- function(Y, cell_types, n_axes = 1L,
                                normalize_blocks = TRUE) {
  block <- build_block_operator(Y, cell_types)
  dec <- eigen(block$matrix, symmetric = TRUE)
  V <- dec$vectors[, seq_len(n_axes), drop = FALSE]
  W <- setNames(lapply(seq_along(cell_types), function(i) {
    x <- V[block$blocks[[i]], , drop = FALSE]
    if (normalize_blocks) {
      for (k in seq_len(ncol(x))) x[, k] <- x[, k] / sqrt(sum(x[, k]^2))
    }
    x
  }), cell_types)
  list(weights = W, values = dec$values[seq_len(n_axes)], raw_vectors = V,
       blocks = block$blocks)
}

polar_factor <- function(G) {
  dec <- svd(G, nu = ncol(G), nv = ncol(G))
  dec$u %*% t(dec$v)
}

joint_stiefel <- function(Y, W, cell_types = names(W), max_iter = 2000L,
                          tol = 1e-10) {
  old_obj <- joint_objective(Y, W, cell_types)
  for (iter in seq_len(max_iter)) {
    for (i in cell_types) {
      gradient <- Reduce("+", lapply(setdiff(cell_types, i), function(j) {
        Y[[i]][[j]] %*% W[[j]]
      }))
      W[[i]] <- polar_factor(gradient)
    }
    new_obj <- joint_objective(Y, W, cell_types)
    if (abs(new_obj - old_obj) <= tol * max(1, abs(old_obj))) break
    old_obj <- new_obj
  }

  # A common rotation leaves the trace objective unchanged. Rotate all blocks
  # so the symmetrized cross-axis score matrix is diagonal and axes are ordered.
  M <- matrix(0, ncol(W[[1L]]), ncol(W[[1L]]))
  pairs <- combn(cell_types, 2L)
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1L, k]
    j <- pairs[2L, k]
    M <- M + crossprod(W[[i]], Y[[i]][[j]] %*% W[[j]])
  }
  rotation <- eigen((M + t(M)) / 2, symmetric = TRUE)$vectors
  W <- lapply(W, function(x) x %*% rotation)
  list(weights = W, iterations = iter,
       objective = joint_objective(Y, W, cell_types))
}

cat("Preparing the real colon D3 PC-space operators...\n")
dat <- readRDS("data-raw/vignette_data/copro_colon_d3.rds")
obj <- newCoProSingle(
  normalizedData = dat$normalizedData,
  locationData = dat$locationData,
  metaData = dat$metaData,
  cellTypes = dat$cellTypes
)
cell_types_3 <- c("Epithelial", "Fibroblast", "Immune")
obj <- subsetData(obj, cellTypesOfInterest = cell_types_3)
obj <- computePCA(obj, nPCA = n_pca, center = TRUE, scale. = TRUE)
obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
obj <- computeKernelMatrix(obj, sigmaValues = sigma, verbose = FALSE)
prepared <- CoPro:::.prepareDataMatrices(
  object = obj, is_multi = FALSE, scalePCs = TRUE, cts = cell_types_3
)
X3 <- prepared$PCmats

Y_time <- time_reps(function() {
  CoPro:::compute_Y_resi(X3, obj@kernelMatrices, sigma, cell_types_3)
}, reps = 10L)
Y3 <- Y_time$value

# -------------------------------------------------------------------------
# 1. Exact two-type reduction and all-axis SVD
# -------------------------------------------------------------------------

cell_types_2 <- cell_types_3[1:2]
X2 <- X3[cell_types_2]
Y2 <- lapply(Y3[cell_types_2], function(row) row[cell_types_2])

direct2 <- two_type_svd(Y2, cell_types_2, n_cc, engine = "base")
current2 <- current_axes_from_Y(Y2, X2, cell_types_2, n_cc, "rank1")

two_equiv <- do.call(rbind, lapply(seq_len(n_cc), function(k) {
  data.frame(
    axis = k,
    cell_type = cell_types_2,
    abs_cosine = vapply(cell_types_2, function(ct) {
      abs(sum(current2[[ct]][, k] * direct2$weights[[ct]][, k]))
    }, numeric(1)),
    current_objective = objective_axes(Y2, current2, cell_types_2)[k],
    svd_singular_value = direct2$values[k]
  )
}))
write.csv(two_equiv, file.path(out_dir, "real_two_type_equivalence.csv"),
          row.names = FALSE)

bench_two <- list(
  prechange_coordinate_end_to_end = time_reps(function() {
    legacy_two_type_end_to_end(
      X2, obj@kernelMatrices, sigma, cell_types_2, n_cc
    )
  }, reps = 10L),
  exported_two_step_exact_svd = time_reps(function() {
    first <- CoPro:::optimize_bilinear(
      X2, obj@kernelMatrices, sigma, max_iter = max_iter, tol = tol
    )
    CoPro:::optimize_bilinear_n(
      X2, obj@kernelMatrices, sigma, first, cell_types_2,
      nCC = n_cc, max_iter = max_iter, tol = tol
    )
  }, reps = 10L),
  form_Y_once_plus_base_svd = time_reps(function() {
    Y <- CoPro:::compute_Y_resi(X2, obj@kernelMatrices, sigma, cell_types_2)
    two_type_svd(Y, cell_types_2, n_cc, engine = "base")
  }, reps = 10L),
  form_Y_once_plus_irlba = time_reps(function() {
    Y <- CoPro:::compute_Y_resi(X2, obj@kernelMatrices, sigma, cell_types_2)
    two_type_svd(Y, cell_types_2, n_cc, engine = "irlba")
  }, reps = 10L),
  solve_only_base_svd = time_reps(function() {
    two_type_svd(Y2, cell_types_2, n_cc, engine = "base")
  }, reps = 30L),
  solve_only_lapack_svd = time_reps(function() {
    two_type_svd(Y2, cell_types_2, n_cc, engine = "lapack")
  }, reps = 30L),
  solve_only_irlba = time_reps(function() {
    two_type_svd(Y2, cell_types_2, n_cc, engine = "irlba")
  }, reps = 30L),
  solve_only_randomized = time_reps(function() {
    two_type_svd(Y2, cell_types_2, n_cc, engine = "randomized")
  }, reps = 30L)
)
two_timing <- do.call(rbind, lapply(names(bench_two), function(method) {
  x <- bench_two[[method]]
  data.frame(method = method, median_sec = x$median_sec,
             mean_sec = x$mean_sec, min_sec = x$min_sec)
}))
two_timing$speedup_vs_prechange <-
  two_timing$median_sec[
    two_timing$method == "prechange_coordinate_end_to_end"
  ] /
  two_timing$median_sec
write.csv(two_timing, file.path(out_dir, "real_two_type_timings.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 2. Synthetic crossover for exact, Krylov, and randomized SVD
# -------------------------------------------------------------------------

synthetic_results <- list()
result_index <- 1L
for (p in c(40L, 100L, 250L, 500L)) {
  U <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  V <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  singular_values <- seq_len(p)^(-0.6)
  A <- U %*% (singular_values * t(V))
  Y <- list(A = list(A = NULL, B = A), B = list(A = t(A), B = NULL))
  truth <- list(values = singular_values[seq_len(n_cc)], u = U[, seq_len(n_cc)],
                v = V[, seq_len(n_cc)])
  reps <- if (p <= 100L) 20L else if (p <= 250L) 10L else 5L

  for (engine in c("base", "lapack", "irlba", "randomized")) {
    timed <- time_reps(function() {
      set.seed(1000L + p)
      two_type_svd(Y, c("A", "B"), n_cc, engine = engine)
    }, reps = reps)
    fit <- timed$value
    Uhat <- fit$weights$A
    Vhat <- fit$weights$B
    singular_rel_error <- max(abs(fit$values - truth$values) / truth$values)
    u_subspace_error <- norm(Uhat %*% t(Uhat) -
                               truth$u %*% t(truth$u), "F") / sqrt(2 * n_cc)
    v_subspace_error <- norm(Vhat %*% t(Vhat) -
                               truth$v %*% t(truth$v), "F") / sqrt(2 * n_cc)
    synthetic_results[[result_index]] <- data.frame(
      p = p, n_axes = n_cc, spectrum = "power_-0.6", method = engine,
      median_sec = timed$median_sec, mean_sec = timed$mean_sec,
      singular_value_max_rel_error = singular_rel_error,
      max_subspace_error = max(u_subspace_error, v_subspace_error)
    )
    result_index <- result_index + 1L
  }
}
synthetic_svd <- do.call(rbind, synthetic_results)
write.csv(synthetic_svd, file.path(out_dir, "synthetic_svd_benchmark.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 3. Multi-type block eigensolver relaxation
# -------------------------------------------------------------------------

current3_first_time <- time_reps(function() {
  init <- CoPro:::initialize_weights_svd(X3, cell_types_3)
  CoPro:::bilinear_w_from_Y_resi(
    init, Y3, n_pca, max_iter = max_iter, tol = tol
  )
}, reps = 10L)
current3_first <- current3_first_time$value

block3_time <- time_reps(function() {
  block_eigen_weights(Y3, cell_types_3, n_axes = 1L, normalize_blocks = TRUE)
}, reps = 10L)
block3 <- block3_time$value

refined_block_time <- time_reps(function() {
  CoPro:::bilinear_w_from_Y_resi(
    block3$weights, Y3, n_pca, max_iter = max_iter, tol = tol
  )
}, reps = 10L)
refined_block <- refined_block_time$value

raw_block_norms <- vapply(seq_along(cell_types_3), function(i) {
  norm(block3$raw_vectors[block3$blocks[[i]], , drop = FALSE], "F")
}, numeric(1))

multi_first <- data.frame(
  method = c("CoPro_coordinate", "block_eigen_then_block_normalize",
             "block_eigen_then_CoPro_refine"),
  objective = c(objective_one(Y3, current3_first),
                objective_one(Y3, block3$weights),
                objective_one(Y3, refined_block)),
  kkt_residual = c(kkt_residual(Y3, current3_first),
                   kkt_residual(Y3, block3$weights),
                   kkt_residual(Y3, refined_block)),
  solver_median_sec = c(current3_first_time$median_sec,
                        block3_time$median_sec,
                        block3_time$median_sec + refined_block_time$median_sec),
  raw_block_norm_min = c(NA, min(raw_block_norms), min(raw_block_norms)),
  raw_block_norm_max = c(NA, max(raw_block_norms), max(raw_block_norms)),
  raw_block_norm_cv = c(NA, sd(raw_block_norms) / mean(raw_block_norms),
                        sd(raw_block_norms) / mean(raw_block_norms))
)
write.csv(multi_first, file.path(out_dir, "real_three_type_block_relaxation.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 4. Higher axes for three types
# -------------------------------------------------------------------------

rank1_time <- time_reps(function() {
  current_axes_from_Y(Y3, X3, cell_types_3, n_cc, "rank1")
}, reps = 5L)
rank1_axes <- rank1_time$value

projection_time <- time_reps(function() {
  current_axes_from_Y(Y3, X3, cell_types_3, n_cc, "projection")
}, reps = 5L)
projection_axes <- projection_time$value

block_init <- block_eigen_weights(Y3, cell_types_3, n_axes = n_cc,
                                  normalize_blocks = FALSE)$weights
block_init <- lapply(block_init, polar_factor)
joint_time <- time_reps(function() {
  joint_stiefel(Y3, block_init, cell_types_3)
}, reps = 5L)
joint_axes <- joint_time$value$weights

higher_summary <- data.frame(
  method = c("sequential_rank1_current", "sequential_projection",
             "joint_product_Stiefel"),
  total_original_objective = c(joint_objective(Y3, rank1_axes),
                               joint_objective(Y3, projection_axes),
                               joint_objective(Y3, joint_axes)),
  cc1_original_objective = c(objective_axes(Y3, rank1_axes)[1L],
                             objective_axes(Y3, projection_axes)[1L],
                             objective_axes(Y3, joint_axes)[1L]),
  max_orthogonality_error = c(max_orthogonality_error(rank1_axes),
                              max_orthogonality_error(projection_axes),
                              max_orthogonality_error(joint_axes)),
  median_sec = c(rank1_time$median_sec, projection_time$median_sec,
                 joint_time$median_sec)
)
write.csv(higher_summary, file.path(out_dir, "real_three_type_higher_axes.csv"),
          row.names = FALSE)

higher_by_axis <- do.call(rbind, lapply(list(
  sequential_rank1_current = rank1_axes,
  sequential_projection = projection_axes,
  joint_product_Stiefel = joint_axes
), function(W) {
  data.frame(axis = seq_len(n_cc), original_objective = objective_axes(Y3, W))
}))
higher_by_axis$method <- rep(c("sequential_rank1_current",
                               "sequential_projection",
                               "joint_product_Stiefel"), each = n_cc)
write.csv(higher_by_axis, file.path(out_dir, "real_three_type_objective_by_axis.csv"),
          row.names = FALSE)

higher_similarity <- do.call(rbind, lapply(seq_len(n_cc), function(k) {
  data.frame(
    axis = k,
    cell_type = cell_types_3,
    rank1_vs_projection_abs_cosine = vapply(cell_types_3, function(ct) {
      abs(sum(rank1_axes[[ct]][, k] * projection_axes[[ct]][, k]))
    }, numeric(1)),
    projection_vs_joint_max_subspace_cosine = vapply(cell_types_3, function(ct) {
      max(abs(crossprod(
        projection_axes[[ct]][, k, drop = FALSE], joint_axes[[ct]]
      )))
    }, numeric(1))
  )
}))
write.csv(higher_similarity,
          file.path(out_dir, "real_three_type_higher_axis_similarity.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# 5. Exact equivalence of sample stacking and operator aggregation
# -------------------------------------------------------------------------

set.seed(20260721)
sample_ids <- paste0("sample", 1:4)
pair_samples <- lapply(seq_along(sample_ids), function(s) {
  n_a <- 25L + 3L * s
  n_b <- 20L + 2L * s
  X_a <- matrix(rnorm(n_a * n_pca), n_a, n_pca)
  X_b <- matrix(rnorm(n_b * n_pca), n_b, n_pca)
  K <- matrix(rnorm(n_a * n_b), n_a, n_b)
  list(X_a = X_a, X_b = X_b, K = K,
       Y = crossprod(X_a, K %*% X_b))
})
Y_sum <- Reduce("+", lapply(pair_samples, `[[`, "Y"))
X_a_stack <- do.call(rbind, lapply(pair_samples, `[[`, "X_a"))
X_b_stack <- do.call(rbind, lapply(pair_samples, `[[`, "X_b"))
K_block <- as(bdiag(lapply(pair_samples, `[[`, "K")), "dgCMatrix")
Y_block <- crossprod(X_a_stack, K_block %*% X_b_stack)

# Deflation commutes with aggregation when weights are shared across samples.
sample_Y <- lapply(pair_samples, function(x) {
  list(A = list(A = NULL, B = x$Y),
       B = list(A = t(x$Y), B = NULL))
})
sample_aggregate <- aggregate_Y(sample_Y, c("A", "B"))
sample_w <- two_type_svd(sample_aggregate, c("A", "B"), 1L)$weights
deflated_then_sum <- aggregate_Y(lapply(sample_Y, function(Y) {
  CoPro:::apply_deflation(Y, sample_w, 1L, c("A", "B"),
                          deflation = "projection")
}), c("A", "B"))
sum_then_deflated <- CoPro:::apply_deflation(
  sample_aggregate, sample_w, 1L, c("A", "B"), deflation = "projection"
)

sample_equivalence <- data.frame(
  comparison = c("block_diagonal_kernel_vs_sum_Y",
                 "deflate_each_then_sum_vs_sum_then_deflate"),
  max_absolute_difference = c(max(abs(Y_sum - Y_block)),
                              max(abs(deflated_then_sum$A$B -
                                        sum_then_deflated$A$B)))
)
write.csv(sample_equivalence,
          file.path(out_dir, "sample_aggregation_equivalence.csv"),
          row.names = FALSE)

# Cost of per-sample rather than aggregate deflation in PC space.
synthetic_sample_cost <- list()
for (n_samples in c(1L, 3L, 10L, 30L)) {
  Y_all <- rep(sample_Y, length.out = n_samples)
  Y_agg <- aggregate_Y(Y_all, c("A", "B"))
  w <- two_type_svd(Y_agg, c("A", "B"), 1L)$weights
  per_sample <- time_reps(function() {
    aggregate_Y(lapply(Y_all, function(Y) {
      CoPro:::apply_deflation(Y, w, 1L, c("A", "B"))
    }), c("A", "B"))
  }, reps = 20L)
  aggregate_first <- time_reps(function() {
    CoPro:::apply_deflation(Y_agg, w, 1L, c("A", "B"))
  }, reps = 20L)
  synthetic_sample_cost[[length(synthetic_sample_cost) + 1L]] <- data.frame(
    n_samples = n_samples,
    per_sample_then_sum_sec = per_sample$median_sec,
    aggregate_then_deflate_sec = aggregate_first$median_sec,
    speedup = per_sample$median_sec / aggregate_first$median_sec
  )
}
sample_cost <- do.call(rbind, synthetic_sample_cost)
write.csv(sample_cost, file.path(out_dir, "sample_aggregation_deflation_timing.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# Machine-readable session information and a compact text summary
# -------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

summary_lines <- c(
  "skrCCA SVD/Krylov/randomized-SVD evaluation",
  sprintf("Real data: colon D3; %d PCs; sigma = %g; %d axes", n_pca, sigma, n_cc),
  "",
  "Two-type exactness (current coordinate method versus exact SVD):",
  capture.output(print(two_equiv)),
  "",
  "Two-type timing:",
  capture.output(print(two_timing)),
  "",
  "Synthetic SVD timing and accuracy:",
  capture.output(print(synthetic_svd)),
  "",
  "Three-type block-eigen relaxation:",
  capture.output(print(multi_first)),
  "",
  "Three-type higher axes:",
  capture.output(print(higher_summary)),
  "",
  "Three-type higher-axis similarity:",
  capture.output(print(higher_similarity)),
  "",
  "Sample aggregation exactness:",
  capture.output(print(sample_equivalence)),
  "",
  "Sample aggregation deflation timing:",
  capture.output(print(sample_cost)),
  "",
  sprintf("Three-type Y formation median: %.6f sec", Y_time$median_sec)
)
writeLines(summary_lines, file.path(out_dir, "summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
