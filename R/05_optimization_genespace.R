#' Gene-Space Average Per-Slide Canonical Correlation Optimization
#'
#' These functions implement the Proposal 1b algorithm: batch-robust CCA
#' in gene space using average per-slide canonical correlation.
#' Each slide's contribution is normalized by its own score variance,
#' preventing batch axes from inflating the objective.
#'
#' @name genespace_optimization
#' @keywords internal
NULL

# ============================================================================
# Internal Helpers
# ============================================================================

.new_genespace_self_operator <- function(Z) {
  structure(
    list(Z = Z, scale = 1 / nrow(Z)),
    class = "CoProGeneSelfOperator"
  )
}

.new_genespace_cross_operator <- function(Z_i, K, Z_j) {
  structure(
    list(
      Z_i = Z_i, K = K, Z_j = Z_j,
      scale = 1 / sqrt(nrow(Z_i) * nrow(Z_j)),
      transposed = FALSE
    ),
    class = "CoProGeneCrossOperator"
  )
}

.transpose_genespace_cross_operator <- function(x) {
  x$transposed <- !isTRUE(x$transposed)
  x
}

.genespace_n_genes <- function(x) {
  if (inherits(x, "CoProGeneSelfOperator")) return(ncol(x$Z))
  nrow(x)
}

.genespace_self_quad <- function(x, w) {
  if (inherits(x, "CoProGeneSelfOperator")) {
    score <- x$Z %*% w
    return(as.numeric(crossprod(score)) * x$scale)
  }
  as.numeric(crossprod(w, x %*% w))
}

.genespace_cross_mult <- function(x, w) {
  if (!inherits(x, "CoProGeneCrossOperator")) return(x %*% w)
  if (!isTRUE(x$transposed)) {
    return(x$scale * crossprod(
      x$Z_i,
      .float32KernelMatMult(x$K, x$Z_j %*% w)
    ))
  }
  x$scale * crossprod(
    x$Z_j,
    .float32KernelMatMult(t(x$K), x$Z_i %*% w)
  )
}

.genespace_cross_bilinear <- function(w_i, x, w_j) {
  as.numeric(crossprod(w_i, .genespace_cross_mult(x, w_j)))
}

#' Retrieve cross-covariance matrix handling directional asymmetry
#' @param C_cross_s Cross-covariance list for one slide
#' @param ct_i First cell type
#' @param ct_j Second cell type
#' @return A cross-covariance matrix or matrix-free operator for
#'   C_{ct_i, ct_j}.
#' @noRd
.get_C_cross <- function(C_cross_s, ct_i, ct_j) {
  key_forward <- paste0(ct_i, "-", ct_j)
  key_reverse <- paste0(ct_j, "-", ct_i)

  if (key_forward %in% names(C_cross_s)) {
    return(C_cross_s[[key_forward]])
  } else if (key_reverse %in% names(C_cross_s)) {
    reverse <- C_cross_s[[key_reverse]]
    if (inherits(reverse, "CoProGeneCrossOperator")) {
      return(.transpose_genespace_cross_operator(reverse))
    }
    return(t(reverse))
  } else {
    stop(paste("Cross-covariance not found for pair:", ct_i, "-", ct_j))
  }
}

#' Compute per-slide score standard deviations for all cell types
#'
#' Because `.prepareGeneSpaceData()` centers every (slide, cell type) block,
#' `mean(Z w) = 0` exactly and this is the per-slide standard deviation of the
#' scores. It is the *unsmoothed* second moment -- no kernel enters -- so the
#' objective's ratio has a kernel-smoothed numerator over an unsmoothed
#' denominator, and is bounded by the kernel's largest singular value rather
#' than by 1. `.sumcorSigma()` uses the same convention in PC space.
#'
#' @param w_list Named list of weight vectors (each G x 1 matrix)
#' @param C_self_slide List of per-slide self-covariance:
#'   \code{C_self_slide[[slide]][[ct]]}
#' @param slides Character vector of slide IDs
#' @param cell_types Character vector of cell type names
#' @param objective `"sumcor"` returns the per-slide scales; `"sumcov"` returns
#'   1 everywhere, which turns the objective below into the plain sum of
#'   kernel-smoothed cross-covariances.
#' @return Named list: \code{sigma_all[[slide]][[ct]]} = scalar sigma value (floored at 1e-12)
#' @noRd
.compute_per_slide_sigma <- function(w_list, C_self_slide, slides, cell_types,
                                     objective = "sumcor") {
  sigma_all <- setNames(vector("list", length(slides)), slides)
  unit <- identical(objective, "sumcov")
  for (s in slides) {
    sigma_all[[s]] <- setNames(vector("list", length(cell_types)), cell_types)
    for (ct in cell_types) {
      if (unit) {
        sigma_all[[s]][[ct]] <- 1
        next
      }
      w <- w_list[[ct]]
      val <- .genespace_self_quad(C_self_slide[[s]][[ct]], w)
      sigma_all[[s]][[ct]] <- max(sqrt(max(val, 0)), 1e-12)
    }
  }
  sigma_all
}

#' Refresh the per-slide scales for a single cell type
#'
#' The Gauss-Seidel sweep updates one block at a time and later blocks in the
#' same sweep must divide by current scales, so only the block just updated is
#' recomputed rather than the whole structure.
#' @noRd
.refresh_slide_sigma <- function(sigma_all, w_list, C_self_slide, slides, ct,
                                 objective = "sumcor") {
  if (identical(objective, "sumcov")) return(sigma_all)
  for (s in slides) {
    val <- .genespace_self_quad(C_self_slide[[s]][[ct]], w_list[[ct]])
    sigma_all[[s]][[ct]] <- max(sqrt(max(val, 0)), 1e-12)
  }
  sigma_all
}

#' Compute the P1b objective value for monitoring
#' @param w_list Named list of weight vectors
#' @param C_self_slide Per-slide self-covariance matrices
#' @param C_cross_slide Per-slide cross-covariance matrices
#' @param slides Slide IDs
#' @param cell_types Cell type names
#' @param objective `"sumcor"` (per-slide self-normalized) or `"sumcov"` (plain
#'   sum of kernel-smoothed cross-covariances).
#' @return Scalar objective value
#' @importFrom utils combn
#' @noRd
.compute_p1b_objective <- function(w_list, C_self_slide, C_cross_slide,
                                   slides, cell_types,
                                   objective = "sumcor") {
  S <- length(slides)
  sigma_all <- .compute_per_slide_sigma(w_list, C_self_slide, slides,
                                        cell_types, objective)

  obj <- 0
  pairs <- combn(cell_types, 2, simplify = FALSE)

  for (s in slides) {
    for (pair in pairs) {
      ct_i <- pair[1]
      ct_j <- pair[2]
      C_ij <- .get_C_cross(C_cross_slide[[s]], ct_i, ct_j)
      rho_s <- .genespace_cross_bilinear(
        w_list[[ct_i]], C_ij, w_list[[ct_j]]
      ) /
        (sigma_all[[s]][[ct_i]] * sigma_all[[s]][[ct_j]])
      obj <- obj + rho_s
    }
  }
  obj / S
}

# ============================================================================
# Exported Optimization Functions
# ============================================================================

#' Gene-space average per-slide CCA — first component
#'
#' Power iteration to find the first canonical component that maximizes
#' the average per-slide canonical correlation across all slides. Each slide's
#' contribution is self-normalized by its own score standard deviation.
#'
#' @param C_self_slide Named list of per-slide self-covariance matrices or
#'   matrix-free operators.
#' @param C_cross_slide Named list of per-slide cross-covariance matrices or
#'   matrix-free operators.
#' @param slides Character vector of slide IDs.
#' @param cell_types Character vector of cell type names.
#' @param max_iter Maximum iterations (default 3000). Must be >= 1.
#' @param tol Convergence tolerance on max weight change (default 1e-6).
#' @param verbose Print progress every 500 iterations (default TRUE).
#' @param sweep Which block sweep to use.
#'
#'   `"gauss-seidel"` (default) updates each cell type using the blocks already
#'   updated in the current sweep. With `sigma` held fixed this makes each block
#'   update the exact maximizer over that block -- the objective is linear in
#'   `w_i` once the others are fixed -- so `w_i' g_i = \|g_i\| >= 0` always and
#'   the iteration cannot come to rest on a negative-objective solution.
#'
#'   `"jacobi"` updates every block from the previous iterate. All blocks then
#'   move off stale values, which is not coordinate ascent: the iteration can
#'   settle on the *negative* singular pair as a period-2 sign orbit, which
#'   `check_convergence()`'s sign-tolerant test accepts as converged. That is
#'   what the post-hoc sign flip below exists to repair, and it is applied only
#'   on this path. The flip is a valid repair for two cell types, where flipping
#'   one block negates the single pairwise term; it is **not** valid for three
#'   or more, where flipping one block negates only the terms touching it and
#'   can lower the objective. Retained so results computed before
#'   `"gauss-seidel"` became the default can be reproduced exactly.
#' @param objective `"sumcor"` (default) divides each slide's cross term by that
#'   slide's own score scales. `"sumcov"` fixes every scale at 1, giving the
#'   plain sum of kernel-smoothed cross-covariances -- the gene-space
#'   counterpart of [runSkrCCA()]'s default.
#'
#' @return Named list of weight vectors, one per cell type (each a G x 1 matrix).
#' @importFrom stats rnorm
#' @keywords internal
#' @export
optimize_genespace_avg_corr <- function(C_self_slide, C_cross_slide,
                                        slides, cell_types,
                                        max_iter = 3000, tol = 1e-6,
                                        verbose = TRUE,
                                        sweep = c("gauss-seidel", "jacobi"),
                                        objective = c("sumcor", "sumcov")) {
  sweep <- match.arg(sweep)
  objective <- match.arg(objective)
  if (length(cell_types) < 2) {
    stop("Gene-space CCA requires at least 2 cell types. Found: ",
         paste(cell_types, collapse = ", "))
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1) {
    stop("max_iter must be a positive integer.")
  }

  S <- length(slides)
  n_genes <- .genespace_n_genes(
    C_self_slide[[slides[1]]][[cell_types[1]]]
  )

  # Initialize with random unit vectors
  w_list <- setNames(
    lapply(cell_types, function(ct) {
      v <- matrix(rnorm(n_genes), ncol = 1)
      v / sqrt(sum(v^2))
    }),
    cell_types
  )

  for (iter in seq_len(max_iter)) {
    w_list_old <- w_list

    # Compute per-slide sigmas
    sigma_all <- .compute_per_slide_sigma(w_list, C_self_slide, slides,
                                          cell_types, objective)

    # Update each cell type. The update below is the gradient of f_avg w.r.t.
    # w_i evaluated with sigma_i, sigma_j held FIXED at their current values
    # (frozen-sigma surrogate). The full gradient would also include
    # a -rho * C_ii * w_i / sigma_i^2 correction term from differentiating
    # 1/sigma_i; omitting it makes this an ALS-style alternating maximization,
    # not exact coordinate ascent. Standard for generalized power methods
    # (NIPALS treats denominators as fixed within a sweep).
    #
    # Which iterate the cross term reads is the sweep choice: `w_source` is the
    # in-progress `w_list` under Gauss-Seidel and the frozen `w_list_old` under
    # Jacobi. See the `sweep` argument for why that decides whether a sign
    # repair is needed at all.
    for (ct_i in cell_types) {
      update <- matrix(0, nrow = n_genes, ncol = 1)
      w_source <- if (identical(sweep, "gauss-seidel")) w_list else w_list_old

      for (s in slides) {
        sig_i <- sigma_all[[s]][[ct_i]]

        for (ct_j in cell_types) {
          if (ct_j == ct_i) next
          sig_j <- sigma_all[[s]][[ct_j]]
          C_ij <- .get_C_cross(C_cross_slide[[s]], ct_i, ct_j)
          update <- update + (1 / sig_i) *
            .genespace_cross_mult(C_ij, w_source[[ct_j]] / sig_j)
        }
      }
      update <- update / S

      # Normalize. Zero-norm means the cross-covariance with all other cell
      # types vanished for w_i (degenerate). Keeping the previous iterate
      # rather than overwriting with random noise; we warn so the user knows.
      norm_val <- sqrt(sum(update^2))
      if (norm_val > 0) {
        w_list[[ct_i]] <- update / norm_val
        if (identical(sweep, "gauss-seidel")) {
          sigma_all <- .refresh_slide_sigma(sigma_all, w_list, C_self_slide,
                                            slides, ct_i, objective)
        }
      } else {
        warning(sprintf(
          "Zero gradient norm for cell type '%s' at iter %d; keeping previous weight.",
          ct_i, iter
        ))
      }
    }

    # Check convergence
    max_diff <- check_convergence(w_list, w_list_old, cell_types)

    if (verbose && (iter %% 500 == 0 || iter == 1)) {
      obj <- .compute_p1b_objective(w_list, C_self_slide, C_cross_slide,
                                    slides, cell_types, objective)
      message(sprintf("  Iter %d: max_diff = %.2e, objective = %.4f", iter, max_diff, obj))
    }

    if (max_diff <= tol) {
      if (verbose) {
        obj <- .compute_p1b_objective(w_list, C_self_slide, C_cross_slide,
                                      slides, cell_types, objective)
        message(sprintf("  Converged at iteration %d (max_diff = %.2e, objective = %.4f)",
                        iter, max_diff, obj))
      }
      break
    }
  }

  if (iter == max_iter && max_diff > tol) {
    warning(sprintf("Did not converge after %d iterations (max_diff = %.2e)", max_iter, max_diff))
  }

  # Ensure matrix format
  for (ct in cell_types) {
    if (!is.matrix(w_list[[ct]])) w_list[[ct]] <- matrix(w_list[[ct]], ncol = 1)
  }

  if (identical(sweep, "jacobi")) {
    # Legacy sign repair. Only correct for two cell types; see the `sweep`
    # argument. Gauss-Seidel needs none, because it cannot converge here.
    obj <- .compute_p1b_objective(w_list, C_self_slide, C_cross_slide,
                                  slides, cell_types, objective)
    if (obj < 0) {
      w_list[[cell_types[1]]] <- -w_list[[cell_types[1]]]
      if (length(cell_types) > 2) {
        warning(sprintf(
          paste0("sweep = \"jacobi\" converged to a negative objective (%.6f) ",
                 "with %d cell types. Flipping one block negates only the ",
                 "pairs touching it, so this repair can lower the objective. ",
                 "Use sweep = \"gauss-seidel\"."),
          obj, length(cell_types)
        ), call. = FALSE)
      }
    }
  }

  w_list
}

#' Gene-space average per-slide CCA — subsequent components
#'
#' Computes components 2 through nCC using Gram-Schmidt deflation in weight
#' space. After computing the gradient update for each cell type, the update
#' is orthogonalized against all previous CC directions before normalizing.
#'
#' @param C_self_slide Per-slide self-covariance matrices (same as first component).
#' @param C_cross_slide Per-slide cross-covariance matrices.
#' @param slides Slide IDs.
#' @param cell_types Cell type names.
#' @param w_list Named list of weight matrices from previous components.
#'   Each entry is a G x k matrix where k = number of components already computed.
#' @param nCC Total number of components desired (must be > existing components).
#' @param max_iter Maximum iterations per component.
#' @param tol Convergence tolerance.
#' @param verbose Print progress.
#' @inheritParams optimize_genespace_avg_corr
#'
#' @return Named list of weight matrices, each G x nCC.
#' @keywords internal
#' @export
optimize_genespace_avg_corr_n <- function(C_self_slide, C_cross_slide,
                                          slides, cell_types,
                                          w_list, nCC = 2,
                                          max_iter = 3000, tol = 1e-6,
                                          verbose = TRUE,
                                          sweep = c("gauss-seidel", "jacobi"),
                                          objective = c("sumcor", "sumcov")) {
  sweep <- match.arg(sweep)
  objective <- match.arg(objective)
  if (!is.numeric(max_iter) || length(max_iter) != 1 || max_iter < 1) {
    stop("max_iter must be a positive integer.")
  }
  S <- length(slides)
  n_genes <- .genespace_n_genes(
    C_self_slide[[slides[1]]][[cell_types[1]]]
  )
  k_start <- ncol(w_list[[cell_types[1]]])

  if (nCC <= k_start) {
    stop(sprintf("nCC (%d) must be greater than existing components (%d)", nCC, k_start))
  }

  for (cc in (k_start + 1):nCC) {
    if (verbose) message(sprintf("  Finding CC %d ...", cc))

    # Initialize current component with random unit vector
    w_current <- setNames(
      lapply(cell_types, function(ct) {
        v <- matrix(rnorm(n_genes), ncol = 1)
        v / sqrt(sum(v^2))
      }),
      cell_types
    )

    for (iter in seq_len(max_iter)) {
      w_current_old <- w_current

      # Compute per-slide sigmas using current weights
      sigma_all <- .compute_per_slide_sigma(w_current, C_self_slide, slides,
                                            cell_types, objective)

      # Update each cell type. Same frozen-sigma surrogate and same sweep
      # choice as the first component; see optimize_genespace_avg_corr().
      for (ct_i in cell_types) {
        update <- matrix(0, nrow = n_genes, ncol = 1)
        w_source <- if (identical(sweep, "gauss-seidel")) {
          w_current
        } else {
          w_current_old
        }

        for (s in slides) {
          sig_i <- sigma_all[[s]][[ct_i]]
          for (ct_j in cell_types) {
            if (ct_j == ct_i) next
            sig_j <- sigma_all[[s]][[ct_j]]
            C_ij <- .get_C_cross(C_cross_slide[[s]], ct_i, ct_j)
            update <- update + (1 / sig_i) *
              .genespace_cross_mult(C_ij, w_source[[ct_j]] / sig_j)
          }
        }
        update <- update / S

        # Gram-Schmidt: project out all previous CC directions
        for (prev_cc in seq_len(cc - 1)) {
          prev_w <- w_list[[ct_i]][, prev_cc, drop = FALSE]
          proj <- as.numeric(t(update) %*% prev_w)
          update <- update - proj * prev_w
        }

        # Normalize. Zero-norm here means deflation has exhausted the signal
        # subspace for this cell type — the random init from the start of
        # this cc would otherwise be silently appended as a "canonical
        # component". Warn and keep the deflated update at zero so the
        # caller can detect the degenerate component (weight is all-zero).
        norm_val <- sqrt(sum(update^2))
        if (norm_val > 0) {
          w_current[[ct_i]] <- update / norm_val
          if (identical(sweep, "gauss-seidel")) {
            sigma_all <- .refresh_slide_sigma(sigma_all, w_current, C_self_slide,
                                              slides, ct_i, objective)
          }
        } else {
          warning(sprintf(
            "Zero norm after Gram-Schmidt deflation for cell type '%s' at CC %d; signal subspace likely exhausted.",
            ct_i, cc
          ))
          w_current[[ct_i]] <- matrix(0, nrow = n_genes, ncol = 1)
        }
      }

      # Check convergence
      max_diff <- check_convergence(w_current, w_current_old, cell_types)

      if (verbose && (iter %% 500 == 0 || iter == 1)) {
        message(sprintf("    Iter %d: max_diff = %.2e", iter, max_diff))
      }

      if (max_diff <= tol) {
        if (verbose) {
          message(sprintf("    CC %d converged at iteration %d", cc, iter))
        }
        break
      }
    }

    if (iter == max_iter && max_diff > tol) {
      warning(sprintf("CC %d did not converge after %d iterations (max_diff = %.2e)",
                      cc, max_iter, max_diff))
    }

    if (identical(sweep, "jacobi")) {
      # Legacy sign repair; see optimize_genespace_avg_corr(). Not applied under
      # Gauss-Seidel, which cannot converge to a negative objective.
      obj_cc <- .compute_p1b_objective(w_current, C_self_slide, C_cross_slide,
                                       slides, cell_types, objective)
      if (obj_cc < 0) {
        w_current[[cell_types[1]]] <- -w_current[[cell_types[1]]]
        if (length(cell_types) > 2) {
          warning(sprintf(
            paste0("sweep = \"jacobi\" gave CC %d a negative objective (%.6f) ",
                   "with %d cell types; the one-block flip can lower it. ",
                   "Use sweep = \"gauss-seidel\"."),
            cc, obj_cc, length(cell_types)
          ), call. = FALSE)
        }
      }
    }

    # Append this component to w_list
    for (ct in cell_types) {
      w_list[[ct]] <- cbind(w_list[[ct]], w_current[[ct]])
    }
  }

  w_list
}
