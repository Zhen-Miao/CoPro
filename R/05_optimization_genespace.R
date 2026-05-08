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

#' Retrieve cross-covariance matrix handling directional asymmetry
#' @param C_cross_s Cross-covariance list for one slide
#' @param ct_i First cell type
#' @param ct_j Second cell type
#' @return G x G cross-covariance matrix C_{ct_i, ct_j}
#' @noRd
.get_C_cross <- function(C_cross_s, ct_i, ct_j) {
  key_forward <- paste0(ct_i, "-", ct_j)
  key_reverse <- paste0(ct_j, "-", ct_i)

  if (key_forward %in% names(C_cross_s)) {
    return(C_cross_s[[key_forward]])
  } else if (key_reverse %in% names(C_cross_s)) {
    return(t(C_cross_s[[key_reverse]]))
  } else {
    stop(paste("Cross-covariance not found for pair:", ct_i, "-", ct_j))
  }
}

#' Compute per-slide score standard deviations for all cell types
#' @param w_list Named list of weight vectors (each G x 1 matrix)
#' @param C_self_slide List of per-slide self-covariance: C_self_slide[[slide]][[ct]]
#' @param slides Character vector of slide IDs
#' @param cell_types Character vector of cell type names
#' @return Named list: sigma_all[[slide]][[ct]] = scalar sigma value (floored at 1e-12)
#' @noRd
.compute_per_slide_sigma <- function(w_list, C_self_slide, slides, cell_types) {
  sigma_all <- setNames(vector("list", length(slides)), slides)
  for (s in slides) {
    sigma_all[[s]] <- setNames(vector("list", length(cell_types)), cell_types)
    for (ct in cell_types) {
      w <- w_list[[ct]]
      val <- as.numeric(t(w) %*% C_self_slide[[s]][[ct]] %*% w)
      sigma_all[[s]][[ct]] <- max(sqrt(max(val, 0)), 1e-12)
    }
  }
  sigma_all
}

#' Compute the P1b objective value for monitoring
#' @param w_list Named list of weight vectors
#' @param C_self_slide Per-slide self-covariance matrices
#' @param C_cross_slide Per-slide cross-covariance matrices
#' @param slides Slide IDs
#' @param cell_types Cell type names
#' @return Scalar objective value
#' @noRd
.compute_p1b_objective <- function(w_list, C_self_slide, C_cross_slide,
                                   slides, cell_types) {
  S <- length(slides)
  sigma_all <- .compute_per_slide_sigma(w_list, C_self_slide, slides, cell_types)

  obj <- 0
  pairs <- combn(cell_types, 2, simplify = FALSE)

  for (s in slides) {
    for (pair in pairs) {
      ct_i <- pair[1]
      ct_j <- pair[2]
      C_ij <- .get_C_cross(C_cross_slide[[s]], ct_i, ct_j)
      rho_s <- as.numeric(t(w_list[[ct_i]]) %*% C_ij %*% w_list[[ct_j]]) /
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
#' @param C_self_slide Named list of per-slide self-covariance matrices.
#'   Structure: \code{C_self_slide[[slide]][[cell_type]]} = G x G matrix.
#' @param C_cross_slide Named list of per-slide cross-covariance matrices.
#'   Structure: \code{C_cross_slide[[slide]][["ctA-ctB"]]} = G x G matrix.
#' @param slides Character vector of slide IDs.
#' @param cell_types Character vector of cell type names.
#' @param max_iter Maximum iterations (default 3000).
#' @param tol Convergence tolerance on max weight change (default 1e-6).
#' @param verbose Print progress every 500 iterations (default TRUE).
#'
#' @return Named list of weight vectors, one per cell type (each a G x 1 matrix).
#' @export
optimize_genespace_avg_corr <- function(C_self_slide, C_cross_slide,
                                        slides, cell_types,
                                        max_iter = 3000, tol = 1e-6,
                                        verbose = TRUE) {
  S <- length(slides)
  n_genes <- nrow(C_self_slide[[slides[1]]][[cell_types[1]]])

  # Initialize with random unit vectors
  set.seed(NULL)
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
    sigma_all <- .compute_per_slide_sigma(w_list, C_self_slide, slides, cell_types)

    # Update each cell type
    for (ct_i in cell_types) {
      update <- matrix(0, nrow = n_genes, ncol = 1)

      for (s in slides) {
        sig_i <- sigma_all[[s]][[ct_i]]

        for (ct_j in cell_types) {
          if (ct_j == ct_i) next
          sig_j <- sigma_all[[s]][[ct_j]]
          C_ij <- .get_C_cross(C_cross_slide[[s]], ct_i, ct_j)
          update <- update + (1 / sig_i) * C_ij %*% (w_list[[ct_j]] / sig_j)
        }
      }
      update <- update / S

      # Normalize
      norm_val <- sqrt(sum(update^2))
      if (norm_val > 0) {
        w_list[[ct_i]] <- update / norm_val
      }
    }

    # Check convergence
    max_diff <- check_convergence(w_list, w_list_old, cell_types)

    if (verbose && (iter %% 500 == 0 || iter == 1)) {
      obj <- .compute_p1b_objective(w_list, C_self_slide, C_cross_slide,
                                    slides, cell_types)
      message(sprintf("  Iter %d: max_diff = %.2e, objective = %.4f", iter, max_diff, obj))
    }

    if (max_diff <= tol) {
      if (verbose) {
        obj <- .compute_p1b_objective(w_list, C_self_slide, C_cross_slide,
                                      slides, cell_types)
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
#' @param nCC Total number of components desired (must be > ncol(w_list[[1]])).
#' @param max_iter Maximum iterations per component.
#' @param tol Convergence tolerance.
#' @param verbose Print progress.
#'
#' @return Named list of weight matrices, each G x nCC.
#' @export
optimize_genespace_avg_corr_n <- function(C_self_slide, C_cross_slide,
                                          slides, cell_types,
                                          w_list, nCC = 2,
                                          max_iter = 3000, tol = 1e-6,
                                          verbose = TRUE) {
  S <- length(slides)
  n_genes <- nrow(C_self_slide[[slides[1]]][[cell_types[1]]])
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
      sigma_all <- .compute_per_slide_sigma(w_current, C_self_slide, slides, cell_types)

      # Update each cell type
      for (ct_i in cell_types) {
        update <- matrix(0, nrow = n_genes, ncol = 1)

        for (s in slides) {
          sig_i <- sigma_all[[s]][[ct_i]]
          for (ct_j in cell_types) {
            if (ct_j == ct_i) next
            sig_j <- sigma_all[[s]][[ct_j]]
            C_ij <- .get_C_cross(C_cross_slide[[s]], ct_i, ct_j)
            update <- update + (1 / sig_i) * C_ij %*% (w_current[[ct_j]] / sig_j)
          }
        }
        update <- update / S

        # Gram-Schmidt: project out all previous CC directions
        for (prev_cc in seq_len(cc - 1)) {
          prev_w <- w_list[[ct_i]][, prev_cc, drop = FALSE]
          proj <- as.numeric(t(update) %*% prev_w)
          update <- update - proj * prev_w
        }

        # Normalize
        norm_val <- sqrt(sum(update^2))
        if (norm_val > 0) {
          w_current[[ct_i]] <- update / norm_val
        }
      }

      # Check convergence
      max_diff <- check_convergence(w_current, w_current_old, cell_types)

      if (verbose && (iter %% 500 == 0)) {
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

    # Append this component to w_list
    for (ct in cell_types) {
      w_list[[ct]] <- cbind(w_list[[ct]], w_current[[ct]])
    }
  }

  w_list
}
