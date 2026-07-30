#' PCA-space SUMCOR optimization for CoPro
#'
#' The SUMCOV objective that `optimize_bilinear*()` maximizes,
#' \eqn{\sum_{i<j} w_i' (\sum_s X_i^{(s)'} K_{ij}^{(s)} X_j^{(s)}) w_j} subject
#' to \eqn{\|w_i\| = 1}, factors exactly as a *slide-weighted* SUMCOR:
#'
#' \deqn{f_{cov}(w) = \sum_{i<j} \sum_s \sqrt{n_i^{(s)} n_j^{(s)}}\,
#'   \sigma_i^{(s)} \sigma_j^{(s)}\, \rho_{ij}^{(s)}}
#'
#' with \eqn{\rho_{ij}^{(s)} = u_i' K u_j / (\|u_i\| \|u_j\|)} and
#' \eqn{u_i^{(s)} = X_i^{(s)} w_i}. The norm constraint pins the *pooled*
#' variance, so per-slide variances stay free and a slide with inflated variance
#' along the canonical direction gets a proportionally larger vote. That is the
#' batch-domination mode SUMCOR removes.
#'
#' Two consequences shape this file. First, the two "slide weight" factors are
#' separable: keeping \eqn{\sqrt{n_i n_j}} preserves "larger slides count for
#' more" while dropping \eqn{\sigma_i \sigma_j} removes the batch-scale
#' sensitivity, which is what `slideWeight = "size"` does. Second, for a single
#' slide the two objectives are the *same problem*: with whitened PCs
#' \eqn{X_i'X_i = (n_i-1) I}, so \eqn{\|w_i\| = 1} is \eqn{Var(X_i w_i) = 1},
#' and SUMCOR -- being scale invariant in each \eqn{w_i} -- attains its maximum
#' on that constraint set, where every denominator is 1 and it reduces to
#' SUMCOV. `optimize_sumcor_pca()` therefore short-circuits to the exact
#' SUMCOV solvers when there is one slide, rather than iterating to the same
#' answer.
#'
#' @name sumcor_optimization
#' @keywords internal
NULL

# Per-slide score scale below which we refuse to divide. Reached when the
# weight lands (nearly) in the null space of a slide's Gram matrix, which is
# possible whenever that slide has fewer cells than PCs. Matches the floor the
# gene-space optimizer uses.
.SUMCOR_SIGMA_FLOOR <- 1e-12

#' Validate and normalize the slide-weight choice
#' @noRd
.resolveSlideWeight <- function(slideWeight) {
  match.arg(slideWeight, c("size", "equal", "covariance"))
}

#' Move the per-slide operators into whitened coordinates
#'
#' Under `scalePCs = FALSE` the PC scores keep their original scale and the CCA
#' constraint is \eqn{w' D w = 1} with \eqn{D = diag(sdev^2)}. Rather than carry
#' that metric through the iteration, the operators are transformed once by
#' \eqn{D^{-1/2}} -- which is exactly the change of variables
#' \eqn{\tilde{X} = X D^{-1/2}}, \eqn{\tilde{w} = D^{1/2} w}, i.e. the
#' `scalePCs = TRUE` matrices -- the whole iteration runs metric-free, and the
#' weights are mapped back at the end by `.unwhitenWeights()`.
#'
#' This is not just tidier, it is *necessary* for `scalePCs` to remain a pure
#' reparametrization when there are more than two cell types. The objective and
#' the constraint set are functions of the canonical scores alone, so the two
#' parametrizations pose the same problem; but the >2-type warm start goes
#' through `initialize_next_component()`, whose singular vectors of
#' \eqn{D_i^{-1/2} Y_{ij} D_j^{-1/2}} are *not* the mapped singular vectors of
#' \eqn{Y_{ij}}. Starting the two runs at non-corresponding points let them
#' converge to different local maxima of a non-concave objective. Whitening up
#' front makes the two runs the identical sequence of operations.
#'
#' @param ops Structure from `.computeSlideOperators()`.
#' @param sdev2_list Named list of squared sdev vectors, or `NULL` for a no-op.
#' @return `ops` with `Y` and `G` transformed.
#' @noRd
.whitenSlideOperators <- function(ops, sdev2_list) {
  if (is.null(sdev2_list)) return(ops)
  inv_sd <- lapply(sdev2_list, function(d) 1 / sqrt(d))

  for (s in ops$slides) {
    for (ct_i in ops$cell_types) {
      ops$G[[s]][[ct_i]] <- sweep(
        sweep(ops$G[[s]][[ct_i]], 1L, inv_sd[[ct_i]], "*"),
        2L, inv_sd[[ct_i]], "*"
      )
      for (ct_j in ops$cell_types) {
        Y <- ops$Y[[s]][[ct_i]][[ct_j]]
        if (is.null(Y)) next
        ops$Y[[s]][[ct_i]][[ct_j]] <- sweep(
          sweep(Y, 1L, inv_sd[[ct_i]], "*"), 2L, inv_sd[[ct_j]], "*"
        )
      }
    }
  }
  ops
}

#' Map whitened weights back to the original PC parametrization
#'
#' \eqn{w = D^{-1/2} \tilde{w}}. A whitened weight with \eqn{\|\tilde{w}\| = 1}
#' maps to one satisfying \eqn{w' D w = 1}, so the CCA constraint holds in the
#' returned parametrization without any further normalization.
#'
#' @param w_list Named list of whitened weight matrices.
#' @param sdev2_list Squared sdev vectors, or `NULL` for a no-op.
#' @noRd
.unwhitenWeights <- function(w_list, sdev2_list) {
  if (is.null(sdev2_list)) return(w_list)
  for (ct in names(w_list)) {
    w_list[[ct]] <- w_list[[ct]] / sqrt(sdev2_list[[ct]])
  }
  w_list
}

#' Map original-parametrization weights into whitened coordinates
#' @noRd
.whitenWeights <- function(w_list, sdev2_list) {
  if (is.null(sdev2_list)) return(w_list)
  for (ct in names(w_list)) {
    w_list[[ct]] <- w_list[[ct]] * sqrt(sdev2_list[[ct]])
  }
  w_list
}

#' Per-slide, per-cell-type score scales for one weight set
#'
#' \eqn{\sigma_i^{(s)} = \sqrt{w_i' G_i^{(s)} w_i} = \|X_i^{(s)} w_i\|} with
#' \eqn{G_i^{(s)} = X_i^{(s)'} X_i^{(s)}}. This is the plain (unsmoothed)
#' per-slide second moment, matching the gene-space convention in
#' `.genespace_self_quad()`: the numerator is kernel smoothed, the denominator
#' is not, so `rho` is bounded by the kernel's largest singular value rather
#' than by 1. Whether the denominator should also be smoothed is a separate
#' question and deliberately not changed here.
#'
#' No division by \eqn{n_i^{(s)}} happens, so this is a norm rather than a
#' root-mean-square. The ratio
#' \eqn{\rho_{ij}^{(s)} = w_i' Y_{ij}^{(s)} w_j / (\sigma_i \sigma_j)} is then
#' cell-count invariant with no further bookkeeping; see `.sumcorSlideWeight()`.
#'
#' @param w_list Named list of single-column weight matrices.
#' @param ops Structure from `.computeSlideOperators()`.
#' @return `sigma[[slide]][[cellType]]`, a scalar per entry, floored.
#' @noRd
.sumcorSigma <- function(w_list, ops) {
  setNames(lapply(ops$slides, function(s) {
    setNames(lapply(ops$cell_types, function(ct) {
      w <- w_list[[ct]]
      val <- as.numeric(crossprod(w, ops$G[[s]][[ct]] %*% w))
      max(sqrt(max(val, 0)), .SUMCOR_SIGMA_FLOOR)
    }), ops$cell_types)
  }), ops$slides)
}

#' Slide weight for one (slide, pair) term
#'
#' `"equal"` gives Kettenring SUMCOR, matching the gene-space objective.
#' `"size"` weights each slide by \eqn{\sqrt{n_i^{(s)} n_j^{(s)}}}, so larger
#' slides count for more without letting per-slide variance back in.
#' `"covariance"` uses \eqn{\sigma_i^{(s)} \sigma_j^{(s)}}, which cancels the
#' denominator and reproduces SUMCOV exactly; it exists so the equivalence can
#' be asserted in a test and is not offered as an analysis choice.
#'
#' Note the cell-count bookkeeping. `.sumcorSigma()` returns
#' \eqn{\|X_i^{(s)} w_i\|}, which already carries a factor \eqn{\sqrt{n_i}}
#' relative to a root-mean-square, so `rho` is cell-count invariant on its own
#' and `"size"` is the factor that *reintroduces* size. That also means
#' `"covariance"` needs no `n` factor: \eqn{\sigma_i \sigma_j \rho_{ij} =
#' w_i' Y_{ij} w_j} is the SUMCOV term exactly.
#' @noRd
.sumcorSlideWeight <- function(slideWeight, ops, s, ct_i, ct_j, sigma_s) {
  switch(
    slideWeight,
    equal = 1,
    size = sqrt(as.numeric(ops$n[[s]][[ct_i]]) * as.numeric(ops$n[[s]][[ct_j]])),
    covariance = sigma_s[[ct_i]] * sigma_s[[ct_j]]
  )
}

#' Does the slide weight stay fixed as the weights move?
#'
#' `"equal"` and `"size"` depend only on the data, so dividing the objective by
#' the total weight is a positive constant and leaves the maximizer alone.
#' `"covariance"` depends on `w`, so dividing by it would turn the objective
#' into a ratio and it would no longer be SUMCOV -- exactly the property the
#' oracle needs. Leave that one unnormalized.
#' @noRd
.sumcorWeightIsConstant <- function(slideWeight) {
  !identical(slideWeight, "covariance")
}

#' The SUMCOR objective value
#'
#' Averaged over the total slide weight so the value is comparable across
#' slide counts and across `slideWeight` settings. The average is a positive
#' constant, so it does not move the maximizer.
#'
#' For one cell type the `i == j` term is used with \eqn{\sigma_i^2} in the
#' denominator, making the objective a weighted sum of Rayleigh quotients --
#' the within-cell-type multi-slide case.
#'
#' @param w_list Named list of single-column weight matrices.
#' @param ops Structure from `.computeSlideOperators()`.
#' @param slideWeight One of `"equal"`, `"size"`, `"covariance"`.
#' @param sigma_all Optional precomputed scales from `.sumcorSigma()`.
#' @return A single numeric value.
#' @noRd
.sumcorObjective <- function(w_list, ops, slideWeight, sigma_all = NULL) {
  if (is.null(sigma_all)) sigma_all <- .sumcorSigma(w_list, ops)
  cell_types <- ops$cell_types
  within <- length(cell_types) == 1L

  pairs <- if (within) {
    list(c(cell_types[[1L]], cell_types[[1L]]))
  } else {
    combn(cell_types, 2, simplify = FALSE)
  }

  total <- 0
  weight_total <- 0
  for (s in ops$slides) {
    for (pair in pairs) {
      ct_i <- pair[[1L]]
      ct_j <- pair[[2L]]
      m <- .sumcorSlideWeight(slideWeight, ops, s, ct_i, ct_j, sigma_all[[s]])
      num <- as.numeric(
        crossprod(w_list[[ct_i]], ops$Y[[s]][[ct_i]][[ct_j]] %*% w_list[[ct_j]])
      )
      total <- total +
        m * num / (sigma_all[[s]][[ct_i]] * sigma_all[[s]][[ct_j]])
      weight_total <- weight_total + m
    }
  }
  if (!.sumcorWeightIsConstant(slideWeight)) return(total)
  if (weight_total <= 0) return(0)
  total / weight_total
}

#' Metric-aware Gram-Schmidt against previously computed axes
#'
#' Deflation happens in weight space rather than operator space because the
#' SUMCOR operator depends on `w` through `sigma`, so there is no fixed operator
#' to deflate. The metric is chosen so `scalePCs` stays a pure
#' reparametrization, mirroring the oblique projection `apply_deflation()` uses:
#' the identity metric for whitened PCs, and `D = diag(sdev^2)` otherwise.
#'
#' Under `scalePCs = TRUE` the pooled Gram matrix is \eqn{(n_i - 1) I}, so
#' identity-metric orthogonality of the weights *is* uncorrelatedness of the
#' pooled canonical scores -- the classical CCA requirement. (In gene space
#' \eqn{C_{ii} \neq I}, so its Gram-Schmidt orthogonalizes weights but not
#' scores.)
#'
#' @param v Candidate direction, single-column matrix.
#' @param prev Matrix of previously accepted axes as columns, or `NULL`.
#' @param d Squared sdev vector, or `NULL` for the identity metric.
#' @return `v` with the span of `prev` removed under the chosen metric.
#' @noRd
.sumcorOrthogonalize <- function(v, prev, d = NULL) {
  if (is.null(prev) || ncol(prev) == 0L) return(v)
  inner <- function(a, b) {
    if (is.null(d)) sum(a * b) else sum(a * d * b)
  }
  v <- as.numeric(v)
  for (k in seq_len(ncol(prev))) {
    p <- as.numeric(prev[, k])
    denom <- inner(p, p)
    if (denom <= 0) next
    v <- v - p * (inner(p, v) / denom)
  }
  matrix(v, ncol = 1L)
}

#' One Gauss-Seidel sweep of the frozen-sigma SUMCOR iteration
#'
#' Blocks are updated in place and later blocks in the same sweep read the
#' already-updated earlier ones -- Gauss-Seidel, not Jacobi. For SUMCOV this
#' makes each block update an exact coordinate maximization (the objective is
#' linear in `w_i` with the other blocks fixed), so the objective cannot
#' decrease and the iteration cannot rest at a negative-objective point. That
#' is why no sign-repair step appears anywhere in this file. Under SUMCOR the
#' same sweep maximizes a frozen-`sigma` surrogate rather than the objective
#' itself, so the caller additionally guards on the true objective.
#'
#' `sigma` for the block just updated is refreshed immediately, so subsequent
#' blocks in the same sweep divide by current scales.
#'
#' @param w_list Weights at the start of the sweep; updated and returned.
#' @param ops Structure from `.computeSlideOperators()`.
#' @param slideWeight Slide-weight choice.
#' @param sdev2_list Optional diagonal CCA metrics per cell type.
#' @param prev_axes Optional named list of previously accepted axis matrices,
#'   for Gram-Schmidt deflation.
#' @return The updated `w_list`.
#' @noRd
.sumcorSweep <- function(w_list, ops, slideWeight, sdev2_list = NULL,
                         prev_axes = NULL) {
  cell_types <- ops$cell_types
  within <- length(cell_types) == 1L
  sigma_all <- .sumcorSigma(w_list, ops)

  for (ct_i in cell_types) {
    update <- matrix(0, nrow = nrow(w_list[[ct_i]]), ncol = 1L)

    for (s in ops$slides) {
      sig_i <- sigma_all[[s]][[ct_i]]

      if (within) {
        m <- .sumcorSlideWeight(slideWeight, ops, s, ct_i, ct_i, sigma_all[[s]])
        update <- update +
          (m / (sig_i * sig_i)) * (ops$Y[[s]][[ct_i]][[ct_i]] %*% w_list[[ct_i]])
        next
      }

      for (ct_j in cell_types) {
        if (ct_j == ct_i) next
        m <- .sumcorSlideWeight(slideWeight, ops, s, ct_i, ct_j, sigma_all[[s]])
        update <- update + (m / sig_i) *
          (ops$Y[[s]][[ct_i]][[ct_j]] %*% (w_list[[ct_j]] / sigma_all[[s]][[ct_j]]))
      }
    }

    d_i <- if (is.null(sdev2_list)) NULL else sdev2_list[[ct_i]]

    # Deflate before normalizing: the axis must be free of earlier directions
    # under the same metric the normalization uses.
    if (!is.null(prev_axes) && !is.null(prev_axes[[ct_i]])) {
      update <- .sumcorOrthogonalize(update, prev_axes[[ct_i]], d_i)
    }

    if (sqrt(sum(as.numeric(update)^2)) <= 0) {
      # Every cross term vanished for this block, or deflation exhausted the
      # subspace. Keep the previous iterate rather than injecting a direction
      # the objective did not choose; the caller reports non-convergence.
      next
    }

    w_list[[ct_i]] <- normalize_gradient_weighted(update, d_i)
    if (!is.null(prev_axes) && !is.null(prev_axes[[ct_i]])) {
      # Normalization can reintroduce a component along the deflated span when
      # the metric is non-identity, so re-project and re-normalize.
      w_list[[ct_i]] <- normalize_vec_weighted(
        .sumcorOrthogonalize(w_list[[ct_i]], prev_axes[[ct_i]], d_i), d_i
      )
    }
    # Refresh only the block just updated -- the others are untouched, so a
    # full .sumcorSigma() recompute would repeat work for every block.
    for (s in ops$slides) {
      val <- as.numeric(
        crossprod(w_list[[ct_i]], ops$G[[s]][[ct_i]] %*% w_list[[ct_i]])
      )
      sigma_all[[s]][[ct_i]] <- max(sqrt(max(val, 0)), .SUMCOR_SIGMA_FLOOR)
    }
  }

  w_list
}

#' Run the SUMCOR iteration for one axis
#'
#' Gauss-Seidel sweeps with three guarantees standing in for the sign-repair
#' step the gene-space Jacobi iteration needs:
#'
#' * the caller supplies a deterministic warm start (the SUMCOV solution), so
#'   there is no RNG dependence and no random initial direction;
#' * the best-objective iterate seen is tracked and returned, so the result can
#'   never be worse on the SUMCOR objective than the warm start;
#' * a sweep that lowers the objective is retried blended back toward the
#'   previous iterate, halving the step until it improves or the step underflows.
#'
#' @param w_init Named list of single-column starting weights.
#' @param ops Structure from `.computeSlideOperators()`.
#' @param slideWeight Slide-weight choice.
#' @param sdev2_list Optional diagonal CCA metrics.
#' @param prev_axes Optional previously accepted axes for deflation.
#' @param max_iter,tol Iteration controls. `tol` applies to both the sign-strict
#'   maximum weight change and the relative objective change.
#' @param verbose Report convergence.
#' @param label Text used in messages.
#' @return A list with `w_list` and `objective`.
#' @noRd
.sumcorIterate <- function(w_init, ops, slideWeight, sdev2_list = NULL,
                           prev_axes = NULL, max_iter = 200, tol = 1e-6,
                           verbose = FALSE, label = "CC") {
  cell_types <- ops$cell_types

  # A deflated axis must start inside the deflated subspace, or the first sweep
  # measures an objective that includes directions already claimed.
  w_list <- w_init
  if (!is.null(prev_axes)) {
    for (ct in cell_types) {
      if (is.null(prev_axes[[ct]])) next
      d_ct <- if (is.null(sdev2_list)) NULL else sdev2_list[[ct]]
      w_list[[ct]] <- normalize_vec_weighted(
        .sumcorOrthogonalize(w_list[[ct]], prev_axes[[ct]], d_ct), d_ct
      )
    }
  }

  obj <- .sumcorObjective(w_list, ops, slideWeight)
  best_w <- w_list
  best_obj <- obj
  converged <- FALSE
  iter <- 0L

  while (iter < max_iter) {
    iter <- iter + 1L
    w_old <- w_list
    obj_old <- obj

    candidate <- .sumcorSweep(w_list, ops, slideWeight, sdev2_list, prev_axes)
    obj_new <- .sumcorObjective(candidate, ops, slideWeight)

    # Backtrack when the frozen-sigma sweep overshot. Blending toward the
    # previous iterate is the same convex-blend device `step_size` provides in
    # bilinear_w_from_Y_resi(); here it is applied only on a decrease.
    step <- 1
    while (obj_new < obj_old && step > 1e-4) {
      step <- step / 2
      blended <- setNames(lapply(cell_types, function(ct) {
        d_ct <- if (is.null(sdev2_list)) NULL else sdev2_list[[ct]]
        mixed <- (1 - step) * w_old[[ct]] + step * candidate[[ct]]
        if (!is.null(prev_axes) && !is.null(prev_axes[[ct]])) {
          mixed <- .sumcorOrthogonalize(mixed, prev_axes[[ct]], d_ct)
        }
        normalize_vec_weighted(mixed, d_ct)
      }), cell_types)
      candidate <- blended
      obj_new <- .sumcorObjective(candidate, ops, slideWeight)
    }

    if (obj_new < obj_old) {
      # Even a tiny step fails to improve: we are at a fixed point of the
      # surrogate that the true objective disagrees with. Stop and keep the
      # best iterate.
      converged <- TRUE
      break
    }

    w_list <- candidate
    obj <- obj_new

    if (obj > best_obj) {
      best_obj <- obj
      best_w <- w_list
    }

    # Sign-strict weight change: unlike the SUMCOV path, a sign flip here is a
    # genuine change of solution, not a harmless power-iteration ambiguity, so
    # check_convergence()'s flip tolerance must not be used.
    weight_diff <- max(vapply(cell_types, function(ct) {
      max(abs(w_list[[ct]] - w_old[[ct]]))
    }, numeric(1)))
    obj_diff <- abs(obj - obj_old) / max(abs(obj_old), 1e-12)

    if (weight_diff <= tol || obj_diff <= tol) {
      converged <- TRUE
      if (verbose) {
        message(sprintf(
          "  %s converged at sweep %d (max weight change = %.2e, objective = %.6f)",
          label, iter, weight_diff, obj
        ))
      }
      break
    }
  }

  if (!converged) {
    warning(sprintf(
      "SUMCOR %s did not converge in %d sweeps; returning the best iterate (objective = %.6f).",
      label, max_iter, best_obj
    ), call. = FALSE)
  }

  if (best_obj < 0) {
    warning(sprintf(
      paste0("SUMCOR %s converged to a negative objective (%.6f). The cell ",
             "types are anti-associated at this sigma, or the axes above it ",
             "have exhausted the positively associated subspace."),
      label, best_obj
    ), call. = FALSE)
  }

  list(w_list = best_w, objective = best_obj)
}

#' PCA-space SUMCOR optimization -- first component
#'
#' Maximizes the per-slide self-normalized objective
#' \deqn{f(w) = \frac{\sum_s \sum_{i<j} m_{ij}^{(s)}\, w_i' Y_{ij}^{(s)} w_j /
#'   (\sigma_i^{(s)} \sigma_j^{(s)})}{\sum_s \sum_{i<j} m_{ij}^{(s)}}}
#' where \eqn{\sigma_i^{(s)} = \sqrt{w_i' X_i^{(s)'} X_i^{(s)} w_i}} and
#' \eqn{m_{ij}^{(s)}} is 1 for `slideWeight = "equal"` or
#' \eqn{\sqrt{n_i^{(s)} n_j^{(s)}}} for `"size"`.
#'
#' With one slide this is the same problem as SUMCOV (see the file header), so
#' the exact SUMCOV solvers are used directly and no iteration is run.
#'
#' @param X_list_all `X_list_all[[slide]][[cellType]]` cell-by-PC matrices.
#' @param flat_kernels Flat list of kernel matrices.
#' @param sigma Kernel bandwidth (numeric scalar).
#' @param slides Slide IDs.
#' @param cell_types Cell types to optimize over.
#' @param slideWeight Per-slide weighting: `"size"` (default) or `"equal"`.
#' @param sdev2_list Optional named list of squared standard deviations per cell
#'   type, supplied when `scalePCs = FALSE`.
#' @param max_iter Maximum Gauss-Seidel sweeps.
#' @param tol Convergence tolerance on weight change and relative objective
#'   change.
#' @param n_cores Cores for the per-slide kernel products.
#' @param verbose Report progress.
#' @param ops Optional precomputed `.computeSlideOperators()` structure.
#' @return Named list of single-column weight matrices, with attributes
#'   `"objective"` and `"slideWeight"`.
#' @keywords internal
#' @export
optimize_sumcor_pca <- function(X_list_all, flat_kernels, sigma, slides,
                                cell_types, slideWeight = "size",
                                sdev2_list = NULL, max_iter = 200,
                                tol = 1e-6, n_cores = 1, verbose = FALSE,
                                ops = NULL) {
  slideWeight <- .resolveSlideWeight(slideWeight)

  if (is.null(ops)) {
    ops <- .computeSlideOperators(X_list_all, flat_kernels, sigma, slides,
                                  cell_types, n_cores)
  }

  # Everything below runs metric-free in whitened coordinates; see
  # .whitenSlideOperators() for why this is required and not merely tidier.
  ops_w <- .whitenSlideOperators(ops, sdev2_list)
  warm <- .sumcorWarmStart(ops_w, cell_types, NULL, nCC = 1L)

  # One slide: SUMCOR and SUMCOV are the same optimization problem, and the
  # SUMCOV route solves it exactly. Iterating would only add tolerance error.
  if (length(ops$slides) == 1L) {
    obj_val <- .sumcorObjective(
      lapply(warm, function(m) m[, 1L, drop = FALSE]), ops_w, slideWeight
    )
    result <- .unwhitenWeights(warm, sdev2_list)
    attr(result, "objective") <- obj_val
    attr(result, "slideWeight") <- slideWeight
    return(result)
  }

  fit <- .sumcorIterate(
    w_init = setNames(lapply(cell_types, function(ct) {
      warm[[ct]][, 1L, drop = FALSE]
    }), cell_types),
    ops = ops_w, slideWeight = slideWeight, sdev2_list = NULL,
    max_iter = max_iter, tol = tol, verbose = verbose, label = "CC 1"
  )

  result <- .unwhitenWeights(fit$w_list, sdev2_list)
  attr(result, "objective") <- fit$objective
  attr(result, "slideWeight") <- slideWeight
  result
}

#' PCA-space SUMCOR optimization -- subsequent components
#'
#' Computes components `k_start + 1` through `nCC` by weight-space Gram-Schmidt
#' deflation in the metric implied by `sdev2_list`, so `scalePCs` remains a pure
#' reparametrization.
#'
#' @inheritParams optimize_sumcor_pca
#' @param w_list Weight matrices holding the components already computed.
#' @param nCC Total number of components wanted.
#' @return Named list of `nPC x nCC` weight matrices, with an `"objectives"`
#'   attribute holding the per-axis objective values.
#' @keywords internal
#' @export
optimize_sumcor_pca_n <- function(X_list_all, flat_kernels, sigma, slides,
                                  cell_types, w_list, nCC = 2,
                                  slideWeight = "size", sdev2_list = NULL,
                                  max_iter = 200, tol = 1e-6, n_cores = 1,
                                  verbose = FALSE, ops = NULL) {
  slideWeight <- .resolveSlideWeight(slideWeight)

  if (is.null(ops)) {
    ops <- .computeSlideOperators(X_list_all, flat_kernels, sigma, slides,
                                  cell_types, n_cores)
  }

  k_start <- ncol(w_list[[cell_types[[1L]]]])
  if (nCC <= k_start) {
    stop(sprintf("nCC (%d) must be greater than the %d component(s) already computed",
                 nCC, k_start))
  }
  max_axes <- min(vapply(cell_types, function(ct) {
    nrow(w_list[[ct]])
  }, integer(1)))
  if (nCC > max_axes) {
    stop(sprintf("nCC (%d) cannot exceed the number of PCs (%d)", nCC, max_axes))
  }

  ops_w <- .whitenSlideOperators(ops, sdev2_list)

  # One slide: all axes come from the exact SUMCOV solvers, for the same reason
  # the first component does.
  if (length(ops$slides) == 1L) {
    result <- .unwhitenWeights(
      .sumcorWarmStart(ops_w, cell_types, NULL, nCC = nCC), sdev2_list
    )
    attr(result, "slideWeight") <- slideWeight
    return(result)
  }

  # The supplied axes arrive in the caller's parametrization; deflation happens
  # in whitened coordinates alongside everything else.
  w_list_w <- .whitenWeights(w_list, sdev2_list)
  objectives <- rep(NA_real_, nCC)

  for (cc in seq(k_start + 1L, nCC)) {
    prev_axes <- setNames(lapply(cell_types, function(ct) {
      w_list_w[[ct]][, seq_len(cc - 1L), drop = FALSE]
    }), cell_types)

    # Warm-start each deflated axis from the corresponding SUMCOV axis, keeping
    # the whole routine deterministic.
    warm_all <- .sumcorWarmStart(ops_w, cell_types, NULL, nCC = cc)
    w_init <- setNames(lapply(cell_types, function(ct) {
      warm_all[[ct]][, cc, drop = FALSE]
    }), cell_types)

    fit <- .sumcorIterate(
      w_init = w_init, ops = ops_w, slideWeight = slideWeight,
      sdev2_list = NULL, prev_axes = prev_axes,
      max_iter = max_iter, tol = tol, verbose = verbose,
      label = sprintf("CC %d", cc)
    )
    objectives[cc] <- fit$objective

    for (ct in cell_types) {
      w_list_w[[ct]] <- cbind(w_list_w[[ct]], fit$w_list[[ct]])
    }
  }

  w_list <- .unwhitenWeights(w_list_w, sdev2_list)
  attr(w_list, "objectives") <- objectives
  attr(w_list, "slideWeight") <- slideWeight
  w_list
}

#' Deterministic SUMCOV warm start for the SUMCOR iteration
#'
#' Uses the exact SUMCOV solvers on the aggregated operators. This replaces the
#' random initialization the gene-space optimizer uses, which made the sign of
#' every axis and its value at the iteration tolerance depend on the RNG stream.
#'
#' @param ops Structure from `.computeSlideOperators()`.
#' @param cell_types Cell types being optimized.
#' @param sdev2_list Optional diagonal CCA metrics.
#' @param nCC Number of axes wanted.
#' @return Named list of `nPC x nCC` weight matrices.
#' @noRd
.sumcorWarmStart <- function(ops, cell_types, sdev2_list, nCC) {
  Y_aggregate <- .aggregateSlideOperators(ops)

  if (length(cell_types) == 1L) {
    return(solve_one_type_eigen(Y_aggregate, cell_types, nCC = nCC,
                                sdev2_list = sdev2_list))
  }
  if (length(cell_types) == 2L) {
    return(solve_two_type_svd(Y_aggregate, cell_types, nCC = nCC,
                              sdev2_list = sdev2_list))
  }

  feature_counts <- setNames(
    vapply(cell_types, function(ct) ncol(ops$G[[ops$slides[[1L]]]][[ct]]),
           integer(1)),
    cell_types
  )
  # A warm start that stops at the iteration cap is fine by construction: the
  # SUMCOR iteration refines it and returns the best-objective iterate, so a
  # non-converged start costs accuracy in neither. Suppressing keeps that
  # implementation detail out of the user's warning stream.
  w1 <- suppressWarnings(bilinear_w_from_Y_resi(
    w_list_new = initialize_next_component(Y_aggregate, cell_types),
    Y_resi = Y_aggregate, n_features = feature_counts,
    max_iter = 200, tol = 1e-6, step_size = 1, sdev2_list = sdev2_list
  ))
  if (nCC == 1L) return(w1)

  # For >2 cell types there is no exact all-axis solver, so reuse the SUMCOV
  # sequential route: deflate the aggregate operator and refine.
  Y_resi <- Y_aggregate
  w_list <- w1
  for (qq in seq_len(nCC - 1L)) {
    Y_resi <- apply_deflation(Y_resi, w_list, qq, cell_types, sdev2_list,
                              deflation = "projection")
    w_next <- suppressWarnings(bilinear_w_from_Y_resi(
      w_list_new = initialize_next_component(Y_resi, cell_types),
      Y_resi = Y_resi, n_features = feature_counts,
      max_iter = 200, tol = 1e-6, step_size = 1, sdev2_list = sdev2_list
    ))
    for (ct in cell_types) {
      w_list[[ct]] <- cbind(w_list[[ct]], w_next[[ct]])
    }
  }
  w_list
}
