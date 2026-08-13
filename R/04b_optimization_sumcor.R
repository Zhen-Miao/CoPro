#' PCA-space SUMCOR optimization for CoPro
#'
#' The SUMCOV objective that `optimize_bilinear*()` maximizes,
#' \eqn{\sum_{i<j} w_i' (\sum_s X_i^{(s)'} K_{ij}^{(s)} X_j^{(s)}) w_j} subject
#' to \eqn{\|w_i\| = 1}, factors exactly as a *slide-weighted* SUMCOR:
#'
#' \deqn{f_{cov}(w) = \sum_{i<j} \sum_s
#'   \sigma_i^{(s)} \sigma_j^{(s)}\, \rho_{ij}^{(s)}}
#'
#' with \eqn{\rho_{ij}^{(s)} = u_i' K u_j / (\|u_i\| \|u_j\|)} and
#' \eqn{u_i^{(s)} = X_i^{(s)} w_i}. There is no \eqn{\sqrt{n_i n_j}} factor
#' here: \eqn{\sigma} is the norm \eqn{\|X_i^{(s)} w_i\|}, not a
#' root-mean-square, so \eqn{\sigma_i \sigma_j \rho_{ij} = w_i' Y_{ij} w_j} is
#' the SUMCOV term already (see `.sumcorSigma()` and `.sumcorSlideWeight()`).
#' The norm constraint pins the *pooled* variance, so per-slide variances stay
#' free and a slide with inflated variance along the canonical direction gets a
#' proportionally larger vote. That is the batch-domination mode SUMCOR removes.
#'
#' Two consequences shape this file. First, the slide weight and the scale
#' factor are separable: `slideWeight = "size"` reintroduces
#' \eqn{\sqrt{n_i n_j}} so larger slides count for more, while dropping
#' \eqn{\sigma_i \sigma_j} removes the batch-scale sensitivity. Second, for a
#' single slide the two objectives often -- but not always -- pose the same
#' problem. With whitened PCs \eqn{X_i'X_i = (n_i-1) I}, so on
#' \eqn{\|w_i\| = 1} the denominators are \eqn{\sigma_i = \sqrt{n_i - 1}} and
#'
#' \deqn{f_{equal}(w) = \sum_{i<j} \frac{w_i' Y_{ij} w_j}{
#'   \sqrt{(n_i-1)(n_j-1)}}}
#'
#' differs from \eqn{f_{cov}} by a *per-pair* constant. A constant that varies
#' by pair leaves the argmax alone only when there is a single pair, or when
#' every \eqn{n_i} is equal. So the reduction to SUMCOV is exact for one or two
#' cell types at any cell counts, and for three or more only when the counts
#' are equal; `.sumcorReducesToSumcov()` is that test. Outside it the criteria
#' have genuinely different maximizers, and the full-gradient optimizer runs
#' rather than short-circuiting. Under `slideWeight = "size"` the mismatch is
#' \eqn{1 + O(1/n)} and usually immaterial; it is `"equal"` -- strict
#' Kettenring SUMCOR and the multi-slide default -- where it can be material.
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

#' The slide token that single-slide operator structures are keyed by
#'
#' Single-slide kernels carry no slide name, so the operator structures invent
#' one. Shared so the two builders cannot drift apart.
#' @noRd
.SINGLE_SLIDE_TOKEN <- ".single_slide"

#' Enumerate the (i, j) terms the objective sums over
#'
#' One cell type gives the within-type term `(i, i)`; otherwise every unordered
#' cross-type pair. Shared by the objective, its gradient, and the reduction
#' test so the three cannot disagree about what is being summed.
#' @noRd
.sumcorPairs <- function(cell_types) {
  if (length(cell_types) == 1L) {
    return(list(c(cell_types[[1L]], cell_types[[1L]])))
  }
  combn(cell_types, 2L, simplify = FALSE)
}

#' Validate the damping factor
#'
#' Same admissible range as the SUMCOV power iteration, so a `step_size` that
#' is valid for one objective is valid for the other.
#' @noRd
.validateStepSize <- function(step_size) {
  if (!is.numeric(step_size) || length(step_size) != 1L ||
      !is.finite(step_size) || step_size <= 0 || step_size > 1) {
    stop("step_size must be a single numeric value in (0, 1]")
  }
  invisible(step_size)
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

#' Do two weight sets share their first `k` axes?
#'
#' Compared by direction only: an axis is defined up to sign, and the two sets
#' can carry different norms because one may have been supplied by the caller
#' rather than produced by a normalizing solver. A zero column can never agree,
#' so it reports `FALSE` rather than dividing by zero.
#'
#' @param a,b Named lists of weight matrices in the same parametrization.
#' @param k Number of leading columns to compare.
#' @param cell_types Cell types to compare over.
#' @param tol Tolerance on `1 - |cos angle|`.
#' @return `TRUE` when every compared column is collinear.
#' @noRd
.axesAgree <- function(a, b, k, cell_types, tol = 1e-8) {
  if (k < 1L) return(TRUE)
  for (ct in cell_types) {
    if (is.null(a[[ct]]) || is.null(b[[ct]])) return(FALSE)
    if (ncol(a[[ct]]) < k || ncol(b[[ct]]) < k) return(FALSE)
    if (nrow(a[[ct]]) != nrow(b[[ct]])) return(FALSE)
    for (j in seq_len(k)) {
      u <- a[[ct]][, j]
      v <- b[[ct]][, j]
      nu <- sqrt(sum(u^2))
      nv <- sqrt(sum(v^2))
      if (nu <= 0 || nv <= 0) return(FALSE)
      if (abs(abs(sum(u * v)) / (nu * nv) - 1) > tol) return(FALSE)
    }
  }
  TRUE
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

#' Does a one-slide SUMCOR problem have the same maximizer as SUMCOV?
#'
#' The reduction requires more than merely having one remaining slide: every
#' per-type Gram matrix must be a scalar identity, \eqn{G_i = c_i I}, so its
#' denominator is constant on \eqn{\|w_i\|=1}. This is true when PCA was fit to
#' that one slide, but is generally false when a multi-slide PCA was fit first
#' and all but one slide were later filtered. Checking the matrices prevents
#' that latter case from taking an invalid shortcut.
#'
#' Once the denominators are constant, SUMCOR is SUMCOV with per-pair weight
#' \eqn{m_{ij}/\sqrt{c_i c_j}}. The maximizers coincide for arbitrary data only
#' when there is at most one pair or all those constants agree.
#' @param ops Structure from `.computeSlideOperators()`, one slide.
#' @param slideWeight One of `"equal"`, `"size"`, `"covariance"`.
#' @return `TRUE` when the SUMCOV solvers give the SUMCOR maximizer exactly.
#' @noRd
.sumcorReducesToSumcov <- function(ops, slideWeight) {
  if (identical(slideWeight, "covariance")) return(TRUE)
  if (length(ops$slides) != 1L) return(FALSE)

  s <- ops$slides[[1L]]
  gram_constants <- vapply(ops$cell_types, function(ct) {
    G <- ops$G[[s]][[ct]]
    c_i <- mean(diag(G))
    target <- diag(c_i, nrow(G))
    scale <- max(1, max(abs(G)), abs(c_i))
    if (max(abs(G - target)) > 1e-8 * scale) return(NA_real_)
    c_i
  }, numeric(1))
  if (anyNA(gram_constants) || any(gram_constants <= 0)) return(FALSE)

  pairs <- .sumcorPairs(ops$cell_types)
  if (length(pairs) <= 1L) return(TRUE)

  pair_constants <- vapply(pairs, function(pair) {
    ct_i <- pair[[1L]]
    ct_j <- pair[[2L]]
    m <- switch(
      slideWeight,
      equal = 1,
      size = sqrt(as.numeric(ops$n[[s]][[ct_i]]) *
                    as.numeric(ops$n[[s]][[ct_j]]))
    )
    m / sqrt(gram_constants[[ct_i]] * gram_constants[[ct_j]])
  }, numeric(1))
  max(pair_constants) - min(pair_constants) <=
    1e-8 * max(1, max(abs(pair_constants)))
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
  pairs <- .sumcorPairs(cell_types)

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

#' Wrap a single-slide PC problem in the per-slide operator structure
#'
#' Single-slide kernels omit a slide token from their flat names, so they
#' cannot be passed through `.computeSlideOperators()`, which always performs
#' a slide-keyed lookup. This small adapter lets an explicit one-slide SUMCOR
#' request use the same exact optimizer as a multi-slide request.
#' @noRd
.computeSingleSlideOperators <- function(X_list, flat_kernels, sigma,
                                         cell_types) {
  slide <- .SINGLE_SLIDE_TOKEN
  Y <- compute_Y_resi(
    X_list, flat_kernels, sigma, cell_types, slide = NULL
  )
  G <- setNames(lapply(cell_types, function(ct) {
    crossprod(X_list[[ct]])
  }), cell_types)
  n <- setNames(vapply(cell_types, function(ct) {
    nrow(X_list[[ct]])
  }, integer(1)), cell_types)
  list(
    Y = setNames(list(Y), slide),
    G = setNames(list(G), slide),
    n = setNames(list(n), slide),
    slides = slide,
    cell_types = cell_types
  )
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

#' Exact Euclidean gradient of the PCA-space SUMCOR objective
#'
#' For one cross-cell-type term
#' \deqn{\rho_{ij}^{(s)} = \frac{w_i'Y_{ij}^{(s)}w_j}{
#'   \sigma_i^{(s)}\sigma_j^{(s)}}}
#' the derivative with respect to \eqn{w_i} is
#' \deqn{\frac{Y_{ij}^{(s)}w_j}{\sigma_i^{(s)}\sigma_j^{(s)}} -
#'   \rho_{ij}^{(s)}\frac{G_i^{(s)}w_i}{(\sigma_i^{(s)})^2}.}
#' The second term is the derivative of the denominator. Omitting it gives the
#' old frozen-scale heuristic, whose fixed points need not be stationary points
#' of SUMCOR. This routine differentiates exactly the value returned by
#' `.sumcorObjective()`, including its constant slide-weight normalization.
#'
#' At the numerical sigma floor the denominator is locally constant, so its
#' derivative is omitted. Adequacy filtering normally keeps the optimizer away
#' from that non-differentiable boundary.
#'
#' `slideWeight = "covariance"` is handled algebraically: the denominators
#' cancel, leaving the exact SUMCOV gradient.
#'
#' @param w_list Named list of single-column weights.
#' @param ops Structure from `.computeSlideOperators()`.
#' @param slideWeight Slide-weight choice.
#' @param sigma_all Optional precomputed scales.
#' @return Named list of gradient columns.
#' @noRd
.sumcorGradient <- function(w_list, ops, slideWeight, sigma_all = NULL) {
  cell_types <- ops$cell_types
  within <- length(cell_types) == 1L
  covariance <- identical(slideWeight, "covariance")
  if (is.null(sigma_all) && !covariance) {
    sigma_all <- .sumcorSigma(w_list, ops)
  }

  gradient <- setNames(lapply(cell_types, function(ct) {
    matrix(0, nrow = nrow(w_list[[ct]]), ncol = 1L)
  }), cell_types)
  weight_total <- 0

  pairs <- .sumcorPairs(cell_types)

  for (s in ops$slides) {
    for (pair in pairs) {
      ct_i <- pair[[1L]]
      ct_j <- pair[[2L]]
      w_i <- w_list[[ct_i]]
      w_j <- w_list[[ct_j]]
      Y <- ops$Y[[s]][[ct_i]][[ct_j]]

      if (covariance) {
        if (within) {
          gradient[[ct_i]] <- gradient[[ct_i]] + (Y + t(Y)) %*% w_i
        } else {
          gradient[[ct_i]] <- gradient[[ct_i]] + Y %*% w_j
          gradient[[ct_j]] <- gradient[[ct_j]] + t(Y) %*% w_i
        }
        next
      }

      sigma_s <- sigma_all[[s]]
      sig_i <- sigma_s[[ct_i]]
      sig_j <- sigma_s[[ct_j]]
      m <- .sumcorSlideWeight(slideWeight, ops, s, ct_i, ct_j, sigma_s)
      weight_total <- weight_total + m

      G_i_w <- ops$G[[s]][[ct_i]] %*% w_i
      q_i <- max(as.numeric(crossprod(w_i, G_i_w)), 0)
      num <- as.numeric(crossprod(w_i, Y %*% w_j))

      if (within) {
        term <- m * (Y + t(Y)) %*% w_i / (sig_i * sig_i)
        if (q_i > .SUMCOR_SIGMA_FLOOR^2) {
          term <- term - m * 2 * num * G_i_w / (sig_i^4)
        }
        gradient[[ct_i]] <- gradient[[ct_i]] + term
        next
      }

      rho <- num / (sig_i * sig_j)
      term_i <- m * (Y %*% w_j) / (sig_i * sig_j)
      if (q_i > .SUMCOR_SIGMA_FLOOR^2) {
        term_i <- term_i - m * rho * G_i_w / (sig_i * sig_i)
      }
      gradient[[ct_i]] <- gradient[[ct_i]] + term_i

      G_j_w <- ops$G[[s]][[ct_j]] %*% w_j
      q_j <- max(as.numeric(crossprod(w_j, G_j_w)), 0)
      term_j <- m * (t(Y) %*% w_i) / (sig_i * sig_j)
      if (q_j > .SUMCOR_SIGMA_FLOOR^2) {
        term_j <- term_j - m * rho * G_j_w / (sig_j * sig_j)
      }
      gradient[[ct_j]] <- gradient[[ct_j]] + term_j
    }
  }

  if (!covariance && weight_total > 0) {
    gradient <- lapply(gradient, `/`, weight_total)
  }
  gradient
}

#' Project an objective gradient onto the feasible tangent space
#'
#' The optimizer itself always runs in whitened coordinates. The feasible set
#' for each block is therefore the unit sphere intersected with the orthogonal
#' complement of the previously accepted axes.
#' @noRd
.sumcorTangentGradient <- function(gradient, w_list, prev_axes = NULL) {
  setNames(lapply(names(w_list), function(ct) {
    g <- gradient[[ct]]
    if (!is.null(prev_axes) && !is.null(prev_axes[[ct]])) {
      g <- .sumcorOrthogonalize(g, prev_axes[[ct]], NULL)
    }
    w <- w_list[[ct]]
    g <- g - w * as.numeric(crossprod(w, g))
    # The two projections commute when `w` and the previous axes are
    # orthogonal. Repeating the first one suppresses accumulated roundoff.
    if (!is.null(prev_axes) && !is.null(prev_axes[[ct]])) {
      g <- .sumcorOrthogonalize(g, prev_axes[[ct]], NULL)
    }
    matrix(as.numeric(g), ncol = 1L)
  }), names(w_list))
}

#' Retract a tangent step back to the feasible product of spheres
#' @noRd
.sumcorRetract <- function(w_list, direction, step, prev_axes = NULL) {
  setNames(lapply(names(w_list), function(ct) {
    candidate <- w_list[[ct]] + step * direction[[ct]]
    if (!is.null(prev_axes) && !is.null(prev_axes[[ct]])) {
      candidate <- .sumcorOrthogonalize(candidate, prev_axes[[ct]], NULL)
    }
    norm <- sqrt(sum(as.numeric(candidate)^2))
    if (!is.finite(norm) || norm <= .SUMCOR_SIGMA_FLOOR) return(NULL)
    matrix(as.numeric(candidate) / norm, ncol = 1L)
  }), names(w_list))
}

#' Run exact projected-gradient SUMCOR optimization for one axis
#'
#' The SUMCOV solution supplies a deterministic warm start. At every iteration
#' the full SUMCOR gradient is projected onto the feasible tangent space and a
#' backtracking Armijo search chooses a step on the product of spheres. Accepted
#' steps are therefore monotone in the *actual* objective, and convergence is
#' based on the norm of its projected gradient (the constrained first-order
#' stationarity condition), not on a frozen-denominator fixed point.
#'
#' Damping (`step_size < 1`) is the same knob the SUMCOV power iteration
#' exposes, expressed in this geometry. There, a damped update is
#' \eqn{\mathrm{normalize}((1-\alpha) w + \alpha w^{+})}. Here \eqn{w} and the
#' retracted candidate \eqn{w^{+} = R_w(t g)} both lie on the great circle
#' through \eqn{w} in direction \eqn{g}, so that blend is *itself* a retraction:
#' \deqn{\mathrm{normalize}\bigl((1-\alpha)w + \alpha R_w(tg)\bigr) = R_w(\tau g),
#'   \qquad \tau = \frac{\alpha t}{(1-\alpha)c + \alpha},
#'   \qquad c = \sqrt{1 + t^2\|g\|^2}.}
#' Damping toward the previous iterate is therefore exactly a shorter step along
#' the same arc, which is why `step_size` is expressed as a step length here
#' rather than as a blend applied after the fact. Applying it to the *trial*
#' step keeps the Armijo test on the point actually returned; blending after a
#' certified step would move off it and forfeit monotonicity.
#'
#' A damped run proposes `step_size` every iteration instead of the adaptive
#' hint, rather than scaling that hint. Scaling would accomplish nothing: the
#' hint is rebuilt from the accepted step each iteration, so a constant factor
#' cancels out of the next proposal and the damping would silently apply to the
#' first iteration only. Replacing the adaptivity with a fixed step is also the
#' closer analogue of the power iteration's fixed \eqn{\alpha}. `step_size = 1`
#' leaves the adaptive iteration bit-for-bit unchanged.
#'
#' @param w_init Named list of single-column starting weights.
#' @param ops Structure from `.computeSlideOperators()`.
#' @param slideWeight Slide-weight choice.
#' @param sdev2_list Optional diagonal CCA metrics.
#' @param prev_axes Optional previously accepted axes for deflation.
#' @param max_iter,tol Iteration controls. `tol` is the maximum projected
#'   gradient norm allowed at convergence.
#' @param step_size Damping factor in (0, 1]. Below 1, every iteration proposes
#'   this fixed step in place of the adaptive hint; at 1 the adaptive iteration
#'   is unchanged. Armijo still runs either way, so damped runs stay monotone.
#' @param verbose Report convergence.
#' @param label Text used in messages.
#' @return A list with `w_list`, `objective`, `gradient_norm`, `iterations`,
#'   and the monotone `objective_trace`.
#' @noRd
.sumcorIterate <- function(w_init, ops, slideWeight, sdev2_list = NULL,
                           prev_axes = NULL, max_iter = 200, tol = 1e-6,
                           step_size = 1, verbose = FALSE, label = "CC") {
  cell_types <- ops$cell_types

  # Carry a non-identity CCA metric only through this change of coordinates.
  # The recursion then optimizes and projects in ordinary Euclidean geometry.
  if (!is.null(sdev2_list)) {
    fit <- .sumcorIterate(
      w_init = .whitenWeights(w_init, sdev2_list),
      ops = .whitenSlideOperators(ops, sdev2_list),
      slideWeight = slideWeight,
      sdev2_list = NULL,
      prev_axes = if (is.null(prev_axes)) NULL else {
        .whitenWeights(prev_axes, sdev2_list)
      },
      max_iter = max_iter, tol = tol, step_size = step_size,
      verbose = verbose, label = label
    )
    fit$w_list <- .unwhitenWeights(fit$w_list, sdev2_list)
    return(fit)
  }

  # Put an arbitrary supplied start exactly on the feasible manifold.
  w_list <- setNames(lapply(cell_types, function(ct) {
    w <- w_init[[ct]]
    if (!is.null(prev_axes) && !is.null(prev_axes[[ct]])) {
      w <- .sumcorOrthogonalize(w, prev_axes[[ct]], NULL)
    }
    norm <- sqrt(sum(as.numeric(w)^2))
    if (!is.finite(norm) || norm <= .SUMCOR_SIGMA_FLOOR) {
      stop("SUMCOR ", label, " has no nonzero feasible warm start for ", ct,
           "; reduce nCC or increase nPCA.")
    }
    matrix(as.numeric(w) / norm, ncol = 1L)
  }), cell_types)

  obj <- .sumcorObjective(w_list, ops, slideWeight)
  objective_trace <- obj
  converged <- FALSE
  iter <- 0L
  gradient_norm <- Inf
  step_hint <- 1
  damped <- step_size < 1
  stalled <- FALSE

  while (iter < max_iter) {
    gradient <- .sumcorGradient(w_list, ops, slideWeight)
    tangent <- .sumcorTangentGradient(gradient, w_list, prev_axes)
    block_norms <- vapply(tangent, function(g) {
      sqrt(sum(as.numeric(g)^2))
    }, numeric(1))
    gradient_norm <- max(block_norms)

    if (!is.finite(gradient_norm)) {
      stop("Non-finite projected gradient in SUMCOR ", label, ".")
    }
    if (gradient_norm <= tol) {
      converged <- TRUE
      break
    }

    slope <- sum(block_norms^2)
    # Damped runs propose a fixed step every iteration; undamped runs propose
    # the adaptive hint. A constant multiplier on an adaptive hint would be
    # absorbed by the hint's own growth within one iteration, so damping has to
    # replace the adaptivity rather than scale it.
    step <- if (damped) step_size else step_hint
    accepted <- FALSE
    candidate <- NULL
    obj_new <- NA_real_
    for (line_search in seq_len(60L)) {
      candidate <- .sumcorRetract(w_list, tangent, step, prev_axes)
      if (!any(vapply(candidate, is.null, logical(1)))) {
        obj_new <- .sumcorObjective(candidate, ops, slideWeight)
        armijo_rhs <- obj + 1e-4 * step * slope
        numerical_slack <- 100 * .Machine$double.eps * max(1, abs(obj))
        if (is.finite(obj_new) && obj_new + numerical_slack >= armijo_rhs) {
          accepted <- TRUE
          break
        }
      }
      step <- step / 2
    }

    if (!accepted) {
      stalled <- TRUE
      break
    }

    w_list <- candidate
    obj <- obj_new
    objective_trace <- c(objective_trace, obj)
    iter <- iter + 1L
    # Reuse successful curvature information without committing to a fixed
    # learning rate; an oversize proposal is harmless because it is searched.
    # Under damping the step is fixed by request, so there is no hint to grow.
    if (!damped) step_hint <- min(step * 1.5, 1e6)
  }

  # Report stationarity at the returned point, not at the previous iterate.
  final_gradient <- .sumcorTangentGradient(
    .sumcorGradient(w_list, ops, slideWeight), w_list, prev_axes
  )
  gradient_norm <- max(vapply(final_gradient, function(g) {
    sqrt(sum(as.numeric(g)^2))
  }, numeric(1)))
  converged <- converged || gradient_norm <= tol

  if (converged && verbose) {
    message(sprintf(
      "  %s converged after %d iteration(s) (projected gradient = %.2e, objective = %.6f)",
      label, iter, gradient_norm, obj
    ))
  }
  if (!converged) {
    warning(sprintf(
      paste0("SUMCOR %s %s after %d iteration(s); returning the last monotone ",
             "iterate (objective = %.6f, projected gradient = %.2e)."),
      label, if (stalled) "line search stalled" else "did not converge",
      iter, obj, gradient_norm
    ), call. = FALSE)
  }

  if (obj < 0) {
    warning(sprintf(
      paste0("SUMCOR %s converged to a negative objective (%.6f). The cell ",
             "types are anti-associated at this sigma, or the axes above it ",
             "have exhausted the positively associated subspace."),
      label, obj
    ), call. = FALSE)
  }

  list(
    w_list = w_list,
    objective = obj,
    gradient_norm = gradient_norm,
    iterations = iter,
    objective_trace = objective_trace
  )
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
#' When one slide has coinciding per-pair constants (see the file header), the
#' existing SUMCOV result is used directly and no SUMCOR iteration is run.
#'
#' @param X_list_all `X_list_all[[slide]][[cellType]]` cell-by-PC matrices.
#' @param flat_kernels Flat list of kernel matrices.
#' @param sigma Kernel bandwidth (numeric scalar).
#' @param slides Slide IDs.
#' @param cell_types Cell types to optimize over.
#' @param slideWeight Per-slide weighting: `"equal"` (default) or `"size"`.
#' @param sdev2_list Optional named list of squared standard deviations per cell
#'   type, supplied when `scalePCs = FALSE`.
#' @param max_iter Maximum projected-gradient iterations.
#' @param tol Convergence tolerance on the projected-gradient norm.
#' @param step_size Damping factor in (0, 1]. Default 1 (undamped). Values below
#'   1 shorten every trial step, the same stabilization the SUMCOV damped power
#'   iteration provides; the Armijo test still runs at the damped step, so the
#'   iteration stays monotone. Also damps the SUMCOV warm start.
#' @param n_cores Cores for the per-slide kernel products.
#' @param verbose Report progress.
#' @param ops Optional precomputed `.computeSlideOperators()` structure.
#' @return Named list of single-column weight matrices, with attributes
#'   `"objective"` and `"slideWeight"`.
#' @keywords internal
#' @export
optimize_sumcor_pca <- function(X_list_all, flat_kernels, sigma, slides,
                                cell_types, slideWeight = "equal",
                                sdev2_list = NULL, max_iter = 200,
                                tol = 1e-6, step_size = 1, n_cores = 1,
                                verbose = FALSE, ops = NULL) {
  slideWeight <- .resolveSlideWeight(slideWeight)
  .validateOptimizerParams(max_iter, tol, step_size)

  if (is.null(ops)) {
    ops <- .computeSlideOperators(X_list_all, flat_kernels, sigma, slides,
                                  cell_types, n_cores)
  }

  # Everything below runs metric-free in whitened coordinates; see
  # .whitenSlideOperators() for why this is required and not merely tidier.
  ops_w <- .whitenSlideOperators(ops, sdev2_list)
  warm <- .sumcorWarmStart(
    ops_w, cell_types, NULL, nCC = 1L, step_size = step_size,
    max_iter = max_iter, tol = tol
  )

  # One slide, and the per-pair constants coincide: SUMCOR and SUMCOV are then
  # the same optimization problem, so the SUMCOV route already targets the
  # requested maximizer. When they do not coincide (three
  # or more cell types with unequal counts) the criteria differ and the
  # iteration below is run instead -- see `.sumcorReducesToSumcov()`.
  if (length(ops$slides) == 1L && .sumcorReducesToSumcov(ops_w, slideWeight)) {
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
    max_iter = max_iter, tol = tol, step_size = step_size,
    verbose = verbose, label = "CC 1"
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
                                  slideWeight = "equal", sdev2_list = NULL,
                                  max_iter = 200, tol = 1e-6, step_size = 1,
                                  n_cores = 1, verbose = FALSE, ops = NULL) {
  slideWeight <- .resolveSlideWeight(slideWeight)
  .validateOptimizerParams(max_iter, tol, step_size)

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

  # The supplied axes arrive in the caller's parametrization; deflation happens
  # in whitened coordinates alongside everything else.
  w_list_w <- .whitenWeights(w_list, sdev2_list)

  # One slide with coinciding per-pair constants: all axes come from the exact
  # SUMCOV solvers, for the same reason the first component does.
  #
  # Only when the supplied axes are the ones those solvers would themselves
  # have produced. `runSkrCCA(transferred_weight_1 = ...)` hands in a first
  # axis that came from somewhere else entirely, and the later axes are then
  # only meaningful conditioned on it; returning the solver's own axes would
  # silently discard the transferred direction. The SUMCOV route already keeps
  # the sequential path in exactly that case.
  if (length(ops$slides) == 1L && .sumcorReducesToSumcov(ops_w, slideWeight)) {
    warm <- .sumcorWarmStart(
      ops_w, cell_types, NULL, nCC = nCC, step_size = step_size,
      max_iter = max_iter, tol = tol
    )
    if (.axesAgree(warm, w_list_w, k_start, cell_types)) {
      result <- .unwhitenWeights(warm, sdev2_list)
      attr(result, "slideWeight") <- slideWeight
      return(result)
    }
  }

  objectives <- rep(NA_real_, nCC)

  for (cc in seq(k_start + 1L, nCC)) {
    prev_axes <- setNames(lapply(cell_types, function(ct) {
      w_list_w[[ct]][, seq_len(cc - 1L), drop = FALSE]
    }), cell_types)

    # Warm-start each deflated axis from the corresponding SUMCOV axis, keeping
    # the whole routine deterministic.
    warm_all <- .sumcorWarmStart(
      ops_w, cell_types, NULL, nCC = cc, step_size = step_size,
      max_iter = max_iter, tol = tol
    )
    w_init <- setNames(lapply(cell_types, function(ct) {
      warm_all[[ct]][, cc, drop = FALSE]
    }), cell_types)

    fit <- .sumcorIterate(
      w_init = w_init, ops = ops_w, slideWeight = slideWeight,
      sdev2_list = NULL, prev_axes = prev_axes,
      max_iter = max_iter, tol = tol, step_size = step_size,
      verbose = verbose, label = sprintf("CC %d", cc)
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
#' @param step_size Damping for the >2-type SUMCOV power iteration, applied to
#'   every axis and not only the first: in the one-slide reducible shortcut this
#'   warm start *is* the returned result, so CC2+ have to honor the requested
#'   damping too, exactly as `optimize_bilinear_n()` does under SUMCOV. The one-
#'   and two-type routes are exact decompositions, so damping does not apply
#'   there.
#' @return Named list of `nPC x nCC` weight matrices.
#' @noRd
.sumcorWarmStart <- function(ops, cell_types, sdev2_list, nCC, step_size = 1,
                             max_iter = 200, tol = 1e-6) {
  .validateOptimizerParams(max_iter, tol, step_size)
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
    max_iter = max_iter, tol = tol, step_size = step_size,
    sdev2_list = sdev2_list
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
      max_iter = max_iter, tol = tol, step_size = step_size,
      sdev2_list = sdev2_list
    ))
    for (ct in cell_types) {
      w_list[[ct]] <- cbind(w_list[[ct]], w_next[[ct]])
    }
  }
  w_list
}
