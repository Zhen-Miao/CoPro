# Diagnostics use the whitened PC coordinates in which SUMCOR is optimized.
# They record numerical behavior without changing the objective or the fit.

.newSumcorDiagnostics <- function(ops) {
  conditioning <- do.call(rbind, lapply(ops$slides, function(s) {
    do.call(rbind, lapply(ops$cell_types, function(ct) {
      G <- ops$G[[s]][[ct]]
      values <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
      largest <- max(values)
      smallest <- min(values)
      threshold <- nrow(G) * .Machine$double.eps * max(largest, 0)
      rank <- sum(values > threshold)
      data.frame(
        slide = s, cell_type = ct, n_cells = ops$n[[s]][[ct]],
        n_features = nrow(G), rank = rank,
        min_eigenvalue = smallest, max_eigenvalue = largest,
        eigenvalue_tolerance = threshold,
        condition_number = if (rank == nrow(G)) largest / smallest else Inf,
        positive_definite = rank == nrow(G), stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(conditioning) <- NULL
  list(
    components = NULL, conditioning = conditioning, score_norms = NULL,
    objective_traces = list(), coordinate_system = "whitened_pc",
    sigma_floor = .SUMCOR_SIGMA_FLOOR
  )
}

# Append diagnostics for one returned axis. A supplied axis is measured but
# never labeled as converged by this optimizer. Direct SUMCOV reductions are
# checked against the SUMCOR residual, including the >2-type iterative case.
.recordSumcorAxis <- function(diagnostics, w_list, ops, slideWeight,
                              component, tol, solver, fit = NULL) {
  w <- lapply(w_list, function(m) m[, component, drop = FALSE])
  prev <- if (component > 1L) {
    lapply(w_list, function(m) m[, seq_len(component - 1L), drop = FALSE])
  } else NULL
  score_norms <- do.call(rbind, lapply(ops$slides, function(s) {
    do.call(rbind, lapply(ops$cell_types, function(ct) {
      value <- as.numeric(crossprod(w[[ct]], ops$G[[s]][[ct]] %*% w[[ct]]))
      norm <- sqrt(max(value, 0))
      data.frame(component = component, slide = s, cell_type = ct,
                 score_norm = norm, floor_active = norm <= .SUMCOR_SIGMA_FLOOR,
                 stringsAsFactors = FALSE)
    }))
  }))
  if (is.null(fit)) {
    gradient <- .sumcorTangentGradient(
      .sumcorGradient(w, ops, slideWeight), w, prev
    )
    residual <- max(vapply(gradient, function(g) sqrt(sum(g^2)), numeric(1)))
    objective <- .sumcorObjective(w, ops, slideWeight)
    supplied <- identical(solver, "supplied")
    converged <- if (supplied) NA else is.finite(residual) && residual <= tol
    status <- if (supplied) "supplied" else if (converged) {
      "gradient_tolerance"
    } else "residual_above_tolerance"
    iterations <- if (!supplied && length(ops$cell_types) <= 2L) 0L else NA_integer_
    floor_encountered <- any(score_norms$floor_active)
    trace <- numeric(0)
  } else {
    residual <- fit$gradient_norm
    objective <- fit$objective
    converged <- fit$converged
    status <- fit$stop_reason
    iterations <- fit$iterations
    floor_encountered <- fit$floor_encountered
    trace <- fit$objective_trace
  }
  diagnostics$components <- rbind(diagnostics$components, data.frame(
    component = component, solver = solver, status = status,
    converged = converged, objective = objective, gradient_norm = residual,
    tolerance = tol, iterations = iterations, floor_encountered = floor_encountered,
    stringsAsFactors = FALSE
  ))
  diagnostics$score_norms <- rbind(diagnostics$score_norms, score_norms)
  diagnostics$objective_traces[[paste0("CC", component)]] <- trace
  diagnostics
}

# Used when the algebraic reduction returns all axes at once, or the caller
# supplies axes without a diagnostic record from this fit.
.sumcorDiagnosticsForAxes <- function(w_list, ops, slideWeight, tol, solver) {
  diagnostics <- .newSumcorDiagnostics(ops)
  nCC <- ncol(w_list[[1L]])
  solver <- rep_len(solver, nCC)
  for (cc in seq_len(nCC)) {
    diagnostics <- .recordSumcorAxis(
      diagnostics, w_list, ops, slideWeight, cc, tol, solver[cc]
    )
  }
  diagnostics
}

#' Inspect saved PC-space SUMCOR solver diagnostics
#'
#' Reads the numerical diagnostics saved by [runSkrCCA()] with
#' `space = "pca", objective = "sumcor"`. No optimization is rerun.
#'
#' @param object A fitted `CoProSingle` or `CoProMulti` object.
#' @param sigma Optional single numeric bandwidth. If omitted, return a named
#'   list of diagnostics for the stored PCA-space bandwidths.
#'
#' @return With `sigma`, a list containing:
#' \describe{
#'   \item{`components`}{One row per axis: solver, stopping status, convergence,
#'     final objective, maximum block projected-gradient norm, tolerance,
#'     iterations, and whether a score denominator reached its floor.}
#'   \item{`conditioning`}{Per-slide, per-cell-type Gram eigenvalues, numerical
#'     rank, condition number, and numerical positive-definiteness. Rank counts
#'     eigenvalues above `n_features * .Machine$double.eps * max_eigenvalue`;
#'     that threshold is also returned.}
#'   \item{`score_norms`}{Unfloored score norms at each returned axis, with
#'     `floor_active` indicating norms at or below `sigma_floor`.}
#'   \item{`objective_traces`}{Initial and accepted objective values for each
#'     full-gradient axis. Empty for supplied axes and SUMCOV reductions.}
#'   \item{`coordinate_system`, `sigma_floor`}{The diagnostic coordinate system
#'     (`"whitened_pc"`) and numerical score-norm floor. This floor is unrelated
#'     to the spatial bandwidth named `sigma`.}
#' }
#' Without `sigma`, these lists are named by the stored `sigma_<value>` keys.
#' A PCA-space fit without saved diagnostics (older objects or SUMCOV)
#' yields `NULL` for that bandwidth; an object without PCA-space fits yields
#' an empty list. Gene-space fitting creates no diagnostic record, but any
#' earlier PCA-space results retained in the same object remain accessible.
#' A requested bandwidth without a PCA-space fit is an error.
#'
#' @details
#' Residuals and conditioning are measured in the whitened coordinates used
#' internally, including when the fit used `scalePCs = FALSE`. For later axes,
#' the gradient is projected off the previously accepted axes as well as the
#' current direction. A small residual is a first-order condition, not a
#' global-optimum or biological-recovery certificate.
#'
#' `solver = "full_gradient"` records the actual stopping reason:
#' `"gradient_tolerance"`, `"max_iter"`, or `"line_search_stalled"`.
#' `floor_encountered` includes the initial point and all evaluated line-search
#' trials, including rejected trials. The warm-start solver is not monitored.
#' `solver = "sumcov_reduction"` means that the ratio problem reduced
#' algebraically to SUMCOV. Its returned residual is checked explicitly;
#' iteration counts are zero for direct one/two-type solutions and unavailable
#' (`NA`) for the iterative multiblock reduction. Only the returned point is
#' checked for denominator flooring in that route.
#' `solver = "supplied"` has `converged = NA`: the axis was not optimized here.
#'
#' The smooth stationarity argument requires positive marginal metrics on the
#' retained feature space and no active denominator floor. Numerical full rank
#' alone does not establish good conditioning or the statistical assumptions
#' of a recovery theorem. Inspect the smallest eigenvalues, condition numbers,
#' score norms, and residuals together. These diagnostics neither add a ridge
#' nor filter additional slides or features.
#'
#' @seealso [runSkrCCA()], [getCCAObjective()]
#' @family scores-and-correlation
#' @export
#' @examples
#' \dontrun{
#' obj <- computePCA(obj, nPCA = 15, center_per_slide = TRUE)
#' obj <- runSkrCCA(obj, space = "pca", objective = "sumcor")
#' diagnostics <- getCCADiagnostics(obj)
#' diagnostics[[1]]$components
#' diagnostics[[1]]$conditioning
#' }
getCCADiagnostics <- function(object, sigma = NULL) {
  if (!inherits(object, "CoPro")) {
    stop("object must be a CoPro object.", call. = FALSE)
  }
  keys <- grep("^sigma_", names(object@skrCCAOut), value = TRUE)
  if (!is.null(sigma)) {
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma)) {
      stop("sigma must be a single finite numeric bandwidth.", call. = FALSE)
    }
    key <- .sigmaName(sigma)
    if (!key %in% keys) {
      stop("No PCA-space fit stored for sigma = ", sigma, ".", call. = FALSE)
    }
    return(attr(object@skrCCAOut[[key]], "ccaDiagnostics", exact = TRUE))
  }
  setNames(lapply(keys, function(key) {
    attr(object@skrCCAOut[[key]], "ccaDiagnostics", exact = TRUE)
  }), keys)
}
