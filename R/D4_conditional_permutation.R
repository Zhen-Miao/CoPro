#' Fit a single conditional canonical axis (internal kernel)
#'
#' Internal numeric kernel for the sequential step-down permutation test. Given
#' (possibly permuted) PC matrices and the FIXED observed lower-axis weight
#' directions, it returns the leading canonical axis of the residual after those
#' lower directions have been deflated, together with that axis' normalized
#' correlation.
#'
#' @details
#' For axis `k = 1` (no lower directions; `k_minus_1 = 0`) this is exactly the
#' first-component optimization used by [runSkrCCAPermu_FairSigma()], so the
#' `k = 1` conditional test reproduces the fair-sigma CC1 test bit-for-bit.
#'
#' For `k >= 2` it deflates the observed CC1..CC(k-1) directions from the
#' cross-product `Y = t(X_i) K_ij X_j` in feature (PC) space, using the SAME
#' fixed observed directions on every permutation, then optimizes the leading
#' residual component. Deflating each permutation by the same observed
#' directions (rather than by that permutation's own leading axis) is what makes
#' the higher-axis null exchangeable with the observed statistic and removes the
#' anti-conservative bias of the naive per-axis permutation p-value.
#'
#' Under whitened PCs (`scalePCs = TRUE`) the fixed directions are removed with
#' the full projection `(I - u u^T) Y (I - v v^T)`. The simpler rank-one
#' subtraction is not used here: it agrees with the projection only when the
#' fixed directions are singular vectors of the current `Y`, which is generally
#' false after permutation. For `scalePCs = FALSE`, the corresponding weighted
#' oblique projection is used.
#'
#' The expensive `compute_Y_resi()` can be computed once per (permutation,
#' sigma) and reused across axes by passing it as `Y_resi`; deflation does not
#' mutate the supplied structure (copy-on-modify), so the same `Y_resi` is safe
#' to reuse for every `k`.
#'
#' @param PCmats Named list of (possibly permuted) PC matrices, one per cell type.
#' @param flat_kernels Flat kernel list from the CoPro object.
#' @param sigma Kernel bandwidth (numeric).
#' @param cts Cell types of interest.
#' @param W_lower Named list of observed weight matrices whose first
#'   `k_minus_1` columns are the fixed deflation directions; ignored when
#'   `k_minus_1 = 0`.
#' @param k_minus_1 Number of lower axes to deflate (0 for the first axis).
#' @param Y_resi Optional precomputed `compute_Y_resi()` structure for this
#'   (PCmats, sigma); recomputed when `NULL` and `k_minus_1 >= 1`.
#' @param kernel_info Optional precomputed kernel and normalizer information
#'   from `.get_ncorr_kernel_info()`.
#' @param sdev2_list Optional diagonal CCA metric used when `scalePCs = FALSE`.
#' @param grams Optional per-cell-type Gram matrices forwarded to
#'   `.compute_ncorr_quick()`.
#' @param maxIter,tol Optimization controls.
#'
#' @return List with `w` (named list of 1-column weight matrices for axis k) and
#'   `ncorr` (its normalized correlation).
#' @importFrom utils capture.output
#' @keywords internal
.fitConditionalAxis <- function(PCmats, flat_kernels, sigma, cts,
                                W_lower = NULL, k_minus_1 = 0,
                                Y_resi = NULL, kernel_info = NULL,
                                sdev2_list = NULL, grams = NULL,
                                maxIter = 200, tol = 1e-5) {
  if (k_minus_1 <= 0) {
    # First axis: identical to optimize_bilinear(). When Y_resi was already
    # computed by the caller, reuse it instead of repeating X' K X.
    if (is.null(Y_resi)) {
      invisible(utils::capture.output(
        w_k <- optimize_bilinear(
          X_list = PCmats, flat_kernels = flat_kernels, sigma = sigma,
          max_iter = maxIter, tol = tol, sdev2_list = sdev2_list
        )
      ))
      names(w_k) <- cts
    } else {
      initial_weights <- if (length(cts) > 2L) {
        initialize_weights_svd(PCmats, cts)
      } else {
        NULL
      }
      invisible(utils::capture.output(
        w_k <- .solveSumcovFirstFromY(
          Y_resi = Y_resi,
          initial_weights = initial_weights,
          cell_types = cts,
          feature_counts = stats::setNames(
            vapply(PCmats[cts], ncol, integer(1)), cts
          ),
          max_iter = maxIter,
          tol = tol,
          sdev2_list = sdev2_list
        )
      ))
    }
  } else {
    # Conditional axis: deflate the fixed observed lower directions from Y,
    # then optimize the leading residual component.
    if (is.null(Y_resi)) {
      Y_resi <- compute_Y_resi(PCmats, flat_kernels, sigma, cts, slide = NULL)
    }
    Yk <- Y_resi
    for (qq in seq_len(k_minus_1)) {
      # Fixed observed directions are not singular vectors of a permuted Y.
      # Full (weighted when needed) projection is therefore required; rank-one
      # subtraction only agrees with projection on the observed operator.
      Yk <- apply_deflation(
        Yk, W_lower, qq, cts, sdev2_list = sdev2_list,
        deflation = "projection"
      )
    }
    initial_weights <- if (length(cts) > 2L) {
      initialize_next_component(Yk, cts)
    } else {
      NULL
    }
    invisible(utils::capture.output(
      w_k <- .solveSumcovFirstFromY(
        Y_resi = Yk,
        initial_weights = initial_weights,
        cell_types = cts,
        feature_counts = stats::setNames(
          vapply(PCmats[cts], ncol, integer(1)), cts
        ),
        max_iter = maxIter,
        tol = tol,
        sdev2_list = sdev2_list
      )
    ))
  }

  # Ensure single-column matrices
  for (ct in cts) {
    w_k[[ct]] <- matrix(w_k[[ct]][, 1], ncol = 1)
  }

  ncorr <- .compute_ncorr_quick(
    PCmats, w_k, flat_kernels, sigma, cts,
    kernel_info = kernel_info, Y_resi = Y_resi, grams = grams
  )
  list(w = w_k, ncorr = ncorr)
}


#' Conditional (sequential step-down) permutation test across canonical axes
#'
#' Performs a permutation test that controls BOTH the sigma-selection
#' multiplicity (via a fair-sigma max-statistic, as in
#' [runSkrCCAPermu_FairSigma()]) AND the canonical-axis multiplicity (via a
#' sequential conditional step-down test). This is the statistically correct
#' way to assign p-values to CC1, CC2, ... jointly, and it supersedes computing
#' a separate per-axis p-value from a full multi-component permutation, which is
#' anti-conservative for axes `k >= 2` when CC1 is strong.
#'
#' @details
#' ## Why a conditional test
#'
#' The observed CC2 is obtained by deflating the observed (real) CC1 direction
#' and optimizing the residual. A naive permutation p-value compares this to the
#' CC2 of fully re-optimized permutations, where each permutation deflates its
#' OWN leading axis. Those two statistics are produced by different operators
#' and are not exchangeable under the null, which biases the CC2 p-value
#' downward (too many false positives). The conditional test instead deflates
#' every permutation by the SAME fixed observed CC1..CC(k-1) directions, so the
#' observed and null axis-`k` statistics share one operator and are exchangeable.
#'
#' ## The statistic and step-down rule
#'
#' For axis `k`, the statistic is the fair-sigma maximum over the candidate
#' bandwidths of the residual normalized correlation after deflating the fixed
#' observed CC1..CC(k-1) directions. The observed value is read from the stored
#' `normalizedCorrelation` (so it matches [getNormCorr()] exactly); each
#' permutation re-optimizes the residual leading axis on permuted PCs. The raw
#' p-value uses the Phipson & Smyth (2010) estimator, which counts the observed
#' configuration as one admissible permutation so the p-value is never exactly
#' zero:
#'
#' \deqn{p_{\mathrm{raw}}(k) = \frac{1 + \#\{\mathrm{perm}_k \ge \mathrm{obs}_k\}}{1 + m},}
#'
#' with Monte-Carlo floor `1 / (m + 1)` (reported as `mc_floor`). Closed
#' step-down control of the family-wise error rate across ordered axes uses
#'
#' \deqn{p_{\mathrm{stepdown}}(k) = \max_{j \le k} p_{\mathrm{raw}}(j),}
#'
#' and testing stops at the first axis with `p_stepdown > alpha`; that axis and
#' all later ones are declared non-significant. No Bonferroni factor is needed.
#' This is the closed/fixed-sequence test of Marcus, Peritz & Gabriel (1976) and
#' the permutation "test of canonical axes" of Legendre, Oksanen & ter Braak
#' (2011).
#'
#' ## Relationship to data residualization
#'
#' The fixed directions are removed using the full projection
#' `(I - u u^T) Y (I - v v^T)` for whitened PCs, or its weighted oblique
#' counterpart for unscaled PCs. Rank-one subtraction is deliberately avoided
#' because the observed directions are not generally singular vectors of a
#' permuted operator. The fair-sigma maximum over the bandwidth family is the
#' Westfall-Young (1993) maxT procedure.
#'
#' @param object A CoPro object with `runSkrCCA()` and
#'   `computeNormalizedCorrelation()` already run.
#' @param nPermu Number of permutations (default 999).
#' @param sigma_values Candidate bandwidths for the fair-sigma maximum. `NULL`
#'   uses all `object@sigmaValues` that have kernels and weights.
#' @param permu_method Permutation null: "bin" (default), "global", "pc", or
#'   "toroidal". See [runSkrCCAPermu()].
#' @param permu_which Which cell types to permute: "second_only" (default),
#'   "both", or "first_only".
#' @param num_bins_x,num_bins_y Bin grid for `permu_method = "bin"`. Default
#'   `NULL` is sigma-aware (see [.sigmaAwareBins()]); the same grid is shared
#'   across the sigma sweep. Because it is sized from the observed selected
#'   bandwidth, the default grid is mildly data-adaptive; pass integers to
#'   remove that second-order circularity.
#' @param match_quantile Whether to use quantile matching for bin permutation.
#' @param alpha Family-wise significance level for the step-down rule (default 0.05).
#' @param maxIter,tol Optimization controls passed to the axis optimizer.
#' @param verbose Whether to print progress and a summary (default TRUE).
#' @param n_cores Number of PSOCK workers. Each worker holds the PCA and kernel
#'   inputs, so choose this with available memory in mind.
#' @inheritParams runSkrCCAPermu
#'
#' @return The CoPro object with results stored in the `@conditionalPermu` slot,
#'   a list whose `per_axis` element is a data frame of `CC_index`,
#'   `observed_stat`, `observed_sigma`, `p_raw`, `p_stepdown`, `mc_floor`, and
#'   `significant`, plus the full null matrices for diagnostics. Use
#'   [calculate_pvalue_stepdown()] to read the per-axis table.
#'
#' @references
#' Phipson B, Smyth GK (2010). Permutation P-values should never be zero.
#' \emph{Stat Appl Genet Mol Biol} 9:Article39.
#' Legendre P, Oksanen J, ter Braak CJF (2011). Testing the significance of
#' canonical axes in redundancy analysis. \emph{Methods Ecol Evol} 2:269-277.
#' Westfall PH, Young SS (1993). \emph{Resampling-based Multiple Testing}.
#'
#' @seealso [runSkrCCAPermu_FairSigma()] (the CC1 / sigma-multiplicity case),
#'   [calculate_pvalue_stepdown()]
#'
#' @examples
#' \dontrun{
#' br <- runSkrCCA(br, scalePCs = TRUE, nCC = 3)
#' br <- computeNormalizedCorrelation(br)
#' br <- runSkrCCAPermu_Conditional(br, nPermu = 200, permu_method = "bin")
#' calculate_pvalue_stepdown(br)
#' }
#'
#' @importFrom stats median
#' @export
runSkrCCAPermu_Conditional <- function(object,
                                       nPermu = 999,
                                       sigma_values = NULL,
                                       permu_method = "bin",
                                       permu_which = "second_only",
                                       num_bins_x = NULL,
                                       num_bins_y = NULL,
                                       match_quantile = FALSE,
                                       alpha = 0.05,
                                       maxIter = 200,
                                       tol = 1e-5,
                                       verbose = TRUE,
                                       n_cores = 1,
                                       factorize =
                                         .defaultFactorizePermutation(),
                                       compactPermutation =
                                         .defaultCompactPermutation()) {

  ## ---- validation ----
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }
  .rejectCellPermutationForMulti(object)
  if (length(object@skrCCAOut) == 0) {
    stop("Please run runSkrCCA() first")
  }
  if (length(object@normalizedCorrelation) == 0) {
    stop("Please run computeNormalizedCorrelation() first")
  }
  nPermu <- .checkPositiveScalarInt(nPermu, "nPermu")
  if (nPermu < 2L) {
    stop("nPermu must be at least 2; one draw produces a degenerate null.")
  }
  if (!(permu_method %in% c("bin", "global", "pc", "toroidal"))) {
    stop("permu_method must be 'bin', 'global', 'pc', or 'toroidal'.")
  }
  if (!(permu_which %in% c("second_only", "both", "first_only"))) {
    stop("permu_which must be 'second_only', 'both', or 'first_only'.")
  }

  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) {
    stop("cellTypesOfInterest is empty.")
  }
  if (length(cts) > 2L) {
    stop("runSkrCCAPermu_Conditional currently supports one predeclared cell-type ",
         "pair. Subset to two cell types; higher-axis and multi-pair selection ",
         "must not be combined implicitly.")
  }
  scalePCs <- object@scalePCs
  nCC <- object@nCC
  sdev2_list <- .permutationSdev2(object, cts)
  ## Restricted to at most two cell types above, where sumcor and sumcov are
  ## the same problem. Verified rather than assumed.
  .resolvePermutationObjective(
    object, cts, scalePCs, supports_sumcor = FALSE, verbose = verbose
  )
  if (length(nCC) == 0 || nCC < 1) {
    stop("nCC must be >= 1; run runSkrCCA() with the desired number of axes.")
  }
  if (length(cts) == 1 && verbose) {
    message("Single cell type: permu_which is ignored and the one type is ",
            "permuted, giving a within-type null that relabels scores against ",
            "their own locations.")
  }

  ## ---- candidate sigma values (fair-sigma sweep) ----
  if (is.null(sigma_values)) {
    sigma_values <- object@sigmaValues
  }
  sigma_values <- sigma_values[sigma_values %in% object@sigmaValues]
  if (length(sigma_values) == 0) {
    stop("No valid sigma values available for the conditional permutation test.")
  }
  sigma_names <- .sigmaName(sigma_values)

  ## Keep only sigma values that have observed skrCCA weights.
  obs_W <- object@skrCCAOut[sigma_names]
  keep <- !vapply(obs_W, is.null, logical(1))
  if (!any(keep)) {
    stop("No skrCCA weights found for the requested sigma values.")
  }
  sigma_values <- sigma_values[keep]
  sigma_names <- sigma_names[keep]
  obs_W <- obs_W[keep]

  ## ---- resolve sigma-aware bins (shared across the sigma sweep and axes) ----
  if (permu_method == "bin" && (is.null(num_bins_x) || is.null(num_bins_y))) {
    bins <- .sigmaAwareBins(object, sigma = object@sigmaValueChoice,
                            verbose = verbose)
    num_bins_x <- bins$num_bins_x
    num_bins_y <- bins$num_bins_y
  }

  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

  if (verbose) {
    cat("Running CONDITIONAL (step-down) permutation test\n")
    cat("================================================\n")
    cat(sprintf("Axes (nCC): %d | sigma values: %d | permutations: %d\n",
                nCC, length(sigma_values), nPermu))
    cat(sprintf("permu_method: %s | permu_which: %s | alpha: %g\n\n",
                permu_method, permu_which, alpha))
  }

  ## ---- observed fair-sigma statistic per axis ----
  ## Read from the stored normalized correlation so it matches getNormCorr()
  ## exactly; the observed CC_k was produced by deflating the observed
  ## CC1..CC(k-1), which is precisely the conditional statistic.
  ncorr_obs <- object@normalizedCorrelation
  obs_stat <- rep(-Inf, nCC)
  obs_sigma <- rep(sigma_values[1], nCC)
  for (si in seq_along(sigma_values)) {
    df <- ncorr_obs[[sigma_names[si]]]
    if (is.null(df)) next
    for (k in seq_len(nCC)) {
      vk <- df$normalizedCorrelation[df$CC_index == k]
      if (length(vk) == 0) next
      vk <- max(vk, na.rm = TRUE)        # over pairs (single pair in manuscript)
      if (is.finite(vk) && vk > obs_stat[k]) {
        obs_stat[k] <- vk
        obs_sigma[k] <- sigma_values[si]
      }
    }
  }

  ## ---- generate permutations once (shared across axes & sigma) ----
  cell_permu <- .getCellPermu(object = object, permu_method = permu_method,
                              nPermu = nPermu, cts = cts,
                              permu_which = permu_which,
                              num_bins_x = num_bins_x, num_bins_y = num_bins_y,
                              match_quantile = match_quantile,
                              compactPermutation = compactPermutation)

  perm_stat <- matrix(NA_real_, nrow = nPermu, ncol = nCC)
  perm_sigma <- matrix(sigma_values[1], nrow = nPermu, ncol = nCC)
  n_failed <- 0L
  failure_messages <- character()

  # These values are kernel/sigma-specific, not permutation/axis-specific.
  # Precomputing them removes nPermu * nCC repeated normalizer products.
  kernel_info <- lapply(
    sigma_values,
    function(sigma) .get_ncorr_kernel_info(
      object@kernelMatrices, sigma, cts,
      normalizer_cache = attr(object, "kernelNormalizerCache", exact = TRUE)
    )
  )

  if (verbose) cat("Running permutations...\n")

  # Fixed cell types let the kernel be applied to their PC matrix once per
  # sigma rather than once per (draw, sigma). See R/D0_permutation_plan.R.
  fixed <- .fixedPermutationTypes(cell_permu, cts)
  plans <- lapply(sigma_values, function(sigma) .buildYPlan(
    PCmats = PCmats, flat_kernels = object@kernelMatrices, sigma = sigma,
    cts = cts, fixed = fixed, factorize = factorize
  ))
  grams <- .permutationGrams(PCmats, cell_permu, cts, factorize = factorize)

  worker <- .makeConditionalWorker(
    PCmats = PCmats, plans = plans, cts = cts, nCC = nCC,
    sigma_values = sigma_values, sigma_names = sigma_names, obs_W = obs_W,
    sdev2_list = sdev2_list, kernel_info = lapply(kernel_info, .slimKernelInfo),
    grams = grams, maxIter = maxIter, tol = tol
  )

  permutation_results <- .runPermutationDraws(
    cell_permu = cell_permu, cts = cts, nPermu = nPermu, worker = worker,
    n_cores = n_cores, verbose = verbose
  )
  for (tt in seq_len(nPermu)) {
    res <- permutation_results[[tt]]
    if (is.null(res) || !is.null(res$error)) {
      n_failed <- n_failed + 1L
      failure_messages <- c(
        failure_messages,
        if (is.null(res)) "worker returned NULL" else res$error
      )
    } else {
      perm_stat[tt, ] <- res$stat
      perm_sigma[tt, ] <- res$sigma
    }
  }
  if (verbose) cat(sprintf("  Completed %d permutations\n", nPermu))

  if (n_failed > 0) {
    error_counts <- sort(table(failure_messages), decreasing = TRUE)
    warning(sprintf(paste0("%d of %d permutations failed to optimize and were ",
                          "dropped from the null. P-values use the remaining ",
                          "permutations. Error counts: %s."), n_failed, nPermu,
                    paste(paste0(names(error_counts), " (",
                                 as.integer(error_counts), ")"),
                          collapse = "; ")))
  } else {
    error_counts <- integer()
  }

  ## ---- step-down p-values ----
  ## Phipson & Smyth (2010): the observed configuration is itself one admissible
  ## permutation (the identity), so a valid Monte-Carlo p-value adds 1 to both
  ## the numerator and the denominator. This guarantees p > 0 (a p-value of
  ## exactly 0 is invalid) and the smallest resolvable value is the Monte-Carlo
  ## floor 1 / (n_eff + 1).
  p_raw <- numeric(nCC)
  mc_floor <- numeric(nCC)
  for (k in seq_len(nCC)) {
    col_k <- perm_stat[!is.na(perm_stat[, k]), k]
    n_eff <- length(col_k)
    if (n_eff == 0) {
      p_raw[k] <- NA_real_
      mc_floor[k] <- NA_real_
    } else {
      p_raw[k] <- (1 + sum(col_k >= obs_stat[k])) / (1 + n_eff)
      mc_floor[k] <- 1 / (1 + n_eff)
    }
  }
  p_stepdown <- cummax(ifelse(is.na(p_raw), 1, p_raw))  # closed step-down

  ## stop at first non-significant axis
  first_ns <- which(p_stepdown > alpha)
  n_sig <- if (length(first_ns) == 0) nCC else (first_ns[1] - 1L)
  significant <- seq_len(nCC) <= n_sig

  per_axis <- data.frame(
    CC_index = seq_len(nCC),
    observed_stat = obs_stat,
    observed_sigma = obs_sigma,
    p_raw = p_raw,
    p_stepdown = p_stepdown,
    mc_floor = mc_floor,
    significant = significant,
    stringsAsFactors = FALSE
  )

  num_bins_out <- c(
    x = if (is.null(num_bins_x)) NA_integer_ else as.integer(num_bins_x),
    y = if (is.null(num_bins_y)) NA_integer_ else as.integer(num_bins_y)
  )

  object@conditionalPermu <- list(
    per_axis = per_axis,
    n_significant_axes = n_sig,
    alpha = alpha,
    mc_floor = mc_floor,
    obs_stats = obs_stat,
    obs_sigma = obs_sigma,
    perm_stats = perm_stat,
    perm_sigma = perm_sigma,
    sigma_values = sigma_values,
    nPermu = as.integer(nPermu),
    n_failed = n_failed,
    error_counts = error_counts,
    permu_method = permu_method,
    permu_which = permu_which,
    num_bins = num_bins_out,
    scalePCs = scalePCs,
    deflation = "projection",
    statistic = "max_over_sigma"
  )
  attr(object, "permutationProvenance") <- list(
    method = "conditional_stepdown",
    sigma_values = as.numeric(sigma_values),
    sigma_aggregation = "max",
    pair_aggregation = "max",
    sigma_predeclared = FALSE,
    selection_adjusted = TRUE,
    scalePCs = scalePCs
  )
  object@nPermu <- as.integer(nPermu)

  if (verbose) {
    cat("\n=== Conditional step-down test complete ===\n")
    for (k in seq_len(nCC)) {
      cat(sprintf(paste0("  CC%d: obs = %.4f (sigma = %g)  p_raw = %.4f  ",
                        "p_stepdown = %.4f  %s\n"),
                  k, obs_stat[k], obs_sigma[k], p_raw[k], p_stepdown[k],
                  if (significant[k]) "significant" else "not significant"))
    }
    cat(sprintf("  -> %d significant canonical axis/axes at alpha = %g\n",
                n_sig, alpha))
    cat(sprintf("  (Phipson-Smyth p-values; Monte-Carlo floor = %.4g with %d permutations)\n",
                1 / (nPermu + 1), nPermu))
  }

  return(object)
}


#' Read the step-down per-axis p-value table
#'
#' Thin reader for the conditional step-down permutation test produced by
#' [runSkrCCAPermu_Conditional()].
#'
#' @param object A CoPro object with `@conditionalPermu` populated.
#'
#' @return A data frame with one row per canonical axis: `CC_index`,
#'   `observed_stat`, `observed_sigma`, `p_raw`, `p_stepdown`, and
#'   `significant`. The number of significant axes, `alpha`, and `nPermu` are
#'   attached as attributes.
#'
#' @seealso [runSkrCCAPermu_Conditional()]
#'
#' @examples
#' \dontrun{
#' br <- runSkrCCAPermu_Conditional(br, nPermu = 200)
#' calculate_pvalue_stepdown(br)
#' }
#'
#' @export
calculate_pvalue_stepdown <- function(object) {
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }
  if (length(object@conditionalPermu) == 0) {
    stop("Run runSkrCCAPermu_Conditional() first.")
  }
  cp <- object@conditionalPermu
  out <- cp$per_axis
  attr(out, "n_significant_axes") <- cp$n_significant_axes
  attr(out, "alpha") <- cp$alpha
  attr(out, "nPermu") <- cp$nPermu
  out
}
