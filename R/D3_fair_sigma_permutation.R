# =============================================================================
# Selection-fair bandwidth inference
# -----------------------------------------------------------------------------
# runSkrCCAPermu_FairSigma() re-optimizes every permutation draw across the
# full bandwidth grid and keeps each draw's best statistic (a Westfall-Young
# maxT over sigma), so the null faces the same bandwidth selection as the
# observed fit. Reuses the same draw within a sweep so the null maximum
# inherits the observed statistic's cross-scale correlation. See
# D0_permutation_plan.R for the architecture overview.
# =============================================================================

.get_ncorr_kernel_info <- function(flat_kernels, sigma, cts,
                                   normalizer_cache = NULL,
                                   slide = NULL) {
  if (length(cts) == 1) {
    ct1 <- ct2 <- cts[1]
  } else {
    ct1 <- cts[1]
    ct2 <- cts[2]
  }

  K <- get_kernel_matrix_flat(flat_kernels, sigma, ct1, ct2, slide = slide)
  Rx <- tryCatch(
    get_kernel_matrix_flat(flat_kernels, sigma, ct1, ct1, slide = slide),
    error = function(e) NULL
  )
  Ry <- tryCatch(
    get_kernel_matrix_flat(flat_kernels, sigma, ct2, ct2, slide = slide),
    error = function(e) NULL
  )

  cache_key <- .kernelNormalizerKey(sigma, ct1, ct2, slide = slide)
  norm_K12 <- .readKernelNormalizer(
    normalizer_cache, cache_key, K, Rx, Ry
  )
  if (is.null(norm_K12)) norm_K12 <- .whitenedFrobNorm(K, Rx, Ry)
  list(K = K, norm_K12 = norm_K12)
}


#' Run Permutation Test with Fair Sigma Selection
#'
#' Performs permutation testing where BOTH observed and permuted data
#' get to optimize sigma selection. This is the statistically correct
#' approach that addresses inflated Type I error caused by sigma selection
#' being applied only to observed data.
#'
#' @details
#' ## The Sigma Selection Problem
#'
#' In standard CoPro analysis, the observed data gets to choose the best
#' sigma (the one maximizing normalized correlation). However, permutation
#' data uses this SAME sigma, which may not be optimal for permuted data.
#' This asymmetry can inflate Type I error.
#'
#' ## The Solution
#'
#' This function runs CCA at EACH sigma value for EACH permutation, then
#' selects the best sigma for that permutation. Both observed and permuted
#' data thus have equal opportunity to optimize sigma selection.
#'
#' ## Computational Cost
#'
#' This is more computationally expensive (nPermu * nSigma CCA runs instead
#' of nPermu runs), but provides statistically correct p-values.
#'
#' @param object A CoPro object with CCA already computed via `runSkrCCA()`
#'   and normalized correlation computed via `computeNormalizedCorrelation()`
#' @param nPermu Number of permutations to run (default: 999)
#' @param sigma_values Vector of sigma values to test. If NULL, uses all
#'   sigma values from the original analysis (object@@sigmaValues)
#' @param permu_method Method of permutation: "bin", "global", "pc", or "toroidal"
#' @param permu_which Which cell types to permute: "second_only", "both",
#'   "first_only". Ignored with a single cell type, which is always permuted.
#' @param num_bins_x Number of bins in x for bin-wise permutation. Default `NULL`
#'   sizes the grid from the observed best bandwidth (`sigmaValueChoice`) via
#'   [.sigmaAwareBins()]; the same grid (and hence the same permutation) is
#'   shared across the sigma sweep. This makes the grid mildly data-adaptive;
#'   pass an integer to remove that second-order circularity.
#' @param num_bins_y Number of bins in y for bin-wise permutation. Default `NULL`
#'   (sigma-aware, as for `num_bins_x`).
#' @param match_quantile Whether to use quantile matching for bin permutation
#' @param maxIter Maximum iterations for CCA optimization
#' @param tol Convergence tolerance
#' @param n_cores Number of PSOCK workers. Each worker holds the PCA and kernel
#'   inputs, so choose this with available memory in mind.
#' @param verbose Whether to print progress messages
#' @inheritParams runSkrCCAPermu
#' @seealso [runSkrCCAPermu_Conditional()] for a sequential step-down test
#'   across canonical axes (the correct treatment when `nCC > 1`).
#'
#' @return CoPro object with fair permutation results stored in:
#'   \itemize{
#'     \item @@skrCCAPermuOut: Best weights for each permutation
#'     \item @@normalizedCorrelationPermu: Best ncorr for each permutation
#'     \item @@fairSigmaPermu: List with sigma selected for each permutation
#'   }
#'
#' @examples
#' \dontrun{
#' # After running standard CoPro analysis
#' br <- runSkrCCA(br, scalePCs = TRUE)
#' br <- computeNormalizedCorrelation(br)
#'
#' # Run fair sigma permutation test
#' br <- runSkrCCAPermu_FairSigma(br, nPermu = 100,
#'                                 permu_method = "toroidal")
#'
#' # Calculate p-value
#' result <- calculate_pvalue(br)
#' }
#'
#' @export
runSkrCCAPermu_FairSigma <- function(object,
                                     nPermu = 999,
                                     sigma_values = NULL,
                                     permu_method = "bin",
                                     permu_which = "second_only",
                                     num_bins_x = NULL,
                                     num_bins_y = NULL,
                                     match_quantile = FALSE,
                                     maxIter = 200,
                                     tol = 1e-5,
                                     n_cores = 1,
                                     verbose = TRUE,
                                     factorize = .defaultFactorizePermutation(),
                                     compactPermutation =
                                       .defaultCompactPermutation()) {

  ## Input validation
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

  # Use all sigma values from original analysis if not specified
  if (is.null(sigma_values)) {
    sigma_values <- object@sigmaValues
  }

  if (length(sigma_values) == 0) {
    stop("No sigma values specified and none found in object")
  }

  # Filter sigma values to only include those that exist in kernel matrices
 # (some sigma values may have been dropped during kernel computation)
  # Use object@sigmaValues which is already filtered during computeKernelMatrix
  available_sigmas <- object@sigmaValues

  if (length(available_sigmas) > 0) {
    original_sigmas <- sigma_values
    sigma_values <- sigma_values[sigma_values %in% available_sigmas]

    if (length(sigma_values) == 0) {
      stop("None of the specified sigma values exist in kernel matrices. ",
           "Available sigmas: ", paste(available_sigmas, collapse = ", "))
    }

    if (length(sigma_values) < length(original_sigmas) && verbose) {
      dropped <- setdiff(original_sigmas, sigma_values)
      cat(paste("Note: Sigma values", paste(dropped, collapse = ", "),
                "not found in kernel matrices and will be skipped.\n"))
    }
  }

  cts <- object@cellTypesOfInterest
  if (length(cts) > 2L) {
    stop("runSkrCCAPermu_FairSigma currently supports one predeclared cell-type ",
         "pair. Subset to two cell types and adjust the resulting pair-level ",
         "p-values across the declared family.")
  }
  if (length(cts) == 1L && verbose) {
    message("Single cell type: permu_which is ignored and the one type is ",
            "permuted, giving a within-type null that relabels scores against ",
            "their own locations.")
  }

  scalePCs <- object@scalePCs
  nCC <- object@nCC
  sdev2_list <- .permutationSdev2(object, cts)
  ## Restricted to at most two cell types above, where sumcor and sumcov are
  ## the same problem. Verified rather than assumed.
  .resolvePermutationObjective(
    object, cts, scalePCs, supports_sumcor = FALSE, verbose = verbose
  )

  if (nCC > 1) {
    warning("Fair sigma permutation tests only the first canonical axis (CC1). ",
            "For axes 2..nCC, use runSkrCCAPermu_Conditional(), which performs ",
            "a sequential step-down test that controls canonical-axis ",
            "multiplicity. Using CC1 here.")
  }

  ## Resolve the bin grid. When num_bins_x / num_bins_y are NULL (the default),
  ## size the patch grid from the observed best bandwidth (sigmaValueChoice).
  ## The SAME permutation is reused across the sigma sweep so that the
  ## max-over-sigma null mirrors the positive cross-scale correlation of the
  ## observed statistic; re-binning per sigma would inflate the null maximum
  ## and make the test conservative.
  if (permu_method == "bin" && (is.null(num_bins_x) || is.null(num_bins_y))) {
    bins <- .sigmaAwareBins(object, sigma = object@sigmaValueChoice,
                            verbose = verbose)
    num_bins_x <- bins$num_bins_x
    num_bins_y <- bins$num_bins_y
  }

  # Get PC matrices
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

  if (verbose) {
    cat("Running FAIR SIGMA permutation test\n")
    cat("===================================\n")
    cat(paste("Testing", length(sigma_values), "sigma values per permutation\n"))
    cat(paste("Total CCA runs:", nPermu * length(sigma_values), "\n"))
    cat(paste("permu_method:", permu_method, "\n"))
    cat(paste("permu_which:", permu_which, "\n\n"))
  }

  # Generate permutation indices
  cell_permu <- .getCellPermu(
    object = object,
    permu_method = permu_method,
    nPermu = nPermu,
    cts = cts,
    permu_which = permu_which,
    num_bins_x = num_bins_x,
    num_bins_y = num_bins_y,
    match_quantile = match_quantile,
    compactPermutation = compactPermutation
  )

  # Kernel matrices and their whitened-Frobenius normalizers are invariant to
  # PC permutations. Cache them once per sigma instead of repeating the
  # cubic normalizer calculation nPermu times.
  kernel_info <- lapply(
    sigma_values,
    function(sigma) .get_ncorr_kernel_info(
      object@kernelMatrices, sigma, cts,
      normalizer_cache = attr(object, "kernelNormalizerCache", exact = TRUE)
    )
  )

  # Fixed cell types let the kernel be applied to their PC matrix once per
  # sigma rather than once per draw. See R/D0_permutation_plan.R.
  fixed <- .fixedPermutationTypes(cell_permu, cts)
  plans <- lapply(sigma_values, function(sigma) .buildYPlan(
    PCmats = PCmats, flat_kernels = object@kernelMatrices, sigma = sigma,
    cts = cts, fixed = fixed, factorize = factorize
  ))
  grams <- .permutationGrams(PCmats, cell_permu, cts, factorize = factorize)

  if (verbose) {
    cat("Running permutations...\n")
  }

  worker <- .makeFairSigmaWorker(
    PCmats = PCmats, plans = plans, cts = cts, sigma_values = sigma_values,
    sdev2_list = sdev2_list, kernel_info = lapply(kernel_info, .slimKernelInfo),
    grams = grams, maxIter = maxIter, tol = tol
  )

  permu_results <- .runPermutationDraws(
    cell_permu = cell_permu, cts = cts, nPermu = nPermu, worker = worker,
    n_cores = n_cores, verbose = verbose
  )
  valid_draw <- vapply(permu_results, function(x) {
    is.finite(x$ncorr) && !is.null(x$weights)
  }, logical(1))
  error_messages <- unlist(lapply(permu_results, `[[`, "errors"),
                           use.names = FALSE)
  n_failed <- sum(!valid_draw)
  if (n_failed > 0L || length(error_messages) > 0L) {
    breakdown <- if (length(error_messages) == 0L) "" else {
      counts <- sort(table(error_messages), decreasing = TRUE)
      paste0(" Error counts: ", paste(
        paste0(names(counts), " (", as.integer(counts), ")"),
        collapse = "; "), ".")
    }
    summary <- if (n_failed > 0L) {
      paste0(n_failed, " of ", nPermu,
             " fair-sigma permutation draws had no finite fit and were dropped.")
    } else {
      paste0("Some sigma-level fits failed, but all ", nPermu,
             " fair-sigma permutation draws retained a finite fit.")
    }
    warning(summary, breakdown, call. = FALSE)
  }
  permu_results <- permu_results[valid_draw]
  if (length(permu_results) == 0L) {
    stop("Every fair-sigma permutation draw failed to produce a finite fit.")
  }
  nPermu_effective <- length(permu_results)
  permu_ncorrs <- vapply(permu_results, `[[`, numeric(1), "ncorr")
  permu_sigmas <- vapply(permu_results, `[[`, numeric(1), "sigma")
  if (verbose) cat("  Completed", nPermu, "permutations (",
                   nPermu_effective, " usable)\n", sep = "")

  # Store results in object
  object@skrCCAPermuOut <- lapply(permu_results, function(x) x$weights)
  names(object@skrCCAPermuOut) <- paste0("permu_", seq_len(nPermu_effective))

  # Create normalized correlation results structure
  pair_cell_types <- .permutationPairTypes(cts)
  correlation_value <- vector("list", nPermu_effective)
  for (tt in seq_len(nPermu_effective)) {
    correlation_value[[tt]] <- data.frame(
      sigmaValues = permu_results[[tt]]$sigma,
      cellType1 = pair_cell_types[1, 1],
      cellType2 = pair_cell_types[2, 1],
      CC_index = 1,
      normalizedCorrelation = permu_results[[tt]]$ncorr,
      stringsAsFactors = FALSE
    )
  }
  names(correlation_value) <- paste0("permu_", seq_len(nPermu_effective))
  object@normalizedCorrelationPermu <- correlation_value
  object@nPermu <- as.integer(nPermu_effective)

  # Get observed best sigma for comparison
  observed_best_sigma <- object@sigmaValueChoice

  # Calculate sigma difference statistics
  sigma_differs <- permu_sigmas != observed_best_sigma
  prop_sigma_differs <- mean(sigma_differs)
  n_sigma_differs <- sum(sigma_differs)

  # Store fair sigma info as attribute (for diagnostics)
  # Note: This is stored as an attribute since there's no dedicated slot
  attr(object, "fairSigmaPermu") <- list(
    sigma_selected = permu_sigmas,
    sigma_values_tested = sigma_values,
    observed_best_sigma = observed_best_sigma,
    sigma_differs = sigma_differs,
    prop_sigma_differs = prop_sigma_differs,
    n_sigma_differs = n_sigma_differs,
    n_requested = as.integer(nPermu),
    n_failed = as.integer(n_failed),
    error_counts = sort(table(error_messages), decreasing = TRUE)
  )
  attr(object, "permutationProvenance") <- list(
    method = "fair_sigma",
    sigma_values = as.numeric(sigma_values),
    sigma_aggregation = "max",
    pair_aggregation = "max",
    sigma_predeclared = FALSE,
    selection_adjusted = TRUE,
    scalePCs = scalePCs
  )
  # Bind this provenance to the null slot as well (see calculate_pvalue()).
  object@normalizedCorrelationPermu <- `attr<-`(
    object@normalizedCorrelationPermu, "provenance",
    attr(object, "permutationProvenance", exact = TRUE)
  )

  # Get observed value for comparison
  observed_values <- getNormCorr(object)$normalizedCorrelation
  observed_ncorr <- if (any(is.finite(observed_values))) {
    max(observed_values[is.finite(observed_values)])
  } else {
    NA_real_
  }

  if (verbose) {
    cat("\n=== Fair Sigma Permutation Complete ===\n")
    cat(paste("Observed best ncorr:", round(observed_ncorr, 4), "\n"))
    cat(paste("Permutation mean:", round(mean(permu_ncorrs), 4), "\n"))
    cat(paste("Permutation SD:", round(sd(permu_ncorrs), 4), "\n"))
    cat(paste("P-value (Phipson-Smyth):",
              round((1 + sum(permu_ncorrs >= observed_ncorr)) /
                      (1 + nPermu_effective), 4),
              paste0("(MC floor ", round(1 / (1 + nPermu_effective), 4), ")"), "\n"))
    cat("\nSigma selection in permutations:\n")
    cat(paste("  Observed best sigma:", observed_best_sigma, "\n"))
    cat(paste("  Sigma values used:", paste(unique(permu_sigmas), collapse = ", "), "\n"))
    cat(paste("  Most common:", names(sort(table(permu_sigmas), decreasing = TRUE))[1], "\n"))
    cat(paste("  Sigma differs from observed:", n_sigma_differs, "/", nPermu_effective,
              "(", round(prop_sigma_differs * 100, 1), "%)\n"))
  }

  return(object)
}
