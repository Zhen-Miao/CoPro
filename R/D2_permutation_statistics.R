

#' Compute Normalized Correlation for Permutation Results
#'
#' Calculates the normalized correlation for each permutation, enabling
#' p-value calculation by comparing observed values to the null distribution.
#'
#' @details
#' The normalized correlation for each permutation is calculated using the
#' same formula as the observed data:
#'
#' \deqn{NC = \frac{s_1^T K_{12} s_2}{||s_1|| \cdot ||s_2|| \cdot ||\tilde K_c||_F}}
#'
#' where \eqn{s_1} and \eqn{s_2} are cell scores, \eqn{K_{12}} is the kernel
#' matrix, and \eqn{||\tilde K_c||_F} is the whitened-Frobenius norm
#' \eqn{||R_x^{1/2} K_c R_y^{1/2}||_F}.
#'
#' @param object A `CoPro` object with permutation results from `runSkrCCAPermu()`
#' @param tol Tolerance for approximate SVD calculation (default: 1e-4)
#' @param verbose Whether to report progress.
#' @inheritParams runSkrCCAPermu
#'
#' @return The `CoPro` object with permutation normalized correlations
#'   stored in `@normalizedCorrelationPermu`
#'
#' @examples
#' \dontrun{
#' # After running permutation testing
#' br <- computeNormalizedCorrelationPermu(br)
#'
#' # Extract permutation values and calculate p-value
#' permu_values <- sapply(br@normalizedCorrelationPermu,
#'                        function(x) x$normalizedCorrelation[1])
#' observed <- max(getNormCorr(br)$normalizedCorrelation)
#'
#' # One-sided p-value (testing if observed > permutation)
#' p_value <- (1 + sum(permu_values >= observed)) /
#'   (1 + length(permu_values))
#' }
#'
#' @export
computeNormalizedCorrelationPermu <- function(
    object, tol = 1e-4,
    factorize = .defaultFactorizePermutation(), verbose = TRUE) {

  ## Input validation
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }
  .rejectCellPermutationForMulti(object)

  if (length(object@skrCCAPermuOut) == 0) {
    stop(paste("skrCCAPermuOut is not available.",
               "Please run runSkrCCAPermu() first."))
  }

  ## Get parameters
  cts <- if (length(object@cellTypesOfInterest) > 0L) {
    object@cellTypesOfInterest
  } else {
    unique(object@cellTypesSub)
  }
  nPermu <- object@nPermu

  if (length(object@scalePCs) == 0) {
    stop("object@scalePCs not specified")
  }
  scalePCs <- object@scalePCs

  provenance <- .permutationProvenance(object)
  sigmaValueChoice <- provenance$sigma_values
  if (length(sigmaValueChoice) != 1L ||
      !identical(provenance$sigma_aggregation, "fixed")) {
    stop("computeNormalizedCorrelationPermu() requires fixed-sigma permutation ",
         "output. Fair-sigma and conditional methods already store their ",
         "normalized null statistics.")
  }
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
  nCC <- object@nCC

  pair_cell_types <- .permutationPairTypes(cts)

  ## Initialize output
  correlation_value <- vector("list", length = nPermu)
  permu_names <- paste("permu", 1:nPermu, sep = "_")
  names(correlation_value) <- permu_names
  s_name <- .sigmaName(sigmaValueChoice)

  ## Calculate whitened-Frobenius normalizers (only need to do this once)
  if (verbose) cat("Calculating whitened-Frobenius normalizers...\n")
  norm_K12 <- setNames(vector(mode = "list", length = 1), s_name)
  norm_K12[[s_name]] <- setNames(vector(mode = "list", length = length(cts)), cts)
  normalizer_cache <- attr(object, "kernelNormalizerCache", exact = TRUE)

  for (i in cts) {
    norm_K12[[s_name]][[i]] <- setNames(
      vector(mode = "list", length = length(cts)), cts)
  }

  for (pp in seq_len(ncol(pair_cell_types))) {
    cellType1 <- pair_cell_types[1, pp]
    cellType2 <- pair_cell_types[2, pp]
    K <- getKernelMatrix(object, sigma = sigmaValueChoice,
                         cellType1 = cellType1, cellType2 = cellType2,
                         verbose = FALSE, materialize = FALSE)
    ## matched-sigma within-type kernels as whitening operators
    Rx <- tryCatch(getKernelMatrix(object, sigma = sigmaValueChoice,
                     cellType1 = cellType1, cellType2 = cellType1,
                     verbose = FALSE, materialize = FALSE),
                   error = function(e) NULL)
    Ry <- tryCatch(getKernelMatrix(object, sigma = sigmaValueChoice,
                     cellType1 = cellType2, cellType2 = cellType2,
                     verbose = FALSE, materialize = FALSE),
                   error = function(e) NULL)
    cache_key <- .kernelNormalizerKey(
      sigmaValueChoice, cellType1, cellType2, slide = NULL
    )
    nrm <- .readKernelNormalizer(normalizer_cache, cache_key, K, Rx, Ry)
    if (is.null(nrm)) {
      nrm <- .whitenedFrobNorm(K, Rx, Ry)
      normalizer_cache <- .cacheKernelNormalizer(
        normalizer_cache, cache_key, K, Rx, Ry, nrm
      )
    }
    norm_K12[[s_name]][[cellType1]][[cellType2]] <- nrm
  }
  if (verbose) cat("Whitened-Frobenius normalizers calculated.\n\n")

  ## The null scores use the same permutations and the same kernel as the fits
  ## that produced them, so the same factorization applies: everything below is
  ## evaluated from the small PC-space operator Y and the per-type Gram
  ## matrices. See R/D0_permutation_plan.R.
  plan <- .buildYPlan(
    PCmats = PCmats,
    flat_kernels = object@kernelMatrices,
    sigma = sigmaValueChoice,
    cts = cts,
    fixed = .fixedPermutationTypes(object@cellPermu, cts),
    factorize = factorize
  )
  grams <- .permutationGrams(PCmats, object@cellPermu, cts,
                             factorize = factorize)

  ## Calculate normalized correlation for each permutation
  if (verbose) cat("Computing normalized correlations for permutations...\n")

  for (tt in seq_len(nPermu)) {
    t <- permu_names[tt]
    correlation_value[[t]] <- data.frame(
      sigmaValues = sigmaValueChoice,
      cellType1 = rep(pair_cell_types[1, ], times = nCC),
      cellType2 = rep(pair_cell_types[2, ], times = nCC),
      CC_index = rep(x = 1:nCC, each = ncol(pair_cell_types)),
      normalizedCorrelation = numeric(length = ncol(pair_cell_types) * nCC),
      stringsAsFactors = FALSE
    )

    PCmats_permuted <- .applyPermutationSpec(
      PCmats, .permutationDrawSpec(object@cellPermu, cts, tt)
    )
    Y_resi <- .yResiFromPlan(plan, PCmats_permuted)

    W_permu <- stats::setNames(lapply(cts, function(ct) {
      object@skrCCAPermuOut[[t]][[ct]][, seq_len(nCC), drop = FALSE]
    }), cts)

    ## ||X w||^2 = w' (X' X) w, and a row permutation leaves X' X unchanged.
    ## A "pc" permutation does not, so those types keep the direct product.
    score_norms <- stats::setNames(lapply(cts, function(ct) {
      W <- W_permu[[ct]]
      G <- grams[[ct]]
      if (is.null(G)) {
        sqrt(colSums((PCmats_permuted[[ct]] %*% W)^2))
      } else {
        sqrt(colSums(W * (G %*% W)))
      }
    }), cts)

    for (pp in seq_len(ncol(pair_cell_types))) {
      cellType1 <- pair_cell_types[1, pp]
      cellType2 <- pair_cell_types[2, pp]
      norm_K12_sel <- norm_K12[[s_name]][[cellType1]][[cellType2]]

      ## s_1' K s_2 for every component at once = diag(W_1' Y W_2).
      numerators <- colSums(
        W_permu[[cellType1]] *
          (Y_resi[[cellType1]][[cellType2]] %*% W_permu[[cellType2]])
      )
      denominators <- score_norms[[cellType1]] * score_norms[[cellType2]] *
        norm_K12_sel
      out_idx <- pp + (seq_len(nCC) - 1L) * ncol(pair_cell_types)
      correlation_value[[t]]$normalizedCorrelation[out_idx] <-
        as.numeric(numerators / denominators)
    }

    # Progress indicator
    if (verbose && (tt %% 20 == 0 || tt == nPermu)) {
      cat(paste("  Completed", tt, "of", nPermu, "permutations\n"))
    }
  }

  ## Store results. Bind the provenance describing this null onto the null
  ## itself, so a later runSkrCCAPermu_Conditional()/_FairSigma() call that
  ## overwrites the shared object-level provenance attribute cannot silently
  ## re-label it when calculate_pvalue() reads it back.
  object@normalizedCorrelationPermu <- `attr<-`(
    correlation_value, "provenance", .permutationProvenance(object)
  )
  attr(object, "kernelNormalizerCache") <- normalizer_cache

  if (verbose) cat("\nNormalized correlation computation complete.\n")

  return(object)
}


#' Compute Normalized Correlation from Ground Truth Scores
#'
#' Computes the normalized correlation using user-provided cell scores
#' (e.g., ground truth scores from simulation) and the pre-calculated
#' kernel matrix from a CoPro object.
#'
#' @details
#' This function is useful for comparing:
#' - Ground truth normalized correlation (from simulated scores)
#' - CoPro's optimized normalized correlation (from CCA)
#' - Permutation distribution
#'
#' Under null simulation (no cross-type coordination), the ground truth
#' normalized correlation should be close to 0, while CoPro's optimized
#' value will be higher because it searches for the maximum correlation.
#'
#' @param object A CoPro object with kernel matrices computed
#' @param scores_ct1 Numeric vector of scores for cell type 1.
#'   Must be in the same order as cells in the CoPro object.
#' @param scores_ct2 Numeric vector of scores for cell type 2.
#'   Must be in the same order as cells in the CoPro object.
#' @param cellType1 Name of cell type 1
#' @param cellType2 Name of cell type 2
#' @param sigma Sigma value for kernel matrix (default: uses sigmaValueChoice)
#' @param tol Tolerance for SVD calculation (default: 1e-4)
#'
#' @return Normalized correlation value
#'
#' @examples
#' \dontrun{
#' # After running CoPro analysis on simulated data
#' # Get ground truth scores from metadata
#' meta <- br@metaDataSub
#' gt_scores_A <- meta$smoothed_score[meta$cell_type == "A"]
#' gt_scores_B <- meta$smoothed_score[meta$cell_type == "B"]
#'
#' # Compute ground truth normalized correlation
#' gt_ncorr <- compute_ground_truth_ncorr(
#'   object = br,
#'   scores_ct1 = gt_scores_A,
#'   scores_ct2 = gt_scores_B,
#'   cellType1 = "A",
#'   cellType2 = "B"
#' )
#' }
#'
#' @export
compute_ground_truth_ncorr <- function(object,
                                       scores_ct1,
                                       scores_ct2,
                                       cellType1,
                                       cellType2,
                                       sigma = NULL,
                                       tol = 1e-4) {

  ## Input validation
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  if (length(object@kernelMatrices) == 0) {
    stop("Kernel matrices not computed. Run computeKernelMatrix() first.")
  }

  ## Get sigma value
  if (is.null(sigma)) {
    sigma <- object@sigmaValueChoice
    if (is.null(sigma) || length(sigma) == 0) {
      stop("sigma not specified and sigmaValueChoice not set in object")
    }
  }

  if (!(sigma %in% object@sigmaValues)) {
    stop(paste("sigma", sigma, "not found in object@sigmaValues"))
  }

  ## Get cell type indices
  cts <- object@cellTypesOfInterest
  if (!(cellType1 %in% cts) || !(cellType2 %in% cts)) {
    stop("cellType1 and cellType2 must be in cellTypesOfInterest")
  }

  ## Validate score lengths
  n_ct1 <- sum(object@cellTypesSub == cellType1)
  n_ct2 <- sum(object@cellTypesSub == cellType2)

  if (length(scores_ct1) != n_ct1) {
    stop(paste("scores_ct1 length", length(scores_ct1),
               "does not match number of", cellType1, "cells:", n_ct1))
  }
  if (length(scores_ct2) != n_ct2) {
    stop(paste("scores_ct2 length", length(scores_ct2),
               "does not match number of", cellType2, "cells:", n_ct2))
  }

  ## Get kernel matrix
  K <- getKernelMatrix(object, sigma = sigma,
                       cellType1 = cellType1, cellType2 = cellType2,
                       verbose = FALSE, materialize = FALSE)

  ## Normalize scores (subtract mean, divide by norm)
  normalize_vec <- function(x) {
    x <- x - mean(x)
    norm_x <- sqrt(sum(x^2))
    if (norm_x > 1e-10) {
      x <- x / norm_x
    }
    return(x)
  }

  s1 <- normalize_vec(scores_ct1)
  s2 <- normalize_vec(scores_ct2)

  ## Whitened-Frobenius normalizer (matched-sigma within-type kernels)
  Rx <- tryCatch(getKernelMatrix(object, sigma = sigma,
                   cellType1 = cellType1, cellType2 = cellType1,
                   verbose = FALSE, materialize = FALSE),
                 error = function(e) NULL)
  Ry <- tryCatch(getKernelMatrix(object, sigma = sigma,
                   cellType1 = cellType2, cellType2 = cellType2,
                   verbose = FALSE, materialize = FALSE),
                 error = function(e) NULL)
  norm_K12 <- .whitenedFrobNorm(K, Rx, Ry)

  ## Calculate normalized correlation
  norm_corr <- .kernelXKY(s1, K, s2)[1, 1] / norm_K12

  return(norm_corr)
}


#' Calculate P-value from Permutation Results
#'
#' Helper function to calculate p-value from permutation testing results.
#'
#' @param object A CoPro object with permutation results
#' @param cc_index Which canonical correlation component to use (default: 1)
#' @param alternative Direction of test: `"greater"` (default) or `"less"`.
#'   A two-sided test is not defined here because the null distribution of the
#'   optimized, max-aggregated normalized-correlation statistic is not symmetric
#'   about zero; use `"greater"` for evidence of co-progression.
#'
#' @return List with the Phipson & Smyth (2010) permutation p-value (`p_value`,
#'   never exactly zero), the Monte-Carlo floor `mc_floor = 1 / (n_permu + 1)`,
#'   the observed value, and the permutation distribution. With more than one
#'   cell-type pair, the observed and each permutation statistic are both the
#'   maximum normalized correlation over pairs for the requested axis.
#' @references Phipson B, Smyth GK (2010). Permutation P-values should never be
#'   zero. \emph{Stat Appl Genet Mol Biol} 9:Article39.
#'
#' @examples
#' \dontrun{
#' result <- calculate_pvalue(br, cc_index = 1, alternative = "greater")
#' print(paste("P-value:", result$p_value))
#' }
#'
#' @export
calculate_pvalue <- function(object, cc_index = 1, alternative = "greater") {

  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }
  .rejectCellPermutationForMulti(object)
  if (identical(alternative, "two.sided")) {
    stop(
      'alternative = "two.sided" is not defined because the optimized, ',
      "max-aggregated normalized-correlation null is not symmetric about zero; use ",
      'alternative = "greater" for evidence of co-progression.'
    )
  }
  alternative <- match.arg(alternative, c("greater", "less"))

  if (length(object@normalizedCorrelationPermu) == 0) {
    stop("Run computeNormalizedCorrelationPermu() first")
  }

  # Prefer the provenance bound to the null itself; the shared object-level
  # attribute can be overwritten by a later runSkrCCAPermu_Conditional() call,
  # which would otherwise re-label this null and change the reported p-value.
  provenance <- attr(object@normalizedCorrelationPermu, "provenance",
                     exact = TRUE)
  if (is.null(provenance)) provenance <- .permutationProvenance(object)

  # Get the observed value for exactly the sigma family represented by the
  # null. This prevents a fixed-sigma null from being compared with an
  # observed max over all fitted sigmas.
  ncorr <- getNormCorr(object)
  ncorr_cc <- ncorr[ncorr$CC_index == cc_index, ]
  sigma_col <- intersect(c("sigmaValues", "sigmaValue"), names(ncorr_cc))
  if (length(provenance$sigma_values) > 0L && length(sigma_col) > 0L) {
    ncorr_cc <- ncorr_cc[
      ncorr_cc[[sigma_col[1L]]] %in% provenance$sigma_values, , drop = FALSE
    ]
  }
  if (nrow(ncorr_cc) == 0L ||
      any(!is.finite(ncorr_cc$normalizedCorrelation))) {
    stop("Observed normalized correlation is missing or non-finite for cc_index.")
  }
  observed <- max(ncorr_cc$normalizedCorrelation)

  # Use the same max-over-pairs aggregation as the observed statistic. The old
  # path selected the first null pair while maximizing the observed rows, which
  # made the two sides different statistics for three or more cell types.
  permu_values <- vapply(object@normalizedCorrelationPermu, function(x) {
    values <- x$normalizedCorrelation[x$CC_index == cc_index]
    if (length(values) == 0L || any(!is.finite(values))) {
      stop("A permutation normalized correlation is missing or non-finite for cc_index.")
    }
    max(values)
  }, numeric(1))

  # Calculate p-value using the Phipson & Smyth (2010) estimator: the observed
  # configuration is itself one admissible permutation, so 1 is added to both
  # numerator and denominator. This is always strictly positive (a permutation
  # p-value of exactly 0 is invalid); the smallest resolvable value is the
  # Monte-Carlo floor 1 / (m + 1).
  m <- length(permu_values)
  if (alternative == "greater") {
    p_value <- (1 + sum(permu_values >= observed)) / (1 + m)
  } else {
    p_value <- (1 + sum(permu_values <= observed)) / (1 + m)
  }

  result <- list(
    p_value = p_value,
    mc_floor = 1 / (1 + m),
    observed = observed,
    permu_mean = mean(permu_values),
    permu_sd = sd(permu_values),
    permu_values = permu_values,
    n_permu = m,
    alternative = alternative,
    pair_aggregation = "max",
    sigma_aggregation = provenance$sigma_aggregation,
    sigma_values = provenance$sigma_values,
    sigma_predeclared = provenance$sigma_predeclared,
    selection_adjusted = if (is.null(provenance$selection_adjusted)) {
      NA
    } else {
      provenance$selection_adjusted
    }
  )

  if (identical(result$selection_adjusted, FALSE)) {
    warning(
      "This is a conditional fixed-sigma p-value at a sigma selected from the ",
      "same data; it is not adjusted for sigma selection. Use ",
      "runSkrCCAPermu_FairSigma() or pass a predeclared sigma to ",
      "runSkrCCAPermu()."
    )
  }

  return(result)
}


#' Quick Normalized Correlation Computation (Internal)
#'
#' Helper function to quickly compute normalized correlation for a given
#' set of PC matrices and weights.
#'
#' @param PCmats Named list of PC matrices
#' @param w_list Named list of weight vectors (matrices with 1 column)
#' @param flat_kernels Flat list of kernel matrices from CoPro object
#' @param sigma Sigma value for kernel selection
#' @param cts Cell types of interest
#' @param tol Tolerance for SVD computation
#' @param kernel_info Optional precomputed list containing `K` and `norm_K12`
#'   for this sigma.
#' @param Y_resi Optional precomputed PC-space operator from
#'   `compute_Y_resi()`. When supplied, it is used for the numerator so the
#'   kernel-vector product is not repeated.
#' @param grams Optional named list of per-cell-type Gram matrices
#'   `crossprod(X)` used for the score norms in the denominator. Valid whenever
#'   the permutation is a bijection on cells; see `.permutationGrams()`.
#'
#' @return Numeric value of normalized correlation
#' @keywords internal
.compute_ncorr_quick <- function(PCmats, w_list, flat_kernels, sigma, cts,
                                 tol = 1e-4, kernel_info = NULL,
                                 Y_resi = NULL, grams = NULL) {
  if (length(cts) > 2L) {
    stop("Quick permutation scoring supports one predeclared cell-type pair; ",
         "subset the CoPro object to one pair or use a symmetric multi-pair path.")
  }
  if (length(cts) == 1) {
    ct1 <- ct2 <- cts[1]
  } else {
    ct1 <- cts[1]
    ct2 <- cts[2]
  }

  A <- PCmats[[ct1]]
  B <- PCmats[[ct2]]
  w1 <- w_list[[ct1]][, 1, drop = FALSE]
  w2 <- w_list[[ct2]][, 1, drop = FALSE]

  if (is.null(kernel_info)) {
    kernel_info <- .get_ncorr_kernel_info(flat_kernels, sigma, cts)
  }
  K <- kernel_info$K
  norm_K12 <- kernel_info$norm_K12

  # Normalized correlation
  Y12 <- if (!is.null(Y_resi) && !is.null(Y_resi[[ct1]])) {
    Y_resi[[ct1]][[ct2]]
  } else {
    NULL
  }
  numerator <- if (!is.null(Y12)) {
    as.numeric(crossprod(w1, Y12 %*% w2))
  } else {
    .kernelXKY(A %*% w1, K, B %*% w2)[1, 1]
  }

  # ||X w||^2 = w' (X' X) w. A row permutation leaves X' X unchanged, so the
  # caller can supply the Gram matrices once and skip an n-by-nPC product per
  # draw. Types with no Gram (or no cache at all) use the direct form.
  score_norm <- function(ct, X, w) {
    G <- if (is.null(grams)) NULL else grams[[ct]]
    if (is.null(G)) return(sqrt(sum((X %*% w)^2)))
    sqrt(as.numeric(crossprod(w, G %*% w)))
  }
  denominator <- score_norm(ct1, A, w1) * score_norm(ct2, B, w2) * norm_K12

  return(numerator / denominator)
}


# Precompute permutation scoring objects that depend only on sigma. The
# whitened-Frobenius normalizer is especially expensive for large kernels and
# must not be recomputed for every permutation or canonical axis.
