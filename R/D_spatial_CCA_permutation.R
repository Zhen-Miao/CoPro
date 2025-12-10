

#' Generate Cell Permutation Indices
#'
#' Internal function to generate permutation indices for each cell type.
#'
#' @details
#' Permutation strategy:
#' - The FIRST cell type is kept FIXED (identity permutation)
#' - All OTHER cell types are permuted (either globally or bin-wise)
#'
#' This is the standard approach for permutation testing: we test whether
#' the relationship between cell types is stronger than expected by chance,
#' while keeping one cell type as reference.
#'
#' @param object A CoPro object
#' @param permu_method "global" or "bin"
#' @param nPermu Number of permutations
#' @param cts Cell types to permute
#' @param num_bins_x Number of bins in x for bin-wise permutation
#' @param num_bins_y Number of bins in y for bin-wise permutation
#'
#' @return List of permutation matrices, one per cell type
#' @keywords internal
.getCellPermu <- function(object, permu_method, nPermu, cts,
                          num_bins_x = 10, num_bins_y = 10) {

  cell_permu <- setNames(vector("list", length = length(cts)), cts)

  if (permu_method == "global") {
    # Global permutation: simple random shuffling
    for (i in cts) {
      n_cell <- sum(object@cellTypesSub == i)
      if (i == cts[1]) {
        # First cell type is kept fixed
        cell_permu[[i]] <- replicate(nPermu, 1:n_cell)
      } else {
        # Other cell types are randomly permuted
        cell_permu[[i]] <- replicate(nPermu,
                                     sample.int(n = n_cell, replace = FALSE))
      }
    }
  } else if (permu_method == "bin") {
    # Bin-wise permutation: preserves local spatial structure

    location_full <- object@locationDataSub
    location_full$"cell_ID" <- rownames(location_full)
    location_full$x_bin <- cut(location_full$x, breaks = num_bins_x, labels = FALSE)
    location_full$y_bin <- cut(location_full$y, breaks = num_bins_y, labels = FALSE)

    # Handle NA bins
    na_x <- is.na(location_full$x_bin)
    na_y <- is.na(location_full$y_bin)
    if (any(na_x)) {
      location_full$x_bin[na_x] <- 1
    }
    if (any(na_y)) {
      location_full$y_bin[na_y] <- 1
    }

    for (i in cts) {
      n_cell <- sum(object@cellTypesSub == i)
      if (i == cts[1]) {
        # First cell type is kept fixed
        cell_permu[[i]] <- replicate(nPermu, 1:n_cell)
      } else {
        # Other cell types are permuted bin-wise
        cell_loc <- location_full[object@cellTypesSub == i, ]
        cell_permu[[i]] <- matrix(ncol = nPermu, nrow = nrow(cell_loc))

        for (j in seq_len(nPermu)) {
          cell_loc_resample <- resample_spatial(location_data = cell_loc,
                                                num_bins_x = num_bins_x,
                                                num_bins_y = num_bins_y)
          cell_permu[[i]][, j] <- match(cell_loc_resample$"cell_ID",
                                        cell_loc$"cell_ID")
        }
      }
    }
  } else {
    stop(paste("permu_method must be 'global' or 'bin'.",
               "Got:", permu_method))
  }

  return(cell_permu)
}


#' Run Spatial CCA with Permutation Testing
#'
#' Performs permutation testing to assess the significance of spatial
#' co-progression detected by CoPro. This generates a null distribution
#' by permuting cell assignments while optionally preserving spatial structure.
#'
#' @details
#' ## Permutation Strategy
#'
#' The function supports two permutation methods:
#'
#' **"global"**: Simple random shuffling of cells. This breaks ALL spatial
#' structure and tests against a null of complete spatial randomness.
#'
#' **"bin"** (recommended): Bin-wise shuffling that preserves local spatial
#' structure. This tests against a null where cells have spatial autocorrelation
#' within their type, but no coordination across types.
#'
#' ## Important Notes
#'
#' - The FIRST cell type is always kept FIXED (not permuted)
#' - All other cell types are permuted
#' - For two cell types A and B, only B is permuted while A stays fixed
#' - This is equivalent to testing: "Is the A-B correlation stronger than
#'   expected if B cells were randomly rearranged (within spatial constraints)?"
#'
#' @param object A `CoPro` object with CCA already computed via `runSkrCCA()`
#' @param tol Tolerance for CCA optimization convergence (default: 1e-5)
#' @param nPermu Number of permutations to run (default: 20).
#'   Increase to 100+ for publication-quality p-values.
#' @param maxIter Maximum iterations for CCA optimization (default: 200)
#' @param permu_method Method of permutation:
#'   \itemize{
#'     \item "bin" (default): Bin-wise permutation preserving local spatial structure
#'     \item "global": Simple random permutation breaking all spatial structure
#'   }
#' @param num_bins_x Number of bins in x direction for bin-wise permutation (default: 10).
#'   Use `diagnose_bin_distribution()` to choose appropriate values.
#' @param num_bins_y Number of bins in y direction for bin-wise permutation (default: 10)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return CoPro object with permutation results stored in `@skrCCAPermuOut`
#'
#' @seealso
#' \code{\link{computeNormalizedCorrelationPermu}} to compute normalized
#' correlation from permutation results
#' \code{\link{diagnose_bin_distribution}} to check bin distribution
#'
#' @examples
#' \dontrun{
#' # After running standard CoPro analysis
#' br <- runSkrCCA(br, scalePCs = TRUE)
#' br <- computeNormalizedCorrelation(br)
#'
#' # Run permutation testing
#' br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "bin")
#' br <- computeNormalizedCorrelationPermu(br)
#'
#' # Calculate p-value
#' observed <- max(getNormCorr(br)$normalizedCorrelation)
#' permu_values <- sapply(br@normalizedCorrelationPermu,
#'                        function(x) x$normalizedCorrelation[1])
#' p_value <- mean(permu_values >= observed)
#' }
#'
#' @importFrom stats setNames
#' @export
runSkrCCAPermu <- function(object, tol = 1e-5, nPermu = 20,
                           maxIter = 200, permu_method = "bin",
                           num_bins_x = 10, num_bins_y = 10, verbose = TRUE) {

  ## Input validation
  if (!is(object, "CoPro")) {
    stop("Input object must be a CoPro object")
  }

  if (!(permu_method %in% c("bin", "global"))) {
    stop("permu_method must be 'bin' or 'global'.")
  }

  if (nPermu < 10) {
    warning("nPermu < 10 may give unreliable p-values. Consider nPermu >= 100.")
  }

  ## Store nPermu
  object@nPermu <- as.integer(nPermu)

  ## Check that runSkrCCA() has been run
  if (length(object@skrCCAOut) == 0) {
    stop("Please run runSkrCCA() before runSkrCCAPermu()")
  }
  nCC <- object@nCC
  scalePCs <- object@scalePCs

  ## Get cell types
  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    warning("No cell types of interest specified, using all cell types.")
    cts <- unique(object@cellTypesSub)
  }

  ## Get sigma value
  sigmaValueChoice <- object@sigmaValueChoice
  if (is.null(sigmaValueChoice) || length(sigmaValueChoice) == 0) {
    stop(paste("sigmaValueChoice not set.",
               "Please run computeNormalizedCorrelation() first."))
  }

  if (!(sigmaValueChoice %in% object@sigmaValues)) {
    stop("sigmaValueChoice does not exist in the list of sigmaValues")
  }

  s_name <- paste("sigma", sigmaValueChoice, sep = "_")

  ## Print bin diagnostics for bin-wise permutation
  if (verbose && permu_method == "bin") {
    cat("Bin-wise permutation settings:\n")
    cat(paste("  num_bins_x:", num_bins_x, "\n"))
    cat(paste("  num_bins_y:", num_bins_y, "\n"))
    cat(paste("  Total bins:", num_bins_x * num_bins_y, "\n"))
    cat("Note: First cell type is kept FIXED, others are permuted.\n\n")
  }

  ## Initialize output
  cca_permu_out <- vector("list", length = nPermu)
  permu_names <- paste("permu", 1:nPermu, sep = "_")
  names(cca_permu_out) <- permu_names

  ## Step 1: Generate cell permutations
  if (verbose) {
    cat("Generating cell permutations...\n")
  }
  cell_permu <- .getCellPermu(object = object, permu_method = permu_method,
                              nPermu = nPermu, cts = cts,
                              num_bins_x = num_bins_x, num_bins_y = num_bins_y)
  if (verbose) {
    cat("Cell permutation indices generated.\n\n")
  }
  object@cellPermu <- cell_permu

  ## Get PCA matrices
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
  PCmats2 <- PCmats

  ## Step 2: Run CCA for each permutation
  if (verbose) {
    cat(paste("Running CCA optimization for", nPermu, "permutations...\n"))
  }

  for (tt in seq_len(nPermu)) {
    t <- permu_names[tt]

    # Apply permutation to PC matrices
    for (i in names(PCmats)) {
      PCmats2[[i]] <- PCmats[[i]][cell_permu[[i]][, tt], ]
    }

    # Run CCA
    cca_result <- optimize_bilinear(
      X_list = PCmats2,
      flat_kernels = object@kernelMatrices,
      sigma = sigmaValueChoice,
      max_iter = maxIter,
      tol = tol
    )
    names(cca_result) <- cts

    if (nCC == 1) {
      cca_permu_out[[t]] <- cca_result
    } else {
      cca_result_n <- optimize_bilinear_n(
        X_list = PCmats2,
        flat_kernels = object@kernelMatrices,
        sigma = sigmaValueChoice,
        w_list = cca_result,
        cellTypesOfInterest = cts,
        nCC = nCC,
        max_iter = maxIter,
        tol = tol
      )
      cca_permu_out[[t]] <- cca_result_n
    }

    # Progress indicator
    if (verbose && (tt %% 10 == 0 || tt == nPermu)) {
      cat(paste("  Completed", tt, "of", nPermu, "permutations\n"))
    }
  }

  object@skrCCAPermuOut <- cca_permu_out

  if (verbose) {
    cat("\nPermutation testing complete.\n")
    cat("Run computeNormalizedCorrelationPermu() to compute p-values.\n")
  }

  return(object)
}


#' Compute Normalized Correlation for Permutation Results
#'
#' Calculates the normalized correlation for each permutation, enabling
#' p-value calculation by comparing observed values to the null distribution.
#'
#' @details
#' The normalized correlation for each permutation is calculated using the
#' same formula as the observed data:
#'
#' \deqn{NC = \frac{s_1^T K_{12} s_2}{||s_1|| \cdot ||s_2|| \cdot ||K_{12}||_2}}
#'
#' where \eqn{s_1} and \eqn{s_2} are cell scores, \eqn{K_{12}} is the kernel
#' matrix, and \eqn{||K_{12}||_2} is the spectral norm.
#'
#' @param object A `CoPro` object with permutation results from `runSkrCCAPermu()`
#' @param tol Tolerance for approximate SVD calculation (default: 1e-4)
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
#' p_value <- mean(permu_values >= observed)
#' }
#'
#' @export
computeNormalizedCorrelationPermu <- function(object, tol = 1e-4) {

  ## Input validation
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  if (length(object@skrCCAPermuOut) == 0) {
    stop(paste("skrCCAPermuOut is not available.",
               "Please run runSkrCCAPermu() first."))
  }

  ## Get parameters
  cts <- object@cellTypesOfInterest
  nPermu <- object@nPermu

  if (length(object@scalePCs) == 0) {
    stop("object@scalePCs not specified")
  }
  scalePCs <- object@scalePCs

  sigmaValueChoice <- object@sigmaValueChoice
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
  nCC <- object@nCC

  pair_cell_types <- combn(cts, 2)

  ## Initialize output
  correlation_value <- vector("list", length = nPermu)
  permu_names <- paste("permu", 1:nPermu, sep = "_")
  names(correlation_value) <- permu_names
  s_name <- paste("sigma", sigmaValueChoice, sep = "_")

  ## Calculate spectral norms (only need to do this once)
  cat("Calculating spectral norms...\n")
  norm_K12 <- setNames(vector(mode = "list", length = 1), s_name)
  norm_K12[[s_name]] <- setNames(vector(mode = "list", length = length(cts)), cts)

  for (i in cts) {
    norm_K12[[s_name]][[i]] <- setNames(
      vector(mode = "list", length = length(cts)), cts)
  }

  for (pp in seq_len(ncol(pair_cell_types))) {
    cellType1 <- pair_cell_types[1, pp]
    cellType2 <- pair_cell_types[2, pp]
    K <- getKernelMatrix(object, sigma = sigmaValueChoice,
                         cellType1 = cellType1, cellType2 = cellType2,
                         verbose = FALSE)
    svd_result <- irlba::irlba(K, nv = 1, tol = tol)
    norm_K12[[s_name]][[cellType1]][[cellType2]] <- svd_result$d[1]
  }
  cat("Spectral norms calculated.\n\n")

  ## Calculate normalized correlation for each permutation
  cat("Computing normalized correlations for permutations...\n")

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

    for (pp in seq_len(ncol(pair_cell_types))) {
      for (cc_index in seq_len(nCC)) {
        cellType1 <- pair_cell_types[1, pp]
        cellType2 <- pair_cell_types[2, pp]

        w_1 <- object@skrCCAPermuOut[[t]][[cellType1]][, cc_index, drop = FALSE]
        w_2 <- object@skrCCAPermuOut[[t]][[cellType2]][, cc_index, drop = FALSE]

        A <- PCmats[[cellType1]][object@cellPermu[[cellType1]][, tt], ]
        B <- PCmats[[cellType2]][object@cellPermu[[cellType2]][, tt], ]

        A_w1 <- A %*% w_1
        B_w2 <- B %*% w_2

        K <- getKernelMatrix(object, sigma = sigmaValueChoice,
                             cellType1 = cellType1, cellType2 = cellType2,
                             verbose = FALSE)
        norm_K12_sel <- norm_K12[[s_name]][[cellType1]][[cellType2]]

        # Calculate normalized correlation
        correlation_value[[t]]$"normalizedCorrelation"[
          pp + (cc_index - 1) * ncol(pair_cell_types)] <-
          (t(A_w1) %*% K %*% B_w2) /
          (sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K12_sel)
      }
    }

    # Progress indicator
    if (tt %% 20 == 0 || tt == nPermu) {
      cat(paste("  Completed", tt, "of", nPermu, "permutations\n"))
    }
  }

  ## Store results
  object@normalizedCorrelationPermu <- correlation_value

  cat("\nNormalized correlation computation complete.\n")

  return(object)
}


#' Calculate P-value from Permutation Results
#'
#' Helper function to calculate p-value from permutation testing results.
#'
#' @param object A CoPro object with permutation results
#' @param cc_index Which canonical correlation component to use (default: 1)
#' @param alternative Direction of test: "greater" (default), "less", or "two.sided"
#'
#' @return List with p-value, observed value, and permutation distribution
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

  if (length(object@normalizedCorrelationPermu) == 0) {
    stop("Run computeNormalizedCorrelationPermu() first")
  }

  # Get observed value
  ncorr <- getNormCorr(object)
  ncorr_cc <- ncorr[ncorr$CC_index == cc_index, ]
  observed <- max(ncorr_cc$normalizedCorrelation)

  # Get permutation values
  permu_values <- sapply(object@normalizedCorrelationPermu, function(x) {
    x$normalizedCorrelation[x$CC_index == cc_index][1]
  })

  # Calculate p-value
  if (alternative == "greater") {
    p_value <- mean(permu_values >= observed)
  } else if (alternative == "less") {
    p_value <- mean(permu_values <= observed)
  } else if (alternative == "two.sided") {
    p_value <- 2 * min(mean(permu_values >= observed),
                       mean(permu_values <= observed))
    p_value <- min(p_value, 1)  # Cap at 1
  } else {
    stop("alternative must be 'greater', 'less', or 'two.sided'")
  }

  result <- list(
    p_value = p_value,
    observed = observed,
    permu_mean = mean(permu_values),
    permu_sd = sd(permu_values),
    permu_values = permu_values,
    n_permu = length(permu_values),
    alternative = alternative
  )

  return(result)
}
