

#' Generate Cell Permutation Indices
#'
#' Internal function to generate permutation indices for each cell type.
#'
#' @details
#' Permutation strategy is controlled by `permu_which`:
#' - "second_only" (default): Keep first cell type FIXED, permute others
#' - "both": Permute ALL cell types independently
#' - "first_only": Keep second+ cell types FIXED, permute only first
#'
#' The default ("second_only") is the standard approach for permutation testing:
#' we test whether the relationship between cell types is stronger than expected
#' by chance, while keeping one cell type as reference.
#'
#' @param object A CoPro object
#' @param permu_method "global", "bin", "pc", or "toroidal"
#' @param nPermu Number of permutations
#' @param cts Cell types to permute
#' @param permu_which Which cell types to permute: "second_only", "both", or "first_only"
#' @param num_bins_x Number of bins in x for bin-wise permutation
#' @param num_bins_y Number of bins in y for bin-wise permutation
#' @param match_quantile Logical. If TRUE and permu_method="bin", matches cells
#'   between tiles based on their relative (quantile) positions to better preserve
#'   within-tile spatial structure. Default: FALSE.
#'
#' @return List of permutation matrices (for "global"/"bin"/"toroidal") or
#'   list of permuted PC matrices (for "pc"), one per cell type
#' @keywords internal
.getCellPermu <- function(object, permu_method, nPermu, cts,
                          permu_which = "second_only",
                          num_bins_x = 10, num_bins_y = 10,
                          match_quantile = FALSE) {

  cell_permu <- setNames(vector("list", length = length(cts)), cts)

  # Determine which cell types should be permuted
  should_permute <- function(ct_name, ct_index) {
    if (permu_which == "second_only") {
      return(ct_index > 1)  # Permute all except first
    } else if (permu_which == "first_only") {
      return(ct_index == 1)  # Permute only first
    } else if (permu_which == "both") {
      return(TRUE)  # Permute all
    } else {
      stop(paste("permu_which must be 'second_only', 'first_only', or 'both'.",
                 "Got:", permu_which))
    }
  }

  if (permu_method == "global") {
    # Global permutation: simple random shuffling
    for (idx in seq_along(cts)) {
      ct <- cts[idx]
      n_cell <- sum(object@cellTypesSub == ct)

      if (should_permute(ct, idx)) {
        # Permute this cell type
        cell_permu[[ct]] <- replicate(nPermu,
                                      sample.int(n = n_cell, replace = FALSE))
      } else {
        # Keep this cell type fixed (identity permutation)
        cell_permu[[ct]] <- replicate(nPermu, 1:n_cell)
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

    for (idx in seq_along(cts)) {
      ct <- cts[idx]
      n_cell <- sum(object@cellTypesSub == ct)

      if (should_permute(ct, idx)) {
        # Permute this cell type bin-wise
        cell_loc <- location_full[object@cellTypesSub == ct, ]
        cell_permu[[ct]] <- matrix(ncol = nPermu, nrow = nrow(cell_loc))

        for (j in seq_len(nPermu)) {
          cell_loc_resample <- resample_spatial(location_data = cell_loc,
                                                num_bins_x = num_bins_x,
                                                num_bins_y = num_bins_y,
                                                match_quantile = match_quantile)
          cell_permu[[ct]][, j] <- match(cell_loc_resample$"cell_ID",
                                         cell_loc$"cell_ID")
        }
      } else {
        # Keep this cell type fixed (identity permutation)
        cell_permu[[ct]] <- replicate(nPermu, 1:n_cell)
      }
    }
  } else if (permu_method == "toroidal") {
    # Toroidal shift permutation: perfectly preserves spatial autocorrelation
    # by shifting coordinates in a wrap-around manner
    location_full <- object@locationDataSub
    location_full$cell_ID <- rownames(location_full)

    for (idx in seq_along(cts)) {
      ct <- cts[idx]
      n_cell <- sum(object@cellTypesSub == ct)

      if (should_permute(ct, idx)) {
        # Apply toroidal shift to this cell type
        cell_loc <- location_full[object@cellTypesSub == ct, ]
        cell_permu[[ct]] <- generate_toroidal_permutations(cell_loc, nPermu)
      } else {
        # Keep this cell type fixed (identity permutation)
        cell_permu[[ct]] <- replicate(nPermu, 1:n_cell)
      }
    }
  } else if (permu_method == "pc") {
    # PC-space permutation (like DIALOGUE): shuffle values within each PC column
    # This breaks cell-to-cell correlation while preserving marginal PC distributions
    # We store random seeds for reproducibility in computeNormalizedCorrelationPermu()

    for (idx in seq_along(cts)) {
      ct <- cts[idx]
      n_cell <- sum(object@cellTypesSub == ct)

      if (should_permute(ct, idx)) {
        # Store seeds for reproducible PC permutation
        # Format: list with "type" = "pc_permute" and "seeds" = vector of seeds
        seeds <- sample.int(.Machine$integer.max, nPermu)
        cell_permu[[ct]] <- list(
          type = "pc_permute",
          seeds = seeds,
          n_cell = n_cell
        )
      } else {
        # Keep this cell type fixed (identity permutation)
        cell_permu[[ct]] <- replicate(nPermu, 1:n_cell)
      }
    }
  } else {
    stop(paste("permu_method must be 'global', 'bin', 'toroidal', or 'pc'.",
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
#' ## Permutation Methods
#'
#' The function supports three permutation methods:
#'
#' **"global"**: Simple random shuffling of cells. This breaks ALL spatial
#' structure and tests against a null of complete spatial randomness.
#'
#' **"bin"** (default): Bin-wise shuffling that preserves local spatial
#' structure. This tests against a null where cells have spatial autocorrelation
#' within their type, but no coordination across types.
#'
#' **"pc"**: PC-space permutation (like DIALOGUE). Shuffles values within each
#' PC dimension across cells, breaking cell-to-cell correlation while preserving
#' the marginal distribution of each PC. This is the same approach used by
#' DIALOGUE's internal significance testing. Use this to compare with DIALOGUE's
#' conservative behavior.
#'
#' **"toroidal"**: Toroidal (wrap-around) shift permutation. Shifts all cells'
#' coordinates by a random amount, wrapping at boundaries. This PERFECTLY
#' preserves spatial autocorrelation within each cell type because relative
#' positions are unchanged. Best choice for reducing FPR when spatial
#' autocorrelation is strong.
#'
#' ## Which Cell Types to Permute
#'
#' The `permu_which` parameter controls which cell types are permuted:
#'
#' **"second_only"** (default): Keep the first cell type FIXED, permute all others.
#' This is the standard approach for permutation testing. For two cell types A and B,
#' only B is permuted while A stays fixed.
#'
#' **"both"**: Permute ALL cell types independently. Both A and B are shuffled
#' with different random permutations. This breaks MORE structure and may
#' lead to higher FPR. Use "second_only" for better FPR control.
#'
#' **"first_only"**: Keep second+ cell types FIXED, permute only the first.
#' Useful if you want to test from the opposite direction.
#'
#' ## Controlling False Positive Rate
#'
#' High FPR under null simulations typically means the permutation is **breaking
#' too much spatial structure**, making the null distribution too low, so that
#' even random observed values appear significant.
#'
#' To **reduce FPR**, you need to **better preserve spatial autocorrelation**:
#'
#' 1. **Use MORE bins** (not fewer): `num_bins_x = 15, num_bins_y = 15` preserves
#'    more local structure than the default 10x10.
#'
#' 2. **Enable quantile matching**: `match_quantile = TRUE` matches cells by their
#'    relative positions within tiles, better preserving within-tile structure.
#'
#' 3. **Only permute one cell type**: `permu_which = "second_only"` (default)
#'    keeps one cell type fixed, preserving more structure than "both".
#'
#' 4. **Avoid global permutation**: `permu_method = "global"` breaks ALL spatial
#'    structure and will likely have high FPR.
#'
#' The key insight: if permutation doesn't adequately preserve within-type spatial
#' autocorrelation, the null distribution will be too low, leading to inflated
#' significance.
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
#'     \item "pc": PC-space permutation (like DIALOGUE) - shuffles values within each
#'       PC dimension, breaking cell correlation while preserving PC distributions
#'     \item "toroidal": Toroidal shift permutation - perfectly preserves spatial
#'       autocorrelation by shifting coordinates in wrap-around manner
#'   }
#' @param permu_which Which cell types to permute:
#'   \itemize{
#'     \item "second_only" (default): Keep first cell type fixed, permute others
#'     \item "both": Permute all cell types independently (more conservative)
#'     \item "first_only": Keep others fixed, permute only the first cell type
#'   }
#' @param num_bins_x Number of bins in x direction for bin-wise permutation (default: 10).
#'   Use `diagnose_bin_distribution()` to choose appropriate values.
#'   **More bins = better preserve local structure = lower FPR.**
#' @param num_bins_y Number of bins in y direction for bin-wise permutation (default: 10).
#'   **More bins = better preserve local structure = lower FPR.**
#' @param match_quantile Logical. If TRUE and `permu_method = "bin"`, matches cells
#'   between tiles based on their relative (quantile) x/y positions. This better
#'   preserves within-tile spatial autocorrelation structure. Default: FALSE.
#'   **Setting TRUE helps reduce FPR by better preserving spatial structure.**
#' @param conservative Logical. If TRUE, automatically uses settings that better
#'   preserve spatial autocorrelation to reduce false positive rate:
#'   `permu_which = "second_only"`, `num_bins_x = 15`, `num_bins_y = 15`,
#'   `match_quantile = TRUE`. Default: FALSE.
#' @param n_cores Number of cores for parallel computation (default: 1).
#'   Set to higher values to speed up permutation testing. Use `parallel::detectCores()`
#'   to find available cores. Parallelization uses the `parallel` package.
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
#' # Standard permutation (only permute second cell type)
#' br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "bin",
#'                      permu_which = "second_only")
#'
#' # Conservative permutation (lower FPR - better preserves spatial structure)
#' br <- runSkrCCAPermu(br, nPermu = 100, conservative = TRUE)
#'
#' # Manual conservative settings: more bins + quantile matching
#' br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "bin",
#'                      permu_which = "second_only",
#'                      num_bins_x = 15, num_bins_y = 15,
#'                      match_quantile = TRUE)
#'
#' # PC-space permutation (like DIALOGUE) - shuffles within PC dimensions
#' br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "pc",
#'                      permu_which = "second_only")
#'
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
                           permu_which = "second_only",
                           num_bins_x = 10, num_bins_y = 10,
                           match_quantile = FALSE,
                           conservative = FALSE,
                           n_cores = 1, verbose = TRUE) {

  ## Input validation
  if (!is(object, "CoPro")) {
    stop("Input object must be a CoPro object")
  }

  ## Apply conservative settings if requested
  ## Conservative = better preserve spatial autocorrelation = lower FPR
  if (conservative) {
    if (verbose) {
      cat("Using CONSERVATIVE mode for lower false positive rate:\n")
      cat("  -> permu_which = 'second_only' (only permute one cell type)\n")
      cat("  -> num_bins_x = 15, num_bins_y = 15 (more bins = better preserve local structure)\n")
      cat("  -> match_quantile = TRUE (preserve within-tile spatial positions)\n\n")
    }
    permu_which <- "second_only"
    num_bins_x <- 15
    num_bins_y <- 15
    match_quantile <- TRUE
  }

  if (!(permu_method %in% c("bin", "global", "pc", "toroidal"))) {
    stop("permu_method must be 'bin', 'global', 'pc', or 'toroidal'.")
  }

  if (!(permu_which %in% c("second_only", "both", "first_only"))) {
    stop("permu_which must be 'second_only', 'both', or 'first_only'.")
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

  ## Print permutation settings
  if (verbose) {
    cat("Permutation settings:\n")
    cat(paste("  permu_method:", permu_method, "\n"))
    cat(paste("  permu_which:", permu_which, "\n"))

    if (permu_which == "second_only") {
      cat(paste("    -> Cell type", cts[1], "is FIXED, others are permuted\n"))
    } else if (permu_which == "first_only") {
      cat(paste("    -> Only cell type", cts[1], "is permuted, others are FIXED\n"))
    } else if (permu_which == "both") {
      cat("    -> ALL cell types are permuted independently\n")
    }

    if (permu_method == "bin") {
      cat(paste("  num_bins_x:", num_bins_x, "\n"))
      cat(paste("  num_bins_y:", num_bins_y, "\n"))
      cat(paste("  Total bins:", num_bins_x * num_bins_y, "\n"))
      cat(paste("  match_quantile:", match_quantile, "\n"))
      if (match_quantile) {
        cat("    -> Cells matched by within-tile quantile position\n")
      }
    } else if (permu_method == "toroidal") {
      cat("    -> Toroidal (wrap-around) shift permutation\n")
      cat("    -> PERFECTLY preserves spatial autocorrelation\n")
      cat("    -> Best for reducing FPR when spatial structure is strong\n")
    } else if (permu_method == "pc") {
      cat("    -> PC-space permutation (like DIALOGUE)\n")
      cat("    -> Shuffles values within each PC dimension across cells\n")
      cat("    -> Breaks cell correlation, preserves PC marginal distributions\n")
    }
    cat("\n")
  }

  ## Initialize output
  cca_permu_out <- vector("list", length = nPermu)
  permu_names <- paste("permu", 1:nPermu, sep = "_")
  names(cca_permu_out) <- permu_names

  ## Step 1: Generate cell permutations
  if (verbose) {
    if (permu_method == "pc") {
      cat("Setting up PC-space permutation (permutation done on-the-fly)...\n")
    } else {
      cat("Generating cell permutation indices...\n")
    }
  }
  cell_permu <- .getCellPermu(object = object, permu_method = permu_method,
                              nPermu = nPermu, cts = cts,
                              permu_which = permu_which,
                              num_bins_x = num_bins_x, num_bins_y = num_bins_y,
                              match_quantile = match_quantile)
  if (verbose) {
    if (permu_method == "pc") {
      cat("PC-space permutation configured.\n\n")
    } else {
      cat("Cell permutation indices generated.\n\n")
    }
  }
  object@cellPermu <- cell_permu

  ## Get PCA matrices
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
  PCmats2 <- PCmats

  ## Step 2: Run CCA for each permutation
  if (verbose) {
    cat(paste("Running CCA optimization for", nPermu, "permutations"))
    if (n_cores > 1) {
      cat(paste(" using", n_cores, "cores...\n"))
    } else {
      cat("...\n")
    }
  }


  # Helper function for PC-space permutation (like DIALOGUE)
  # Shuffles values within each PC column across cells
  # Uses seed for reproducibility
  permute_pc_matrix <- function(pc_mat, seed = NULL) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    # Apply permutation within each column (PC dimension)
    # This is exactly what DIALOGUE does: apply(X1, 2, function(x) sample(x, length(x)))
    permuted <- apply(pc_mat, 2, function(x) sample(x, length(x)))
    rownames(permuted) <- rownames(pc_mat)
    return(permuted)
  }

  # Define worker function for single permutation
  run_single_permu <- function(tt) {
    PCmats_local <- PCmats

    for (i in names(PCmats_local)) {
      if (is.list(cell_permu[[i]]) && !is.null(cell_permu[[i]]$type) && 
          cell_permu[[i]]$type == "pc_permute") {
        # PC-space permutation: shuffle within each PC column using stored seed
        seed <- cell_permu[[i]]$seeds[tt]
        PCmats_local[[i]] <- permute_pc_matrix(PCmats[[i]], seed = seed)
      } else {
        # Standard index-based permutation (global or bin)
        PCmats_local[[i]] <- PCmats[[i]][cell_permu[[i]][, tt], ]
      }
    }

    cca_result <- optimize_bilinear(
      X_list = PCmats_local,
      flat_kernels = object@kernelMatrices,
      sigma = sigmaValueChoice,
      max_iter = maxIter,
      tol = tol
    )
    names(cca_result) <- cts

    if (nCC == 1) {
      return(cca_result)
    } else {
      cca_result_n <- optimize_bilinear_n(
        X_list = PCmats_local,
        flat_kernels = object@kernelMatrices,
        sigma = sigmaValueChoice,
        w_list = cca_result,
        cellTypesOfInterest = cts,
        nCC = nCC,
        max_iter = maxIter,
        tol = tol
      )
      return(cca_result_n)
    }
  }

  if (n_cores > 1) {
    # Parallel execution using PSOCK cluster
    # Note: Requires CoPro to be INSTALLED (not just loaded with devtools::load_all())
    # Fork-based parallelism (mclapply) causes segfaults with BLAS/irlba on macOS

    # Check if CoPro is installed
    if (!requireNamespace("CoPro", quietly = TRUE)) {
      warning("Parallel execution requires CoPro to be installed (not just loaded).\n",
              "Falling back to sequential execution.\n",
              "To enable parallelism, run: devtools::install() first.")
      n_cores <- 1
    }
  }

  if (n_cores > 1) {
    if (verbose) {
      cat("  Setting up parallel cluster with", n_cores, "workers...\n")
    }

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Load CoPro package in workers
    parallel::clusterEvalQ(cl, library(CoPro))

    # Export required objects to workers
    flat_kernels_local <- object@kernelMatrices
    parallel::clusterExport(cl, c("PCmats", "cell_permu", "cts", "nCC",
                                  "sigmaValueChoice", "maxIter", "tol",
                                  "flat_kernels_local"),
                            envir = environment())

    # Helper function for PC-space permutation (for parallel workers)
    permute_pc_matrix_parallel <- function(pc_mat, seed = NULL) {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      permuted <- apply(pc_mat, 2, function(x) sample(x, length(x)))
      rownames(permuted) <- rownames(pc_mat)
      return(permuted)
    }

    # Export the helper function
    parallel::clusterExport(cl, "permute_pc_matrix_parallel", envir = environment())

    # Worker function
    worker_fn <- function(tt) {
      PCmats_local <- PCmats

      for (i in names(PCmats_local)) {
        if (is.list(cell_permu[[i]]) && !is.null(cell_permu[[i]]$type) && 
            cell_permu[[i]]$type == "pc_permute") {
          # PC-space permutation: shuffle within each PC column using stored seed
          seed <- cell_permu[[i]]$seeds[tt]
          PCmats_local[[i]] <- permute_pc_matrix_parallel(PCmats[[i]], seed = seed)
        } else {
          # Standard index-based permutation (global or bin)
          PCmats_local[[i]] <- PCmats[[i]][cell_permu[[i]][, tt], ]
        }
      }

      cca_result <- CoPro::optimize_bilinear(
        X_list = PCmats_local,
        flat_kernels = flat_kernels_local,
        sigma = sigmaValueChoice,
        max_iter = maxIter,
        tol = tol
      )
      names(cca_result) <- cts

      if (nCC == 1) {
        return(cca_result)
      } else {
        cca_result_n <- CoPro::optimize_bilinear_n(
          X_list = PCmats_local,
          flat_kernels = flat_kernels_local,
          sigma = sigmaValueChoice,
          w_list = cca_result,
          cellTypesOfInterest = cts,
          nCC = nCC,
          max_iter = maxIter,
          tol = tol
        )
        return(cca_result_n)
      }
    }

    cca_results_list <- parallel::parLapply(cl, seq_len(nPermu), worker_fn)

    # Convert list to named list
    names(cca_results_list) <- permu_names
    cca_permu_out <- cca_results_list

    if (verbose) {
      cat("  Parallel computation complete.\n")
    }

  } else {
    # Sequential execution (original behavior)
    for (tt in seq_len(nPermu)) {
      t <- permu_names[tt]
      cca_permu_out[[t]] <- run_single_permu(tt)

      # Progress indicator
      if (verbose && (tt %% 10 == 0 || tt == nPermu)) {
        cat(paste("  Completed", tt, "of", nPermu, "permutations\n"))
      }
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

  ## Helper function for PC-space permutation (recreate using stored seed)
  permute_pc_matrix_ncorr <- function(pc_mat, seed) {
    set.seed(seed)
    permuted <- apply(pc_mat, 2, function(x) sample(x, length(x)))
    rownames(permuted) <- rownames(pc_mat)
    return(permuted)
  }

  ## Helper to get permuted PC matrix for a cell type
  get_permuted_pcmat <- function(ct, tt) {
    cell_permu_ct <- object@cellPermu[[ct]]

    if (is.list(cell_permu_ct) && !is.null(cell_permu_ct$type) &&
        cell_permu_ct$type == "pc_permute") {
      # PC-space permutation: recreate using stored seed
      seed <- cell_permu_ct$seeds[tt]
      return(permute_pc_matrix_ncorr(PCmats[[ct]], seed))
    } else {
      # Standard index-based permutation (global or bin)
      return(PCmats[[ct]][cell_permu_ct[, tt], ])
    }
  }

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

        A <- get_permuted_pcmat(cellType1, tt)
        B <- get_permuted_pcmat(cellType2, tt)

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
                       verbose = FALSE)

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

  ## Calculate spectral norm of kernel matrix
  svd_result <- irlba::irlba(K, nv = 1, tol = tol)
  norm_K12 <- svd_result$d[1]

  ## Calculate normalized correlation
  norm_corr <- as.numeric(t(s1) %*% K %*% s2) / norm_K12

  return(norm_corr)
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
#'
#' @return Numeric value of normalized correlation
#' @keywords internal
.compute_ncorr_quick <- function(PCmats, w_list, flat_kernels, sigma, cts,
                                 tol = 1e-4) {
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

  A_w1 <- A %*% w1
  B_w2 <- B %*% w2

  # Get kernel matrix
  K <- get_kernel_matrix_flat(flat_kernels, sigma, ct1, ct2, slide = NULL)

  # Compute spectral norm
  svd_result <- irlba::irlba(K, nv = 1, tol = tol)
  norm_K12 <- svd_result$d[1]

  # Normalized correlation
  numerator <- as.numeric(t(A_w1) %*% K %*% B_w2)
  denominator <- sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K12

  return(numerator / denominator)
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
#' @param nPermu Number of permutations to run (default: 100)
#' @param sigma_values Vector of sigma values to test. If NULL, uses all
#'   sigma values from the original analysis (object@@sigmaValues)
#' @param permu_method Method of permutation: "bin", "global", "pc", or "toroidal"
#' @param permu_which Which cell types to permute: "second_only", "both", "first_only"
#' @param num_bins_x Number of bins in x for bin-wise permutation
#' @param num_bins_y Number of bins in y for bin-wise permutation
#' @param match_quantile Whether to use quantile matching for bin permutation
#' @param maxIter Maximum iterations for CCA optimization
#' @param tol Convergence tolerance
#' @param n_cores Number of cores for parallel computation (not yet implemented)
#' @param verbose Whether to print progress messages
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
                                     nPermu = 100,
                                     sigma_values = NULL,
                                     permu_method = "bin",
                                     permu_which = "second_only",
                                     num_bins_x = 10,
                                     num_bins_y = 10,
                                     match_quantile = FALSE,
                                     maxIter = 200,
                                     tol = 1e-5,
                                     n_cores = 1,
                                     verbose = TRUE) {

  ## Input validation
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  if (length(object@skrCCAOut) == 0) {
    stop("Please run runSkrCCA() first")
  }

  if (length(object@normalizedCorrelation) == 0) {
    stop("Please run computeNormalizedCorrelation() first")
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
  scalePCs <- object@scalePCs
  nCC <- object@nCC

  if (nCC > 1) {
    warning("Fair sigma permutation currently only supports nCC = 1. Using first CC.")
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
    match_quantile = match_quantile
  )

  # Store results
  permu_results <- vector("list", nPermu)
  permu_ncorrs <- numeric(nPermu)
  permu_sigmas <- numeric(nPermu)

  # Helper function for PC-space permutation
  permute_pc_matrix <- function(pc_mat, seed) {
    set.seed(seed)
    permuted <- apply(pc_mat, 2, function(x) sample(x, length(x)))
    rownames(permuted) <- rownames(pc_mat)
    return(permuted)
  }

  if (verbose) {
    cat("Running permutations...\n")
  }

  for (tt in seq_len(nPermu)) {
    # Apply permutation to get permuted PC matrices
    PCmats_local <- PCmats

    for (i in names(PCmats_local)) {
      if (is.list(cell_permu[[i]]) && !is.null(cell_permu[[i]]$type) &&
          cell_permu[[i]]$type == "pc_permute") {
        # PC-space permutation
        seed <- cell_permu[[i]]$seeds[tt]
        PCmats_local[[i]] <- permute_pc_matrix(PCmats[[i]], seed)
      } else {
        # Standard index-based permutation (global, bin, or toroidal)
        PCmats_local[[i]] <- PCmats[[i]][cell_permu[[i]][, tt], ]
      }
    }

    # Test ALL sigma values and pick the best
    best_ncorr <- -Inf
    best_sigma <- sigma_values[1]
    best_weights <- NULL

    for (sigma in sigma_values) {
      # Run CCA for this sigma
      cca_result <- optimize_bilinear(
        X_list = PCmats_local,
        flat_kernels = object@kernelMatrices,
        sigma = sigma,
        max_iter = maxIter,
        tol = tol
      )
      names(cca_result) <- cts

      # Compute normalized correlation
      ncorr <- .compute_ncorr_quick(
        PCmats_local, cca_result,
        object@kernelMatrices, sigma, cts
      )

      if (ncorr > best_ncorr) {
        best_ncorr <- ncorr
        best_sigma <- sigma
        best_weights <- cca_result
      }
    }

    permu_results[[tt]] <- list(
      weights = best_weights,
      sigma = best_sigma,
      ncorr = best_ncorr
    )
    permu_ncorrs[tt] <- best_ncorr
    permu_sigmas[tt] <- best_sigma

    # Progress indicator
    if (verbose && (tt %% 10 == 0 || tt == nPermu)) {
      cat(paste("  Completed", tt, "of", nPermu, "permutations\n"))
    }
  }

  # Store results in object
  object@skrCCAPermuOut <- lapply(permu_results, function(x) x$weights)
  names(object@skrCCAPermuOut) <- paste0("permu_", seq_len(nPermu))

  # Create normalized correlation results structure
  pair_cell_types <- combn(cts, 2)
  correlation_value <- vector("list", nPermu)
  for (tt in seq_len(nPermu)) {
    correlation_value[[tt]] <- data.frame(
      sigmaValues = permu_results[[tt]]$sigma,
      cellType1 = pair_cell_types[1, 1],
      cellType2 = pair_cell_types[2, 1],
      CC_index = 1,
      normalizedCorrelation = permu_results[[tt]]$ncorr,
      stringsAsFactors = FALSE
    )
  }
  names(correlation_value) <- paste0("permu_", seq_len(nPermu))
  object@normalizedCorrelationPermu <- correlation_value
  object@nPermu <- as.integer(nPermu)

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
    n_sigma_differs = n_sigma_differs
  )

  # Get observed value for comparison
  observed_ncorr <- max(getNormCorr(object)$normalizedCorrelation)

  if (verbose) {
    cat("\n=== Fair Sigma Permutation Complete ===\n")
    cat(paste("Observed best ncorr:", round(observed_ncorr, 4), "\n"))
    cat(paste("Permutation mean:", round(mean(permu_ncorrs), 4), "\n"))
    cat(paste("Permutation SD:", round(sd(permu_ncorrs), 4), "\n"))
    cat(paste("P-value:", round(mean(permu_ncorrs >= observed_ncorr), 4), "\n"))
    cat("\nSigma selection in permutations:\n")
    cat(paste("  Observed best sigma:", observed_best_sigma, "\n"))
    cat(paste("  Sigma values used:", paste(unique(permu_sigmas), collapse = ", "), "\n"))
    cat(paste("  Most common:", names(sort(table(permu_sigmas), decreasing = TRUE))[1], "\n"))
    cat(paste("  Sigma differs from observed:", n_sigma_differs, "/", nPermu,
              "(", round(prop_sigma_differs * 100, 1), "%)\n"))
  }

  return(object)
}
