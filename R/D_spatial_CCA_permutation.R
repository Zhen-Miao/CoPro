

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
#' structure and tests against a null of complete spatial randomness. It is a
#' *deliberately-broken reference*: by destroying within-type autocorrelation it
#' inflates the effective sample size (Clifford-Richardson-Hemon 1989) and is
#' therefore anti-conservative by construction. Keep it for calibration, not as
#' a default.
#'
#' **"bin"** (default): Bin-wise shuffling that preserves local spatial
#' structure. This tests against a null where cells have spatial autocorrelation
#' within their type, but no coordination across types. By default the patch
#' grid is sized from the bandwidth (see [.sigmaAwareBins()]). This is a
#' restricted / approximate permutation that is valid under stationarity; under
#' autocorrelation the within-type cells are not exchangeable, which is exactly
#' the point (Anderson & ter Braak 2003).
#'
#' **"pc"**: PC-space permutation (like DIALOGUE). Shuffles values within each
#' PC dimension across cells, breaking cell-to-cell correlation while preserving
#' the marginal distribution of each PC. Like "global" this is a
#' *deliberately-broken reference* that destroys within-type autocorrelation;
#' use it to reproduce DIALOGUE-style complete-spatial-randomness behaviour and
#' to demonstrate FPR inflation, not as a default.
#'
#' **"toroidal"**: Toroidal (wrap-around) shift permutation. Shifts all cells'
#' coordinates by a random amount, wrapping at boundaries, then re-ranks to a
#' permutation. This *approximately* preserves within-type spatial
#' autocorrelation (exactly only on a regular lattice; see
#' [generate_toroidal_permutations()]) and assumes spatial stationarity and
#' periodic wrap-around, which is biologically false at tissue edges. Benchmark
#' it against the sigma-aware "bin" null rather than preferring it a priori.
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
#'     \item "toroidal": Toroidal shift permutation - approximately preserves
#'       spatial autocorrelation by shifting coordinates in a wrap-around manner
#'       (assumes stationarity/periodicity; see [generate_toroidal_permutations()])
#'   }
#' @param permu_which Which cell types to permute:
#'   \itemize{
#'     \item "second_only" (default): Keep first cell type fixed, permute others
#'     \item "both": Permute all cell types independently (more conservative)
#'     \item "first_only": Keep others fixed, permute only the first cell type
#'   }
#' @param num_bins_x Number of bins in x direction for bin-wise permutation.
#'   Default `NULL` sizes the grid automatically from the kernel bandwidth so
#'   each patch is ~2*sigma wide on the normalized distance scale (see
#'   [.sigmaAwareBins()]). Pass an explicit integer to override.
#'   Use `diagnose_bin_distribution()` to inspect a chosen grid.
#'   **More bins = better preserve local structure = lower FPR.**
#' @param num_bins_y Number of bins in y direction for bin-wise permutation.
#'   Default `NULL` (sigma-aware, as for `num_bins_x`). Pass an integer to override.
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
                           num_bins_x = NULL, num_bins_y = NULL,
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

  ## Resolve the bin grid for bin-wise permutation.
  ## When num_bins_x / num_bins_y are NULL (the default), size the patch grid
  ## from the kernel bandwidth so each patch is ~2*sigma wide on the normalized
  ## distance scale: large enough to preserve within-type spatial
  ## autocorrelation, small enough that shuffling whole patches still breaks
  ## cross-type coordination. Explicit integer values override this.
  if (permu_method == "bin" && (is.null(num_bins_x) || is.null(num_bins_y))) {
    bins <- .sigmaAwareBins(object, sigma = sigmaValueChoice, verbose = verbose)
    num_bins_x <- bins$num_bins_x
    num_bins_y <- bins$num_bins_y
  }

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
      cat("    -> Approximately preserves spatial autocorrelation (assumes stationarity/periodicity)\n")
      cat("    -> Benchmark against the sigma-aware 'bin' null rather than preferring a priori\n")
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
#' \deqn{NC = \frac{s_1^T K_{12} s_2}{||s_1|| \cdot ||s_2|| \cdot ||\tilde K_c||_F}}
#'
#' where \eqn{s_1} and \eqn{s_2} are cell scores, \eqn{K_{12}} is the kernel
#' matrix, and \eqn{||\tilde K_c||_F} is the whitened-Frobenius norm
#' \eqn{||R_x^{1/2} K_c R_y^{1/2}||_F}.
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

  ## Calculate whitened-Frobenius normalizers (only need to do this once)
  cat("Calculating whitened-Frobenius normalizers...\n")
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
    ## matched-sigma within-type kernels as whitening operators
    Rx <- tryCatch(getKernelMatrix(object, sigma = sigmaValueChoice,
                     cellType1 = cellType1, cellType2 = cellType1,
                     verbose = FALSE), error = function(e) NULL)
    Ry <- tryCatch(getKernelMatrix(object, sigma = sigmaValueChoice,
                     cellType1 = cellType2, cellType2 = cellType2,
                     verbose = FALSE), error = function(e) NULL)
    norm_K12[[s_name]][[cellType1]][[cellType2]] <- .whitenedFrobNorm(K, Rx, Ry)
  }
  cat("Whitened-Frobenius normalizers calculated.\n\n")

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

  ## Whitened-Frobenius normalizer (matched-sigma within-type kernels)
  Rx <- tryCatch(getKernelMatrix(object, sigma = sigma,
                   cellType1 = cellType1, cellType2 = cellType1,
                   verbose = FALSE), error = function(e) NULL)
  Ry <- tryCatch(getKernelMatrix(object, sigma = sigma,
                   cellType1 = cellType2, cellType2 = cellType2,
                   verbose = FALSE), error = function(e) NULL)
  norm_K12 <- .whitenedFrobNorm(K, Rx, Ry)

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
#' @return List with the Phipson & Smyth (2010) permutation p-value (`p_value`,
#'   never exactly zero), the Monte-Carlo floor `mc_floor = 1 / (n_permu + 1)`,
#'   the observed value, and the permutation distribution.
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

  # Calculate p-value using the Phipson & Smyth (2010) estimator: the observed
  # configuration is itself one admissible permutation, so 1 is added to both
  # numerator and denominator. This is always strictly positive (a permutation
  # p-value of exactly 0 is invalid); the smallest resolvable value is the
  # Monte-Carlo floor 1 / (m + 1).
  m <- length(permu_values)
  if (alternative == "greater") {
    p_value <- (1 + sum(permu_values >= observed)) / (1 + m)
  } else if (alternative == "less") {
    p_value <- (1 + sum(permu_values <= observed)) / (1 + m)
  } else if (alternative == "two.sided") {
    p_greater <- (1 + sum(permu_values >= observed)) / (1 + m)
    p_less <- (1 + sum(permu_values <= observed)) / (1 + m)
    p_value <- min(2 * min(p_greater, p_less), 1)  # Cap at 1
  } else {
    stop("alternative must be 'greater', 'less', or 'two.sided'")
  }

  result <- list(
    p_value = p_value,
    mc_floor = 1 / (1 + m),
    observed = observed,
    permu_mean = mean(permu_values),
    permu_sd = sd(permu_values),
    permu_values = permu_values,
    n_permu = m,
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

  # Whitened-Frobenius normalizer (matched-sigma within-type kernels)
  Rx <- tryCatch(get_kernel_matrix_flat(flat_kernels, sigma, ct1, ct1, slide = NULL),
                 error = function(e) NULL)
  Ry <- tryCatch(get_kernel_matrix_flat(flat_kernels, sigma, ct2, ct2, slide = NULL),
                 error = function(e) NULL)
  norm_K12 <- .whitenedFrobNorm(K, Rx, Ry)

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
#' @param num_bins_x Number of bins in x for bin-wise permutation. Default `NULL`
#'   sizes the grid from the observed best bandwidth (`sigmaValueChoice`) via
#'   [.sigmaAwareBins()]; the same grid (and hence the same permutation) is
#'   shared across the sigma sweep. Pass an integer to override.
#' @param num_bins_y Number of bins in y for bin-wise permutation. Default `NULL`
#'   (sigma-aware, as for `num_bins_x`).
#' @param match_quantile Whether to use quantile matching for bin permutation
#' @param maxIter Maximum iterations for CCA optimization
#' @param tol Convergence tolerance
#' @param n_cores Number of cores for parallel computation (not yet implemented)
#' @param verbose Whether to print progress messages
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
                                     nPermu = 100,
                                     sigma_values = NULL,
                                     permu_method = "bin",
                                     permu_which = "second_only",
                                     num_bins_x = NULL,
                                     num_bins_y = NULL,
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
    cat(paste("P-value (Phipson-Smyth):",
              round((1 + sum(permu_ncorrs >= observed_ncorr)) / (1 + nPermu), 4),
              paste0("(MC floor ", round(1 / (1 + nPermu), 4), ")"), "\n"))
    cat("\nSigma selection in permutations:\n")
    cat(paste("  Observed best sigma:", observed_best_sigma, "\n"))
    cat(paste("  Sigma values used:", paste(unique(permu_sigmas), collapse = ", "), "\n"))
    cat(paste("  Most common:", names(sort(table(permu_sigmas), decreasing = TRUE))[1], "\n"))
    cat(paste("  Sigma differs from observed:", n_sigma_differs, "/", nPermu,
              "(", round(prop_sigma_differs * 100, 1), "%)\n"))
  }

  return(object)
}


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
#' Under whitened PCs (`scalePCs = TRUE`, so `X^T X = c I`) this Y-space
#' deflation is algebraically identical to the Freedman-Lane / ter Braak
#' residualization that removes the observed lower canonical variates from the
#' data before recomputing `Y`, because
#' `(I - u u^T) Y (I - v v^T) = Y - (u^T Y v)\, u v^T`.
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
#' @param maxIter,tol Optimization controls.
#'
#' @return List with `w` (named list of 1-column weight matrices for axis k) and
#'   `ncorr` (its normalized correlation).
#' @importFrom utils capture.output
#' @keywords internal
.fitConditionalAxis <- function(PCmats, flat_kernels, sigma, cts,
                                W_lower = NULL, k_minus_1 = 0,
                                Y_resi = NULL, maxIter = 200, tol = 1e-5) {
  if (k_minus_1 <= 0) {
    # First axis: identical to the fair-sigma CC1 optimization (X-based).
    invisible(utils::capture.output(
      w_k <- optimize_bilinear(
        X_list = PCmats, flat_kernels = flat_kernels, sigma = sigma,
        max_iter = maxIter, tol = tol
      )
    ))
    names(w_k) <- cts
  } else {
    # Conditional axis: deflate the fixed observed lower directions from Y,
    # then optimize the leading residual component.
    if (is.null(Y_resi)) {
      Y_resi <- compute_Y_resi(PCmats, flat_kernels, sigma, cts, slide = NULL)
    }
    Yk <- Y_resi
    for (qq in seq_len(k_minus_1)) {
      Yk <- apply_deflation(Yk, W_lower, qq, cts)
    }
    w_new <- initialize_next_component(Yk, cts)
    invisible(utils::capture.output(
      w_k <- bilinear_w_from_Y_resi(
        w_list_new = w_new, Y_resi = Yk,
        n_features = ncol(PCmats[[cts[1]]]), max_iter = maxIter, tol = tol
      )
    ))
  }

  # Ensure single-column matrices
  for (ct in cts) {
    w_k[[ct]] <- matrix(w_k[[ct]][, 1], ncol = 1)
  }

  ncorr <- .compute_ncorr_quick(PCmats, w_k, flat_kernels, sigma, cts)
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
#' Deflating `Y` by the fixed observed weight directions is algebraically
#' identical to the Freedman-Lane / ter Braak residualization that removes the
#' observed lower canonical *variates* from the data before recomputing the
#' cross-product, whenever the PCs are whitened (`scalePCs = TRUE`, so
#' `X^T X = c I`): `(I - u u^T) Y (I - v v^T) = Y - (u^T Y v)\, u v^T`. The
#' fair-sigma maximum over the bandwidth family is the Westfall-Young (1993)
#' maxT procedure.
#'
#' @param object A CoPro object with `runSkrCCA()` and
#'   `computeNormalizedCorrelation()` already run.
#' @param nPermu Number of permutations (default 100).
#' @param sigma_values Candidate bandwidths for the fair-sigma maximum. `NULL`
#'   uses all `object@sigmaValues` that have kernels and weights.
#' @param permu_method Permutation null: "bin" (default), "global", "pc", or
#'   "toroidal". See [runSkrCCAPermu()].
#' @param permu_which Which cell types to permute: "second_only" (default),
#'   "both", or "first_only".
#' @param num_bins_x,num_bins_y Bin grid for `permu_method = "bin"`. Default
#'   `NULL` is sigma-aware (see [.sigmaAwareBins()]); the same grid is shared
#'   across the sigma sweep. Pass integers to override.
#' @param match_quantile Whether to use quantile matching for bin permutation.
#' @param alpha Family-wise significance level for the step-down rule (default 0.05).
#' @param maxIter,tol Optimization controls passed to the axis optimizer.
#' @param verbose Whether to print progress and a summary (default TRUE).
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
                                       nPermu = 100,
                                       sigma_values = NULL,
                                       permu_method = "bin",
                                       permu_which = "second_only",
                                       num_bins_x = NULL,
                                       num_bins_y = NULL,
                                       match_quantile = FALSE,
                                       alpha = 0.05,
                                       maxIter = 200,
                                       tol = 1e-5,
                                       verbose = TRUE) {

  ## ---- validation ----
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }
  if (length(object@skrCCAOut) == 0) {
    stop("Please run runSkrCCA() first")
  }
  if (length(object@normalizedCorrelation) == 0) {
    stop("Please run computeNormalizedCorrelation() first")
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
  scalePCs <- object@scalePCs
  nCC <- object@nCC
  if (length(nCC) == 0 || nCC < 1) {
    stop("nCC must be >= 1; run runSkrCCA() with the desired number of axes.")
  }
  if (length(cts) == 1 && permu_which == "second_only") {
    warning("With a single cell type, permu_which = 'second_only' gives the ",
            "identity permutation (no shuffling). Use 'both' or 'first_only' ",
            "for a within-type permutation test.")
  }

  ## ---- candidate sigma values (fair-sigma sweep) ----
  if (is.null(sigma_values)) {
    sigma_values <- object@sigmaValues
  }
  sigma_values <- sigma_values[sigma_values %in% object@sigmaValues]
  if (length(sigma_values) == 0) {
    stop("No valid sigma values available for the conditional permutation test.")
  }
  sigma_names <- paste("sigma", sigma_values, sep = "_")

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
                              match_quantile = match_quantile)

  permute_pc_matrix <- function(pc_mat, seed) {
    set.seed(seed)
    permuted <- apply(pc_mat, 2, function(x) sample(x, length(x)))
    rownames(permuted) <- rownames(pc_mat)
    permuted
  }

  perm_stat <- matrix(NA_real_, nrow = nPermu, ncol = nCC)
  perm_sigma <- matrix(sigma_values[1], nrow = nPermu, ncol = nCC)
  n_failed <- 0L

  if (verbose) cat("Running permutations...\n")

  for (tt in seq_len(nPermu)) {
    ## build permuted PC matrices
    PCmats_local <- PCmats
    for (ct in names(PCmats_local)) {
      cp <- cell_permu[[ct]]
      if (is.list(cp) && !is.null(cp$type) && cp$type == "pc_permute") {
        PCmats_local[[ct]] <- permute_pc_matrix(PCmats[[ct]], cp$seeds[tt])
      } else {
        PCmats_local[[ct]] <- PCmats[[ct]][cp[, tt], ]
      }
    }

    res <- tryCatch({
      stat_k <- rep(-Inf, nCC)
      sig_k <- rep(sigma_values[1], nCC)
      for (si in seq_along(sigma_values)) {
        s <- sigma_values[si]
        sname <- sigma_names[si]
        ## compute Y once per (permutation, sigma); reuse across axes >= 2
        Y0 <- if (nCC >= 2) {
          compute_Y_resi(PCmats_local, object@kernelMatrices, s, cts,
                         slide = NULL)
        } else {
          NULL
        }
        for (k in seq_len(nCC)) {
          fit <- .fitConditionalAxis(
            PCmats = PCmats_local, flat_kernels = object@kernelMatrices,
            sigma = s, cts = cts, W_lower = obs_W[[sname]], k_minus_1 = k - 1,
            Y_resi = Y0, maxIter = maxIter, tol = tol
          )
          if (is.finite(fit$ncorr) && fit$ncorr > stat_k[k]) {
            stat_k[k] <- fit$ncorr
            sig_k[k] <- s
          }
        }
      }
      list(stat = stat_k, sigma = sig_k)
    }, error = function(e) NULL)

    if (is.null(res)) {
      n_failed <- n_failed + 1L
    } else {
      perm_stat[tt, ] <- res$stat
      perm_sigma[tt, ] <- res$sigma
    }

    if (verbose && (tt %% 10 == 0 || tt == nPermu)) {
      cat(sprintf("  Completed %d of %d permutations\n", tt, nPermu))
    }
  }

  if (n_failed > 0) {
    warning(sprintf(paste0("%d of %d permutations failed to optimize and were ",
                          "dropped from the null. P-values use the remaining ",
                          "permutations."), n_failed, nPermu))
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
    permu_method = permu_method,
    permu_which = permu_which,
    num_bins = num_bins_out
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
