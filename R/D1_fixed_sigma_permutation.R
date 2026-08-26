

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
#' @param compactPermutation Logical. Store the permuted side as one seed per
#'   draw instead of an explicit index matrix.
#'
#' @return List of permutation matrices (for "global"/"bin"/"toroidal") or
#'   list of permuted PC matrices (for "pc"), one per cell type. Held-fixed
#'   cell types are stored as the compact identity marker built by
#'   `.identityPermutation()` rather than as an explicit index matrix.
#' @keywords internal
.getCellPermu <- function(object, permu_method, nPermu, cts,
                          permu_which = "second_only",
                          num_bins_x = 10, num_bins_y = 10,
                          match_quantile = FALSE,
                          compactPermutation = .defaultCompactPermutation()) {

  cell_permu <- setNames(vector("list", length = length(cts)), cts)
  compact <- isTRUE(compactPermutation)

  # Determine which cell types should be permuted
  should_permute <- function(ct_name, ct_index) {
    # With one cell type there is no other type to hold fixed, so `permu_which`
    # has nothing to choose between and the only meaningful null relabels that
    # type's scores against its own locations. The literal "second_only" rule
    # (`ct_index > 1`) returns FALSE here, which makes every draw the identity:
    # the null then equals the observed statistic and the p-value is 1 by
    # construction, with nothing in the output to say so.
    if (length(cts) == 1L) {
      return(TRUE)
    }
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
        cell_permu[[ct]] <- if (compact) {
          list(type = "global_seed", n_cell = n_cell,
               seeds = sample.int(.Machine$integer.max, nPermu))
        } else {
          replicate(nPermu, sample.int(n = n_cell, replace = FALSE))
        }
      } else {
        # Keep this cell type fixed (identity permutation)
        cell_permu[[ct]] <- .identityPermutation(n_cell)
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

        # Bins, memberships, and neighboring-bin lookups are invariant across
        # permutations. Preparing them once avoids repeatedly scanning the full
        # data frame for every occupied bin in every permutation.
        prepared <- .prepareSpatialResampling(
          location_data = cell_loc,
          num_bins_x = num_bins_x,
          num_bins_y = num_bins_y
        )

        if (compact) {
          # `prepared` is O(n); the draws it generates are O(n * nPermu). Keep
          # the former and regenerate the latter from a per-draw seed.
          cell_permu[[ct]] <- list(
            type = "bin_seed", n_cell = n_cell,
            seeds = sample.int(.Machine$integer.max, nPermu),
            prepared = prepared, match_quantile = match_quantile
          )
        } else {
          cell_permu[[ct]] <- matrix(ncol = nPermu, nrow = nrow(cell_loc))
          for (j in seq_len(nPermu)) {
            cell_permu[[ct]][, j] <- .drawSpatialPermutation(
              prepared, match_quantile = match_quantile
            )
          }
        }
      } else {
        # Keep this cell type fixed (identity permutation)
        cell_permu[[ct]] <- .identityPermutation(n_cell)
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
        cell_permu[[ct]] <- .identityPermutation(n_cell)
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
        cell_permu[[ct]] <- .identityPermutation(n_cell)
      }
    }
  } else {
    stop(paste("permu_method must be 'global', 'bin', 'toroidal', or 'pc'.",
               "Got:", permu_method))
  }

  return(cell_permu)
}

.permutationSdev2 <- function(object, cts) {
  if (isTRUE(object@scalePCs)) return(NULL)
  stats::setNames(lapply(cts, function(ct) {
    sdev <- object@pcaGlobal[[ct]]$sdev
    if (is.null(sdev) || length(sdev) == 0L) {
      stop("PCA standard deviations are missing for cell type '", ct, "'.")
    }
    sdev^2
  }), cts)
}

.permutationProvenance <- function(object) {
  provenance <- attr(object, "permutationProvenance", exact = TRUE)
  if (!is.null(provenance)) return(provenance)

  # Backward-compatible inference for objects created before provenance was
  # recorded. Every stored null row still identifies the sigma it tested.
  null_sigmas <- unique(unlist(lapply(object@normalizedCorrelationPermu, function(x) {
    sigma_col <- intersect(c("sigmaValues", "sigmaValue"), names(x))
    if (length(sigma_col) == 0L) return(numeric())
    as.numeric(x[[sigma_col[1L]]])
  })))
  list(
    method = "legacy",
    sigma_values = null_sigmas,
    sigma_aggregation = if (length(null_sigmas) > 1L) "max" else "fixed",
    pair_aggregation = "max",
    sigma_predeclared = NA
  )
}

.rejectCellPermutationForMulti <- function(object) {
  if (inherits(object, "CoProMulti")) {
    stop(
      "Cell-level permutation is not replicate-level inference for CoProMulti. ",
      "Use runSlideLevelInference() to obtain leave-one-replicate-out effects, ",
      "a replicate sign-flip p-value, and a replicate bootstrap interval."
    )
  }
  invisible(NULL)
}

#' Decide which criterion the permutation null must be built with
#'
#' A permutation p-value is only meaningful when the null draws are optimized by
#' the same criterion as the observed statistic. Every caller of this function
#' has already passed `.rejectCellPermutationForMulti()`, so the object has a
#' single slide and the SUMCOR reduction test applies in full:
#'
#' * Recorded as `sumcov` (or unrecorded, which pre-dates the attribute) --
#'   nothing to do.
#' * Recorded as `sumcor` where `.sumcorReducesToSumcov()` holds -- the fitted
#'   weights *are* the SUMCOV maximizer, so the existing SUMCOV null is already
#'   the matching null. This covers one or two cell types at any counts, and
#'   three or more at equal counts, i.e. nearly every real call.
#' * Recorded as `sumcor` where it does not hold -- three or more cell types at
#'   unequal counts. Here the criteria genuinely differ and the null draws are
#'   re-optimized under SUMCOR.
#'
#' A *bijective* within-slide label permutation permutes the rows of \eqn{X_i},
#' which leaves \eqn{G_i = X_i' X_i} unchanged, so the per-slide scales SUMCOR
#' divides by can be built once and reused. That covers `"global"` and
#' `"toroidal"` but **not** the default `"bin"` null, which resamples cells
#' rather than permuting them, nor `"pc"`, which shuffles each PC column
#' independently. The Grams returned here are therefore the unpermuted baseline;
#' `runSkrCCAPermu()` runs them through `.permutationGrams()` to find out which
#' a draw may inherit, and `.drawGrams()` recomputes the rest per draw. The `Y`
#' operator-reuse factorization is untouched either way.
#'
#' @param object A single-slide `CoPro` object with `@skrCCAOut` populated.
#' @param cts Cell types of interest.
#' @param scalePCs Whether the PC scores are whitened, matching the caller.
#' @param supports_sumcor Whether the calling test can re-optimize its draws
#'   under SUMCOR. `FALSE` for the fair-sigma and conditional tests, which are
#'   restricted to at most two cell types and so can only ever reach the
#'   reducible case; the check is defensive against that restriction being
#'   relaxed without revisiting how their null is built.
#' @return `list(objective, slideWeight, grams, n_cells)`. `grams`/`n_cells` are
#'   `NULL` unless a SUMCOR null is actually required.
#' @noRd
.resolvePermutationObjective <- function(object, cts, scalePCs = TRUE,
                                         supports_sumcor = TRUE,
                                         verbose = TRUE) {
  sumcov <- list(objective = "sumcov", slideWeight = NULL,
                 grams = NULL, n_cells = NULL)
  record <- attr(object@skrCCAOut, "ccaObjective")
  if (is.null(record) || !identical(record$objective, "sumcor")) return(sumcov)

  # Everything below reads @pcaGlobal and re-optimizes with the PC-space
  # solvers. A gene-space fit lives in a different feature space entirely, so
  # its weights cannot be tested here whatever criterion they were fitted
  # under. Unreachable today -- runGeneSpaceCCA() requires CoProMulti and every
  # caller has already refused CoProMulti -- but the two guards are independent
  # and this one should not depend on that.
  if (identical(record$space, "gene")) {
    stop(
      "These weights were fitted in gene space (space = \"gene\"), but the ",
      "cell-level permutation routines build their null from the PCA-space ",
      "operators. Comparing the two would mix feature spaces. Use ",
      "runSlideLevelInference() for replicate-level inference, or re-fit with ",
      "runSkrCCA(space = \"pca\") before testing."
    )
  }

  slideWeight <- record$slideWeight
  if (is.null(slideWeight)) slideWeight <- "equal"

  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)[cts]
  grams <- stats::setNames(lapply(cts, function(ct) {
    crossprod(PCmats[[ct]])
  }), cts)
  n_cells <- stats::setNames(
    vapply(cts, function(ct) nrow(PCmats[[ct]]), integer(1)), cts
  )
  ops <- .singleSlidePermutationOps(NULL, grams, n_cells, cts)

  # The reduction is a property of the whitened operators, so test the same
  # matrices the optimizer would see under this object's scalePCs setting.
  ops_w <- .whitenSlideOperators(ops, .permutationSdev2(object, cts))
  if (.sumcorReducesToSumcov(ops_w, slideWeight)) {
    if (verbose) message(
      "Weights were fitted under objective = \"sumcor\", but with ",
      length(cts), " cell type(s) on one slide that is the same optimization ",
      "problem as \"sumcov\". The existing null is the matching null."
    )
    return(sumcov)
  }

  if (!supports_sumcor) {
    stop(
      "These weights were fitted under objective = \"sumcor\" with ",
      length(cts), " cell types at unequal counts, where sumcor and sumcov ",
      "have different maximizers. This test builds its null with the sumcov ",
      "solvers, and comparing the two would mix criteria. Re-fit with ",
      "runSkrCCA(objective = \"sumcov\"), or use runSkrCCAPermu(), whose null ",
      "is re-optimized under sumcor."
    )
  }

  list(objective = "sumcor", slideWeight = slideWeight,
       grams = grams, n_cells = n_cells)
}

#' Wrap a permuted PC-space operator in the one-slide SUMCOR operator structure
#'
#' `Y_resi` already has the `[[ct_i]][[ct_j]]` shape the SUMCOR routines expect;
#' only the slide level and the permutation-invariant Gram matrices need adding.
#' @noRd
.singleSlidePermutationOps <- function(Y_resi, grams, n_cells, cts) {
  slide <- .SINGLE_SLIDE_TOKEN
  if (is.null(Y_resi)) {
    Y_resi <- stats::setNames(lapply(cts, function(i) {
      stats::setNames(vector("list", length(cts)), cts)
    }), cts)
  }
  list(
    Y = stats::setNames(list(Y_resi), slide),
    G = stats::setNames(list(grams), slide),
    n = stats::setNames(list(n_cells), slide),
    slides = slide,
    cell_types = cts
  )
}

#' Optimize one permuted draw under SUMCOR
#'
#' Only ever reached with three or more cell types, because fewer than that
#' reduces to SUMCOV and never gets here. `prev_axes` carries the conditional
#' test's fixed lower directions: SUMCOR's operator depends on `w` through
#' `sigma`, so there is no operator to deflate and orthogonality is imposed in
#' weight space, exactly as `optimize_sumcor_pca_n()` does.
#' @noRd
.fitSumcorPermutedAxes <- function(Y_resi, grams, n_cells, cts, nCC,
                                   sdev2_list, slideWeight, prev_axes = NULL,
                                   maxIter = 200, tol = 1e-6, step_size = 1) {
  ops <- .singleSlidePermutationOps(Y_resi, grams, n_cells, cts)
  ops_w <- .whitenSlideOperators(ops, sdev2_list)
  prev_w <- if (is.null(prev_axes)) NULL else {
    .whitenWeights(prev_axes, sdev2_list)
  }

  warm <- .sumcorWarmStart(
    ops_w, cts, NULL, nCC = 1L, step_size = step_size,
    max_iter = maxIter, tol = tol
  )
  fit <- suppressWarnings(.sumcorIterate(
    w_init = stats::setNames(lapply(cts, function(ct) {
      warm[[ct]][, 1L, drop = FALSE]
    }), cts),
    ops = ops_w, slideWeight = slideWeight, sdev2_list = NULL,
    prev_axes = prev_w, max_iter = maxIter, tol = tol,
    step_size = step_size, label = "permutation CC 1"
  ))
  w_list <- fit$w_list

  if (nCC > 1L) {
    for (cc in seq(2L, nCC)) {
      axes <- stats::setNames(lapply(cts, function(ct) {
        w_list[[ct]][, seq_len(cc - 1L), drop = FALSE]
      }), cts)
      if (!is.null(prev_w)) {
        axes <- stats::setNames(lapply(cts, function(ct) {
          cbind(prev_w[[ct]], axes[[ct]])
        }), cts)
      }
      warm_cc <- .sumcorWarmStart(
        ops_w, cts, NULL, nCC = cc, step_size = step_size,
        max_iter = maxIter, tol = tol
      )
      fit_cc <- suppressWarnings(.sumcorIterate(
        w_init = stats::setNames(lapply(cts, function(ct) {
          warm_cc[[ct]][, cc, drop = FALSE]
        }), cts),
        ops = ops_w, slideWeight = slideWeight, sdev2_list = NULL,
        prev_axes = axes, max_iter = maxIter, tol = tol,
        step_size = step_size, label = sprintf("permutation CC %d", cc)
      ))
      for (ct in cts) {
        w_list[[ct]] <- cbind(w_list[[ct]], fit_cc$w_list[[ct]])
      }
    }
  }

  .unwhitenWeights(w_list, sdev2_list)
}

# Resolve the library directory that holds an *installed* CoPro, or NULL when
# CoPro is only available via devtools::load_all(). PSOCK workers start fresh R
# sessions and must library(CoPro); under load_all() find.package() points at the
# source tree (which has no Meta/package.rds and cannot be loaded in a worker).
# Callers fall back to sequential execution when this returns NULL.
.installedCoProLibrary <- function() {
  pkg_dir <- tryCatch(find.package("CoPro"), error = function(e) NULL)
  if (is.null(pkg_dir) ||
      !file.exists(file.path(pkg_dir, "Meta", "package.rds"))) {
    return(NULL)
  }
  dirname(pkg_dir)
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
#' With a **single** cell type none of these apply -- there is no other type to
#' hold fixed -- so `permu_which` is ignored and that type is permuted. The
#' resulting null is the within-type one: scores are relabelled against their
#' own locations, breaking the association between a cell's score and where it
#' sits while leaving the spatial configuration and the score distribution
#' intact.
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
#' @param nPermu Number of permutations to run (default: 999), giving a
#'   Monte-Carlo floor of 0.001 under the Phipson--Smyth correction.
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
#'   Ignored when the object has a single cell type: there is nothing to hold
#'   fixed, so that type is permuted and the null is the within-type one.
#' @param num_bins_x Number of bins in x direction for bin-wise permutation.
#'   Default `NULL` sizes the grid automatically from the kernel bandwidth so
#'   each patch is ~2*sigma wide on the normalized distance scale (see
#'   [.sigmaAwareBins()]). When that bandwidth is `object@@sigmaValueChoice`,
#'   the null grid is mildly data-adaptive (a second-order circularity). Pass an
#'   explicit integer, or a predeclared `sigma`, to avoid it.
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
#' @param sigma Optional predeclared sigma value. When supplied, both observed
#'   and null statistics are evaluated at this fixed sigma. When `NULL`, the
#'   selected `object@@sigmaValueChoice` is used conditionally and the returned
#'   p-value is marked as not adjusted for sigma selection; use
#'   [runSkrCCAPermu_FairSigma()] for a max-over-sigma test.
#' @param factorize Apply the fixed-side operator factorization (default: TRUE).
#'   A cell type held fixed across every draw lets its side of `X' K X` be
#'   multiplied by the kernel once instead of once per draw, and lets the score
#'   norms be read off a cached Gram matrix. The identity is exact, so this
#'   changes speed and memory only, never a p-value. Set `FALSE` to route every
#'   pair through the original sparse product when comparing the two paths.
#'   Defaults to `getOption("CoPro.factorizePermutation", TRUE)`, the global
#'   flag this argument replaced.
#' @param compactPermutation Store the permuted side as one seed per draw rather
#'   than an explicit index matrix (default: FALSE). Saves the `n * nPermu * 4`
#'   bytes per permuted cell type that an index matrix costs, but it changes
#'   *which* permutations are drawn, so re-running a saved analysis moves its
#'   p-values within Monte Carlo error. Held-fixed types are always stored
#'   compactly, which changes no number at all. Defaults to
#'   `getOption("CoPro.compactPermutation", FALSE)`, the global flag this
#'   argument replaced.
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
#' p_value <- (1 + sum(permu_values >= observed)) /
#'   (1 + length(permu_values))
#' }
#'
#' @importFrom stats setNames
#' @export
runSkrCCAPermu <- function(object, tol = 1e-5, nPermu = 999,
                           maxIter = 200, permu_method = "bin",
                           permu_which = "second_only",
                           num_bins_x = NULL, num_bins_y = NULL,
                           match_quantile = FALSE,
                           conservative = FALSE,
                           n_cores = 1, verbose = TRUE,
                           sigma = NULL,
                           factorize = .defaultFactorizePermutation(),
                           compactPermutation = .defaultCompactPermutation()) {

  ## Input validation
  if (!is(object, "CoPro")) {
    stop("Input object must be a CoPro object")
  }
  .rejectCellPermutationForMulti(object)

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

  nPermu <- .checkPositiveScalarInt(nPermu, "nPermu")
  if (nPermu < 2L) {
    stop("nPermu must be at least 2; one draw produces a degenerate null.")
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

  ## The null must be optimized by the same criterion as the observed weights.
  permu_objective <- .resolvePermutationObjective(
    object, cts, scalePCs, verbose = verbose
  )

  if (length(cts) == 1L && verbose) {
    message("Single cell type: permu_which is ignored and the one type is ",
            "permuted, giving a within-type null that relabels scores against ",
            "their own locations.")
  }

  ## Get sigma value
  sigma_predeclared <- !is.null(sigma)
  sigmaValueChoice <- if (sigma_predeclared) sigma else object@sigmaValueChoice
  if (is.null(sigmaValueChoice) || length(sigmaValueChoice) == 0) {
    stop(paste("sigmaValueChoice not set.",
               "Please run computeNormalizedCorrelation() first."))
  }

  if (!(sigmaValueChoice %in% object@sigmaValues)) {
    stop("sigmaValueChoice does not exist in the list of sigmaValues")
  }

  s_name <- .sigmaName(sigmaValueChoice)

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
  permu_names <- paste("permu", 1:nPermu, sep = "_")

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
                              match_quantile = match_quantile,
                              compactPermutation = compactPermutation)
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
  sdev2_list <- .permutationSdev2(object, cts)

  ## Step 2: Run CCA for each permutation
  if (verbose) {
    cat(paste("Running CCA optimization for", nPermu, "permutations"))
    if (n_cores > 1) {
      cat(paste(" using", n_cores, "cores...\n"))
    } else {
      cat("...\n")
    }
  }


  # Cell types held fixed across every draw let their side of X' K X be
  # multiplied by the kernel once, instead of once per draw. See
  # R/D0_permutation_plan.R for the identity, its cost, and the escape hatch.
  plan <- .buildYPlan(
    PCmats = PCmats,
    flat_kernels = object@kernelMatrices,
    sigma = sigmaValueChoice,
    cts = cts,
    fixed = .fixedPermutationTypes(cell_permu, cts),
    factorize = factorize
  )

  # A SUMCOR null divides by per-type score scales, so it needs the Gram
  # matrices of the data each draw actually fits. Only a bijection on cells
  # leaves them where they were; the same test that decides which score norms
  # can be factorized decides which Grams a draw may inherit, and the worker
  # recomputes the rest.
  permu_grams <- if (identical(permu_objective$objective, "sumcor")) {
    .permutationGrams(PCmats[cts], cell_permu, cts, factorize = factorize)
  } else {
    NULL
  }

  worker <- .makeSkrCCAPermuWorker(
    PCmats = PCmats, plan = plan, cts = cts, nCC = nCC,
    sdev2_list = sdev2_list, maxIter = maxIter, tol = tol,
    permu_objective = permu_objective, permu_grams = permu_grams
  )

  cca_permu_out <- .runPermutationDraws(
    cell_permu = cell_permu, cts = cts, nPermu = nPermu, worker = worker,
    n_cores = n_cores, verbose = verbose, progress_every = 10L
  )
  names(cca_permu_out) <- permu_names

  object@skrCCAPermuOut <- cca_permu_out
  attr(object, "permutationProvenance") <- list(
    method = "fixed_sigma",
    sigma_values = as.numeric(sigmaValueChoice),
    sigma_aggregation = "fixed",
    pair_aggregation = "max",
    sigma_predeclared = sigma_predeclared,
    selection_adjusted = sigma_predeclared || length(object@sigmaValues) == 1L,
    scalePCs = scalePCs
  )

  if (verbose) {
    cat("\nPermutation testing complete.\n")
    cat("Run computeNormalizedCorrelationPermu() to compute p-values.\n")
  }

  return(object)
}
