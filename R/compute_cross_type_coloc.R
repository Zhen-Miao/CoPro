# Colocalization score for two cell types on a slide
# -------------------------------------------------
# Single-number per-slide statistic based on the inhomogeneous
# bivariate pair correlation function g12(r), standardized against
# a random-labelling null, and averaged over a biologically relevant
# distance band [r_min, r_max].
#
# CES_z = mean_r ( g12_obs(r) - mean_null(r) ) / sd_null(r )
# Positive values => co-localization; negative => segregation.
#
# Requires: spatstat.geom, spatstat.core, spatstat.explore

# Global variables for spatstat functions to avoid linting warnings
utils::globalVariables(c("getSlideList", "density.ppp", "pcfcross.inhom", "rlabel", 
                         "nndist", "nncross", "owin", "ppp", "superimpose", "marks"))

# Note: spatstat packages should be in Suggests and used with :: notation
# Check availability at runtime since these are optional dependencies

# ---------------------------
# Helper: compute g12_inhom
# ---------------------------
.compute_g12_inhom <- function(X, type_i, type_j, r_vec, sigma_intensity, edge_correction = "translation") {
  # Estimate inhomogeneous intensities for each type (using marks)
  # We'll fit separate intensity surfaces via kernel smoothing; sigma_intensity in the same units as coords.
  # Build intensity weights via rhohat? Simpler: use density.ppp per type and evaluate at point locations.

  # Extract subpatterns
  Xi <- X[X$marks == type_i]
  Xj <- X[X$marks == type_j]

  # Intensity estimates on a pixel grid then interpolate automatically inside pcfcross.inhom via 'lambdaI' and 'lambdaJ'
  lambdaI <- spatstat.explore::density.ppp(Xi, sigma = sigma_intensity, at = "pixels", edge = TRUE)
  lambdaJ <- spatstat.explore::density.ppp(Xj, sigma = sigma_intensity, at = "pixels", edge = TRUE)

  # pcfcross.inhom uses intensity pixel images for both types
  g <- spatstat.explore::pcfcross.inhom(X, i = type_i, j = type_j,
                      lambdaI = lambdaI, lambdaJ = lambdaJ,
                      r = r_vec, correction = edge_correction)
  # Return numeric vector of the "trans" (translation) correction if present else "iso"
  valname <- if ("trans" %in% colnames(g)) "trans" else if ("iso" %in% colnames(g)) "iso" else colnames(g)[2]
  as.numeric(g[[valname]])
}

# ----------------------------------------------
# Main: compute standardized colocalization score
# ----------------------------------------------
coloc_score <- function(
  A, B,                    # data.frames with columns x, y (in microns or pixels; see `units`)
  window_range,            # list(xrange = c(xmin,xmax), yrange = c(ymin,ymax)) in same units as A/B
  r_um_range = c(10, 60),  # distance band in microns over which to summarize co-location
  pixel_size_um = 1,       # microns per coordinate unit (set to actual pixel size if A/B are in pixels)
  cell_diam_um = 10,       # typical cell diameter (microns) to set smoothing
  nsim = 199,              # number of random-labelling simulations
  r_step_um = 2,           # step (microns) for r grid
  min_points_per_type = 20 # guard for low counts
) {
  stopifnot(all(c("x","y") %in% names(A)), all(c("x","y") %in% names(B)))

  # Convert inputs to microns (standardized units)
  s <- 1 / pixel_size_um  # Conversion factor: coordinate_units -> microns
  Au <- data.frame(x = A$x * s, y = A$y * s)
  Bu <- data.frame(x = B$x * s, y = B$y * s)

  # Window in microns (same units as converted coordinates)
  win <- spatstat.geom::owin(xrange = window_range$xrange * s, yrange = window_range$yrange * s)

  # Build multitype pattern
  X <- spatstat.geom::superimpose(
    A = spatstat.geom::ppp(Au$x, Au$y, window = win, check = FALSE),
    B = spatstat.geom::ppp(Bu$x, Bu$y, window = win, check = FALSE),
    W = win, check = FALSE
  )
  spatstat.geom::marks(X) <- factor(spatstat.geom::marks(X), levels = c("A","B"))

  nA <- sum(spatstat.geom::marks(X) == "A"); nB <- sum(spatstat.geom::marks(X) == "B")
  if (nA < min_points_per_type || nB < min_points_per_type) {
    warning(sprintf("Too few points (A=%d, B=%d). Returning NA.", nA, nB))
    return(list(
      score = NA_real_,
      r_um = numeric(0),
      z_of_r = numeric(0),
      g_obs = numeric(0),
      g_mean = numeric(0),
      g_sd = numeric(0),
      meta = list(nA = nA, nB = nB, note = "insufficient points")
    ))
  }

  # r grid in microns (after coordinate conversion, both coordinates and r are in microns)
  # spatstat requires first r value to be 0, so we start from 0 and filter later
  r_vec_um <- seq(0, r_um_range[2], by = r_step_um)
  r_vec <- r_vec_um  # No scaling needed - coordinates are now in microns, r should be too

  # Intensity smoothing bandwidth in microns (same units as coordinates after conversion)
  sigma_intensity <- (3 * cell_diam_um)

  # Observed g12_inhom
  g_obs <- .compute_g12_inhom(X, "A", "B", r_vec, sigma_intensity)

  # Simulate random labelling nsim times, compute g each time
  g_mat <- matrix(NA_real_, nrow = length(r_vec), ncol = nsim)
  set.seed(1)
  for (k in seq_len(nsim)) {
    X_sim <- spatstat.random::rlabel(X)  # keeps locations, shuffles marks
    # Ensure mark levels preserved
    spatstat.geom::marks(X_sim) <- factor(spatstat.geom::marks(X_sim), levels = c("A","B"))
    g_mat[, k] <- .compute_g12_inhom(X_sim, "A", "B", r_vec, sigma_intensity)
  }

  g_mean <- rowMeans(g_mat, na.rm = TRUE)
  g_sd   <- apply(g_mat, 1, sd, na.rm = TRUE)

  # Filter to the analysis distance band (exclude r values below r_um_range[1])
  analysis_indices <- r_vec_um >= r_um_range[1] & r_vec_um <= r_um_range[2]
  
  if (sum(analysis_indices) == 0) {
    # No r values in the analysis range
    return(list(
      score = NA_real_,
      r_um = r_vec_um,
      z_of_r = rep(NA_real_, length(r_vec_um)),
      g_obs = g_obs,
      g_mean = g_mean,
      g_sd = g_sd,
      meta = list(nA = nA, nB = nB, note = "no r values in analysis range")
    ))
  }
  
  # Z-scores per r; guard against zero sd
  g_sd[g_sd == 0 | is.na(g_sd)] <- Inf
  z_of_r <- (g_obs - g_mean) / g_sd

  # Single-number score: average z across the analysis band only
  ces_z <- mean(z_of_r[analysis_indices & is.finite(z_of_r)])

  list(
    score = ces_z,              # higher => stronger co-localization
    r_um = r_vec_um,
    z_of_r = z_of_r,
    g_obs = g_obs,
    g_mean = g_mean,
    g_sd = g_sd,
    meta = list(
      nA = nA, nB = nB,
      r_um_range = r_um_range,
      r_step_um = r_step_um,
      pixel_size_um = pixel_size_um,
      sigma_intensity_um = 3 * cell_diam_um,
      nsim = nsim,
      window_range = window_range
    )
  )
}

# ---------------------------
# Convenience: quick wrapper
# ---------------------------
# Example usage:
# A <- data.frame(x = runif(120, 0, 500), y = runif(120, 0, 400))
# B <- data.frame(x = runif(150, 0, 500), y = runif(150, 0, 400))
# window_range <- list(xrange = c(0, 500), yrange = c(0, 400))
# out <- coloc_score(A, B, window_range, r_um_range = c(10, 60),
#                    pixel_size_um = 1, cell_diam_um = 10, nsim = 199)
# out$score  # single-number per slide

# ---------------------------
# Optional: nearest-neighbor fallback (robust when one type is sparse)
# ---------------------------
# Returns proportion of cross-type nearest neighbors averaged both ways.
# Bounded in [0,1]; ~0.5 under random labelling (balanced abundances).

nn_cross_fraction <- function(A, B, window_range, pixel_size_um = 1) {
  s <- 1 / pixel_size_um  # Conversion factor: coordinate_units -> microns
  Au <- data.frame(x = A$x * s, y = A$y * s)
  Bu <- data.frame(x = B$x * s, y = B$y * s)
  win <- spatstat.geom::owin(xrange = window_range$xrange * s, yrange = window_range$yrange * s)
  X <- spatstat.geom::superimpose(
    A = spatstat.geom::ppp(Au$x, Au$y, window = win, check = FALSE),
    B = spatstat.geom::ppp(Bu$x, Bu$y, window = win, check = FALSE),
    W = win, check = FALSE
  )
  spatstat.geom::marks(X) <- factor(spatstat.geom::marks(X), levels = c("A","B"))

  # Compare distance to nearest *cross-type* vs *same-type* neighbor.
  dA_B <- spatstat.geom::nncross(X[X$marks == "A"], X[X$marks == "B"], what = "dist")$dist
  dA_A <- spatstat.geom::nndist(X[X$marks == "A"], what = "dist")
  dB_A <- spatstat.geom::nncross(X[X$marks == "B"], X[X$marks == "A"], what = "dist")$dist
  dB_B <- spatstat.geom::nndist(X[X$marks == "B"], what = "dist")

  pA <- mean(dA_B < dA_A)
  pB <- mean(dB_A < dB_B)
  mean(c(pA, pB))
}

#' Compute Colocalization Scores for All Cell Type Pairs
#'
#' This function computes cross-type colocalization scores for all pairs of cell types
#' in a CoPro object using the inhomogeneous cross pair-correlation function approach.
#' The colocalization score is a normalized measure based on the standardized difference
#' between observed and expected (under random labelling) spatial correlation.
#'
#' @param object A `CoProSingle` or `CoProMulti` object with location data and cell types.
#' @param r_um_range Numeric vector of length 2 specifying the distance band in microns
#'   over which to summarize colocalization (default: c(10, 60)).
#' @param pixel_size_um Numeric scalar specifying microns per coordinate unit. 
#'   This converts your input coordinates to microns for analysis. For example:
#'   - If 1 coordinate unit = 50 microns, set pixel_size_um = 50
#'   - If 1 coordinate unit = 0.325 microns, set pixel_size_um = 0.325
#'   - If coordinates are already in microns, set pixel_size_um = 1 (default)
#' @param cell_diam_um Numeric scalar specifying typical cell diameter in microns, used
#'   for intensity smoothing bandwidth (default: 10).
#' @param nsim Integer specifying number of random-labelling simulations for null
#'   distribution (default: 199).
#' @param r_step_um Numeric scalar specifying step size in microns for the r grid
#'   (default: 2).
#' @param min_points_per_type Integer specifying minimum number of points required per
#'   cell type to compute colocalization (default: 20).
#' @param edge_correction Character specifying edge correction method for pair correlation
#'   function (default: "translation").
#' @param verbose Logical specifying whether to print progress messages (default: TRUE).
#' @param include_self Logical specifying whether to include self-type pairs (same cell type).
#'   Default FALSE since colocalization is typically measured between different cell types.
#'
#' @return For `CoProSingle` objects: A data.frame with columns:
#'   - `cellType1`: First cell type name
#'   - `cellType2`: Second cell type name  
#'   - `colocScore`: Colocalization score (higher = more colocalized)
#'   - `nCells1`: Number of cells of type 1
#'   - `nCells2`: Number of cells of type 2
#'   - `r_um_min`: Minimum distance in analysis band
#'   - `r_um_max`: Maximum distance in analysis band
#'   
#'   For `CoProMulti` objects: Same columns plus:
#'   - `slideID`: Slide identifier
#'
#' @details The colocalization score is computed as:
#'   \deqn{CES_z = mean_r((g12_obs(r) - mean_null(r)) / sd_null(r))}
#'   
#'   Where:
#'   \itemize{
#'   - \eqn{g12_obs(r)} is the observed inhomogeneous cross pair-correlation function
#'   - \eqn{mean_null(r)} and \eqn{sd_null(r)} are the mean and standard deviation from random labelling simulations
#'   - The average is taken over the distance band \eqn{[r_um_range[1], r_um_range[2]]}
#'   }
#'   
#'   Positive scores indicate colocalization, negative scores indicate segregation.
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' coloc_results <- getColocScores(object)
#' 
#' # Custom parameters for high-resolution data
#' coloc_results <- getColocScores(object, 
#'                                r_um_range = c(5, 30),
#'                                pixel_size_um = 0.325,  # Convert from pixels to microns
#'                                cell_diam_um = 8,
#'                                nsim = 99)
#'                                
#' # For multi-slide objects
#' coloc_results <- getColocScores(multi_object)
#' # Results will include slideID column
#' }
#'
#' @export
#' @importFrom utils combn
getColocScores <- function(object,
                          r_um_range = c(10, 60),
                          pixel_size_um = 1,
                          cell_diam_um = 10,
                          nsim = 199,
                          r_step_um = 2,
                          min_points_per_type = 20,
                          edge_correction = "translation",
                          verbose = TRUE,
                          include_self = FALSE) {
  
  # Input validation
  if (!(is(object, "CoProSingle") || is(object, "CoProMulti"))) {
    stop("object must be a CoProSingle or CoProMulti object")
  }
  
  if (length(r_um_range) != 2 || r_um_range[1] >= r_um_range[2]) {
    stop("r_um_range must be a vector of length 2 with r_um_range[1] < r_um_range[2]")
  }
  
  if (pixel_size_um <= 0) {
    stop("pixel_size_um must be positive")
  }
  
  if (cell_diam_um <= 0) {
    stop("cell_diam_um must be positive")
  }
  
  if (nsim < 1) {
    stop("nsim must be at least 1")
  }
  
  if (min_points_per_type < 1) {
    stop("min_points_per_type must be at least 1")
  }
  
  # Check required data
  if (nrow(object@locationDataSub) == 0) {
    stop("No location data found. Ensure the object has been properly initialized with location data.")
  }
  
  if (length(object@cellTypesSub) == 0) {
    stop("No cell type data found. Ensure the object has been properly initialized with cell types.")
  }
  
  if (!all(c("x", "y") %in% colnames(object@locationDataSub))) {
    stop("Location data must contain 'x' and 'y' columns")
  }
  
  # Get cell types
  cts <- unique(object@cellTypesSub)
  if (length(cts) < 2 && !include_self) {
    stop("Need at least 2 cell types for colocalization analysis, or set include_self = TRUE")
  }
  
  is_multi <- is(object, "CoProMulti")
  
  if (is_multi) {
    return(.getColocScoresMulti(object, r_um_range, pixel_size_um, cell_diam_um, 
                               nsim, r_step_um, min_points_per_type, edge_correction,
                               verbose, include_self))
  } else {
    return(.getColocScoresSingle(object, r_um_range, pixel_size_um, cell_diam_um,
                                nsim, r_step_um, min_points_per_type, edge_correction,
                                verbose, include_self))
  }
}

#' Compute colocalization scores for single-slide object
#' @noRd
.getColocScoresSingle <- function(object, r_um_range, pixel_size_um, cell_diam_um,
                                 nsim, r_step_um, min_points_per_type, edge_correction,
                                 verbose, include_self) {
  
  cts <- unique(object@cellTypesSub)
  
  # Generate cell type pairs
  if (include_self) {
    # Include all pairs including self-pairs
    all_pairs <- expand.grid(cellType1 = cts, cellType2 = cts, stringsAsFactors = FALSE)
  } else {
    # Only cross-type pairs
    if (length(cts) < 2) {
      stop("Need at least 2 cell types for cross-type colocalization analysis")
    }
    pair_matrix <- combn(cts, 2)
    all_pairs <- data.frame(
      cellType1 = pair_matrix[1, ],
      cellType2 = pair_matrix[2, ],
      stringsAsFactors = FALSE
    )
  }
  
  # Determine window range from data
  window_range <- list(
    xrange = range(object@locationDataSub$x, na.rm = TRUE),
    yrange = range(object@locationDataSub$y, na.rm = TRUE)
  )
  
  if (verbose) {
    cat("Computing colocalization scores for", nrow(all_pairs), "cell type pairs\n")
    cat("Window range: x =", sprintf("%.2f", window_range$xrange), 
        ", y =", sprintf("%.2f", window_range$yrange), "\n")
    cat("Distance band:", r_um_range[1], "-", r_um_range[2], "microns\n")
  }
  
  # Initialize results data frame
  results <- data.frame(
    cellType1 = character(nrow(all_pairs)),
    cellType2 = character(nrow(all_pairs)),
    colocScore = numeric(nrow(all_pairs)),
    nCells1 = integer(nrow(all_pairs)),
    nCells2 = integer(nrow(all_pairs)),
    r_um_min = numeric(nrow(all_pairs)),
    r_um_max = numeric(nrow(all_pairs)),
    stringsAsFactors = FALSE
  )
  
  # Compute colocalization for each pair
  for (i in seq_len(nrow(all_pairs))) {
    ct1 <- all_pairs$cellType1[i]
    ct2 <- all_pairs$cellType2[i]
    
    if (verbose) {
      cat("Processing pair", i, "of", nrow(all_pairs), ":", ct1, "vs", ct2, "\n")
    }
    
    # Extract cell locations for each type
    cells1_idx <- object@cellTypesSub == ct1
    cells2_idx <- object@cellTypesSub == ct2
    
    A <- object@locationDataSub[cells1_idx, c("x", "y"), drop = FALSE]
    B <- object@locationDataSub[cells2_idx, c("x", "y"), drop = FALSE]
    
    # Compute colocalization score
    coloc_result <- tryCatch({
      coloc_score(A = A, B = B, 
                 window_range = window_range,
                 r_um_range = r_um_range,
                 pixel_size_um = pixel_size_um,
                 cell_diam_um = cell_diam_um,
                 nsim = nsim,
                 r_step_um = r_step_um,
                 min_points_per_type = min_points_per_type)
    }, error = function(e) {
      warning(paste("Error computing colocalization for", ct1, "vs", ct2, ":", e$message))
      list(score = NA_real_, meta = list(nA = nrow(A), nB = nrow(B)))
    })
    
    # Store results
    results$cellType1[i] <- ct1
    results$cellType2[i] <- ct2
    results$colocScore[i] <- coloc_result$score
    results$nCells1[i] <- coloc_result$meta$nA
    results$nCells2[i] <- coloc_result$meta$nB
    results$r_um_min[i] <- r_um_range[1]
    results$r_um_max[i] <- r_um_range[2]
  }
  
  if (verbose) {
    valid_scores <- !is.na(results$colocScore)
    cat("Completed:", sum(valid_scores), "valid scores,", sum(!valid_scores), "failed\n")
    if (sum(valid_scores) > 0) {
      cat("Score range:", sprintf("%.3f", range(results$colocScore, na.rm = TRUE)), "\n")
    }
  }
  
  return(results)
}

#' Compute colocalization scores for multi-slide object
#' @noRd
.getColocScoresMulti <- function(object, r_um_range, pixel_size_um, cell_diam_um,
                                nsim, r_step_um, min_points_per_type, edge_correction,
                                verbose, include_self) {
  
  slides <- getSlideList(object)
  cts <- unique(object@cellTypesSub)
  
  # Generate cell type pairs
  if (include_self) {
    all_pairs <- expand.grid(cellType1 = cts, cellType2 = cts, stringsAsFactors = FALSE)
  } else {
    if (length(cts) < 2) {
      stop("Need at least 2 cell types for cross-type colocalization analysis")
    }
    pair_matrix <- combn(cts, 2)
    all_pairs <- data.frame(
      cellType1 = pair_matrix[1, ],
      cellType2 = pair_matrix[2, ],
      stringsAsFactors = FALSE
    )
  }
  
  if (verbose) {
    cat("Computing colocalization scores for", length(slides), "slides and", 
        nrow(all_pairs), "cell type pairs\n")
  }
  
  # Initialize results list
  all_results <- vector("list", length = length(slides))
  names(all_results) <- slides
  
  # Process each slide
  for (slide_idx in seq_along(slides)) {
    slide_id <- slides[slide_idx]
    
    if (verbose) {
      cat("Processing slide", slide_idx, "of", length(slides), ":", slide_id, "\n")
    }
    
    # Get slide-specific indices
    slide_indices <- .getSlideIndices(object, slide_id)
    
    if (sum(slide_indices) == 0) {
      warning(paste("No cells found for slide:", slide_id))
      next
    }
    
    # Extract slide-specific data
    slide_location <- object@locationDataSub[slide_indices, , drop = FALSE]
    slide_celltypes <- object@cellTypesSub[slide_indices]
    
    # Determine window range for this slide
    window_range <- list(
      xrange = range(slide_location$x, na.rm = TRUE),
      yrange = range(slide_location$y, na.rm = TRUE)
    )
    
    # Initialize slide results
    slide_results <- data.frame(
      slideID = character(nrow(all_pairs)),
      cellType1 = character(nrow(all_pairs)),
      cellType2 = character(nrow(all_pairs)),
      colocScore = numeric(nrow(all_pairs)),
      nCells1 = integer(nrow(all_pairs)),
      nCells2 = integer(nrow(all_pairs)),
      r_um_min = numeric(nrow(all_pairs)),
      r_um_max = numeric(nrow(all_pairs)),
      stringsAsFactors = FALSE
    )
    
    # Compute colocalization for each pair on this slide
    for (i in seq_len(nrow(all_pairs))) {
      ct1 <- all_pairs$cellType1[i]
      ct2 <- all_pairs$cellType2[i]
      
      # Extract cell locations for each type on this slide
      cells1_idx <- slide_celltypes == ct1
      cells2_idx <- slide_celltypes == ct2
      
      if (sum(cells1_idx) == 0 || sum(cells2_idx) == 0) {
        # One or both cell types not present on this slide
        slide_results$slideID[i] <- slide_id
        slide_results$cellType1[i] <- ct1
        slide_results$cellType2[i] <- ct2
        slide_results$colocScore[i] <- NA_real_
        slide_results$nCells1[i] <- sum(cells1_idx)
        slide_results$nCells2[i] <- sum(cells2_idx)
        slide_results$r_um_min[i] <- r_um_range[1]
        slide_results$r_um_max[i] <- r_um_range[2]
        next
      }
      
      A <- slide_location[cells1_idx, c("x", "y"), drop = FALSE]
      B <- slide_location[cells2_idx, c("x", "y"), drop = FALSE]
      
      # Compute colocalization score
      coloc_result <- tryCatch({
        coloc_score(A = A, B = B, 
                   window_range = window_range,
                   r_um_range = r_um_range,
                   pixel_size_um = pixel_size_um,
                   cell_diam_um = cell_diam_um,
                   nsim = nsim,
                   r_step_um = r_step_um,
                   min_points_per_type = min_points_per_type)
      }, error = function(e) {
        if (verbose) {
          warning(paste("Error computing colocalization for", ct1, "vs", ct2, 
                       "on slide", slide_id, ":", e$message))
        }
        list(score = NA_real_, meta = list(nA = nrow(A), nB = nrow(B)))
      })
      
      # Store results
      slide_results$slideID[i] <- slide_id
      slide_results$cellType1[i] <- ct1
      slide_results$cellType2[i] <- ct2
      slide_results$colocScore[i] <- coloc_result$score
      slide_results$nCells1[i] <- coloc_result$meta$nA
      slide_results$nCells2[i] <- coloc_result$meta$nB
      slide_results$r_um_min[i] <- r_um_range[1]
      slide_results$r_um_max[i] <- r_um_range[2]
    }
    
    all_results[[slide_idx]] <- slide_results
  }
  
  # Combine results from all slides
  final_results <- do.call(rbind, all_results)
  rownames(final_results) <- NULL
  
  if (verbose) {
    valid_scores <- !is.na(final_results$colocScore)
    cat("Completed:", sum(valid_scores), "valid scores,", sum(!valid_scores), "failed across all slides\n")
    if (sum(valid_scores) > 0) {
      cat("Score range:", sprintf("%.3f", range(final_results$colocScore, na.rm = TRUE)), "\n")
    }
  }
  
  return(final_results)
}
