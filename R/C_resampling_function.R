
#' Get Neighboring Bins
#'
#' Internal helper function to find neighboring bins (3x3 grid around target).
#'
#' @param bin_coords Data frame with bin_id, x_bin, y_bin columns
#' @param tar_bin_id Target bin ID (string like "3_5")
#' @param num_bins_x Number of bins in x direction
#' @param num_bins_y Number of bins in y direction
#'
#' @return Character vector of neighboring bin IDs
#' @keywords internal
.get_neighbor_bins <- function(bin_coords, tar_bin_id,
                               num_bins_x, num_bins_y) {
  # Get coordinates of target bin
  tar_row <- bin_coords[bin_coords$bin_id == tar_bin_id, ]
  if (nrow(tar_row) == 0) {
    return(character(0))
  }

  x_bin <- tar_row$x_bin[1]
  y_bin <- tar_row$y_bin[1]

  neighbor_coords <- expand.grid(
    x_bin + -1:1,
    y_bin + -1:1
  )
  # Remove out-of-bounds bins

  neighbor_coords <- neighbor_coords[
    neighbor_coords$Var1 >= 1 & neighbor_coords$Var1 <= num_bins_x &
      neighbor_coords$Var2 >= 1 & neighbor_coords$Var2 <= num_bins_y, ]
  neighbor_bins <- paste(neighbor_coords$Var1, neighbor_coords$Var2, sep = "_")
  return(as.character(neighbor_bins))
}


#' Diagnose Bin Distribution
#'
#' Helper function to check how cells are distributed across bins.
#' Useful for choosing appropriate bin numbers.
#'
#' @param location_data Data frame with x, y columns
#' @param num_bins_x Number of bins in x direction
#' @param num_bins_y Number of bins in y direction
#'
#' @return List with bin statistics
#' @export
diagnose_bin_distribution <- function(location_data,
                                      num_bins_x = 10,
                                      num_bins_y = 10) {
  # Create bins
  x_bin <- cut(location_data$x, breaks = num_bins_x, labels = FALSE)
  y_bin <- cut(location_data$y, breaks = num_bins_y, labels = FALSE)
  bin_id <- paste(x_bin, y_bin, sep = "_")

  # Count cells per bin
  bin_counts <- table(bin_id)

  # Summary statistics
  stats <- list(
    total_cells = nrow(location_data),
    total_bins = num_bins_x * num_bins_y,
    occupied_bins = length(bin_counts),
    empty_bins = num_bins_x * num_bins_y - length(bin_counts),
    cells_per_bin = list(
      min = min(bin_counts),
      max = max(bin_counts),
      mean = mean(bin_counts),
      median = median(bin_counts),
      sd = sd(bin_counts)
    ),
    sparse_bins = sum(bin_counts < 5),
    na_cells = sum(is.na(x_bin) | is.na(y_bin))
  )

  # Recommendations
  if (stats$na_cells > 0) {
    warning(paste(stats$na_cells, "cells have NA bin assignments.",
                  "Check if data range is appropriate for bin numbers."))
  }

  if (stats$sparse_bins > stats$occupied_bins * 0.5) {
    message("Many bins have < 5 cells. Consider reducing num_bins_x/num_bins_y.")
  }

  return(stats)
}


#' Recover the raw-to-normalized distance scaling factor
#'
#' Internal helper that recovers the single scalar mapping raw
#' spatial-coordinate distances to the normalized distance scale on which
#' `sigmaValues` are defined. `computeDistance(normalizeDistance = TRUE)`
#' multiplies every stored distance by one constant, so for any entry that was
#' not affected by low-distance truncation the ratio (normalized distance) /
#' (raw Euclidean distance) equals that constant. We probe random off-diagonal
#' entries and take the median ratio, which is robust to the small fraction of
#' truncated entries and avoids copying the full distance matrix.
#'
#' @details
#' This assumes isotropic coordinates (`xDistScale = yDistScale`, the default).
#' Under anisotropic scaling the recovered value is a representative, not exact,
#' factor; the `[sigma, 4*sigma]` guardrail in [.sigmaAwareBins()] will flag a
#' grossly mis-sized grid. When `normalizeDistance = FALSE` the recovered factor
#' is simply the coordinate scale (1 by default), which is also correct because
#' `sigma` is then on the raw coordinate scale.
#'
#' @param object A CoPro object with `@distances` populated.
#' @param n_probe Number of random cell pairs to probe (default 200).
#'
#' @return Numeric scaling factor (normalized = raw * factor), or `NA_real_`
#'   if it cannot be recovered.
#' @keywords internal
.recoverDistanceScaleFactor <- function(object, n_probe = 200L) {
  if (length(object@distances) == 0) {
    return(NA_real_)
  }
  loc <- as.data.frame(object@locationDataSub)
  if (!all(c("x", "y") %in% colnames(loc)) || is.null(rownames(loc))) {
    return(NA_real_)
  }

  D <- object@distances[[1]]
  rn <- rownames(D)
  cn <- colnames(D)
  if (is.null(rn) || is.null(cn) || is.null(dim(D)) ||
      nrow(D) < 1 || ncol(D) < 1) {
    return(NA_real_)
  }

  loc_ids <- rownames(loc)
  ratios <- numeric(0)
  max_tries <- n_probe * 20L
  tries <- 0L
  while (length(ratios) < n_probe && tries < max_tries) {
    tries <- tries + 1L
    i <- sample.int(nrow(D), 1L)
    j <- sample.int(ncol(D), 1L)
    if (identical(rn[i], cn[j])) next         # skip self (within-type diagonal)
    d_norm <- D[i, j]
    if (!is.finite(d_norm) || d_norm <= 0) next
    if (!(rn[i] %in% loc_ids) || !(cn[j] %in% loc_ids)) next
    dx <- loc[rn[i], "x"] - loc[cn[j], "x"]
    dy <- loc[rn[i], "y"] - loc[cn[j], "y"]
    d_raw <- sqrt(dx * dx + dy * dy)
    if (!is.finite(d_raw) || d_raw <= 0) next
    ratios <- c(ratios, d_norm / d_raw)
  }

  if (length(ratios) == 0) {
    return(NA_real_)
  }
  stats::median(ratios)
}


#' Compute sigma-aware bin counts for spatial permutation
#'
#' Internal helper that chooses the number of bins for bin-wise spatial
#' permutation from the kernel bandwidth `sigma` instead of a fixed hard-coded
#' grid. The target patch side is `2 * sigma` on the normalized distance scale,
#' so each patch is large enough to preserve within-type spatial
#' autocorrelation (length scale `~sigma`) while shuffling whole patches still
#' breaks longer-range cross-type coordination. A hard-coded grid (e.g. the
#' historical 10x10) ignores both the tissue extent and `sigma`, which is the
#' root cause of mis-calibrated bin-wise permutation.
#'
#' @param object A CoPro object (uses `@locationDataSub` and `@distances`).
#' @param sigma Kernel bandwidth (numeric, on the normalized distance scale).
#' @param min_bins Minimum number of bins per axis (default 2).
#' @param verbose Whether to message the chosen grid (default TRUE).
#'
#' @return A list with integer `num_bins_x`, `num_bins_y`, and the recovered
#'   `scale_factor`.
#' @keywords internal
.sigmaAwareBins <- function(object, sigma, min_bins = 2L, verbose = TRUE) {
  loc <- as.data.frame(object@locationDataSub)
  ext_x <- diff(range(loc$x))
  ext_y <- diff(range(loc$y))

  sf <- .recoverDistanceScaleFactor(object)
  if (!is.finite(sf) || sf <= 0 || !is.finite(sigma) || sigma <= 0 ||
      !is.finite(ext_x) || !is.finite(ext_y) || ext_x <= 0 || ext_y <= 0) {
    warning("sigma-aware bins: could not recover distance scale; ",
            "falling back to a 10 x 10 grid. Pass num_bins_x/num_bins_y ",
            "explicitly to override.")
    return(list(num_bins_x = 10L, num_bins_y = 10L, scale_factor = NA_real_))
  }

  # Target patch side = 2*sigma on the normalized scale, in raw coordinate units
  patch_raw <- (2 * sigma) / sf
  nbx <- max(min_bins, floor(ext_x / patch_raw))
  nby <- max(min_bins, floor(ext_y / patch_raw))

  # Realized patch sides, mapped back to the normalized scale, for the guardrail
  side_x <- (ext_x / nbx) * sf
  side_y <- (ext_y / nby) * sf
  for (ax in c("x", "y")) {
    s <- if (ax == "x") side_x else side_y
    if (s < sigma) {
      warning(sprintf(paste0("sigma-aware bins: %s patch side (%.3g) < sigma ",
                             "(%.3g); patches may be too small and over-shuffle ",
                             "within-type structure (anti-conservative)."),
                      ax, s, sigma))
    } else if (s > 4 * sigma) {
      warning(sprintf(paste0("sigma-aware bins: %s patch side (%.3g) > 4*sigma ",
                             "(%.3g); patches may be too large to break ",
                             "cross-type coordination (conservative)."),
                      ax, s, 4 * sigma))
    }
  }

  if (verbose) {
    message(sprintf(paste0("sigma-aware bins (sigma = %g): %d x %d grid ",
                          "(patch side ~ %.3g x %.3g on the normalized scale)."),
                    sigma, nbx, nby, side_x, side_y))
  }

  list(num_bins_x = as.integer(nbx), num_bins_y = as.integer(nby),
       scale_factor = sf)
}


#' Match Cells by Within-Tile Quantile Position
#'
#' Internal helper to match cells from target tile to original tile based on
#' their relative (quantile) positions within each tile. This preserves
#' spatial structure better than random sampling.
#'
#' @param orig_points Data frame of original bin cells with x, y columns
#' @param candidate_points Data frame of candidate cells to sample from
#' @param n_points Number of cells to sample
#'
#' @return Integer vector of indices into candidate_points
#' @keywords internal
.match_by_quantile_position <- function(orig_points, candidate_points, n_points) {

  n_candidates <- nrow(candidate_points)

  # Calculate quantile positions for original points within their tile
  if (nrow(orig_points) == 1) {
    orig_x_quant <- 0.5
    orig_y_quant <- 0.5
  } else {
    orig_x_quant <- rank(orig_points$x, ties.method = "average") / (nrow(orig_points) + 1)
    orig_y_quant <- rank(orig_points$y, ties.method = "average") / (nrow(orig_points) + 1)
  }

  # Calculate quantile positions for candidate points
  if (n_candidates == 1) {
    cand_x_quant <- 0.5
    cand_y_quant <- 0.5
  } else {
    cand_x_quant <- rank(candidate_points$x, ties.method = "average") / (n_candidates + 1)
    cand_y_quant <- rank(candidate_points$y, ties.method = "average") / (n_candidates + 1)
  }

  # Fast path: if enough candidates, use vectorized nearest-neighbor matching
  if (n_candidates >= n_points) {
    # Compute full distance matrix (n_points x n_candidates)
    dist_matrix <- outer(orig_x_quant, cand_x_quant, "-")^2 +
                   outer(orig_y_quant, cand_y_quant, "-")^2

    # Greedy assignment: for each original point, find closest unused candidate
    sampled_indices <- integer(n_points)
    for (i in seq_len(n_points)) {
      best_idx <- which.min(dist_matrix[i, ])
      sampled_indices[i] <- best_idx
      dist_matrix[, best_idx] <- Inf  # Mark as used
    }
  } else {
    # Not enough candidates: need to allow reuse
    # Compute distances and allow multiple assignments
    dist_matrix <- outer(orig_x_quant, cand_x_quant, "-")^2 +
                   outer(orig_y_quant, cand_y_quant, "-")^2

    sampled_indices <- integer(n_points)
    used_counts <- integer(n_candidates)
    max_reuse <- ceiling(n_points / n_candidates)

    for (i in seq_len(n_points)) {
      best_idx <- which.min(dist_matrix[i, ])
      sampled_indices[i] <- best_idx
      used_counts[best_idx] <- used_counts[best_idx] + 1

      # If this candidate has been used max times, mark as unavailable
      if (used_counts[best_idx] >= max_reuse) {
        dist_matrix[, best_idx] <- Inf
      }
    }
  }

  return(sampled_indices)
}


#' Spatial Resampling for Permutation Testing
#'
#' Performs bin-wise spatial resampling to preserve local spatial structure
#' while breaking cross-type coordination. This is used for permutation testing
#' in CoPro to generate null distributions that account for spatial autocorrelation.
#'
#' @details
#' The algorithm:
#' 1. Divides the spatial domain into a grid of bins
#' 2. Creates a random mapping (shuffle) between bins
#' 3. For each original bin, samples cells from the mapped target bin
#' 4. If target bin has insufficient cells, expands to neighboring bins
#'
#' When `match_quantile = TRUE`, step 3 uses quantile-based matching instead of
#' random sampling. This matches cells based on their relative x/y positions
#' within each tile (using rank quantiles), better preserving within-tile
#' spatial autocorrelation structure.
#'
#' This preserves spatial autocorrelation within cell types while breaking
#' the cross-type coordination that would indicate true co-progression.
#'
#' @param location_data Data frame with columns: x, y, cell_ID.
#'   Optionally can include pre-computed x_bin, y_bin columns.
#' @param num_bins_x Number of bins in the x direction (default: 10).
#'   Larger values preserve more local structure but may have sparse bins.
#' @param num_bins_y Number of bins in the y direction (default: 10).
#' @param match_quantile Logical. If TRUE, matches cells between original and
#'   target tiles based on their relative (quantile) positions within each tile.
#'   This better preserves within-tile spatial structure. Default: FALSE.
#'
#' @return Data frame with resampled cell assignments. The x, y coordinates
#'   remain fixed, but cell_ID is shuffled according to the bin-wise resampling.
#'
#' @examples
#' \dontrun{
#' # Create example data
#' loc_data <- data.frame(
#'   x = runif(100, 0, 10),
#'   y = runif(100, 0, 10),
#'   cell_ID = paste0("cell_", 1:100)
#' )
#'
#' # Check bin distribution first
#' diagnose_bin_distribution(loc_data, num_bins_x = 5, num_bins_y = 5)
#'
#' # Perform spatial resampling (random within tile)
#' resampled <- resample_spatial(loc_data, num_bins_x = 5, num_bins_y = 5)
#'
#' # Perform spatial resampling (quantile-matched within tile)
#' resampled_matched <- resample_spatial(loc_data, num_bins_x = 5, num_bins_y = 5,
#'                                       match_quantile = TRUE)
#' }
#'
#' @export
resample_spatial <- function(location_data,
                             num_bins_x = 10, num_bins_y = 10,
                             match_quantile = FALSE) {

  # Input validation
 if (is.matrix(location_data) || is.data.frame(location_data)) {
    if (is.matrix(location_data)) {
      location_data <- as.data.frame(location_data)
    }
    if (!all(c("x", "y", "cell_ID") %in% colnames(location_data))) {
      stop("Input location_data must have 'x', 'y', and 'cell_ID' columns.")
    }
  } else {
    stop("Input location_data must be a matrix or data frame.")
  }

  # Validate bin numbers
  if (num_bins_x < 2 || num_bins_y < 2) {
    stop("num_bins_x and num_bins_y must be at least 2.")
  }
  if (num_bins_x > nrow(location_data) || num_bins_y > nrow(location_data)) {
    warning("Number of bins exceeds number of cells. Consider reducing bin numbers.")
  }

  original_cell_loc_order <- paste(location_data$x, location_data$y, sep = "_")
  rownames(location_data) <- original_cell_loc_order

  # Create bins for x and y coordinates if not already present
  if (!all(c("x_bin", "y_bin") %in% colnames(location_data))) {
    location_data$x_bin <- cut(location_data$x, breaks = num_bins_x, labels = FALSE)
    location_data$y_bin <- cut(location_data$y, breaks = num_bins_y, labels = FALSE)
  }

  # Handle NA bins (cells outside the binning range)
  na_bins <- is.na(location_data$x_bin) | is.na(location_data$y_bin)
  if (any(na_bins)) {
    warning(paste(sum(na_bins), "cells have NA bin assignments and will be",
                  "assigned to nearest valid bin."))
    # Assign NA bins to the closest valid bin
    if (any(is.na(location_data$x_bin))) {
      location_data$x_bin[is.na(location_data$x_bin)] <- 1
    }
    if (any(is.na(location_data$y_bin))) {
      location_data$y_bin[is.na(location_data$y_bin)] <- 1
    }
  }

  # Assign a bin ID to each point
  location_data$bin_id <- paste(location_data$x_bin, location_data$y_bin, sep = "_")

  # Get unique bins and shuffle them
  unique_bins <- unique(location_data$bin_id)
  shuffled_bins <- sample(unique_bins)

  # Create a mapping from original bins to shuffled bins
  bin_mapping <- setNames(shuffled_bins, unique_bins)

  # Create a lookup table for bin coordinates
  bin_coords <- unique(location_data[, c("bin_id", "x_bin", "y_bin")])
  rownames(bin_coords) <- bin_coords$bin_id

  # Initialize a list to store resampled location_data
  resampled_list <- vector("list", length(unique_bins))

  # Resample points for each bin
  for (i in seq_along(unique_bins)) {
    orig_bin <- unique_bins[i]
    target_bin <- bin_mapping[[orig_bin]]

    id_ori <- location_data$bin_id == orig_bin
    id_tar <- location_data$bin_id == target_bin

    # Points in the original bin
    orig_points <- location_data[id_ori, ]
    n_points <- sum(id_ori)

    # Candidate points from the target bin
    candidate_bins <- target_bin
    candidate_points <- location_data[id_tar, ]
    n_points_tar <- sum(id_tar)

    # If not enough points, include neighboring bins
    if (n_points_tar < n_points) {
      neighbors <- .get_neighbor_bins(bin_coords = bin_coords,
                                      tar_bin_id = target_bin,
                                      num_bins_x = num_bins_x,
                                      num_bins_y = num_bins_y)
      # Exclude bins already in candidate_bins
      additional_bins <- setdiff(neighbors, candidate_bins)
      for (neighbor_bin in additional_bins) {
        neighbor_idx <- location_data$bin_id == neighbor_bin
        if (any(neighbor_idx)) {
          neighbor_points <- location_data[neighbor_idx, ]
          candidate_points <- rbind(candidate_points, neighbor_points)
          candidate_bins <- c(candidate_bins, neighbor_bin)
        }
        if (nrow(candidate_points) >= n_points) {
          break
        }
      }
    }

    # Sample cells from candidates
    if (match_quantile) {
      # Use quantile-based matching to preserve within-tile structure
      sampled_indices <- .match_by_quantile_position(orig_points, candidate_points, n_points)
    } else {
      # Random sampling (original behavior)
      if (nrow(candidate_points) < n_points) {
        sampled_indices <- sample(nrow(candidate_points), n_points, replace = TRUE)
      } else {
        sampled_indices <- sample(nrow(candidate_points), n_points)
      }
    }

    # Replace the cell_ID in the original points
    orig_points$"cell_ID" <- candidate_points[sampled_indices, "cell_ID"]
    orig_points$"bin_id" <- candidate_points[sampled_indices, "bin_id"]

    # Store the resampled points
    resampled_list[[i]] <- orig_points
  }

  # Combine all resampled points
  resampled_location_data <- do.call(rbind, resampled_list)

  # Ensure the row names are in the original order
  rownames(resampled_location_data) <- paste(resampled_location_data$x,
                                             resampled_location_data$y,
                                             sep = "_")
  resampled_location_data <- resampled_location_data[original_cell_loc_order, ]

  return(resampled_location_data)
}


#' Generate Toroidal Shift Permutation Indices
#'
#' Shifts spatial coordinates in a toroidal (wrap-around) manner to break
#' cross-type coordination while approximately preserving within-type spatial
#' autocorrelation. Useful as one of several autocorrelation-respecting nulls
#' for permutation testing (see also the sigma-aware `bin` null and Moran
#' Spectral Randomization for irregular tissue).
#'
#' @details
#' The toroidal shift works by:
#' 1. Applying a random shift to all coordinates (wrapping at boundaries)
#' 2. Re-ranking cells by their shifted positions and matching back to the
#'    original ordering.
#'
#' Important caveats (the docstring previously over-claimed "perfect"
#' preservation):
#' \itemize{
#'   \item The position matching is by coordinate *rank*, which equals a rigid
#'     torus translation only on a regular lattice. On an irregular point cloud
#'     it is a monotone rearrangement that preserves pairwise distances only
#'     approximately, so within-type autocorrelation is approximately (not
#'     exactly) preserved.
#'   \item The torus-translation test (Harms et al. 2001) assumes spatial
#'     stationarity and periodic wrap-around. Gluing opposite tissue edges
#'     creates artificial neighbours at the seam and is biologically false for
#'     non-rectangular / hole-bearing tissue, and the family of distinct shifts
#'     is effectively small (a coarse null with limited power).
#' }
#' For these reasons toroidal should be benchmarked against, not preferred over,
#' the sigma-aware patch null and graph-based surrogates; let calibration choose.
#'
#' @param location_data Data frame with x, y, cell_ID columns
#' @param n_permu Number of permutations to generate
#' @param seed Optional random seed for reproducibility
#'
#' @return Matrix of permutation indices (n_cells x n_permu). Each column
#'   contains a permutation of row indices that can be used to reorder cells.
#'
#' @examples
#' \dontrun{
#' # Create example location data
#' loc_data <- data.frame(
#'   x = runif(100, 0, 10),
#'   y = runif(100, 0, 10),
#'   cell_ID = paste0("cell_", 1:100)
#' )
#'
#' # Generate 100 toroidal permutations
#' perm_matrix <- generate_toroidal_permutations(loc_data, n_permu = 100)
#'
#' # Apply first permutation
#' permuted_cells <- loc_data$cell_ID[perm_matrix[, 1]]
#' }
#'
#' @importFrom stats runif
#' @export
generate_toroidal_permutations <- function(location_data, n_permu = 100,
                                           seed = NULL) {

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
    on.exit({
      if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv)
      else assign(".Random.seed", old_seed, envir = .GlobalEnv)
    })
    set.seed(seed)
  }

  # Input validation
  if (!is.data.frame(location_data)) {
    stop("location_data must be a data frame")
  }
  if (!all(c("x", "y") %in% colnames(location_data))) {
    stop("location_data must have 'x' and 'y' columns")
  }

  n_cells <- nrow(location_data)

  # Get spatial extent
  x_range <- range(location_data$x)
  y_range <- range(location_data$y)
  x_width <- diff(x_range)
  y_width <- diff(y_range)

  # Handle edge case of zero width

if (x_width < 1e-10 || y_width < 1e-10) {
    warning("Spatial extent is very small. Toroidal shift may not work well.")
    # Return identity permutations
    return(replicate(n_permu, seq_len(n_cells)))
  }

  # Original ordering by position (for matching)
  orig_order <- order(location_data$x, location_data$y)

  perm_matrix <- matrix(NA_integer_, nrow = n_cells, ncol = n_permu)

  for (tt in seq_len(n_permu)) {
    # Random shifts (ensure non-trivial shift)
    shift_x <- runif(1, x_width * 0.1, x_width * 0.9)
    shift_y <- runif(1, y_width * 0.1, y_width * 0.9)

    # Apply toroidal shift (wrap-around)
    new_x <- ((location_data$x - x_range[1] + shift_x) %% x_width) + x_range[1]
    new_y <- ((location_data$y - y_range[1] + shift_y) %% y_width) + y_range[1]

    # Create permutation by matching positions
    # Sort both original and shifted by coordinates
    new_order <- order(new_x, new_y)

    # The cell at original position i gets data from cell at shifted position
    perm_matrix[, tt] <- new_order[order(orig_order)]
  }

  return(perm_matrix)
}
