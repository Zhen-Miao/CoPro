#' computeDistance between pairs of cell types
#'
#' @importFrom fields rdist
#' @importFrom utils combn
#' @importFrom stats setNames quantile lm predict coef
#' @param object A `CoPro` object
#' @param distType Type of distance to compute: "Euclidean2D",
#'  "Euclidean3D", or "Morphology-Aware"
#' @param xDistScale Scale for x distance
#' @param yDistScale Scale for y distance
#' @param zDistScale Scale for z distance
#' @param verbose Whether to print info about the quantile of the distance
#' @param normalizeDistance Whether to normalize distance? The normalization
#'  will make sure that the 0.01% cell-cell distance will become 0.01, thus
#'  ensuring no matter which input scale is used for the distance matrix,
#'  the output will roughly be in mm^3. This ensures that the kernel sizes
#'  from 0.001 to 0.1 will make sense. Default = TRUE
#' @param truncateLowDist Whether to truncate small distances so that the cells
#'  that are nearly overlapping with each other do not have a super small
#'  distance. Default = TRUE.
#' @param knn_k Number of nearest neighbors for KNN graph construction
#'  (used only for Morphology-Aware distance). Default = 10.
#' @param geodesic_threshold Geodesic distance threshold for regression fitting.
#'  Cell pairs with geodesic distance <= this value are used to fit the
#'  regression. Default = 10.
#' @param geodesic_cutoff Geodesic distance value at which to evaluate the
#'  regression for determining the Euclidean distance cutoff. Default = 7.
#' @return `CoPro` object with distance matrix computed
#' @rdname computeDistance
#' @aliases computeDistance,CoProSingle-method
#' @aliases computeDistance,CoProMulti-method
#' @export
#' @details
#' For "Morphology-Aware" distance:
#' This method addresses the issue where cells appear spatially close (small
#' Euclidean distance) but are actually topologically distant due to tissue
#' morphology (e.g., semilunar folds in colon tissue).
#'
#' The algorithm:
#' \enumerate{
#'   \item Compute Euclidean distances (d_E) between all cells
#'   \item Build a K-nearest neighbor (KNN) graph and compute geodesic
#'         distances (d_g) as shortest path distances on this graph
#'   \item Fit a linear regression between d_g and d_E for cell pairs
#'         where d_g <= geodesic_threshold (default: 10)
#'   \item Use the fitted model to predict d_E at d_g = geodesic_cutoff
#'         (default: 7), this becomes cutoff_d_E
#'   \item For cell pairs where d_g > geodesic_threshold AND d_E < cutoff_d_E,
#'         set distance to the maximum distance in the matrix
#' }
#'
#' This requires the igraph package to be installed.
setGeneric(
  "computeDistance",
  function(object, distType =
             c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
           xDistScale = 1, yDistScale = 1,
           zDistScale = 1, normalizeDistance = TRUE, truncateLowDist = TRUE,
           verbose = TRUE, knn_k = 10, geodesic_threshold = 10,
           geodesic_cutoff = 7) standardGeneric("computeDistance")
)

# Helper function to choose cell types
.choose_cts <- function(object) {
  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    stop(paste("no cell type of interest specified,"))
  }
  return(cts)
}

# Helper function to check distance type
.check_dist_type <- function(distType, object) {
  if (distType == "Euclidean2D") {
    if (!all(c("x", "y") %in% colnames(object@locationDataSub))) {
      stop(paste("please make sure x, y are all available to run",
                 "2D Euclidean distance calculation"))
    }
  } else if (distType == "Euclidean3D") {
    if (!all(c("x", "y", "z") %in% colnames(object@locationDataSub))) {
      stop(paste("please make sure x, y, z are all available to run",
                 "3D Euclidean distance calculation"))
    }
  } else if (distType == "Morphology-Aware") {
    # Morphology-aware requires 2D coordinates
    if (!all(c("x", "y") %in% colnames(object@locationDataSub))) {
      stop(paste("please make sure x, y are all available to run",
                 "Morphology-Aware distance calculation"))
    }
    # Check if igraph is available
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop(paste("Package 'igraph' is required for Morphology-Aware distance.",
                 "Please install it with: install.packages('igraph')"))
    }
  }
  return(TRUE)
}

# Helper function to compute KNN graph from coordinate matrix
# Returns a sparse adjacency matrix
#
# Note on memory usage: For very large datasets (>50k cells), computing the
# full distance matrix could be memory-intensive. The current implementation
# computes n_cells x n_cells matrices. For extremely large datasets, consider
# using approximate nearest neighbor methods (e.g., Annoy, HNSW).
.computeKnnGraph <- function(coord_mat, k = 10) {
  n_cells <- nrow(coord_mat)
  
  # Validate k parameter
  if (k >= n_cells) {
    warning(sprintf("k (%d) >= number of cells (%d). Setting k = %d.",
                    k, n_cells, n_cells - 1))
    k <- n_cells - 1
  }
  
  if (k < 1) {
    stop("k must be at least 1 for KNN graph construction.")
  }
  
  # Compute pairwise Euclidean distances
  dist_mat <- fields::rdist(coord_mat)
  
  # For each cell, find the k nearest neighbors
  # We use k+1 because the cell itself will be included
  knn_indices <- t(apply(dist_mat, 1, function(row) {
    order(row)[2:(k + 1)]  # Exclude self (index 1 after sorting)
  }))
  
  # Build sparse adjacency matrix
  # Create edge list
  i_indices <- rep(seq_len(n_cells), each = k)
  j_indices <- as.vector(t(knn_indices))
  
  # Create symmetric adjacency (undirected graph)
  adj_mat <- Matrix::sparseMatrix(
    i = c(i_indices, j_indices),
    j = c(j_indices, i_indices),
    x = 1,
    dims = c(n_cells, n_cells),
    dimnames = list(rownames(coord_mat), rownames(coord_mat))
  )
  
  # Ensure binary adjacency (no duplicates count as 2)
  adj_mat@x[adj_mat@x > 0] <- 1
  
  return(adj_mat)
}

# Helper function to compute geodesic distances from KNN adjacency matrix
# Uses igraph for shortest path computation
#
# Note on disconnected components: If the KNN graph has disconnected components,
# igraph::distances() returns Inf for pairs that are not connected. This is
# handled correctly in .applyMorphologyFilter() by including is.infinite(d_g)
# in the filter mask, treating unreachable cells as topologically distant.
.computeGeodesicDistance <- function(adj_mat) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for geodesic distance computation.")
  }
  
  # Convert adjacency matrix to igraph object
  graph <- igraph::graph_from_adjacency_matrix(
    adj_mat,
    mode = "undirected",
    weighted = NULL
  )
  
  # Compute shortest paths (geodesic distances)
  geodesic_dist <- igraph::distances(graph, mode = "all")
  
  # Preserve row/column names from adjacency matrix
  if (!is.null(rownames(adj_mat))) {
    rownames(geodesic_dist) <- rownames(adj_mat)
    colnames(geodesic_dist) <- colnames(adj_mat)
  }
  
  return(geodesic_dist)
}

# Helper function to apply morphology-aware filtering
# d_E: Euclidean distance matrix
# d_g: Geodesic distance matrix
# geodesic_threshold: threshold for regression fitting (default: 10)
# geodesic_cutoff: geodesic value for cutoff prediction (default: 7)
#
# Note on regression quality: The effectiveness of this filtering depends on
# a reasonably strong linear relationship between geodesic and Euclidean
# distances for nearby cells (d_g <= geodesic_threshold). If R-squared is low,
# the cutoff may not be meaningful. In such cases, consider:
#   - Increasing knn_k to create a denser graph
#   - Checking if the tissue has the expected folded morphology
#   - Using standard Euclidean distance instead
.applyMorphologyFilter <- function(d_E, d_g, geodesic_threshold = 10,
                                   geodesic_cutoff = 7, verbose = TRUE) {
  
  # Validate parameters
  if (geodesic_cutoff > geodesic_threshold) {
    warning(sprintf(paste("geodesic_cutoff (%.1f) > geodesic_threshold (%.1f).",
                          "This means extrapolating beyond the regression data.",
                          "Consider setting geodesic_cutoff <= geodesic_threshold."),
                    geodesic_cutoff, geodesic_threshold))
  }
  
  # Get the maximum distance for replacement
  max_dist <- max(d_E[is.finite(d_E)], na.rm = TRUE)
  
  # Step 1: Select pairs with d_g <= geodesic_threshold for regression
  mask_fit <- d_g <= geodesic_threshold & is.finite(d_g) & is.finite(d_E)
  
  # Extract values for regression
  d_g_fit <- d_g[mask_fit]
  d_E_fit <- d_E[mask_fit]
  
  if (length(d_g_fit) < 10) {
    warning(paste("Too few cell pairs with geodesic distance <=", geodesic_threshold,
                  "for reliable regression. Returning Euclidean distances."))
    return(d_E)
  }
  
  # Step 2: Fit linear regression: d_E ~ d_g
  fit_data <- data.frame(d_g = d_g_fit, d_E = d_E_fit)
  model <- stats::lm(d_E ~ d_g, data = fit_data)
  r_squared <- summary(model)$r.squared
  
  if (verbose) {
    message("Morphology-aware distance regression:")
    message(sprintf("  Intercept: %.4f", coef(model)[1]))
    message(sprintf("  Slope: %.4f", coef(model)[2]))
    message(sprintf("  R-squared: %.4f", r_squared))
  }
  
  # Warn if R-squared is too low (weak relationship)
  if (r_squared < 0.3) {
    warning(sprintf(paste(
      "Low R-squared (%.3f) in morphology-aware regression.",
      "The relationship between geodesic and Euclidean distances is weak.",
      "The distance cutoff may not be meaningful. Consider:",
      "  - Increasing knn_k for a denser graph",
      "  - Checking if the tissue has folded morphology",
      "  - Using 'Euclidean2D' distance type instead"),
      r_squared))
  } else if (r_squared < 0.5) {
    warning(sprintf(paste(
      "Moderate R-squared (%.3f) in morphology-aware regression.",
      "Results may be less reliable. Check verbose output for diagnostics."),
      r_squared))
  }
  
  # Step 3: Predict d_E at geodesic_cutoff
  cutoff_d_E <- stats::predict(model, newdata = data.frame(d_g = geodesic_cutoff))
  
  if (verbose) {
    message(sprintf("  Cutoff d_E (at d_g = %.1f): %.4f", geodesic_cutoff, cutoff_d_E))
  }
  
  # Step 4: Apply filtering
  # Mask for pairs that are:
  # - topologically distant (d_g > geodesic_threshold)
  # - but spatially close (d_E < cutoff_d_E)
  mask_filter <- (d_g > geodesic_threshold | is.infinite(d_g)) & (d_E < cutoff_d_E)
  
  n_filtered <- sum(mask_filter, na.rm = TRUE)
  if (verbose) {
    message(sprintf("  Filtering %d cell pairs (%.2f%% of total)",
                    n_filtered, 100 * n_filtered / length(d_E)))
  }
  
  # Create filtered distance matrix
  d_filtered <- d_E
  d_filtered[mask_filter] <- max_dist
  
  return(d_filtered)
}

# Helper function to validate inputs
.checkInputDistance <- function(object, distType, xDistScale, yDistScale, zDistScale) {
  # Check distance type
  .check_dist_type(distType, object)
  
  # Check cell types
  cts <- .choose_cts(object)
  
  # Check scale parameters
  if (any(c(xDistScale, yDistScale, zDistScale) <= 0)) {
    stop("Distance scales must be positive")
  }
  
  # Additional check for 3D distance type and zDistScale
  if (distType == "Euclidean3D" && zDistScale == 1 && (xDistScale != 1 || yDistScale != 1)) {
    warning("Using default zDistScale = 1 while x/y scales are different. Consider setting zDistScale explicitly.")
  }
  
  return(cts)
}

# Core dispatcher function
.computeDistanceCore <- function(object, distType, xDistScale, yDistScale, zDistScale,
                                normalizeDistance, truncateLowDist, verbose,
                                knn_k = 10, geodesic_threshold = 10,
                                geodesic_cutoff = 7) {
  cts <- .checkInputDistance(object, distType, xDistScale, yDistScale, zDistScale)
  
  # Determine whether to compute pairwise or within-cell-type distances
  if (length(cts) == 1) {
    return(.computeDistanceWithin(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                 normalizeDistance, truncateLowDist, verbose,
                                 knn_k, geodesic_threshold, geodesic_cutoff))
  } else {
    return(.computeDistancePairs(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                normalizeDistance, truncateLowDist, verbose,
                                knn_k, geodesic_threshold, geodesic_cutoff))
  }
}

# Helper function to extract coordinate matrix
.getCoordinateMatrix <- function(object, cellType, distType, xDistScale, yDistScale, zDistScale, slideID = NULL) {
  # Get cell type subset indices
  if (is.null(slideID)) {
    # Single slide case
    ct_indices <- object@cellTypesSub == cellType
  } else {
    # Multi-slide case - filter by both cell type and slide
    ct_indices <- .getSlideCellTypeIndices(object, slide = slideID, cellType = cellType)
  }
  
  # Check if any cells were found
  if (!any(ct_indices)) {
    if (is.null(slideID)) {
      stop(paste("No cells found for cell type:", cellType))
    } else {
      stop(paste("No cells found for cell type:", cellType, "in slide:", slideID))
    }
  }
  
  if (distType == "Euclidean2D" || distType == "Morphology-Aware") {
    # Morphology-Aware uses 2D coordinates as the base
    mat <- cbind(
      object@locationDataSub$x[ct_indices] * xDistScale,
      object@locationDataSub$y[ct_indices] * yDistScale
    )
  } else if (distType == "Euclidean3D") {
    mat <- cbind(
      object@locationDataSub$x[ct_indices] * xDistScale,
      object@locationDataSub$y[ct_indices] * yDistScale,
      object@locationDataSub$z[ct_indices] * zDistScale
    )
  }
  
  # Set rownames to match the cell IDs
  if (!is.null(slideID)) {
    rownames(mat) <- .getSlideCellTypeIDs(object, slide = slideID, cellType = cellType)
  } else {
    rownames(mat) <- rownames(object@locationDataSub)[ct_indices]
  }
  
  return(mat)
}

# Helper function to process distance matrix
.processDistanceMatrix <- function(distances_ij, truncateLowDist, 
                                  percentile_choice = NULL, set_diag_inf = FALSE) {
  # Set diagonal to infinity if requested (for within-cell-type distances)
  if (set_diag_inf) {
    diag(distances_ij) <- Inf
  }
  
  # Handle zero distances
  if (any(distances_ij == 0, na.rm = TRUE)) {
    warning(paste("Zero distances detected, replacing with",
                  "the smallest non-zero distances, please",
                  "consider checking the location of cells",
                  "for potential errors"))
    min_nonzero <- min(distances_ij[distances_ij != 0 & !is.infinite(distances_ij)], na.rm = TRUE)
    distances_ij[distances_ij == 0] <- min_nonzero
  }
  
  # Calculate percentile if not provided
  if (is.null(percentile_choice)) {
    percentile_choice <- min(1e-3, 2/(max(nrow(distances_ij), ncol(distances_ij))))
  }
  
  # Get distance percentile (excluding infinite values)
  finite_distances <- distances_ij[is.finite(distances_ij) & distances_ij != 0]
  if (length(finite_distances) == 0) {
    stop("No finite non-zero distances found")
  }
  
  dist_percentile <- quantile(finite_distances, percentile_choice)
  
  # Truncate low distances if requested
  if (truncateLowDist) {
    distances_ij[distances_ij < dist_percentile & is.finite(distances_ij)] <- dist_percentile
  }
  
  return(list(distances = distances_ij, percentile = dist_percentile))
}

# Function for computing distances between pairs of cell types
.computeDistancePairs <- function(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                 normalizeDistance, truncateLowDist, verbose,
                                 knn_k = 10, geodesic_threshold = 10,
                                 geodesic_cutoff = 7) {
  
  # Initialize flat distance structure
  distances <- list()
  
  pair_cell_types <- combn(cts, 2)
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance) {
    cat("normalizeDistance is set to TRUE, so distance will be",
        "normalized, so that 0.01 percentile distance will be scaled",
        "to 0.01\n")
  }
  
  dist_percentiles <- vector(mode = "numeric", length = ncol(pair_cell_types))
  
  # For Morphology-Aware: compute KNN and geodesic on ALL cells first
  geodesic_all <- NULL
  cell_ids_all <- NULL
  if (distType == "Morphology-Aware") {
    if (verbose) message("Computing KNN graph and geodesic distances for all cells...")
    
    # Get coordinates for ALL cells (not just cell types of interest)
    all_coords <- cbind(
      object@locationDataSub$x * xDistScale,
      object@locationDataSub$y * yDistScale
    )
    rownames(all_coords) <- rownames(object@locationDataSub)
    cell_ids_all <- rownames(all_coords)
    
    # Compute KNN graph on all cells
    knn_adj <- .computeKnnGraph(all_coords, k = knn_k)
    
    # Compute geodesic distances
    geodesic_all <- .computeGeodesicDistance(knn_adj)
    
    if (verbose) message("Geodesic distance computation complete.")
  }
  
  # Calculate distances for each pair
  for (pp in seq_len(ncol(pair_cell_types))) {
    i <- pair_cell_types[1, pp]
    j <- pair_cell_types[2, pp]
    
    # Get coordinate matrices
    mat1 <- .getCoordinateMatrix(object, i, distType, xDistScale, yDistScale, zDistScale)
    mat2 <- .getCoordinateMatrix(object, j, distType, xDistScale, yDistScale, zDistScale)
    
    # Compute Euclidean distance
    distances_ij <- fields::rdist(mat1, mat2)
    
    # Apply morphology-aware filtering if requested
    if (distType == "Morphology-Aware") {
      # Extract geodesic distances for this cell type pair
      idx_i <- match(rownames(mat1), cell_ids_all)
      idx_j <- match(rownames(mat2), cell_ids_all)
      if (any(is.na(idx_i)) || any(is.na(idx_j))) {
        stop("Cell IDs in coordinate matrix do not match geodesic distance matrix. ",
             "Check that cell IDs are consistent across inputs.")
      }
      geodesic_ij <- geodesic_all[idx_i, idx_j, drop = FALSE]

      if (verbose) message(sprintf("Applying morphology filter for pair: %s - %s", i, j))
      
      # Apply filtering
      distances_ij <- .applyMorphologyFilter(
        d_E = distances_ij,
        d_g = geodesic_ij,
        geodesic_threshold = geodesic_threshold,
        geodesic_cutoff = geodesic_cutoff,
        verbose = verbose
      )
    }
    
    # Process distance matrix
    processed <- .processDistanceMatrix(distances_ij, truncateLowDist)
    distances_ij <- processed$distances
    dist_percentiles[pp] <- processed$percentile
    
    # Save the distances using flat structure
    flat_name <- .createDistMatrixName(i, j, slide = NULL)
    distances[[flat_name]] <- distances_ij
    
    if (verbose) {
      if (verbose) print(quantile(distances_ij, na.rm = TRUE))
    }
  }
  
  # Apply normalization if requested
  if (normalizeDistance) {
    min_percentile <- min(dist_percentiles)
    scaling_factor <- 0.01 / min_percentile
    cat("The scaling factor for normalizing distance is", scaling_factor, "\n")
    
    for (pp in seq_len(ncol(pair_cell_types))) {
      i <- pair_cell_types[1, pp]
      j <- pair_cell_types[2, pp]
      flat_name <- .createDistMatrixName(i, j, slide = NULL)
      distances[[flat_name]] <- distances[[flat_name]] * scaling_factor
    }
  }
  
  object@distances <- distances
  return(object)
}

# Function for computing within-cell-type distances
.computeDistanceWithin <- function(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                  normalizeDistance, truncateLowDist, verbose,
                                  knn_k = 10, geodesic_threshold = 10,
                                  geodesic_cutoff = 7) {
  
  # Initialize flat distance structure
  distances <- list()
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance) {
    cat("normalizeDistance is set to TRUE, so distance will be",
        "normalized, so that 0.01 percentile distance will be scaled",
        "to 0.01\n")
  }
  
  # Get coordinate matrix for the single cell type
  mat1 <- .getCoordinateMatrix(object, cts, distType, xDistScale, yDistScale, zDistScale)
  
  # Compute Euclidean distance matrix
  distances_ij <- fields::rdist(mat1)
  
  # Apply morphology-aware filtering if requested
  if (distType == "Morphology-Aware") {
    if (verbose) message("Computing KNN graph and geodesic distances for all cells...")
    
    # Get coordinates for ALL cells (not just cell types of interest)
    all_coords <- cbind(
      object@locationDataSub$x * xDistScale,
      object@locationDataSub$y * yDistScale
    )
    rownames(all_coords) <- rownames(object@locationDataSub)
    cell_ids_all <- rownames(all_coords)
    
    # Compute KNN graph on all cells
    knn_adj <- .computeKnnGraph(all_coords, k = knn_k)
    
    # Compute geodesic distances
    geodesic_all <- .computeGeodesicDistance(knn_adj)
    
    if (verbose) message("Geodesic distance computation complete.")
    
    # Extract geodesic distances for this cell type
    idx <- match(rownames(mat1), cell_ids_all)
    geodesic_ij <- geodesic_all[idx, idx, drop = FALSE]
    
    if (verbose) message(sprintf("Applying morphology filter for cell type: %s", cts))
    
    # Apply filtering
    distances_ij <- .applyMorphologyFilter(
      d_E = distances_ij,
      d_g = geodesic_ij,
      geodesic_threshold = geodesic_threshold,
      geodesic_cutoff = geodesic_cutoff,
      verbose = verbose
    )
  }
  
  # Process distance matrix (set diagonal to Inf)
  processed <- .processDistanceMatrix(distances_ij, truncateLowDist, 
                                     percentile_choice = 1e-4, set_diag_inf = TRUE)
  distances_ij <- processed$distances
  dist_percentile <- processed$percentile
  
  # Save the distances using flat structure
  flat_name <- .createDistMatrixName(cts, cts, slide = NULL)
  distances[[flat_name]] <- distances_ij
  
  if (verbose) {
    if (verbose) print(quantile(distances_ij[is.finite(distances_ij)], na.rm = TRUE))
  }
  
  # Apply normalization if requested
  if (normalizeDistance) {
    scaling_factor <- 0.01 / dist_percentile
    cat("The scaling factor for normalizing distance is", scaling_factor, "\n")
    flat_name <- .createDistMatrixName(cts, cts, slide = NULL)
    distances[[flat_name]] <- distances_ij * scaling_factor
  }
  
  object@distances <- distances
  return(object)
}

#' @rdname computeDistance
#' @aliases computeDistance,CoProSingle-method
#' @export
setMethod("computeDistance", "CoProSingle", function(object, distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
                                                    xDistScale = 1, yDistScale = 1, zDistScale = 1,
                                                    normalizeDistance = TRUE, truncateLowDist = TRUE, verbose = TRUE,
                                                    knn_k = 10, geodesic_threshold = 10, geodesic_cutoff = 7) {
  distType <- match.arg(distType)
  .computeDistanceCore(object, distType, xDistScale, yDistScale, zDistScale,
                      normalizeDistance, truncateLowDist, verbose,
                      knn_k, geodesic_threshold, geodesic_cutoff)
})

#' @rdname computeDistance
#' @aliases computeDistance,CoProMulti-method
#' @export
setMethod("computeDistance", "CoProMulti", function(object, distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
                                                   xDistScale = 1, yDistScale = 1, zDistScale = 1,
                                                   normalizeDistance = TRUE, truncateLowDist = TRUE, verbose = TRUE,
                                                   knn_k = 10, geodesic_threshold = 10, geodesic_cutoff = 7) {
  distType <- match.arg(distType)
  .computeDistanceCoreMulti(object, distType, xDistScale, yDistScale, zDistScale,
                           normalizeDistance, truncateLowDist, verbose,
                           knn_k, geodesic_threshold, geodesic_cutoff)
})

# Core dispatcher for multi-slide objects
.computeDistanceCoreMulti <- function(object, distType, xDistScale, yDistScale, zDistScale,
                                     normalizeDistance, truncateLowDist, verbose,
                                     knn_k = 10, geodesic_threshold = 10,
                                     geodesic_cutoff = 7) {
  cts <- .checkInputDistance(object, distType, xDistScale, yDistScale, zDistScale)
  
  # Determine whether to compute pairwise or within-cell-type distances across slides
  if (length(cts) == 1) {
    return(.computeDistanceMultiWithin(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                      normalizeDistance, truncateLowDist, verbose,
                                      knn_k, geodesic_threshold, geodesic_cutoff))
  } else {
    return(.computeDistanceMultiPairs(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                     normalizeDistance, truncateLowDist, verbose,
                                     knn_k, geodesic_threshold, geodesic_cutoff))
  }
}

# Helper function to get slide-specific location data and cell types
.getSlideData <- function(object, slideID) {
  slide_indices <- which(.getSlideIndices(object, slideID))
  list(
    locationData = object@locationDataSub[slide_indices, , drop = FALSE],
    cellTypes = object@cellTypesSub[slide_indices],
    cellIDs = rownames(object@locationDataSub)[slide_indices]
  )
}

# Helper function to initialize distance structure for multi-slide (flat)
.initializeDistanceStructureMulti <- function(slides, cts) {
  # Return empty list - will be populated with flat names
  return(list())
}

# Helper function to process multi-slide distance normalization
.normalizeDistancesMulti <- function(distances_all, slides, cts, global_min_percentile, 
                                    pair_cell_types = NULL, verbose = TRUE) {
  if (is.infinite(global_min_percentile)) {
    warning("Cannot normalize distances - no valid non-zero distances found across slides.")
    return(distances_all)
  }
  
  scaling_factor <- 0.01 / global_min_percentile
  if (verbose) cat("Global distance scaling factor:", scaling_factor, "\n")
  
  if (is.null(pair_cell_types)) {
    # Single cell type case
    ct <- cts
    for (sID in slides) {
      flat_name <- .createDistMatrixName(ct, ct, slide = sID)
      if (flat_name %in% names(distances_all) && !is.null(distances_all[[flat_name]])) {
        distances_all[[flat_name]] <- distances_all[[flat_name]] * scaling_factor
      }
    }
  } else {
    # Multiple cell types case
    for (sID in slides) {
      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]
        flat_name <- .createDistMatrixName(ct_i, ct_j, slide = sID)
        if (flat_name %in% names(distances_all) && !is.null(distances_all[[flat_name]])) {
          distances_all[[flat_name]] <- distances_all[[flat_name]] * scaling_factor
        }
      }
    }
  }
  
  return(distances_all)
}

# Function for computing within-cell-type distances across multiple slides
.computeDistanceMultiWithin <- function(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                       normalizeDistance, truncateLowDist, verbose,
                                       knn_k = 10, geodesic_threshold = 10,
                                       geodesic_cutoff = 7) {
  
  slides <- getSlideList(object)
  distances_all <- .initializeDistanceStructureMulti(slides, cts)
  global_min_percentile <- Inf
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance) {
    cat("normalizeDistance is set to TRUE, so distance will be",
        "normalized across all slides, so that 0.01 percentile distance will be scaled",
        "to 0.01\n")
  }
  
  for (sID in slides) {
    if (verbose) message(paste("Computing within-cell-type distances for slide:", sID))
    
    # Check if there are enough cells for this cell type in this slide
    slide_ct_count <- .countSlideCellType(object, slide = sID, cellType = cts)
    if (slide_ct_count <= 5) {
      if (verbose) message(paste("Skipping slide", sID, "- insufficient cells of type", cts, "(", slide_ct_count, "cells)"))
      next
    }
    
    # Get coordinate matrix for this slide and cell type
    mat1 <- .getCoordinateMatrix(object, cts, distType, xDistScale, yDistScale, zDistScale, slideID = sID)
    
    # Compute Euclidean distance matrix
    distances_ij <- fields::rdist(mat1)
    
    # Apply morphology-aware filtering if requested
    if (distType == "Morphology-Aware") {
      if (verbose) message("Computing KNN graph and geodesic distances for slide...")
      
      # Get coordinates for ALL cells in this slide
      slide_indices <- which(.getSlideIndices(object, sID))
      all_coords <- cbind(
        object@locationDataSub$x[slide_indices] * xDistScale,
        object@locationDataSub$y[slide_indices] * yDistScale
      )
      rownames(all_coords) <- rownames(object@locationDataSub)[slide_indices]
      cell_ids_all <- rownames(all_coords)
      
      # Compute KNN graph on all cells in this slide
      knn_adj <- .computeKnnGraph(all_coords, k = knn_k)
      
      # Compute geodesic distances
      geodesic_all <- .computeGeodesicDistance(knn_adj)
      
      if (verbose) message("Geodesic distance computation complete.")
      
      # Extract geodesic distances for this cell type
      idx <- match(rownames(mat1), cell_ids_all)
      geodesic_ij <- geodesic_all[idx, idx, drop = FALSE]
      
      if (verbose) message(sprintf("Applying morphology filter for cell type: %s", cts))
      
      # Apply filtering
      distances_ij <- .applyMorphologyFilter(
        d_E = distances_ij,
        d_g = geodesic_ij,
        geodesic_threshold = geodesic_threshold,
        geodesic_cutoff = geodesic_cutoff,
        verbose = verbose
      )
    }
    
    # Process distance matrix
    processed <- .processDistanceMatrix(distances_ij, truncateLowDist, 
                                       percentile_choice = 1e-4, set_diag_inf = TRUE)
    distances_ij <- processed$distances
    dist_percentile <- processed$percentile
    
    if (!is.na(dist_percentile) && is.finite(dist_percentile)) {
      global_min_percentile <- min(global_min_percentile, dist_percentile, na.rm = TRUE)
    }
    
    # Save the distances using flat structure
    flat_name <- .createDistMatrixName(cts, cts, slide = sID)
    distances_all[[flat_name]] <- distances_ij
    
    if (verbose) {
      cat("Slide:", sID, ", Cell type:", cts, "\n")
      if (verbose) print(quantile(distances_ij[is.finite(distances_ij)], na.rm = TRUE))
    }
  }
  
  # Apply normalization if requested
  if (normalizeDistance) {
    distances_all <- .normalizeDistancesMulti(distances_all, slides, cts, global_min_percentile, verbose = verbose)
  }
  
  object@distances <- distances_all
  return(object)
}

# Function for computing pairwise distances across multiple slides
.computeDistanceMultiPairs <- function(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                      normalizeDistance, truncateLowDist, verbose,
                                      knn_k = 10, geodesic_threshold = 10,
                                      geodesic_cutoff = 7) {
  
  slides <- getSlideList(object)
  distances_all <- .initializeDistanceStructureMulti(slides, cts)
  pair_cell_types <- combn(cts, 2)
  global_min_percentile <- Inf
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance) {
    cat("normalizeDistance is set to TRUE, so distance will be",
        "normalized across all slides, so that 0.01 percentile distance will be scaled",
        "to 0.01\n")
  }
  
  for (sID in slides) {
    if (verbose) message(paste("Computing pairwise distances for slide:", sID))
    
    # For Morphology-Aware: compute KNN and geodesic once per slide
    geodesic_all <- NULL
    cell_ids_all <- NULL
    if (distType == "Morphology-Aware") {
      if (verbose) message("Computing KNN graph and geodesic distances for slide...")
      
      # Get coordinates for ALL cells in this slide
      slide_indices <- which(.getSlideIndices(object, sID))
      all_coords <- cbind(
        object@locationDataSub$x[slide_indices] * xDistScale,
        object@locationDataSub$y[slide_indices] * yDistScale
      )
      rownames(all_coords) <- rownames(object@locationDataSub)[slide_indices]
      cell_ids_all <- rownames(all_coords)
      
      # Compute KNN graph on all cells in this slide
      knn_adj <- .computeKnnGraph(all_coords, k = knn_k)
      
      # Compute geodesic distances
      geodesic_all <- .computeGeodesicDistance(knn_adj)
      
      if (verbose) message("Geodesic distance computation complete.")
    }
    
    for (pp in seq_len(ncol(pair_cell_types))) {
      ct_i <- pair_cell_types[1, pp]
      ct_j <- pair_cell_types[2, pp]
      
      # Check cell counts for both cell types in this slide
      ct_i_count <- .countSlideCellType(object, slide = sID, cellType = ct_i)
      ct_j_count <- .countSlideCellType(object, slide = sID, cellType = ct_j)
      
      if (ct_i_count <= 5 || ct_j_count <= 5) {
        if (verbose) message(paste("Skipping pair", ct_i, "-", ct_j, 
                                  "in slide", sID, "(insufficient cells:", ct_i_count, "vs", ct_j_count, ")"))
        next
      }
      
      # Get coordinate matrices using the improved helper function
      mat1 <- .getCoordinateMatrix(object, ct_i, distType, xDistScale, yDistScale, zDistScale, slideID = sID)
      mat2 <- .getCoordinateMatrix(object, ct_j, distType, xDistScale, yDistScale, zDistScale, slideID = sID)
      
      # Compute Euclidean distance matrix
      distances_ij <- fields::rdist(mat1, mat2)
      
      # Apply morphology-aware filtering if requested
      if (distType == "Morphology-Aware") {
        # Extract geodesic distances for this cell type pair
        idx_i <- match(rownames(mat1), cell_ids_all)
        idx_j <- match(rownames(mat2), cell_ids_all)
        if (any(is.na(idx_i)) || any(is.na(idx_j))) {
          stop("Cell IDs in coordinate matrix do not match geodesic distance matrix. ",
               "Check that cell IDs are consistent across inputs.")
        }
        geodesic_ij <- geodesic_all[idx_i, idx_j, drop = FALSE]

        if (verbose) message(sprintf("Applying morphology filter for pair: %s - %s", ct_i, ct_j))
        
        # Apply filtering
        distances_ij <- .applyMorphologyFilter(
          d_E = distances_ij,
          d_g = geodesic_ij,
          geodesic_threshold = geodesic_threshold,
          geodesic_cutoff = geodesic_cutoff,
          verbose = verbose
        )
      }
      
      # Process distance matrix
      processed <- .processDistanceMatrix(distances_ij, truncateLowDist)
      distances_ij <- processed$distances
      dist_percentile <- processed$percentile
      
      if (!is.na(dist_percentile) && is.finite(dist_percentile)) {
        global_min_percentile <- min(global_min_percentile, dist_percentile, na.rm = TRUE)
      }
      
      # Save the distances using flat structure
      flat_name <- .createDistMatrixName(ct_i, ct_j, slide = sID)
      distances_all[[flat_name]] <- distances_ij
      
      if (verbose) {
        cat("Slide:", sID, ", Pair:", ct_i, "-", ct_j, "\n")
        if (verbose) print(quantile(distances_ij, na.rm = TRUE))
      }
    }
  }
  
  # Apply normalization if requested
  if (normalizeDistance) {
    distances_all <- .normalizeDistancesMulti(distances_all, slides, cts, global_min_percentile, 
                                             pair_cell_types, verbose = verbose)
  }
  
  object@distances <- distances_all
  return(object)
}
