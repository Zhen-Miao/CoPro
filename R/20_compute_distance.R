#' computeDistance between pairs of cell types
#'
#' @importFrom fields rdist
#' @importFrom utils combn
#' @importFrom stats setNames quantile
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
#' @return `CoPro` object with distance matrix computed
#' @rdname computeDistance
#' @aliases computeDistance,CoProSingle-method
#' @aliases computeDistance,CoProMulti-method
#' @export
#' @note To-do: add morphology-aware kernel
setGeneric(
  "computeDistance",
  function(object, distType =
             c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
           xDistScale = 1, yDistScale = 1,
           zDistScale = 1, normalizeDistance = TRUE, truncateLowDist = TRUE,
           verbose = TRUE) standardGeneric("computeDistance")
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
    stop("morphology-aware kernel is not available at this moment")
  }
  return(TRUE)
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
                                normalizeDistance, truncateLowDist, verbose) {
  cts <- .checkInputDistance(object, distType, xDistScale, yDistScale, zDistScale)
  
  # Determine whether to compute pairwise or within-cell-type distances
  if (length(cts) == 1) {
    return(.computeDistanceWithin(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                 normalizeDistance, truncateLowDist, verbose))
  } else {
    return(.computeDistancePairs(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                normalizeDistance, truncateLowDist, verbose))
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
    ct_indices <- object@cellTypesSub == cellType & object@metaDataSub$slideID == slideID
  }
  
  # Check if any cells were found
  if (!any(ct_indices)) {
    if (is.null(slideID)) {
      stop(paste("No cells found for cell type:", cellType))
    } else {
      stop(paste("No cells found for cell type:", cellType, "in slide:", slideID))
    }
  }
  
  if (distType == "Euclidean2D") {
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
  } else if (distType == "Morphology-Aware") {
    stop("morphology-aware kernel is not available at this moment")
  }
  
  # Set rownames to match the cell IDs
  if (!is.null(slideID)) {
    slide_cell_ids <- rownames(object@locationDataSub)[object@metaDataSub$slideID == slideID]
    ct_cell_ids <- slide_cell_ids[object@cellTypesSub[object@metaDataSub$slideID == slideID] == cellType]
    rownames(mat) <- ct_cell_ids
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
                                 normalizeDistance, truncateLowDist, verbose) {
  
  # Initialize distance list structure
  distances <- setNames(rep(list(), length = length(cts)), cts)
  for (i in cts) {
    distances[[i]] <- setNames(rep(list(), length = length(cts)), cts)
  }
  
  pair_cell_types <- combn(cts, 2)
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance) {
    cat("normalizeDistance is set to TRUE, so distance will be",
        "normalized, so that 0.01 percentile distance will be scaled",
        "to 0.01\n")
  }
  
  dist_percentiles <- vector(mode = "numeric", length = ncol(pair_cell_types))
  
  # Calculate distances for each pair
  for (pp in seq_len(ncol(pair_cell_types))) {
    i <- pair_cell_types[1, pp]
    j <- pair_cell_types[2, pp]
    
    # Get coordinate matrices
    mat1 <- .getCoordinateMatrix(object, i, distType, xDistScale, yDistScale, zDistScale)
    mat2 <- .getCoordinateMatrix(object, j, distType, xDistScale, yDistScale, zDistScale)
    
    # Compute distance
    distances_ij <- fields::rdist(mat1, mat2)
    
    # Process distance matrix
    processed <- .processDistanceMatrix(distances_ij, truncateLowDist)
    distances_ij <- processed$distances
    dist_percentiles[pp] <- processed$percentile
    
    # Save the distances
    distances[[i]][[j]] <- distances_ij
    
    if (verbose) {
      cat("quantile of the distances between", i, "and", j, "is: \n")
      print(quantile(distances_ij, na.rm = TRUE))
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
      distances[[i]][[j]] <- distances[[i]][[j]] * scaling_factor
    }
  }
  
  object@distances <- distances
  return(object)
}

# Function for computing within-cell-type distances
.computeDistanceWithin <- function(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                  normalizeDistance, truncateLowDist, verbose) {
  
  # Initialize distance list structure
  distances <- setNames(rep(list(), length = length(cts)), cts)
  for (i in cts) {
    distances[[i]] <- setNames(rep(list(), length = length(cts)), cts)
  }
  
  # Notify users if normalizeDistance = TRUE
  if (normalizeDistance) {
    cat("normalizeDistance is set to TRUE, so distance will be",
        "normalized, so that 0.01 percentile distance will be scaled",
        "to 0.01\n")
  }
  
  # Get coordinate matrix for the single cell type
  mat1 <- .getCoordinateMatrix(object, cts, distType, xDistScale, yDistScale, zDistScale)
  
  # Compute distance matrix
  distances_ij <- fields::rdist(mat1)
  
  # Process distance matrix (set diagonal to Inf)
  processed <- .processDistanceMatrix(distances_ij, truncateLowDist, 
                                     percentile_choice = 1e-4, set_diag_inf = TRUE)
  distances_ij <- processed$distances
  dist_percentile <- processed$percentile
  
  # Save the distances
  distances[[cts]][[cts]] <- distances_ij
  
  if (verbose) {
    cat("quantile of the distances within", cts, "is: \n")
    print(quantile(distances_ij[is.finite(distances_ij)], na.rm = TRUE))
  }
  
  # Apply normalization if requested
  if (normalizeDistance) {
    scaling_factor <- 0.01 / dist_percentile
    cat("The scaling factor for normalizing distance is", scaling_factor, "\n")
    distances[[cts]][[cts]] <- distances_ij * scaling_factor
  }
  
  object@distances <- distances
  return(object)
}

#' @rdname computeDistance
#' @aliases computeDistance,CoProSingle-method
#' @export
setMethod("computeDistance", "CoProSingle", function(object, distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
                                                    xDistScale = 1, yDistScale = 1, zDistScale = 1,
                                                    normalizeDistance = TRUE, truncateLowDist = TRUE, verbose = TRUE) {
  distType <- match.arg(distType)
  .computeDistanceCore(object, distType, xDistScale, yDistScale, zDistScale,
                      normalizeDistance, truncateLowDist, verbose)
})

#' @rdname computeDistance
#' @aliases computeDistance,CoProMulti-method
#' @export
setMethod("computeDistance", "CoProMulti", function(object, distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
                                                   xDistScale = 1, yDistScale = 1, zDistScale = 1,
                                                   normalizeDistance = TRUE, truncateLowDist = TRUE, verbose = TRUE) {
  distType <- match.arg(distType)
  .computeDistanceCoreMulti(object, distType, xDistScale, yDistScale, zDistScale,
                           normalizeDistance, truncateLowDist, verbose)
})

# Core dispatcher for multi-slide objects
.computeDistanceCoreMulti <- function(object, distType, xDistScale, yDistScale, zDistScale,
                                     normalizeDistance, truncateLowDist, verbose) {
  cts <- .checkInputDistance(object, distType, xDistScale, yDistScale, zDistScale)
  
  # Determine whether to compute pairwise or within-cell-type distances across slides
  if (length(cts) == 1) {
    return(.computeDistanceMultiWithin(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                      normalizeDistance, truncateLowDist, verbose))
  } else {
    return(.computeDistanceMultiPairs(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                     normalizeDistance, truncateLowDist, verbose))
  }
}

# Helper function to get slide-specific location data and cell types
.getSlideData <- function(object, slideID) {
  slide_indices <- which(object@metaDataSub$slideID == slideID)
  list(
    locationData = object@locationDataSub[slide_indices, , drop = FALSE],
    cellTypes = object@cellTypesSub[slide_indices],
    cellIDs = rownames(object@locationDataSub)[slide_indices]
  )
}

# Helper function to initialize distance structure for multi-slide
.initializeDistanceStructureMulti <- function(slides, cts) {
  distances_all <- setNames(vector("list", length = length(slides)), slides)
  for (sID in slides) {
    distances_slide <- setNames(vector("list", length = length(cts)), cts)
    for (i in cts) {
      distances_slide[[i]] <- setNames(vector("list", length = length(cts)), cts)
    }
    distances_all[[sID]] <- distances_slide
  }
  return(distances_all)
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
      if (!is.null(distances_all[[sID]][[ct]][[ct]])) {
        distances_all[[sID]][[ct]][[ct]] <- distances_all[[sID]][[ct]][[ct]] * scaling_factor
      }
    }
  } else {
    # Multiple cell types case
    for (sID in slides) {
      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]
        if (!is.null(distances_all[[sID]][[ct_i]][[ct_j]])) {
          distances_all[[sID]][[ct_i]][[ct_j]] <- distances_all[[sID]][[ct_i]][[ct_j]] * scaling_factor
        }
      }
    }
  }
  
  return(distances_all)
}

# Function for computing within-cell-type distances across multiple slides
.computeDistanceMultiWithin <- function(object, cts, distType, xDistScale, yDistScale, zDistScale,
                                       normalizeDistance, truncateLowDist, verbose) {
  
  slides <- object@slideList
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
    slide_ct_count <- sum(object@cellTypesSub == cts & object@metaDataSub$slideID == sID)
    if (slide_ct_count <= 5) {
      if (verbose) message(paste("Skipping slide", sID, "- insufficient cells of type", cts, "(", slide_ct_count, "cells)"))
      next
    }
    
    # Get coordinate matrix for this slide and cell type
    mat1 <- .getCoordinateMatrix(object, cts, distType, xDistScale, yDistScale, zDistScale, slideID = sID)
    
    # Compute distance matrix
    distances_ij <- fields::rdist(mat1)
    
    # Process distance matrix
    processed <- .processDistanceMatrix(distances_ij, truncateLowDist, 
                                       percentile_choice = 1e-4, set_diag_inf = TRUE)
    distances_ij <- processed$distances
    dist_percentile <- processed$percentile
    
    if (!is.na(dist_percentile) && is.finite(dist_percentile)) {
      global_min_percentile <- min(global_min_percentile, dist_percentile, na.rm = TRUE)
    }
    
    # Save the distances
    distances_all[[sID]][[cts]][[cts]] <- distances_ij
    
    if (verbose) {
      cat("Slide:", sID, ", Cell type:", cts, "\n")
      print(quantile(distances_ij[is.finite(distances_ij)], na.rm = TRUE))
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
                                      normalizeDistance, truncateLowDist, verbose) {
  
  slides <- object@slideList
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
    
    for (pp in seq_len(ncol(pair_cell_types))) {
      ct_i <- pair_cell_types[1, pp]
      ct_j <- pair_cell_types[2, pp]
      
      # Check cell counts for both cell types in this slide
      ct_i_count <- sum(object@cellTypesSub == ct_i & object@metaDataSub$slideID == sID)
      ct_j_count <- sum(object@cellTypesSub == ct_j & object@metaDataSub$slideID == sID)
      
      if (ct_i_count <= 5 || ct_j_count <= 5) {
        if (verbose) message(paste("Skipping pair", ct_i, "-", ct_j, 
                                  "in slide", sID, "(insufficient cells:", ct_i_count, "vs", ct_j_count, ")"))
        next
      }
      
      # Get coordinate matrices using the improved helper function
      mat1 <- .getCoordinateMatrix(object, ct_i, distType, xDistScale, yDistScale, zDistScale, slideID = sID)
      mat2 <- .getCoordinateMatrix(object, ct_j, distType, xDistScale, yDistScale, zDistScale, slideID = sID)
      
      # Compute distance matrix
      distances_ij <- fields::rdist(mat1, mat2)
      
      # Process distance matrix
      processed <- .processDistanceMatrix(distances_ij, truncateLowDist)
      distances_ij <- processed$distances
      dist_percentile <- processed$percentile
      
      if (!is.na(dist_percentile) && is.finite(dist_percentile)) {
        global_min_percentile <- min(global_min_percentile, dist_percentile, na.rm = TRUE)
      }
      
      # Save the distances
      distances_all[[sID]][[ct_i]][[ct_j]] <- distances_ij
      
      if (verbose) {
        cat("Slide:", sID, ", Pair:", ct_i, "-", ct_j, "\n")
        print(quantile(distances_ij, na.rm = TRUE))
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
