

## 0. rescale the scores to (0,1)
rescale_scores <- function(current_scores){
  scaled_scores <- (current_scores - min(current_scores)) /
    (max(current_scores) - min(current_scores))
  return(scaled_scores)
}


#### 1. simulate points in the space, and assign a score




simulate_smooth_points <- function(n_points,
                                   x_range = c(0, 10),
                                   y_range = c(0, 10),
                                   bandwidth = 1.5,
                                   n_rounds = 1,
                                   score_sd = 0.2,
                                   rescale = FALSE,
                                   seed = NULL,
                                   labels, label_prob = NULL,
                                   label_name = "label"
                                   ) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (is.null(label_prob)) {
    # If no label_prob provided, use equal label_prob
    label_prob <- rep(1/length(labels), length(labels))
  }

  # Additional validation
  if (length(labels) != length(label_prob)) {
    stop("Length of labels must match length of label_prob")
  }
  if (abs(sum(label_prob) - 1) > 1e-10) {
    stop("label_prob must sum to 1")
  }
  if (any(label_prob < 0)) {
    stop("label_prob must be non-negative")
  }

  # 1. Simulate points in 2D space
  x <- runif(n_points, x_range[1], x_range[2])
  y <- runif(n_points, y_range[1], y_range[2])
  points <- data.frame(x = x, y = y)

  # 2. Assign initial random scores between 0 and 1
  points$initial_score <- rnorm(n_points, mean = 0, sd = score_sd)

  # 3. Calculate distance matrix efficiently
  dist_matrix <- as.matrix(dist(points[, c("x", "y")]))

  # 4. Calculate Gaussian kernel matrix
  kernel_matrix <- exp(-0.5 * (dist_matrix / bandwidth)^2)
  kernel_matrix <- kernel_matrix / rowSums(kernel_matrix)

  # 5. Perform spatial smoothing for specified number of rounds
  current_scores <- points$initial_score
  for(i in 1:n_rounds) {
    current_scores <- as.vector(kernel_matrix %*% current_scores)
  }

  # 6. Put the score back to the 0-1 range
  if(rescale){
    current_scores <- rescale_scores(current_scores)
  }

  # Add smoothed scores to the data frame
  points$smoothed_score <- current_scores


  # Sample from labels using provided label_prob
  sampled_labels <- sample(labels,
                           size = n_points,
                           replace = TRUE,
                           prob = label_prob)
  points[[label_name]] = sampled_labels

  return(points)
}




## calculate norm corr for the data
norm_corr_sim <- function(points, label_name = "label",
                          score_name = "smoothed_score",
                          bandwidth, rm_diag = FALSE,
                          ct1, ct2, row_norm_kernel = TRUE){

  if(length(points[[label_name]]) == 0){
    stop("the label does not exist, check label_name input")
  }else{
    lb <- points[[label_name]]
  }

  if(rm_diag && (length(unique(lb)) != 1 || ct1 != ct2 ) ){
    stop("please set rm_diag to TRUE only when there are one cell type")
  }

  # distance matrix
  dist_matrix <- fields::rdist( x1 = as.matrix(points[lb == ct1, c("x", "y")]),
                 x2 = as.matrix(points[lb == ct2, c("x", "y")]))

  # Gaussian kernel matrix
  kernel_matrix <- exp(-0.5 * (dist_matrix / bandwidth)^2)

  if(rm_diag){
    diag(kernel_matrix) <- 0
  }

  if(row_norm_kernel){
    rs_kernel <- rowSums(kernel_matrix)
    nz_ind <- rs_kernel > 1e-4
    kernel_matrix[nz_ind,] <- kernel_matrix[nz_ind,] / rs_kernel[nz_ind]
  }

  s1T = matrix(normalize_vec(points[[score_name]][lb == ct1]), nrow = 1)
  s2 = matrix(normalize_vec(points[[score_name]][lb == ct2]), ncol = 1)

  svd_result <- irlba::irlba(kernel_matrix, nv = 1, tol = 1e-5)
  norm_K12 <- svd_result$d[1]

  norm_corr <- (s1T %*% kernel_matrix %*% s2 ) / norm_K12

  return(norm_corr)
}




generate_prob_vector <- function(length, seed = NULL) {

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate random positive numbers
  raw_vector <- runif(length)

  # Normalize the vector so it sums to 1
  prob_vector <- raw_vector / sum(raw_vector)

  return(prob_vector)
}


## ============================================================================
## Null Simulation Functions
## ============================================================================
## These functions simulate spatial data where cells have spatial autocorrelation
## WITHIN their own cell type, but NO coordination ACROSS cell types.
## This is useful for evaluating the null distribution of co-progression statistics.
## ============================================================================

#' Simulate Spatial Points with Null Hypothesis (No Cross-Type Coordination)
#'
#' This function simulates points in 2D space where each cell type has its own
#' independent smoothed scores. Cells have spatial autocorrelation within their
#' own type, but there is no coordination across cell types.
#'
#' @param n_points Total number of points to simulate
#' @param x_range Range for x coordinates (default: c(0, 10))
#' @param y_range Range for y coordinates (default: c(0, 10))
#' @param bandwidth Bandwidth for Gaussian kernel smoothing
#' @param n_rounds Number of smoothing rounds (default: 1)
#' @param score_sd Standard deviation for initial random scores (default: 0.2)
#' @param rescale Whether to rescale scores to (0,1) (default: FALSE)
#' @param seed Random seed for reproducibility
#' @param labels Vector of cell type labels
#' @param label_prob Probability of each label (default: equal probabilities)
#' @param label_name Name of the label column (default: "label")
#'
#' @return A data.frame with columns: x, y, initial_score, smoothed_score, and cell type
#'
simulate_smooth_points_null <- function(n_points,
                                        x_range = c(0, 10),
                                        y_range = c(0, 10),
                                        bandwidth = 1.5,
                                        n_rounds = 1,
                                        score_sd = 0.2,
                                        rescale = FALSE,
                                        seed = NULL,
                                        labels,
                                        label_prob = NULL,
                                        label_name = "label") {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Input validation
  if (is.null(label_prob)) {
    label_prob <- rep(1/length(labels), length(labels))
  }
  if (length(labels) != length(label_prob)) {
    stop("Length of labels must match length of label_prob")
  }
  if (abs(sum(label_prob) - 1) > 1e-10) {
    stop("label_prob must sum to 1")
  }
  if (any(label_prob < 0)) {
    stop("label_prob must be non-negative")
  }

  # 1. Simulate points in 2D space
  x <- runif(n_points, x_range[1], x_range[2])
  y <- runif(n_points, y_range[1], y_range[2])
  points <- data.frame(x = x, y = y)

  # 2. Assign cell type labels FIRST
  sampled_labels <- sample(labels,
                           size = n_points,
                           replace = TRUE,
                           prob = label_prob)
  points[[label_name]] <- sampled_labels

  # 3. Initialize columns for scores
  points$initial_score <- NA_real_
  points$smoothed_score <- NA_real_

  # 4. For each cell type, assign INDEPENDENT initial scores and smooth WITHIN type

  for (lbl in labels) {
    idx <- which(points[[label_name]] == lbl)
    n_type <- length(idx)

    if (n_type == 0) next

    # Assign independent initial scores for this cell type
    points$initial_score[idx] <- rnorm(n_type, mean = 0, sd = score_sd)

    # Calculate distance matrix for this cell type only
    type_coords <- points[idx, c("x", "y")]
    dist_matrix <- as.matrix(dist(type_coords))

    # Calculate Gaussian kernel matrix
    kernel_matrix <- exp(-0.5 * (dist_matrix / bandwidth)^2)
    kernel_matrix <- kernel_matrix / rowSums(kernel_matrix)

    # Perform spatial smoothing within this cell type
    current_scores <- points$initial_score[idx]
    for (i in 1:n_rounds) {
      current_scores <- as.vector(kernel_matrix %*% current_scores)
    }

    # Rescale if requested
    if (rescale && length(unique(current_scores)) > 1) {
      current_scores <- rescale_scores(current_scores)
    }

    # Store smoothed scores
    points$smoothed_score[idx] <- current_scores
  }

  return(points)
}


#' Calculate Null Distribution of Normalized Correlation via Permutation
#'
#' This function calculates the null distribution by repeatedly simulating
#' data with no cross-type coordination and computing normalized correlation.
#'
#' @param n_sim Number of simulations (default: 100)
#' @param n_points Number of points per simulation
#' @param x_range Range for x coordinates
#' @param y_range Range for y coordinates
#' @param bandwidth_sim Bandwidth for simulation smoothing
#' @param bandwidth_test Bandwidth(s) for computing normalized correlation
#' @param n_rounds Number of smoothing rounds
#' @param score_sd Standard deviation for initial scores
#' @param labels Vector of cell type labels
#' @param label_prob Probability of each label
#' @param label_name Name of the label column
#' @param row_norm_kernel Whether to row-normalize kernel matrix
#' @param seed_base Base seed for reproducibility
#' @param verbose Whether to print progress (default: TRUE)
#'
#' @return A data.frame with simulation results
#'
simulate_null_distribution <- function(n_sim = 100,
                                       n_points,
                                       x_range = c(0, 10),
                                       y_range = c(0, 10),
                                       bandwidth_sim = 1.5,
                                       bandwidth_test = c(0.1, 0.2, 0.4),
                                       n_rounds = 1,
                                       score_sd = 0.2,
                                       labels,
                                       label_prob = NULL,
                                       label_name = "label",
                                       row_norm_kernel = TRUE,
                                       seed_base = NULL,
                                       verbose = TRUE) {

  results <- data.frame()

  for (i in 1:n_sim) {
    if (verbose && i %% 10 == 0) {
      message(paste("Simulation", i, "of", n_sim))
    }

    # Set seed for this iteration
    seed_i <- if (!is.null(seed_base)) seed_base + i else NULL

    # Simulate null data
    sim_data <- simulate_smooth_points_null(
      n_points = n_points,
      x_range = x_range,
      y_range = y_range,
      bandwidth = bandwidth_sim,
      n_rounds = n_rounds,
      score_sd = score_sd,
      seed = seed_i,
      labels = labels,
      label_prob = label_prob,
      label_name = label_name
    )

    # Compute normalized correlation for each bandwidth
    for (bw in bandwidth_test) {
      nc <- tryCatch({
        norm_corr_sim(
          points = sim_data,
          label_name = label_name,
          score_name = "smoothed_score",
          bandwidth = bw,
          ct1 = labels[1],
          ct2 = labels[2],
          row_norm_kernel = row_norm_kernel
        )
      }, error = function(e) NA)

      results <- rbind(results, data.frame(
        sim_id = i,
        bandwidth = bw,
        norm_corr = as.numeric(nc)
      ))
    }
  }

  return(results)
}


match_by_percentile <- function(sim_df, real_cell_score,rm_outlier = TRUE) {
  ## make sure names exist
  if(length(names(real_cell_score)) == 0){
    stop("real_cell_score must have cell indices as names")
  }

  if(rm_outlier){
    sd_score = sd(real_cell_score)
    rs <- real_cell_score
    rs = rs[rs > (-3 * sd_score ) & rs < (3 * sd_score)]
    real_cell_score <- rs
  }

  # Get the number of simulated points
  n_sim <- nrow(sim_df)
  n_real <- length(real_cell_score)

  if(n_sim > n_real) {
    stop(paste("Number of simulated points cannot",
               "exceed number of real cells (after filtering)"))
  }

  # Calculate percentiles for simulated points and real cells
  sim_percentiles <- rank(sim_df$smoothed_score) / n_sim
  real_percentiles <- rank(real_cell_score) / n_real

  # First match: get closest percentile match for each point
  distances <- outer(sim_percentiles, real_percentiles, FUN = "-")
  matched_indices <- apply(abs(distances), 1, which.min)

  # Check for duplicates
  duplicate_matches <- duplicated(matched_indices)

  if(any(duplicate_matches)) {
    # For points with duplicate matches, find alternative matches
    duplicate_points <- which(duplicate_matches)

    # Create mask of already used indices
    used_indices <- unique(matched_indices[!duplicate_matches])

    # For each duplicate, find the nearest unused cell
    for(i in duplicate_points) {
      # Get distances to all cells for this point
      point_distances <- abs(distances[i,])
      # Order by distance, excluding used indices
      available_order <- order(point_distances)
      # Find first unused index
      new_match <- available_order[!available_order %in% used_indices][1]

      # Update matches
      matched_indices[i] <- new_match
      used_indices <- c(used_indices, new_match)
    }
  }

  # Add matched cell names and scores
  sim_df$matched_cell <- names(real_cell_score)[matched_indices]
  sim_df$matched_score <- real_cell_score[matched_indices]

  return(sim_df)
}

