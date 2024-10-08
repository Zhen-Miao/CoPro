
# Function to find neighboring bins
.get_neighbor_bins <- function(location_data, tar_bin_id, num_bins_x, num_bins_y) {
  x_bin <- location_data[tar_bin_id, "x_bin"]
  y_bin <- location_data[tar_bin_id, "y_bin"]
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


resample_spatial <- function(location_data,
                             num_bins_x = 10, num_bins_y = 10) {

  # Ensure the location_data has x and y columns
  if (is.matrix(location_data) || is.data.frame(location_data)) {
    # Convert to location_data frame if it's a matrix
    if (is.matrix(location_data)) {
      location_data <- as.data.frame(location_data)
    }
    # Check for x and y columns
    if (!all(c("x", "y", "cell_ID") %in% colnames(location_data))) {
      stop("Input location_data must have 'x', 'y', and 'cell_ID' columns.")
    }
  } else {
    stop("Input location_data must be a matrix or location_data frame.")
  }

  original_cell_loc_order <- paste(location_data$x, location_data$y, sep = "_")
  rownames(location_data) <- original_cell_loc_order

  # Create bins for x and y coordinates
  if(!all(c("x_bin", "y_bin") %in% colnames(location_data))){
    location_data$x_bin <- cut(location_data$x, breaks = num_bins_x, labels = FALSE)
    location_data$y_bin <- cut(location_data$y, breaks = num_bins_y, labels = FALSE)
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
      neighbors <- .get_neighbor_bins(location_data = location_data,
                                      tar_bin_id = target_bin,
                                      num_bins_x = num_bins_x,
                                      num_bins_y = num_bins_y)
      # Exclude bins already in candidate_bins
      additional_bins <- setdiff(neighbors, candidate_bins)
      for (neighbor_bin in additional_bins) {
        neighbor_points <- location_data[location_data$bin_id == neighbor_bin, ]
        candidate_points <- rbind(candidate_points, neighbor_points)
        candidate_bins <- c(candidate_bins, neighbor_bin)
        if (nrow(candidate_points) >= n_points) {
          break
        }
      }
    }

    # If still not enough points, sample with replacement as a last resort
    if (nrow(candidate_points) < n_points) {
      sampled_indices <- sample(nrow(candidate_points), n_points, replace = TRUE)
    } else {
      sampled_indices <- sample(nrow(candidate_points), n_points)
    }

    # Replace the coordinates in the original points
    orig_points$"cell_ID" <- candidate_points[sampled_indices,"cell_ID"]
    orig_points$"bin_id" <- candidate_points[sampled_indices,"bin_id"]

    # Store the resampled points
    resampled_list[[i]] <- orig_points
  }

  # Combine all resampled points
  resampled_location_data <- do.call(rbind, resampled_list)

  # Ensure the row names are in the original order
  rownames(resampled_location_data) <- paste(resampled_location_data$x,
                                             resampled_location_data$y,
                                             sep = "_")
  resampled_location_data <- resampled_location_data[original_cell_loc_order,]

  # # Remove auxiliary columns
  # resampled_location_data$bin_id <- NULL
  # resampled_location_data$x_bin <- NULL
  # resampled_location_data$y_bin <- NULL

  return(resampled_location_data)
}

#
# # Sample location_d
# set.seed(123)
# location_data <- data.frame(
#   x = runif(1000, 0, 100),
#   y = runif(1000, 0, 100),
#   cell_ID = paste0("c_1", 1:1000)
#
# )
#
#
# # Create bins for x and y coordinates
# location_data$x_bin1 <- cut(location_data$x, breaks = num_bins_x, labels = FALSE)
# location_data$y_bin1 <- cut(location_data$y, breaks = num_bins_y, labels = FALSE)
#
# # Assign a bin ID to each point
# location_data$bin_id1 <- paste(location_data$x_bin1, location_data$y_bin1, sep = "_")
#
#
# # Resample the location_d
# resampled_location_data <- resample_spatial(location_data, num_bins_x = 10, num_bins_y = 10)
#
#
# library(ggplot2)
# ggplot(location_data)+
#   geom_point(aes(x = x, y = y, color = bin_id1))+
#   theme_minimal()
# ggsave("orginal simulated data.pdf")
#
# ggplot(resampled_location_data)+
#   geom_point(aes(x = x, y = y, color = factor(bin_id)))+
#   theme_minimal()
# ggsave("resampled simulated data.pdf")



