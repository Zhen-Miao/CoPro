library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)

#' Check package versions for UMAP compatibility
#'
#' @param verbose Whether to print detailed information
#' @return List of package versions
#' @export
check_umap_packages <- function(verbose = TRUE) {
  
  packages <- c("Matrix", "irlba", "uwot", "RSpectra")
  versions <- list()
  
  for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      versions[[pkg]] <- as.character(packageVersion(pkg))
    } else {
      versions[[pkg]] <- "Not installed"
    }
  }
  
  if (verbose) {
    cat("Package versions for UMAP:\n")
    for (pkg in names(versions)) {
      cat(sprintf("  %s: %s\n", pkg, versions[[pkg]]))
    }
    
    # Check for known problematic combinations
    matrix_version <- versions[["Matrix"]]
    irlba_version <- versions[["irlba"]]
    
    if (matrix_version != "Not installed" && 
        compareVersion(matrix_version, "1.5.0") >= 0) {
      cat("\nNote: Matrix version >= 1.5.0 may cause compatibility issues.\n")
      cat("If you encounter 'as_cholmod_sparse' errors, try:\n")
      cat("  install.packages('Matrix', version = '1.4-1', repos = 'https://cran.r-project.org')\n")
      cat("  # or update all packages:\n")
      cat("  install.packages(c('Matrix', 'irlba', 'uwot'), force = TRUE)\n")
    }
  }
  
  return(invisible(versions))
}


#' Convert CCA components to UMAP coordinates (Seurat-style)
#'
#' @param cca_matrix Matrix of shape (n_cells, n_canonical_components)
#' @param n_components Number of CCA components to use (default: 4)
#' @param n_neighbors UMAP parameter for local neighborhood size (default: 30)
#' @param min_dist UMAP parameter for minimum distance between points (default: 0.3)
#' @param metric Distance metric for UMAP (default: 'euclidean')
#' @param n_epochs Number of epochs for UMAP optimization (default: NULL, auto-determined)
#' @param learning_rate Learning rate for UMAP optimization (default: 1)
#' @param spread Effective scale of embedded points (default: 1)
#' @param negative_sample_rate Negative sampling rate (default: 5)
#' @param a Parameter controlling embedding spread (default: NULL, calculated from spread/min_dist)
#' @param b Parameter controlling embedding spread (default: NULL, calculated from spread/min_dist)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Show progress messages (default: TRUE)
#' @param n_threads Number of threads to use (default: 0 for auto)
#' @param ... Additional parameters to pass to uwot::umap
#'
#' @return Matrix of UMAP coordinates with shape (n_cells, 2)
#' @importFrom uwot umap
#' @export
cca_to_umap <- function(cca_matrix,
                        n_components = 4,
                        n_neighbors = 30,
                        min_dist = 0.3,
                        metric = "euclidean",
                        n_epochs = NULL,
                        learning_rate = 1,
                        spread = 1,
                        negative_sample_rate = 5,
                        a = NULL,
                        b = NULL,
                        seed = 42,
                        verbose = TRUE,
                        n_threads = 0,
                        ...) {
  
  # Check if uwot is installed
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required. Please install it with: install.packages('uwot')")
  }
  
  # Validate input
  if (!is.matrix(cca_matrix) && !is.data.frame(cca_matrix)) {
    stop("cca_matrix must be a matrix or data frame")
  }
  
  # Convert to matrix if data frame
  if (is.data.frame(cca_matrix)) {
    cca_matrix <- as.matrix(cca_matrix)
  }
  
  if (ncol(cca_matrix) < n_components) {
    stop(sprintf("CCA matrix has only %d components, but %d were requested",
                 ncol(cca_matrix), n_components))
  }
  
  # Select first n_components
  cca_subset <- cca_matrix[, 1:n_components, drop = FALSE]
  
  # Print status if verbose
  if (verbose) {
    message(sprintf("Running UMAP on %d cells with %d CCA components", 
                    nrow(cca_subset), n_components))
  }
  
  # Set random seed
  set.seed(seed)
  
  # Determine number of threads
  if (n_threads == 0) {
    n_threads <- max(1, parallel::detectCores() - 1)
  }
  
  # Run UMAP with uwot (Seurat's backend)
  # Try different initialization methods if Matrix package compatibility issues occur
  umap_coords <- NULL
  
  # First try with default initialization (normalized Laplacian)
  if (is.null(umap_coords)) {
    umap_coords <- tryCatch({
      if (verbose) message("Trying UMAP with normalized Laplacian initialization...")
      uwot::umap(
        X = cca_subset,
        n_neighbors = n_neighbors,
        n_components = 2,
        metric = metric,
        n_epochs = n_epochs,
        learning_rate = learning_rate,
        min_dist = min_dist,
        spread = spread,
        negative_sample_rate = negative_sample_rate,
        a = a,
        b = b,
        n_threads = n_threads,
        verbose = verbose,
        n_sgd_threads = n_threads,
        grain_size = 1,
        init = "normlaplacian",
        ...
      )
    }, error = function(e) {
      if (verbose) message("Laplacian initialization failed: ", e$message)
      NULL
    })
  }
  
  # If that fails, try with spectral initialization
  if (is.null(umap_coords)) {
    umap_coords <- tryCatch({
      if (verbose) message("Trying UMAP with spectral initialization...")
      uwot::umap(
        X = cca_subset,
        n_neighbors = n_neighbors,
        n_components = 2,
        metric = metric,
        n_epochs = n_epochs,
        learning_rate = learning_rate,
        min_dist = min_dist,
        spread = spread,
        negative_sample_rate = negative_sample_rate,
        a = a,
        b = b,
        n_threads = n_threads,
        verbose = verbose,
        n_sgd_threads = n_threads,
        grain_size = 1,
        init = "spectral",
        ...
      )
    }, error = function(e) {
      if (verbose) message("Spectral initialization failed: ", e$message)
      NULL
    })
  }
  
  # If that fails, try with random initialization
  if (is.null(umap_coords)) {
    umap_coords <- tryCatch({
      if (verbose) message("Trying UMAP with random initialization...")
      uwot::umap(
        X = cca_subset,
        n_neighbors = n_neighbors,
        n_components = 2,
        metric = metric,
        n_epochs = n_epochs,
        learning_rate = learning_rate,
        min_dist = min_dist,
        spread = spread,
        negative_sample_rate = negative_sample_rate,
        a = a,
        b = b,
        n_threads = n_threads,
        verbose = verbose,
        n_sgd_threads = n_threads,
        grain_size = 1,
        init = "random",
        ...
      )
    }, error = function(e) {
      if (verbose) message("Random initialization failed: ", e$message)
      NULL
    })
  }
  
  # If all initialization methods fail, try with reduced threads and simpler parameters
  if (is.null(umap_coords)) {
    if (verbose) message("Trying UMAP with simplified parameters...")
    umap_coords <- tryCatch({
      uwot::umap(
        X = cca_subset,
        n_neighbors = min(n_neighbors, 15),  # Reduce neighbors
        n_components = 2,
        metric = "euclidean",  # Force euclidean
        n_epochs = if(is.null(n_epochs)) 200 else min(n_epochs, 200),
        learning_rate = 1,
        min_dist = min_dist,
        spread = 1,
        negative_sample_rate = 5,
        n_threads = 1,  # Single thread
        verbose = verbose,
        init = "random"
      )
    }, error = function(e) {
      stop("UMAP failed with all initialization methods. Error: ", e$message,
           "\n\nThis is likely due to Matrix package version incompatibility.",
           "\nTry: install.packages(c('Matrix', 'irlba', 'uwot'), force = TRUE)")
    })
  }
  
  # Add column names
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  # Return coordinates
  return(umap_coords)
}



#' Plot UMAP coordinates with optional labels
#'
#' @param umap_coords UMAP coordinates from cca_to_umap function
#' @param labels Labels for coloring points (e.g., cell types, clusters)
#' @param point_size Size of scatter points (default: 1)
#' @param alpha Transparency of points (default: 0.6)
#' @param palette Color palette for labels (default: "Set3")
#' @param title Plot title
#'
#' @return ggplot object
plot_cca_umap <- function(umap_coords, 
                          labels = NULL, 
                          point_size = 1, 
                          alpha = 0.6, 
                          palette = "Set3", 
                          title = "UMAP of CCA Components") {
  
  # Create data frame
  plot_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2]
  )
  
  # Add labels if provided
  if (!is.null(labels)) {
    plot_data$Label <- as.factor(labels)
    
    # Create plot with labels
    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Label)) +
      geom_point(size = point_size, alpha = alpha) +
      scale_color_brewer(palette = palette) +
      theme_minimal() +
      theme(legend.position = "right") +
      labs(title = title, x = "UMAP 1", y = "UMAP 2")
    
    # Adjust colors if too many categories
    n_categories <- length(unique(labels))
    if (n_categories > 12) {
      p <- p + scale_color_manual(values = viridis(n_categories))
    }
  } else {
    # Create plot without labels
    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2)) +
      geom_point(size = point_size, alpha = alpha, color = "blue") +
      theme_minimal() +
      labs(title = title, x = "UMAP 1", y = "UMAP 2")
  }
  
  return(p)
}


#' Complete pipeline for CCA to UMAP analysis with optional diagnostics
#'
#' @param cca_matrix CCA output matrix
#' @param labels Cell labels for visualization
#' @param n_components Number of CCA components to use
#' @param plot_variance Whether to plot explained variance of CCs
#' @param ... Parameters to pass to UMAP
#'
#' @return List containing UMAP coordinates and plots
analyze_cca_umap <- function(cca_matrix, 
                             labels = NULL, 
                             n_components = 4, 
                             plot_variance = TRUE,
                             ...) {
  
  plots <- list()
  
  # Optionally plot variance explained by CCs
  if (plot_variance && ncol(cca_matrix) >= n_components) {
    # Calculate variance of each CC
    cc_vars <- apply(cca_matrix, 2, var)
    n_show <- min(20, length(cc_vars))
    
    # Variance by CC plot
    var_data <- data.frame(
      CC = 1:n_show,
      Variance = cc_vars[1:n_show]
    )
    
    p1 <- ggplot(var_data[1:min(10, nrow(var_data)), ], 
                 aes(x = CC, y = Variance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(title = "Variance by CC", 
           x = "Canonical Component", 
           y = "Variance")
    
    # Cumulative variance plot
    cum_var <- cumsum(cc_vars) / sum(cc_vars) * 100
    cum_data <- data.frame(
      CC = 1:n_show,
      CumVar = cum_var[1:n_show]
    )
    
    p2 <- ggplot(cum_data, aes(x = CC, y = CumVar)) +
      geom_line() +
      geom_point() +
      geom_vline(xintercept = n_components, 
                 color = "red", 
                 linetype = "dashed") +
      theme_minimal() +
      labs(title = "Cumulative Variance Explained", 
           x = "Number of CCs", 
           y = "Cumulative Variance (%)") +
      annotate("text", 
               x = n_components + 0.5, 
               y = cum_var[n_components] - 5, 
               label = paste0("Selected: ", n_components, " CCs"),
               color = "red",
               hjust = 0)
    
    # Combine plots
    variance_plot <- plot_grid(p1, p2, ncol = 2)
    plots$variance_plot <- variance_plot
    print(variance_plot)
  }
  
  # Compute UMAP
  umap_coords <- cca_to_umap(cca_matrix, 
                              n_components = n_components, 
                              ...)
  
  # Plot UMAP
  umap_plot <- plot_cca_umap(umap_coords, labels = labels)
  plots$umap_plot <- umap_plot
  print(umap_plot)
  
  # Return results
  return(list(
    umap_coords = umap_coords,
    plots = plots
  ))
}


#' Helper function to create density plots for UMAP visualization
#'
#' @param umap_coords UMAP coordinates
#' @param labels Optional labels for faceting
#' @param bins Number of bins for 2D density (default: 50)
#'
#' @return ggplot object
plot_umap_density <- function(umap_coords, 
                              labels = NULL, 
                              bins = 50) {
  
  plot_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2]
  )
  
  if (!is.null(labels)) {
    plot_data$Label <- as.factor(labels)
    
    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2)) +
      geom_density_2d_filled(bins = bins) +
      facet_wrap(~ Label) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = "UMAP Density by Group", 
           x = "UMAP 1", 
           y = "UMAP 2")
  } else {
    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2)) +
      geom_density_2d_filled(bins = bins) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = "UMAP Density", 
           x = "UMAP 1", 
           y = "UMAP 2")
  }
  
  return(p)
}


#' Save UMAP results to files
#'
#' @param umap_coords UMAP coordinates
#' @param output_prefix Prefix for output files
#' @param labels Optional labels to save
#'
#' @return NULL (saves files)
save_umap_results <- function(umap_coords, 
                              output_prefix, 
                              labels = NULL) {
  
  # Save coordinates
  coords_df <- data.frame(
    Cell_ID = 1:nrow(umap_coords),
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2]
  )
  
  if (!is.null(labels)) {
    coords_df$Label <- labels
  }
  
  write.csv(coords_df, 
            file = paste0(output_prefix, "_umap_coords.csv"), 
            row.names = FALSE)
  
  message(paste("UMAP coordinates saved to", 
                paste0(output_prefix, "_umap_coords.csv")))
}


# Example usage:
# # Assuming you have your CCA output matrix
# # cca_output <- your_cca_matrix  # shape: (n_cells, n_canonical_components)
# # cell_labels <- your_cell_labels  # optional
# 
# # Basic usage
# # umap_coords <- cca_to_umap(cca_output, n_components = 4)
# 
# # With visualization and analysis
# # results <- analyze_cca_umap(cca_output, 
# #                             labels = cell_labels, 
# #                             n_components = 4)
# 
# # Access coordinates
# # umap_coords <- results$umap_coords
# 
# # Custom UMAP parameters
# # umap_coords <- cca_to_umap(cca_output, 
# #                            n_components = 4, 
# #                            n_neighbors = 50, 
# #                            min_dist = 0.1, 
# #                            metric = "cosine")
# 
# # Create density plot
# # density_plot <- plot_umap_density(umap_coords, labels = cell_labels)
# # print(density_plot)
# 
# # Save results
# # save_umap_results(umap_coords, "my_analysis", labels = cell_labels)