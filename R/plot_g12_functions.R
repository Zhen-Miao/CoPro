# Global variables for ggplot2 and dplyr to avoid linting warnings
utils::globalVariables(c("r_um", "g12_obs", "pair_label", "g12_lower", "g12_upper", 
                         "slideID", "cellType1", "cellType2", "nCells1", "nCells2"))

#' Plot g_12(r) pair correlation functions for colocalization analysis
#'
#' This function creates plots of the inhomogeneous cross pair-correlation function
#' g_12(r) values across different radii for all cell type pairs. It can generate
#' both individual plots for each pair and a combined plot with all pairs overlaid.
#'
#' @param object A CoProSingle or CoProMulti object with location and cell type data.
#' @param r_um_range Numeric vector of length 2 specifying the distance range in 
#'   microns for analysis (default: c(10, 60)).
#' @param pixel_size_um Numeric value specifying microns per coordinate unit.
#'   This converts your input coordinates to microns for analysis. For example:
#'   - If 1 coordinate unit = 50 microns, set pixel_size_um = 50
#'   - If 1 coordinate unit = 0.325 microns, set pixel_size_um = 0.325
#'   - If coordinates are already in microns, set pixel_size_um = 1 (default)
#' @param cell_diam_um Numeric value specifying typical cell diameter in microns 
#'   for intensity smoothing (default: 10).
#' @param nsim Integer specifying number of random labelling simulations for 
#'   null distribution (default: 199).
#' @param r_step_um Numeric value specifying step size for radius grid in microns 
#'   (default: 2).
#' @param min_points_per_type Integer specifying minimum number of points required 
#'   per cell type (default: 20).
#' @param edge_correction Character specifying edge correction method 
#'   (default: "translation").
#' @param plot_type Character specifying plot type: "combined" for all pairs in 
#'   one plot, "individual" for separate plots per pair, or "both" (default: "combined").
#' @param include_confidence Logical; whether to include confidence bands from 
#'   simulations (default: TRUE).
#' @param confidence_level Numeric; confidence level for bands (default: 0.95).
#' @param colors Character vector of colors for different cell type pairs. If NULL,
#'   uses default color palette.
#' @param line_width Numeric; width of the g_12(r) lines (default: 1.2).
#' @param alpha_bands Numeric; transparency for confidence bands (default: 0.3).
#' @param verbose Logical; whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#'   - `plot`: ggplot object(s) with the g_12(r) plots
#'   - `data`: data.frame with g_12(r) values, confidence intervals, and metadata
#'   - `summary`: summary statistics for each cell type pair
#'
#' @details The function computes the inhomogeneous cross pair-correlation function
#' g_12(r) for each cell type pair and plots the observed values along with 
#' confidence bands derived from random labelling simulations.
#' 
#' Values of g_12(r):
#' - g_12(r) = 1: Random spatial relationship at distance r
#' - g_12(r) > 1: Attraction/colocalization at distance r  
#' - g_12(r) < 1: Repulsion/segregation at distance r
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' g12_plots <- plotG12Functions(object)
#' 
#' # Custom parameters with individual plots
#' g12_plots <- plotG12Functions(object, 
#'                              r_um_range = c(5, 40),
#'                              plot_type = "individual",
#'                              include_confidence = TRUE)
#' 
#' # Access the plot and data
#' print(g12_plots$plot)
#' head(g12_plots$data)
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline facet_wrap 
#'   labs theme_minimal theme element_text scale_color_manual scale_fill_manual
#' @importFrom dplyr bind_rows mutate
#' @importFrom grDevices rainbow
#' @importFrom stats quantile
plotG12Functions <- function(object,
                            r_um_range = c(10, 60),
                            pixel_size_um = 1,
                            cell_diam_um = 10,
                            nsim = 199,
                            r_step_um = 2,
                            min_points_per_type = 20,
                            edge_correction = "translation",
                            plot_type = c("combined", "individual", "both"),
                            include_confidence = TRUE,
                            confidence_level = 0.95,
                            colors = NULL,
                            line_width = 1.2,
                            alpha_bands = 0.3,
                            verbose = TRUE) {
  
  # Input validation
  if (!(is(object, "CoProSingle") || is(object, "CoProMulti"))) {
    stop("object must be a CoProSingle or CoProMulti object")
  }
  
  plot_type <- match.arg(plot_type)
  
  if (length(r_um_range) != 2 || r_um_range[1] >= r_um_range[2]) {
    stop("r_um_range must be a vector of length 2 with r_um_range[1] < r_um_range[2]")
  }
  
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be between 0 and 1")
  }
  
  # Check required data
  if (nrow(object@locationDataSub) == 0) {
    stop("No location data found. Ensure the object has been properly initialized.")
  }
  
  if (length(object@cellTypesSub) == 0) {
    stop("No cell type data found. Ensure the object has been properly initialized.")
  }
  
  # Determine processing approach based on object type
  if (is(object, "CoProSingle")) {
    plot_data <- .plotG12Single(object, r_um_range, pixel_size_um, cell_diam_um,
                               nsim, r_step_um, min_points_per_type, edge_correction,
                               include_confidence, confidence_level, verbose)
  } else {
    plot_data <- .plotG12Multi(object, r_um_range, pixel_size_um, cell_diam_um,
                              nsim, r_step_um, min_points_per_type, edge_correction,
                              include_confidence, confidence_level, verbose)
  }
  
  # Generate plots
  plots <- .generateG12Plots(plot_data, plot_type, colors, line_width, 
                            alpha_bands, include_confidence, r_um_range)
  
  # Generate summary statistics
  summary_stats <- .generateG12Summary(plot_data, r_um_range)
  
  return(list(
    plot = plots,
    data = plot_data,
    summary = summary_stats
  ))
}

#' Process g_12(r) data for single-slide objects
#' @noRd
.plotG12Single <- function(object, r_um_range, pixel_size_um, cell_diam_um,
                          nsim, r_step_um, min_points_per_type, edge_correction,
                          include_confidence, confidence_level, verbose) {
  
  # Get cell types and generate pairs
  cts <- unique(object@cellTypesSub)
  if (length(cts) < 2) {
    stop("Need at least 2 cell types for cross-type g_12(r) analysis")
  }
  
  pair_matrix <- combn(cts, 2)
  all_pairs <- data.frame(
    cellType1 = pair_matrix[1, ],
    cellType2 = pair_matrix[2, ],
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("Computing g_12(r) functions for", nrow(all_pairs), "cell type pairs\n")
  }
  
  # Determine window range
  window_range <- list(
    xrange = range(object@locationDataSub$x, na.rm = TRUE),
    yrange = range(object@locationDataSub$y, na.rm = TRUE)
  )
  
  # Process each pair
  all_results <- vector("list", length = nrow(all_pairs))
  
  for (i in seq_len(nrow(all_pairs))) {
    ct1 <- all_pairs$cellType1[i]
    ct2 <- all_pairs$cellType2[i]
    
    if (verbose) {
      cat("Processing pair", i, "of", nrow(all_pairs), ":", ct1, "vs", ct2, "\n")
    }
    
    # Extract cell locations for each type
    cells1_idx <- object@cellTypesSub == ct1
    cells2_idx <- object@cellTypesSub == ct2
    
    if (sum(cells1_idx) < min_points_per_type || sum(cells2_idx) < min_points_per_type) {
      if (verbose) {
        warning(paste("Skipping pair", ct1, "vs", ct2, 
                     ": insufficient points (", sum(cells1_idx), ",", sum(cells2_idx), ")"))
      }
      next
    }
    
    A <- object@locationDataSub[cells1_idx, c("x", "y"), drop = FALSE]
    B <- object@locationDataSub[cells2_idx, c("x", "y"), drop = FALSE]
    
    # Compute detailed g_12(r) data
    g12_result <- tryCatch({
      .computeDetailedG12(A = A, B = B, window_range = window_range,
                         r_um_range = r_um_range, pixel_size_um = pixel_size_um,
                         cell_diam_um = cell_diam_um, nsim = nsim,
                         r_step_um = r_step_um, include_confidence = include_confidence,
                         confidence_level = confidence_level)
    }, error = function(e) {
      if (verbose) {
        warning(paste("Error computing g_12(r) for", ct1, "vs", ct2, ":", e$message))
      }
      return(NULL)
    })
    
    if (!is.null(g12_result)) {
      g12_result$cellType1 <- ct1
      g12_result$cellType2 <- ct2
      g12_result$pair_label <- paste(ct1, "vs", ct2)
      all_results[[i]] <- g12_result
    }
  }
  
  # Combine results
  valid_results <- all_results[!sapply(all_results, is.null)]
  if (length(valid_results) == 0) {
    stop("No valid g_12(r) results computed")
  }
  
  combined_data <- dplyr::bind_rows(valid_results)
  return(combined_data)
}

#' Process g_12(r) data for multi-slide objects
#' @noRd
.plotG12Multi <- function(object, r_um_range, pixel_size_um, cell_diam_um,
                         nsim, r_step_um, min_points_per_type, edge_correction,
                         include_confidence, confidence_level, verbose) {
  
  slides <- getSlideList(object)
  cts <- unique(object@cellTypesSub)
  
  if (length(cts) < 2) {
    stop("Need at least 2 cell types for cross-type g_12(r) analysis")
  }
  
  pair_matrix <- combn(cts, 2)
  all_pairs <- data.frame(
    cellType1 = pair_matrix[1, ],
    cellType2 = pair_matrix[2, ],
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("Computing g_12(r) functions for", length(slides), "slides and", 
        nrow(all_pairs), "cell type pairs\n")
  }
  
  # Process each slide
  all_slide_results <- vector("list", length = length(slides))
  
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
    
    # Process each pair on this slide
    slide_results <- vector("list", length = nrow(all_pairs))
    
    for (i in seq_len(nrow(all_pairs))) {
      ct1 <- all_pairs$cellType1[i]
      ct2 <- all_pairs$cellType2[i]
      
      # Extract cell locations for each type on this slide
      cells1_idx <- slide_celltypes == ct1
      cells2_idx <- slide_celltypes == ct2
      
      if (sum(cells1_idx) < min_points_per_type || sum(cells2_idx) < min_points_per_type) {
        next
      }
      
      A <- slide_location[cells1_idx, c("x", "y"), drop = FALSE]
      B <- slide_location[cells2_idx, c("x", "y"), drop = FALSE]
      
      # Compute detailed g_12(r) data
      g12_result <- tryCatch({
        .computeDetailedG12(A = A, B = B, window_range = window_range,
                           r_um_range = r_um_range, pixel_size_um = pixel_size_um,
                           cell_diam_um = cell_diam_um, nsim = nsim,
                           r_step_um = r_step_um, include_confidence = include_confidence,
                           confidence_level = confidence_level)
      }, error = function(e) {
        if (verbose) {
          warning(paste("Error computing g_12(r) for", ct1, "vs", ct2, 
                       "on slide", slide_id, ":", e$message))
        }
        return(NULL)
      })
      
      if (!is.null(g12_result)) {
        g12_result$slideID <- slide_id
        g12_result$cellType1 <- ct1
        g12_result$cellType2 <- ct2
        g12_result$pair_label <- paste(ct1, "vs", ct2)
        slide_results[[i]] <- g12_result
      }
    }
    
    valid_slide_results <- slide_results[!sapply(slide_results, is.null)]
    if (length(valid_slide_results) > 0) {
      all_slide_results[[slide_idx]] <- dplyr::bind_rows(valid_slide_results)
    }
  }
  
  # Combine results from all slides
  valid_slide_results <- all_slide_results[!sapply(all_slide_results, is.null)]
  if (length(valid_slide_results) == 0) {
    stop("No valid g_12(r) results computed across all slides")
  }
  
  combined_data <- dplyr::bind_rows(valid_slide_results)
  return(combined_data)
}

#' Compute detailed g_12(r) data with confidence intervals
#' @noRd
.computeDetailedG12 <- function(A, B, window_range, r_um_range, pixel_size_um,
                               cell_diam_um, nsim, r_step_um, include_confidence,
                               confidence_level) {
  
  # Use the same logic as coloc_score but return detailed data
  nA <- nrow(A)
  nB <- nrow(B)
  
  # Convert inputs to microns (standardized units)
  s <- 1 / pixel_size_um  # Conversion factor: coordinate_units -> microns
  Au <- data.frame(x = A$x * s, y = A$y * s)
  Bu <- data.frame(x = B$x * s, y = B$y * s)
  
  # Window in microns (same units as converted coordinates)
  win <- spatstat.geom::owin(xrange = window_range$xrange * s, 
                            yrange = window_range$yrange * s)
  
  # Build multitype pattern
  all_coords <- rbind(
    data.frame(x = Au$x, y = Au$y, type = "A"),
    data.frame(x = Bu$x, y = Bu$y, type = "B")
  )
  
  X <- spatstat.geom::ppp(x = all_coords$x, y = all_coords$y, window = win,
                         marks = factor(all_coords$type, levels = c("A", "B")))
  
  # r grid in microns (after coordinate conversion, both coordinates and r are in microns)
  r_vec_um <- seq(0, r_um_range[2], by = r_step_um)
  r_vec <- r_vec_um  # No scaling needed - coordinates are now in microns, r should be too
  
  # Intensity smoothing bandwidth in microns (same units as coordinates after conversion)
  sigma_intensity <- (3 * cell_diam_um)
  
  # Observed g12_inhom
  g_obs <- .compute_g12_inhom(X, "A", "B", r_vec, sigma_intensity)
  
  # Prepare result data frame
  result_data <- data.frame(
    r_um = r_vec_um,
    g12_obs = g_obs,
    stringsAsFactors = FALSE
  )
  
  # Add confidence intervals if requested
  if (include_confidence && nsim > 0) {
    # Simulate random labelling nsim times
    g_mat <- matrix(NA_real_, nrow = length(r_vec), ncol = nsim)

    for (k in seq_len(nsim)) {
      X_sim <- spatstat.random::rlabel(X)
      spatstat.geom::marks(X_sim) <- factor(spatstat.geom::marks(X_sim), levels = c("A", "B"))
      g_mat[, k] <- .compute_g12_inhom(X_sim, "A", "B", r_vec, sigma_intensity)
    }
    
    # Calculate confidence intervals
    alpha <- 1 - confidence_level
    lower_quantile <- alpha / 2
    upper_quantile <- 1 - alpha / 2
    
    g_lower <- apply(g_mat, 1, quantile, probs = lower_quantile, na.rm = TRUE)
    g_upper <- apply(g_mat, 1, quantile, probs = upper_quantile, na.rm = TRUE)
    g_mean <- rowMeans(g_mat, na.rm = TRUE)
    
    result_data$g12_mean_null <- g_mean
    result_data$g12_lower <- g_lower
    result_data$g12_upper <- g_upper
  }
  
  # Add metadata
  result_data$nCells1 <- nA
  result_data$nCells2 <- nB
  
  return(result_data)
}

#' Generate g_12(r) plots
#' @noRd
.generateG12Plots <- function(plot_data, plot_type, colors, line_width, 
                             alpha_bands, include_confidence, r_um_range) {
  
  # Set up colors
  n_pairs <- length(unique(plot_data$pair_label))
  if (is.null(colors)) {
    # Use a colorblind-friendly palette
    colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7", "#999999")[1:n_pairs]
    if (n_pairs > 8) {
      colors <- rainbow(n_pairs)
    }
  } else if (length(colors) < n_pairs) {
    warning("Not enough colors provided, using default palette")
    colors <- rainbow(n_pairs)
  }
  
  names(colors) <- unique(plot_data$pair_label)
  
  # Base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = r_um, y = g12_obs, color = pair_label))
  
  # Add confidence bands if available
  if (include_confidence && "g12_lower" %in% colnames(plot_data)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = g12_lower, ymax = g12_upper, 
                                              fill = pair_label),
                                 alpha = alpha_bands, color = NA)
  }
  
  # Add observed g_12(r) lines
  p <- p + ggplot2::geom_line(linewidth = line_width)
  
  # Add reference line at g_12(r) = 1 (random expectation)
  p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dashed", 
                              color = "black", alpha = 0.7)
  
  # Customize appearance
  p <- p + ggplot2::scale_color_manual(values = colors, name = "Cell Type Pairs") +
    ggplot2::labs(
      x = "Distance r (um)",
      y = expression(g[12](r)),
      title = "Cross-type Pair Correlation Functions",
      subtitle = paste0("Distance range: ", r_um_range[1], "-", r_um_range[2], " um")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10)
    )
  
  if (include_confidence && "g12_lower" %in% colnames(plot_data)) {
    p <- p + ggplot2::scale_fill_manual(values = colors, name = "Cell Type Pairs")
  }
  
  # Handle different plot types
  if (plot_type == "individual") {
    # Check if we have slides
    if ("slideID" %in% colnames(plot_data)) {
      p <- p + ggplot2::facet_wrap(~ pair_label + slideID, scales = "free_y")
    } else {
      p <- p + ggplot2::facet_wrap(~ pair_label, scales = "free_y")
    }
    p <- p + ggplot2::theme(legend.position = "none")
  } else if (plot_type == "both") {
    # Return both combined and individual plots
    p_individual <- p + ggplot2::facet_wrap(~ pair_label, scales = "free_y") +
      ggplot2::theme(legend.position = "none")
    
    return(list(combined = p, individual = p_individual))
  }
  
  return(p)
}

#' Generate summary statistics for g_12(r) data
#' @noRd
.generateG12Summary <- function(plot_data, r_um_range) {
  
  # Filter to analysis range
  analysis_data <- plot_data[plot_data$r_um >= r_um_range[1] & 
                            plot_data$r_um <= r_um_range[2], ]
  
  # Calculate summary statistics by pair
  if ("slideID" %in% colnames(plot_data)) {
    # Multi-slide case
    grouped_data <- dplyr::group_by(analysis_data, slideID, cellType1, cellType2, pair_label)
    summary_stats <- dplyr::summarise(
      grouped_data,
      mean_g12 = mean(g12_obs, na.rm = TRUE),
      max_g12 = max(g12_obs, na.rm = TRUE),
      min_g12 = min(g12_obs, na.rm = TRUE),
      r_at_max = r_um[which.max(g12_obs)],
      nCells1 = dplyr::first(nCells1),
      nCells2 = dplyr::first(nCells2),
      .groups = "drop"
    )
  } else {
    # Single-slide case
    grouped_data <- dplyr::group_by(analysis_data, cellType1, cellType2, pair_label)
    summary_stats <- dplyr::summarise(
      grouped_data,
      mean_g12 = mean(g12_obs, na.rm = TRUE),
      max_g12 = max(g12_obs, na.rm = TRUE),
      min_g12 = min(g12_obs, na.rm = TRUE),
      r_at_max = r_um[which.max(g12_obs)],
      nCells1 = dplyr::first(nCells1),
      nCells2 = dplyr::first(nCells2),
      .groups = "drop"
    )
  }
  
  return(summary_stats)
}
