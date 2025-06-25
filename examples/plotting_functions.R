#' Plot Superimposed Cell Scores with Different Color Gradients
#'
#' This function creates a superimposed plot of cell scores for multiple cell types
#' using different color gradients. Each cell type is plotted as a separate layer
#' with its own color scale to enable visual comparison of score distributions.
#'
#' @param object A CoPro object containing cell scores
#' @param sigmaValueChoice Sigma value to use for plotting
#' @param ccIndex Canonical component index to plot (default: 1)
#' @param cell_types Vector of cell types to include in the plot (default: all available)
#' @param color_scales Named list of color scales for each cell type. 
#'   Options: "viridis", "red_pink", "blue_purple", "green_yellow", "plasma", "inferno", "magma"
#' @param alpha_values Named list of alpha values for each cell type (default: 0.7)
#' @param point_size Size of points (default: 0.8)
#' @param title Plot title
#' @param quantile_range Quantile range to map colors (default: c(0.01, 0.99))
#'
#' @return A ggplot object with superimposed cell type scores
#' @export
#' @importFrom ggplot2 ggplot geom_point theme_minimal coord_fixed labs
#' @importFrom viridis viridis plasma inferno magma
plotSuperimposedCellScores <- function(object, 
                                       sigmaValueChoice,
                                       ccIndex = 1,
                                       cell_types = NULL,
                                       color_scales = NULL,
                                       alpha_values = NULL,
                                       point_size = 0.8,
                                       title = "Superimposed Cell Type Scores",
                                       quantile_range = c(0.01, 0.99)) {
  
  # Load required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("viridis package is required for this function")
  }
  
  # Get cell types to plot
  if (is.null(cell_types)) {
    if (length(object@cellTypesOfInterest) != 0) {
      cell_types <- object@cellTypesOfInterest
    } else {
      cell_types <- unique(object@cellTypesSub)
    }
  }
  
  if (length(cell_types) < 2) {
    stop("At least two cell types are required for superimposed plotting")
  }
  
  # Set default alpha values
  if (is.null(alpha_values)) {
    alpha_values <- setNames(rep(0.7, length(cell_types)), cell_types)
  }
  
  # Set default color scales
  if (is.null(color_scales)) {
    default_scales <- c("viridis", "red_pink", "blue_purple", "plasma", "inferno", "magma")
    color_scales <- setNames(default_scales[1:length(cell_types)], cell_types)
  }
  
  # Get cell scores data for each cell type
  sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")
  
  # Check if multi-slide or single slide
  if (!isMultiSlide(object)) {
    plot_data_list <- .getSingleSlideDataForSuperimposed(object, sigma_name, ccIndex, cell_types)
  } else {
    plot_data_list <- .getMultiSlideDataForSuperimposed(object, sigma_name, ccIndex, cell_types)
  }
  
  # Define color palettes
  color_palettes <- list(
    viridis = viridis::viridis(100),
    red_pink = colorRampPalette(c("#FFE0E6", "#FF1493", "#8B0000"))(100),
    blue_purple = colorRampPalette(c("#E6F3FF", "#4169E1", "#191970"))(100),
    green_yellow = colorRampPalette(c("#F0FFF0", "#32CD32", "#006400"))(100),
    plasma = viridis::plasma(100),
    inferno = viridis::inferno(100),
    magma = viridis::magma(100)
  )
  
  # Create the base plot
  p <- ggplot2::ggplot()
  
  # Add layers for each cell type
  for (ct in cell_types) {
    ct_data <- plot_data_list[[ct]]
    if (is.null(ct_data) || nrow(ct_data) == 0) {
      warning(paste("No data available for cell type:", ct))
      next
    }
    
    # Normalize scores to 0-1 range based on quantiles
    score_range <- quantile(ct_data$cellScores, probs = quantile_range, na.rm = TRUE)
    ct_data$cellScores_norm <- pmax(0, pmin(1, (ct_data$cellScores - score_range[1]) / 
                                              (score_range[2] - score_range[1])))
    
    # Get color palette for this cell type
    palette <- color_palettes[[color_scales[ct]]]
    if (is.null(palette)) {
      palette <- viridis::viridis(100)
      warning(paste("Unknown color scale for", ct, "using viridis"))
    }
    
    # Map scores to colors
    color_indices <- pmax(1, pmin(100, round(ct_data$cellScores_norm * 99) + 1))
    ct_data$point_colors <- palette[color_indices]
    
    # Add points for this cell type
    p <- p + ggplot2::geom_point(
      data = ct_data,
      ggplot2::aes(x = x, y = y),
      color = ct_data$point_colors,
      size = point_size,
      alpha = alpha_values[ct]
    )
  }
  
  # Add theme and labels
  p <- p + 
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      subtitle = paste("Cell types:", paste(cell_types, collapse = ", ")),
      x = "X coordinate",
      y = "Y coordinate",
      caption = paste("Color scales:", paste(names(color_scales), color_scales, sep = " = ", collapse = "; "))
    )
  
  return(p)
}

#' Create Superimposed Plot with Custom Legends
#'
#' This version creates separate legend plots for each cell type to show the color mapping
#'
#' @param object A CoPro object containing cell scores
#' @param sigmaValueChoice Sigma value to use for plotting
#' @param ccIndex Canonical component index to plot (default: 1)
#' @param cell_types Vector of cell types to include in the plot
#' @param color_scales Named list of color scales for each cell type
#' @param alpha_values Named list of alpha values for each cell type (default: 0.7)
#' @param point_size Size of points (default: 0.8)
#' @param create_legends Whether to create separate legend plots (default: TRUE)
#'
#' @return A list containing the main plot and legend plots (if requested)
#' @export
plotSuperimposedWithLegends <- function(object, 
                                         sigmaValueChoice,
                                         ccIndex = 1,
                                         cell_types = NULL,
                                         color_scales = NULL,
                                         alpha_values = NULL,
                                         point_size = 0.8,
                                         create_legends = TRUE) {
  
  # Get the main plot
  main_plot <- plotSuperimposedCellScores(
    object = object,
    sigmaValueChoice = sigmaValueChoice,
    ccIndex = ccIndex,
    cell_types = cell_types,
    color_scales = color_scales,
    alpha_values = alpha_values,
    point_size = point_size
  )
  
  result <- list(main_plot = main_plot)
  
  if (create_legends) {
    # Get cell types to plot
    if (is.null(cell_types)) {
      if (length(object@cellTypesOfInterest) != 0) {
        cell_types <- object@cellTypesOfInterest
      } else {
        cell_types <- unique(object@cellTypesSub)
      }
    }
    
    # Set default color scales
    if (is.null(color_scales)) {
      default_scales <- c("viridis", "red_pink", "blue_purple", "plasma", "inferno", "magma")
      color_scales <- setNames(default_scales[1:length(cell_types)], cell_types)
    }
    
    # Get data for creating legends
    sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")
    if (!isMultiSlide(object)) {
      plot_data_list <- .getSingleSlideDataForSuperimposed(object, sigma_name, ccIndex, cell_types)
    } else {
      plot_data_list <- .getMultiSlideDataForSuperimposed(object, sigma_name, ccIndex, cell_types)
    }
    
    # Create legend plots
    legend_plots <- list()
    
    for (ct in cell_types) {
      ct_data <- plot_data_list[[ct]]
      if (is.null(ct_data) || nrow(ct_data) == 0) {
        next
      }
      
      # Create a simple legend plot
      legend_data <- data.frame(
        x = 1:100,
        y = rep(1, 100),
        score = seq(min(ct_data$cellScores, na.rm = TRUE), 
                   max(ct_data$cellScores, na.rm = TRUE), 
                   length.out = 100)
      )
      
      if (color_scales[ct] == "viridis") {
        legend_plot <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y, fill = score)) +
          ggplot2::geom_tile() +
          viridis::scale_fill_viridis(name = paste(ct, "Score")) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste(ct, "Legend"))
      } else if (color_scales[ct] == "red_pink") {
        legend_plot <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y, fill = score)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient(low = "#FFE0E6", high = "#8B0000", name = paste(ct, "Score")) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste(ct, "Legend"))
      } else {
        # Default to viridis for other options
        legend_plot <- ggplot2::ggplot(legend_data, ggplot2::aes(x = x, y = y, fill = score)) +
          ggplot2::geom_tile() +
          viridis::scale_fill_viridis(option = color_scales[ct], name = paste(ct, "Score")) +
          ggplot2::theme_void() +
          ggplot2::labs(title = paste(ct, "Legend"))
      }
      
      legend_plots[[ct]] <- legend_plot
    }
    
    result$legend_plots <- legend_plots
  }
  
  return(result)
}

#' Helper function to get single slide data for superimposed plotting
#' @noRd
.getSingleSlideDataForSuperimposed <- function(object, sigma_name, ccIndex, cell_types) {
  plot_data_list <- setNames(vector("list", length(cell_types)), cell_types)
  
  for (ct in cell_types) {
    # Get location data for this cell type
    loc_data <- object@locationDataSub[object@cellTypesSub == ct, , drop = FALSE]
    
    if (nrow(loc_data) == 0) {
      warning(paste("No cells found for cell type", ct))
      next
    }
    
    # Check if cell scores exist for this cell type
    if (!ct %in% names(object@cellScores[[sigma_name]])) {
      warning(paste("Cell scores not found for cell type", ct))
      next
    }
    
    cell_scores_matrix <- object@cellScores[[sigma_name]][[ct]]
    common_cells <- intersect(rownames(loc_data), rownames(cell_scores_matrix))
    
    if (length(common_cells) == 0) {
      warning(paste("No matching cells between location data and cell scores for cell type", ct))
      next
    }
    
    # Subset to common cells
    loc_data <- loc_data[common_cells, , drop = FALSE]
    loc_data$cellScores <- cell_scores_matrix[common_cells, ccIndex]
    loc_data$cellType <- ct
    
    plot_data_list[[ct]] <- loc_data
  }
  
  return(plot_data_list)
}

#' Helper function to get multi-slide data for superimposed plotting
#' @noRd
.getMultiSlideDataForSuperimposed <- function(object, sigma_name, ccIndex, cell_types) {
  all_slides <- unique(getSlideID(object))
  plot_data_list <- setNames(vector("list", length(cell_types)), cell_types)
  
  for (ct in cell_types) {
    ct_data_list <- list()
    
    # Check if cell scores exist for this cell type (now using aggregated structure)
    if (!ct %in% names(object@cellScores[[sigma_name]])) {
      warning(paste("Cell scores not found for cell type", ct))
      next
    }
    
    cell_scores_matrix <- object@cellScores[[sigma_name]][[ct]]
    
    if (is.null(cell_scores_matrix) || nrow(cell_scores_matrix) == 0) {
      warning(paste("Empty cell scores matrix for cell type", ct))
      next
    }
    
    if (ccIndex > ncol(cell_scores_matrix)) {
      stop(paste("ccIndex", ccIndex, "exceeds number of canonical components"))
    }
    
    for (q in all_slides) {
      # Get location data for this slide and cell type
      slide_mask <- object@cellTypesSub == ct & object@metaDataSub$slideID == q
      loc_data <- object@locationDataSub[slide_mask, , drop = FALSE]
      
      if (nrow(loc_data) == 0) {
        next
      }
      
      # Get cell IDs for this slide and cell type
      slide_cell_ids <- rownames(object@metaDataSub)[slide_mask]
      
      # Find common cells between location data and cell scores
      common_cells <- intersect(rownames(loc_data), slide_cell_ids)
      common_cells <- intersect(common_cells, rownames(cell_scores_matrix))
      
      if (length(common_cells) == 0) {
        next
      }
      
      # Subset to common cells
      loc_data <- loc_data[common_cells, , drop = FALSE]
      loc_data$cellScores <- cell_scores_matrix[common_cells, ccIndex]
      loc_data$cellType <- ct
      loc_data$slideID <- q
      
      ct_data_list <- append(ct_data_list, list(loc_data))
    }
    
    if (length(ct_data_list) > 0) {
      plot_data_list[[ct]] <- do.call(rbind, ct_data_list)
    }
  }
  
  return(plot_data_list)
}

#' Create Side-by-Side Comparison Plot
#'
#' This function creates separate subplots for each cell type for easy comparison
#'
#' @param object A CoPro object containing cell scores
#' @param sigmaValueChoice Sigma value to use for plotting
#' @param ccIndex Canonical component index to plot (default: 1)
#' @param cell_types Vector of cell types to include in the plot
#' @param use_viridis Whether to use viridis color scale for all plots (default: TRUE)
#' @param point_size Size of points (default: 0.8)
#' @param ncol Number of columns for faceting (default: 2)
#'
#' @return A ggplot object with faceted cell type scores
#' @export
#' @importFrom ggplot2 ggplot geom_point facet_wrap theme_minimal coord_fixed
plotCellTypeComparison <- function(object,
                                   sigmaValueChoice,
                                   ccIndex = 1,
                                   cell_types = NULL,
                                   use_viridis = TRUE,
                                   point_size = 0.8,
                                   ncol = 2) {
  
  # Get all cell type data
  if (is.null(cell_types)) {
    if (length(object@cellTypesOfInterest) != 0) {
      cell_types <- object@cellTypesOfInterest
    } else {
      cell_types <- unique(object@cellTypesSub)
    }
  }
  
  # Get data for all cell types
  sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")
  
  if (!isMultiSlide(object)) {
    plot_data_list <- .getSingleSlideDataForSuperimposed(object, sigma_name, ccIndex, cell_types)
  } else {
    plot_data_list <- .getMultiSlideDataForSuperimposed(object, sigma_name, ccIndex, cell_types)
  }
  
  # Combine all data
  combined_data <- do.call(rbind, plot_data_list)
  
  if (is.null(combined_data) || nrow(combined_data) == 0) {
    stop("No valid data found for any cell type")
  }
  
  # Create faceted plot
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = x, y = y, color = cellScores)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(~ cellType, ncol = ncol) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Cell Type Score Comparison",
      x = "X coordinate",
      y = "Y coordinate",
      color = "Cell Score"
    )
  
  if (use_viridis) {
    p <- p + viridis::scale_color_viridis(option = "D")
  } else {
    p <- p + ggplot2::scale_color_gradient(low = "lightblue", high = "darkblue")
  }
  
  return(p)
} 