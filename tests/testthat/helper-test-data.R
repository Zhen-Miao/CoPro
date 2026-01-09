#' Helper functions for generating test data
#' 
#' These functions create synthetic spatial transcriptomics data for testing
#' the CoPro package functionality.

#' Generate simple test data for CoPro
#' 
#' Creates a minimal dataset suitable for testing basic CoPro functionality
#' 
#' @param n_cells Total number of cells
#' @param n_genes Number of genes
#' @param n_cell_types Number of cell types
#' @param seed Random seed for reproducibility
#' @return List with normalizedData, locationData, metaData, cellTypes
generate_test_data_single <- function(n_cells = 100, 
                                      n_genes = 50, 
                                      n_cell_types = 2,
                                      seed = 42) {
  set.seed(seed)
  
  # Generate cell type labels
  cell_type_names <- paste0("CellType", LETTERS[1:n_cell_types])
  cell_types <- sample(cell_type_names, n_cells, replace = TRUE)
  
  # Generate cell IDs
  cell_ids <- paste0("cell_", seq_len(n_cells))
  
  # Generate gene names
  gene_names <- paste0("Gene", seq_len(n_genes))
  
  # Generate normalized expression data (log-transformed counts)
  normalizedData <- matrix(
    rnorm(n_cells * n_genes, mean = 2, sd = 1),
    nrow = n_cells,
    ncol = n_genes
  )
  normalizedData[normalizedData < 0] <- 0
  rownames(normalizedData) <- cell_ids
  colnames(normalizedData) <- gene_names
  
  # Generate spatial coordinates in a 10x10 grid
  locationData <- data.frame(
    x = runif(n_cells, 0, 10),
    y = runif(n_cells, 0, 10)
  )
  rownames(locationData) <- cell_ids
  
  # Generate metadata
  metaData <- data.frame(
    cell_id = cell_ids,
    cell_type = cell_types,
    total_counts = rowSums(normalizedData),
    row.names = cell_ids
  )
  
  return(list(
    normalizedData = normalizedData,
    locationData = locationData,
    metaData = metaData,
    cellTypes = cell_types
  ))
}

#' Generate test data for multi-slide CoPro
#' 
#' Creates a dataset with multiple slides for testing CoProMulti functionality
#' 
#' @param n_cells_per_slide Number of cells per slide
#' @param n_slides Number of slides
#' @param n_genes Number of genes
#' @param n_cell_types Number of cell types
#' @param seed Random seed for reproducibility
#' @return List with normalizedData, locationData, metaData, cellTypes, slideID
generate_test_data_multi <- function(n_cells_per_slide = 50,
                                     n_slides = 2,
                                     n_genes = 30,
                                     n_cell_types = 2,
                                     seed = 42) {
  set.seed(seed)
  
  n_cells <- n_cells_per_slide * n_slides
  
  # Generate slide IDs
  slide_names <- paste0("Slide", seq_len(n_slides))
  slideID <- rep(slide_names, each = n_cells_per_slide)
  
  # Generate cell type labels
  cell_type_names <- paste0("CellType", LETTERS[1:n_cell_types])
  cell_types <- sample(cell_type_names, n_cells, replace = TRUE)
  
  # Generate unique cell IDs (with slide prefix to ensure uniqueness)
  cell_ids <- paste0(rep(slide_names, each = n_cells_per_slide), "_cell_", 
                     rep(seq_len(n_cells_per_slide), n_slides))
  
  # Generate gene names
  gene_names <- paste0("Gene", seq_len(n_genes))
  
  # Generate normalized expression data
  normalizedData <- matrix(
    rnorm(n_cells * n_genes, mean = 2, sd = 1),
    nrow = n_cells,
    ncol = n_genes
  )
  normalizedData[normalizedData < 0] <- 0
  rownames(normalizedData) <- cell_ids
  colnames(normalizedData) <- gene_names
  
  # Generate spatial coordinates (each slide has its own coordinate system)
  locationData <- data.frame(
    x = runif(n_cells, 0, 10),
    y = runif(n_cells, 0, 10)
  )
  rownames(locationData) <- cell_ids
  
  # Generate metadata
  metaData <- data.frame(
    cell_id = cell_ids,
    cell_type = cell_types,
    total_counts = rowSums(normalizedData),
    row.names = cell_ids
  )
  
  return(list(
    normalizedData = normalizedData,
    locationData = locationData,
    metaData = metaData,
    cellTypes = cell_types,
    slideID = slideID
  ))
}

#' Generate test data with spatial patterns
#' 
#' Creates data where cell scores have spatial autocorrelation,
#' useful for testing whether CoPro can detect spatial patterns.
#' 
#' @param n_cells Number of cells
#' @param n_genes Number of genes
#' @param bandwidth Spatial smoothing bandwidth
#' @param seed Random seed
#' @return List with test data including spatially correlated expression
generate_test_data_with_pattern <- function(n_cells = 100,
                                            n_genes = 50,
                                            bandwidth = 1.5,
                                            seed = 42) {
  set.seed(seed)
  
  # Generate two cell types
  cell_types <- rep(c("TypeA", "TypeB"), each = n_cells / 2)
  cell_ids <- paste0("cell_", seq_len(n_cells))
  gene_names <- paste0("Gene", seq_len(n_genes))
  
  # Generate spatial coordinates
  x <- runif(n_cells, 0, 10)
  y <- runif(n_cells, 0, 10)
  
  # Create a spatially smooth score for the first 10 genes
  dist_matrix <- as.matrix(dist(cbind(x, y)))
  kernel_matrix <- exp(-0.5 * (dist_matrix / bandwidth)^2)
  kernel_matrix <- kernel_matrix / rowSums(kernel_matrix)
  
  # Generate base expression data
  normalizedData <- matrix(
    rnorm(n_cells * n_genes, mean = 2, sd = 0.5),
    nrow = n_cells,
    ncol = n_genes
  )
  
  # Add spatial pattern to first few genes
  initial_pattern <- rnorm(n_cells, mean = 0, sd = 1)
  spatial_pattern <- as.vector(kernel_matrix %*% initial_pattern)
  spatial_pattern <- (spatial_pattern - mean(spatial_pattern)) / sd(spatial_pattern)
  
  # Add pattern to first 10 genes
  for (g in 1:min(10, n_genes)) {
    normalizedData[, g] <- normalizedData[, g] + spatial_pattern * 0.5
  }
  
  normalizedData[normalizedData < 0] <- 0
  rownames(normalizedData) <- cell_ids
  colnames(normalizedData) <- gene_names
  
  locationData <- data.frame(x = x, y = y)
  rownames(locationData) <- cell_ids
  
  metaData <- data.frame(
    cell_id = cell_ids,
    cell_type = cell_types,
    spatial_score = spatial_pattern,
    row.names = cell_ids
  )
  
  return(list(
    normalizedData = normalizedData,
    locationData = locationData,
    metaData = metaData,
    cellTypes = cell_types,
    true_pattern = spatial_pattern
  ))
}

#' Create a minimal CoPro object for testing
#' 
#' Convenience function that creates and returns a CoProSingle object
#' ready for analysis
#' 
#' @param n_cells Number of cells
#' @param n_genes Number of genes
#' @param n_cell_types Number of cell types
#' @param seed Random seed
#' @return A CoProSingle object
create_test_copro_single <- function(n_cells = 100,
                                     n_genes = 50,
                                     n_cell_types = 2,
                                     seed = 42) {
  test_data <- generate_test_data_single(n_cells, n_genes, n_cell_types, seed)
  
  obj <- newCoProSingle(
    normalizedData = test_data$normalizedData,
    locationData = test_data$locationData,
    metaData = test_data$metaData,
    cellTypes = test_data$cellTypes
  )
  
  return(obj)
}

#' Create a minimal CoProMulti object for testing
#' 
#' @param n_cells_per_slide Number of cells per slide
#' @param n_slides Number of slides
#' @param n_genes Number of genes
#' @param n_cell_types Number of cell types
#' @param seed Random seed
#' @return A CoProMulti object
create_test_copro_multi <- function(n_cells_per_slide = 50,
                                    n_slides = 2,
                                    n_genes = 30,
                                    n_cell_types = 2,
                                    seed = 42) {
  test_data <- generate_test_data_multi(n_cells_per_slide, n_slides, 
                                        n_genes, n_cell_types, seed)
  
  obj <- newCoProMulti(
    normalizedData = test_data$normalizedData,
    locationData = test_data$locationData,
    metaData = test_data$metaData,
    cellTypes = test_data$cellTypes,
    slideID = test_data$slideID
  )
  
  return(obj)
}

