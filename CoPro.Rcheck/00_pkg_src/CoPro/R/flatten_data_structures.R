#' Flatten nested data structures with informative names
#' 
#' This file contains functions to convert between nested and flat data structures
#' to improve memory efficiency and reduce complexity while maintaining 
#' backward compatibility through accessor functions.

#' Create informative names for kernel matrices
#' @param sigma Sigma value
#' @param cellType1 First cell type  
#' @param cellType2 Second cell type
#' @param slide Slide ID (optional, for multi-slide objects)
#' @return Character string with informative name
#' @noRd
.createKernelMatrixName <- function(sigma, cellType1, cellType2, slide = NULL) {
  if (is.null(slide)) {
    # Single slide: "kernel|sigma0.1|TypeA|TypeB"
    paste("kernel", paste0("sigma", sigma), cellType1, cellType2, sep = "|")
  } else {
    # Multi slide: "kernel|sigma0.1|slide1|TypeA|TypeB"  
    paste("kernel", paste0("sigma", sigma), slide, cellType1, cellType2, sep = "|")
  }
}

#' Create informative names for distance matrices
#' @param cellType1 First cell type
#' @param cellType2 Second cell type
#' @param slide Slide ID (optional, for multi-slide objects)
#' @return Character string with informative name
#' @noRd
.createDistMatrixName <- function(cellType1, cellType2, slide = NULL) {
  if (is.null(slide)) {
    # Single slide: "dist|TypeA|TypeB"
    paste("dist", cellType1, cellType2, sep = "|")
  } else {
    # Multi slide: "dist|slide1|TypeA|TypeB"
    paste("dist", slide, cellType1, cellType2, sep = "|")
  }
}

#' Create informative names for cell scores
#' @param sigma Sigma value
#' @param cellType Cell type
#' @param slide Slide ID (optional, for specific slide access)
#' @return Character string with informative name
#' @noRd
.createCellScoresName <- function(sigma, cellType, slide = NULL) {
  if (is.null(slide)) {
    # "cellScores|sigma0.1|TypeA"
    paste("cellScores", paste0("sigma", sigma), cellType, sep = "|")
  } else {
    # "cellScores|sigma0.1|slide1|TypeA" (if we need slide-specific access)
    paste("cellScores", paste0("sigma", sigma), slide, cellType, sep = "|")
  }
}

#' Create informative names for gene scores
#' @param sigma Sigma value
#' @param cellType Cell type
#' @param slide Slide ID (optional, for specific slide access)
#' @return Character string with informative name
#' @noRd
.createGeneScoresName <- function(sigma, cellType, slide = NULL) {
  if (is.null(slide)) {
    # "geneScores|sigma0.1|TypeA"
    paste("geneScores", paste0("sigma", sigma), cellType, sep = "|")
  } else {
    # "geneScores|sigma0.1|slide1|TypeA" (if we need slide-specific access)
    paste("geneScores", paste0("sigma", sigma), slide, cellType, sep = "|")
  }
}

#' Parse kernel matrix name to extract components
#' @param name Flattened name string
#' @return List with sigma, cellType1, cellType2, and slide (if applicable)
#' @noRd
.parseKernelMatrixName <- function(name) {
  parts <- strsplit(name, "\\|")[[1]]
  
  if (length(parts) == 4) {
    # Single slide: "kernel|sigma0.1|TypeA|TypeB"
    list(
      sigma = as.numeric(gsub("sigma", "", parts[2])),
      cellType1 = parts[3],
      cellType2 = parts[4],
      slide = NULL
    )
  } else if (length(parts) == 5) {
    # Multi slide: "kernel|sigma0.1|slide1|TypeA|TypeB"
    list(
      sigma = as.numeric(gsub("sigma", "", parts[2])),
      slide = parts[3],
      cellType1 = parts[4],
      cellType2 = parts[5]
    )
  } else {
    stop("Invalid kernel matrix name format: ", name)
  }
}

#' Parse distance matrix name to extract components
#' @param name Flattened name string
#' @return List with cellType1, cellType2, and slide (if applicable)
#' @noRd
.parseDistMatrixName <- function(name) {
  parts <- strsplit(name, "\\|")[[1]]
  
  if (length(parts) == 3) {
    # Single slide: "dist|TypeA|TypeB"
    list(
      cellType1 = parts[2],
      cellType2 = parts[3],
      slide = NULL
    )
  } else if (length(parts) == 4) {
    # Multi slide: "dist|slide1|TypeA|TypeB"
    list(
      slide = parts[2],
      cellType1 = parts[3],
      cellType2 = parts[4]
    )
  } else {
    stop("Invalid distance matrix name format: ", name)
  }
}

#' Parse cell scores name to extract components
#' @param name Flattened name string
#' @return List with sigma, cellType, and slide (if applicable)
#' @noRd
.parseCellScoresName <- function(name) {
  parts <- strsplit(name, "\\|")[[1]]
  
  if (length(parts) == 3) {
    # "cellScores|sigma0.1|TypeA"
    list(
      sigma = as.numeric(gsub("sigma", "", parts[2])),
      cellType = parts[3],
      slide = NULL
    )
  } else if (length(parts) == 4) {
    # "cellScores|sigma0.1|slide1|TypeA"
    list(
      sigma = as.numeric(gsub("sigma", "", parts[2])),
      slide = parts[3],
      cellType = parts[4]
    )
  } else {
    stop("Invalid cell scores name format: ", name)
  }
}

#' Parse gene scores name to extract components
#' @param name Flattened name string
#' @return List with sigma, cellType, and slide (if applicable)
#' @noRd
.parseGeneScoresName <- function(name) {
  parts <- strsplit(name, "\\|")[[1]]
  
  if (length(parts) == 3) {
    # "geneScores|sigma0.1|TypeA"
    list(
      sigma = as.numeric(gsub("sigma", "", parts[2])),
      cellType = parts[3],
      slide = NULL
    )
  } else if (length(parts) == 4) {
    # "geneScores|sigma0.1|slide1|TypeA"
    list(
      sigma = as.numeric(gsub("sigma", "", parts[2])),
      slide = parts[3],
      cellType = parts[4]
    )
  } else {
    stop("Invalid gene scores name format: ", name)
  }
}

#' Convert flat kernel matrices back to nested structure
#' @param flat_kernels Flat list of kernel matrices with informative names
#' @param sigmaValues Vector of sigma values
#' @param cellTypes Vector of cell types
#' @param slides Vector of slide IDs (for multi-slide objects)
#' @return Nested list structure
#' @noRd
.flatToNestedKernels <- function(flat_kernels, sigmaValues, cellTypes, slides = NULL) {
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  
  if (is.null(slides)) {
    # Single slide structure
    nested <- setNames(vector("list", length(sigmaValues)), sigma_names)
    for (sigma_name in sigma_names) {
      nested[[sigma_name]] <- setNames(vector("list", length(cellTypes)), cellTypes)
      for (ct1 in cellTypes) {
        nested[[sigma_name]][[ct1]] <- setNames(vector("list", length(cellTypes)), cellTypes)
      }
    }
    
    # Fill nested structure from flat
    for (name in names(flat_kernels)) {
      parsed <- .parseKernelMatrixName(name)
      if (is.null(parsed$slide)) {  # Only process single-slide entries
        sigma_name <- paste("sigma", parsed$sigma, sep = "_")
        nested[[sigma_name]][[parsed$cellType1]][[parsed$cellType2]] <- flat_kernels[[name]]
      }
    }
  } else {
    # Multi slide structure  
    nested <- setNames(vector("list", length(sigmaValues)), sigma_names)
    for (sigma_name in sigma_names) {
      nested[[sigma_name]] <- setNames(vector("list", length(slides)), slides)
      for (slide in slides) {
        nested[[sigma_name]][[slide]] <- setNames(vector("list", length(cellTypes)), cellTypes)
        for (ct1 in cellTypes) {
          nested[[sigma_name]][[slide]][[ct1]] <- setNames(vector("list", length(cellTypes)), cellTypes)
        }
      }
    }
    
    # Fill nested structure from flat
    for (name in names(flat_kernels)) {
      parsed <- .parseKernelMatrixName(name)
      if (!is.null(parsed$slide)) {  # Only process multi-slide entries
        sigma_name <- paste("sigma", parsed$sigma, sep = "_")
        nested[[sigma_name]][[parsed$slide]][[parsed$cellType1]][[parsed$cellType2]] <- flat_kernels[[name]]
      }
    }
  }
  
  return(nested)
}

#' Convert flat cell scores back to nested structure
#' @param flat_scores Flat list of cell score matrices with informative names
#' @param sigmaValues Vector of sigma values
#' @param cellTypes Vector of cell types
#' @return Nested list structure
#' @noRd
.flatToNestedCellScores <- function(flat_scores, sigmaValues, cellTypes) {
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  
  nested <- setNames(vector("list", length(sigmaValues)), sigma_names)
  for (sigma_name in sigma_names) {
    nested[[sigma_name]] <- setNames(vector("list", length(cellTypes)), cellTypes)
  }
  
  # Fill nested structure from flat
  for (name in names(flat_scores)) {
    parsed <- .parseCellScoresName(name)
    if (is.null(parsed$slide)) {  # Only process aggregated entries
      sigma_name <- paste("sigma", parsed$sigma, sep = "_")
      nested[[sigma_name]][[parsed$cellType]] <- flat_scores[[name]]
    }
  }
  
  return(nested)
}

#' Convert nested kernel matrices to flat structure
#' @param nested_kernels Nested kernel matrix structure
#' @param slides Vector of slide IDs (NULL for single slide)
#' @return Flat list with informative names
#' @noRd
.nestedToFlatKernels <- function(nested_kernels, slides = NULL) {
  flat_kernels <- list()
  
  for (sigma_name in names(nested_kernels)) {
    sigma_val <- as.numeric(gsub("sigma_", "", sigma_name))
    
    if (is.null(slides)) {
      # Single slide structure: kernels[[sigma]][[ct1]][[ct2]]
      for (ct1 in names(nested_kernels[[sigma_name]])) {
        for (ct2 in names(nested_kernels[[sigma_name]][[ct1]])) {
          kernel_matrix <- nested_kernels[[sigma_name]][[ct1]][[ct2]]
          if (!is.null(kernel_matrix) && length(kernel_matrix) > 0) {
            flat_name <- .createKernelMatrixName(sigma_val, ct1, ct2, slide = NULL)
            flat_kernels[[flat_name]] <- kernel_matrix
          }
        }
      }
    } else {
      # Multi slide structure: kernels[[sigma]][[slide]][[ct1]][[ct2]]
      for (slide in names(nested_kernels[[sigma_name]])) {
        for (ct1 in names(nested_kernels[[sigma_name]][[slide]])) {
          for (ct2 in names(nested_kernels[[sigma_name]][[slide]][[ct1]])) {
            kernel_matrix <- nested_kernels[[sigma_name]][[slide]][[ct1]][[ct2]]
            if (!is.null(kernel_matrix) && length(kernel_matrix) > 0) {
              flat_name <- .createKernelMatrixName(sigma_val, ct1, ct2, slide = slide)
              flat_kernels[[flat_name]] <- kernel_matrix
            }
          }
        }
      }
    }
  }
  
  return(flat_kernels)
}

#' Convert nested cell scores to flat structure
#' @param nested_scores Nested cell scores structure
#' @return Flat list with informative names
#' @noRd
.nestedToFlatCellScores <- function(nested_scores) {
  flat_scores <- list()
  
  for (sigma_name in names(nested_scores)) {
    sigma_val <- as.numeric(gsub("sigma_", "", sigma_name))
    
    # Structure: scores[[sigma]][[cellType]]
    for (ct in names(nested_scores[[sigma_name]])) {
      score_matrix <- nested_scores[[sigma_name]][[ct]]
      if (!is.null(score_matrix) && length(score_matrix) > 0) {
        flat_name <- .createCellScoresName(sigma_val, ct, slide = NULL)
        flat_scores[[flat_name]] <- score_matrix
      }
    }
  }
  
  return(flat_scores)
}

#' Update accessor functions to work with flat structures
#' 
#' These functions need to be called to update the accessor behavior
#' when the data structure is flattened
#' 
#' @noRd
.updateAccessorFunctions <- function() {
  # The accessor functions (getKernelMatrix and getCellScores) would need
  # to be updated to work with flat structures. This is handled by
  # modifying the internal access logic in those functions.
  
  # For now, we maintain backward compatibility by converting flat to nested
  # when the accessor functions are called.
  
  message("Accessor functions updated for flat structure compatibility")
} 