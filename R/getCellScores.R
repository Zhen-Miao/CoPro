#' Get cell scores from CoPro object
#'
#' This function provides a safe way to access cell scores from CoPro objects,
#' using the flat structure for efficient access.
#' It works with both single-slide (CoProSingle) and multi-slide (CoProMulti) objects.
#'
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @param sigma Numeric sigma value for which to retrieve cell scores
#' @param cellType Character string specifying the cell type
#' @param slide Character string specifying slide ID (only for CoProMulti objects, ignored for CoProSingle)
#' @param ccIndex Optional integer specifying which canonical component to extract (1-based). 
#'               If NULL (default), returns the entire matrix
#' @param cells Optional character vector of specific cell IDs to extract. 
#'             If NULL (default), returns all cells
#' @param verbose Logical indicating whether to print informative messages (default: TRUE)
#'
#' @return A matrix of cell scores with cells as rows and canonical components as columns.
#'         If ccIndex is specified, returns a vector of scores for that component.
#'         If cells is specified, returns only those specific cells.
#'
#' @details
#' The function uses a flat structure for efficient access:
#' - Cell scores are stored with names like "cellScores|sigma0.1|TypeA"
#' 
#' For multi-slide objects, cell scores are aggregated across slides and specific
#' slide filtering is done using the metaDataSub slot.
#'
#' @examples
#' \dontrun{
#' # Get all cell scores for a specific sigma and cell type
#' scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA")
#' 
#' # Get scores for a specific canonical component
#' cc1_scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA", ccIndex = 1)
#' 
#' # Get scores for specific cells
#' specific_scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA", 
#'                                  cells = c("cell_1", "cell_2"))
#' 
#' # For multi-slide object, specify slide
#' slide_scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA", 
#'                               slide = "slide1")
#' }
#'
#' @export
#' @rdname getCellScores
setGeneric("getCellScores", function(object, sigma, cellType, slide = NULL, 
                                     ccIndex = NULL, cells = NULL, verbose = TRUE) {
  standardGeneric("getCellScores")
})

#' @rdname getCellScores
#' @export
setMethod("getCellScores", "CoProSingle", function(object, sigma, cellType, slide = NULL, 
                                                    ccIndex = NULL, cells = NULL, verbose = TRUE) {
  # Input validation
  if (length(object@cellScores) == 0) {
    stop("cellScores slot is empty. Please run computeGeneAndCellScores() first.")
  }
  
  if (length(sigma) != 1 || !is.numeric(sigma)) {
    stop("sigma must be a single numeric value")
  }
  
  if (length(cellType) != 1 || !is.character(cellType)) {
    stop("cellType must be a single character string")
  }
  
  if (!is.null(slide) && verbose) {
    message("Note: slide parameter is ignored for CoProSingle objects")
  }
  
  # Check sigma value exists
  if (!sigma %in% object@sigmaValues) {
    stop("Sigma value ", sigma, " not found in object@sigmaValues: ", 
         paste(object@sigmaValues, collapse = ", "))
  }
  
  # Always use flat structure access
  flat_name <- .createCellScoresName(sigma, cellType, slide = NULL)
  if (!flat_name %in% names(object@cellScores)) {
    stop("Cell scores not found for sigma=", sigma, " and cellType='", cellType, "'")
  }
  scores_matrix <- object@cellScores[[flat_name]]
  
  if (is.null(scores_matrix)) {
    stop("Cell scores matrix is NULL for sigma=", sigma, " and cellType='", cellType, "'")
  }
  
  if (verbose) {
    message("Retrieved cell scores for sigma=", sigma, ", cellType='", cellType, 
            "' (", nrow(scores_matrix), " cells x ", ncol(scores_matrix), " components)")
  }
  
  # Apply ccIndex filtering if specified
  if (!is.null(ccIndex)) {
    if (length(ccIndex) != 1 || !is.numeric(ccIndex) || ccIndex < 1) {
      stop("ccIndex must be a single positive integer")
    }
    
    if (ccIndex > ncol(scores_matrix)) {
      stop("ccIndex ", ccIndex, " exceeds number of canonical components (", 
           ncol(scores_matrix), ")")
    }
    
    scores_matrix <- scores_matrix[, ccIndex, drop = TRUE]
    
    if (verbose) {
      message("Extracted canonical component ", ccIndex)
    }
  }
  
  # Apply cell filtering if specified
  if (!is.null(cells)) {
    if (!is.character(cells)) {
      stop("cells must be a character vector of cell IDs")
    }
    
    available_cells <- if (is.null(ccIndex)) rownames(scores_matrix) else names(scores_matrix)
    missing_cells <- setdiff(cells, available_cells)
    
    if (length(missing_cells) > 0) {
      stop("Cells not found: ", paste(missing_cells, collapse = ", "))
    }
    
    scores_matrix <- if (is.null(ccIndex)) {
      scores_matrix[cells, , drop = FALSE]
    } else {
      scores_matrix[cells]
    }
    
    if (verbose) {
      message("Filtered to ", length(cells), " specific cells")
    }
  }
  
  return(scores_matrix)
})

#' @rdname getCellScores
#' @export
setMethod("getCellScores", "CoProMulti", function(object, sigma, cellType, slide = NULL, 
                                                   ccIndex = NULL, cells = NULL, verbose = TRUE) {
  # Input validation
  if (length(object@cellScores) == 0) {
    stop("cellScores slot is empty. Please run computeGeneAndCellScores() first.")
  }
  
  if (length(sigma) != 1 || !is.numeric(sigma)) {
    stop("sigma must be a single numeric value")
  }
  
  if (length(cellType) != 1 || !is.character(cellType)) {
    stop("cellType must be a single character string")
  }
  
  # Check sigma value exists
  if (!sigma %in% object@sigmaValues) {
    stop("Sigma value ", sigma, " not found in object@sigmaValues: ", 
         paste(object@sigmaValues, collapse = ", "))
  }
  
  # Always use flat structure access
  flat_name <- .createCellScoresName(sigma, cellType, slide = NULL)
  if (!flat_name %in% names(object@cellScores)) {
    stop("Cell scores not found for sigma=", sigma, " and cellType='", cellType, "'")
  }
  scores_matrix <- object@cellScores[[flat_name]]
  
  if (is.null(scores_matrix)) {
    stop("Cell scores matrix is NULL for sigma=", sigma, " and cellType='", cellType, "'")
  }
  
  # Apply slide filtering if specified
  if (!is.null(slide)) {
    if (length(slide) != 1 || !is.character(slide)) {
      stop("slide must be a single character string")
    }
    
    if (!slide %in% getSlideList(object)) {
      stop("Slide '", slide, "' not found. Available slides: ", 
           paste(getSlideList(object), collapse = ", "))
    }
    
    # Filter to cells from the specified slide
    slide_cells <- .getSlideCellTypeIDs(object, slide = slide, cellType = cellType)
    
    if (length(slide_cells) == 0) {
      stop("No cells found for slide '", slide, "' and cellType '", cellType, "'")
    }
    
    # Check if slide cells exist in scores matrix
    available_cells <- intersect(slide_cells, rownames(scores_matrix))
    if (length(available_cells) == 0) {
      stop("No matching cells found in scores matrix for slide '", slide, 
           "' and cellType '", cellType, "'")
    }
    
    scores_matrix <- scores_matrix[available_cells, , drop = FALSE]
    
    if (verbose) {
      message("Filtered to slide '", slide, "' (", nrow(scores_matrix), " cells)")
    }
  }
  
  if (verbose && is.null(slide)) {
    message("Retrieved cell scores for sigma=", sigma, ", cellType='", cellType, 
            "' (", nrow(scores_matrix), " cells x ", ncol(scores_matrix), " components, aggregated across slides)")
  } else if (verbose) {
    message("Retrieved cell scores for sigma=", sigma, ", cellType='", cellType, 
            "', slide='", slide, "' (", nrow(scores_matrix), " cells x ", ncol(scores_matrix), " components)")
  }
  
  # Apply ccIndex filtering if specified
  if (!is.null(ccIndex)) {
    if (length(ccIndex) != 1 || !is.numeric(ccIndex) || ccIndex < 1) {
      stop("ccIndex must be a single positive integer")
    }
    
    if (ccIndex > ncol(scores_matrix)) {
      stop("ccIndex ", ccIndex, " exceeds number of canonical components (", 
           ncol(scores_matrix), ")")
    }
    
    scores_matrix <- scores_matrix[, ccIndex, drop = TRUE]
    
    if (verbose) {
      message("Extracted canonical component ", ccIndex)
    }
  }
  
  # Apply cell filtering if specified
  if (!is.null(cells)) {
    if (!is.character(cells)) {
      stop("cells must be a character vector of cell IDs")
    }
    
    available_cells <- if (is.null(ccIndex)) rownames(scores_matrix) else names(scores_matrix)
    missing_cells <- setdiff(cells, available_cells)
    
    if (length(missing_cells) > 0) {
      stop("Cells not found: ", paste(missing_cells, collapse = ", "))
    }
    
    scores_matrix <- if (is.null(ccIndex)) {
      scores_matrix[cells, , drop = FALSE]
    } else {
      scores_matrix[cells]
    }
    
    if (verbose) {
      message("Filtered to ", length(cells), " specific cells")
    }
  }
  
  return(scores_matrix)
}) 