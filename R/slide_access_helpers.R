#' Helper functions for consistent slide access
#' 
#' These functions provide standardized ways to access slide information
#' and filter data by slides, ensuring consistency across the codebase.

#' Get slide indices for filtering data
#' @param object A CoPro object  
#' @param slide Slide ID to filter by
#' @return Logical vector for filtering (TRUE for cells in specified slide)
#' @noRd
.getSlideIndices <- function(object, slide) {
  if (!isMultiSlide(object)) {
    # For single slide objects, return all TRUE
    return(rep(TRUE, length(object@cellTypesSub)))
  }
  
  slide_ids <- getSlideID(object)
  if (!slide %in% getSlideList(object)) {
    stop(paste("Slide", slide, "not found. Available slides:", 
               paste(getSlideList(object), collapse = ", ")))
  }
  
  return(slide_ids == slide)
}

#' Get cells for specific slide and cell type combination
#' @param object A CoPro object
#' @param slide Slide ID (ignored for CoProSingle)
#' @param cellType Cell type to filter by
#' @return Character vector of cell IDs
#' @noRd
.getSlideCellTypeIndices <- function(object, slide = NULL, cellType = NULL) {
  # Start with all cells
  indices <- rep(TRUE, length(object@cellTypesSub))
  
  # Filter by slide if specified and object is multi-slide
  if (!is.null(slide) && isMultiSlide(object)) {
    indices <- indices & .getSlideIndices(object, slide)
  }
  
  # Filter by cell type if specified
  if (!is.null(cellType)) {
    indices <- indices & (object@cellTypesSub == cellType)
  }
  
  return(indices)
}

#' Count cells for specific slide and cell type combination
#' @param object A CoPro object
#' @param slide Slide ID (ignored for CoProSingle)  
#' @param cellType Cell type to count
#' @return Integer count of cells
#' @noRd
.countSlideCellType <- function(object, slide = NULL, cellType = NULL) {
  indices <- .getSlideCellTypeIndices(object, slide, cellType)
  return(sum(indices))
}

#' Get cell IDs for specific slide and cell type combination
#' @param object A CoPro object
#' @param slide Slide ID (ignored for CoProSingle)
#' @param cellType Cell type to filter by  
#' @return Character vector of cell IDs
#' @noRd
.getSlideCellTypeIDs <- function(object, slide = NULL, cellType = NULL) {
  indices <- .getSlideCellTypeIndices(object, slide, cellType)
  return(rownames(object@metaDataSub)[indices])
}

#' Validate slide parameter for functions
#' @param object A CoPro object
#' @param slide Slide ID to validate
#' @param require_slide Whether slide is required for this object type
#' @return TRUE if valid, stops with error if invalid
#' @noRd
.validateSlideParameter <- function(object, slide, require_slide = NULL) {
  if (isMultiSlide(object)) {
    if (is.null(require_slide)) require_slide <- TRUE
    
    if (require_slide && is.null(slide)) {
      stop("slide parameter is required for CoProMulti objects")
    }
    
    if (!is.null(slide) && !slide %in% getSlideList(object)) {
      stop(paste("Slide", slide, "not found. Available slides:", 
                 paste(getSlideList(object), collapse = ", ")))
    }
  } else {
    # For single slide objects, slide parameter should be NULL or ignored
    if (!is.null(slide)) {
      warning("slide parameter ignored for CoProSingle objects")
    }
  }
  
  return(TRUE)
}

#' Check for slide consistency between slideList and metadata
#' @param object A CoProMulti object
#' @return TRUE if consistent, stops with error if inconsistent
#' @noRd
.validateSlideConsistency <- function(object) {
  if (!isMultiSlide(object)) {
    return(TRUE)
  }
  
  # Check metaDataSub first (for subsetted objects), then fall back to metaData (for new objects)
  metadata_to_check <- if (nrow(object@metaDataSub) > 0) {
    object@metaDataSub
  } else {
    object@metaData
  }
  
  # Check that slideList contains all unique values from metadata slideID
  if ("slideID" %in% colnames(metadata_to_check)) {
    expected_slides <- sort(unique(metadata_to_check$slideID))
    actual_slides <- sort(getSlideList(object))
    
    if (!identical(expected_slides, actual_slides)) {
      stop("slideList does not match unique values in metadata$slideID")
    }
  } else {
    # If no slideID in metadata, slideList should be empty
    if (length(getSlideList(object)) > 0) {
      stop("slideList is not empty but metadata lacks slideID column")
    }
  }
  
  return(TRUE)
} 