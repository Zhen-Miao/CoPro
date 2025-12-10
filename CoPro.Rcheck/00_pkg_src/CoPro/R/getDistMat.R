#' Get Distance Matrix
#' 
#' A unified accessor function for distance matrices that uses flat structure
#' for efficient access.
#' 
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @param cellType1 First cell type name
#' @param cellType2 Second cell type name  
#' @param slide Slide ID (required for CoProMulti objects, ignored for CoProSingle)
#' @param returnTranspose If TRUE, forces return of transpose when accessing symmetric matrices
#' @param verbose Whether to print detailed error messages
#'
#' @return Distance matrix as a numeric matrix
#' @export
#' @rdname getDistMat
#' @aliases getDistMat,CoProSingle-method
#' @aliases getDistMat,CoProMulti-method
setGeneric("getDistMat", 
           function(object, cellType1, cellType2, slide = NULL, 
                   returnTranspose = FALSE, verbose = TRUE) 
             standardGeneric("getDistMat"))

#' @rdname getDistMat
#' @export
setMethod("getDistMat", "CoProSingle",
          function(object, cellType1, cellType2, slide = NULL, 
                   returnTranspose = FALSE, verbose = TRUE) {
            .getDistMatSingle(object, cellType1, cellType2, 
                             returnTranspose, verbose)
          })

#' @rdname getDistMat  
#' @export
setMethod("getDistMat", "CoProMulti",
          function(object, cellType1, cellType2, slide = NULL, 
                   returnTranspose = FALSE, verbose = TRUE) {
            .getDistMatMulti(object, cellType1, cellType2, 
                            slide, returnTranspose, verbose)
          })

#' Internal function for single slide distance matrix access
#' @noRd
.getDistMatSingle <- function(object, cellType1, cellType2, 
                             returnTranspose = FALSE, verbose = TRUE) {
  
  # Input validation
  .validateDistMatInputs(object, cellType1, cellType2, verbose)
  
  # Check if distance matrices exist
  if (length(object@distances) == 0) {
    stop("No distance matrices found. Run computeDistance() first.")
  }
  
  # Always use flat structure access
  return(.getDistMatFlat(object@distances, cellType1, cellType2, 
                        slide = NULL, returnTranspose, verbose))
}

#' Internal function for multi slide distance matrix access
#' @noRd
.getDistMatMulti <- function(object, cellType1, cellType2, slide,
                            returnTranspose = FALSE, verbose = TRUE) {
  
  # Input validation
  .validateDistMatInputs(object, cellType1, cellType2, verbose)
  
  if (is.null(slide)) {
    stop("slide parameter is required for CoProMulti objects")
  }
  
  if (!slide %in% getSlideList(object)) {
    stop(paste("Slide", slide, "not found. Available slides:", 
               paste(getSlideList(object), collapse = ", ")))
  }
  
  # Check if distance matrices exist
  if (length(object@distances) == 0) {
    stop("No distance matrices found. Run computeDistance() first.")
  }
  
  # Always use flat structure access
  return(.getDistMatFlat(object@distances, cellType1, cellType2, 
                        slide = slide, returnTranspose, verbose))
}

#' Input validation for distance matrix access
#' @noRd
.validateDistMatInputs <- function(object, cellType1, cellType2, verbose) {
  
  # Check if cell types are valid
  if (!is.character(cellType1) || !is.character(cellType2)) {
    stop("cellType1 and cellType2 must be character strings")
  }
  
  if (length(cellType1) != 1 || length(cellType2) != 1) {
    stop("cellType1 and cellType2 must be single character strings")
  }
  
  # Check if cell types exist in object
  if (length(object@cellTypesOfInterest) > 0) {
    if (!cellType1 %in% object@cellTypesOfInterest) {
      stop(paste("cellType1", cellType1, "not in cellTypesOfInterest:", 
                 paste(object@cellTypesOfInterest, collapse = ", ")))
    }
    if (!cellType2 %in% object@cellTypesOfInterest) {
      stop(paste("cellType2", cellType2, "not in cellTypesOfInterest:", 
                 paste(object@cellTypesOfInterest, collapse = ", ")))
    }
  }
}

#' Throw standardized distance matrix not found error
#' @noRd
.throwDistMatNotFoundError <- function(cellType1, cellType2, slide = NULL) {
  if (is.null(slide)) {
    stop(paste("Distance matrix not found for", cellType1, "->", cellType2, 
               ". Check if computeDistance() was run properly."))
  } else {
    stop(paste("Distance matrix not found for", cellType1, "->", cellType2, 
               "on slide", slide, ". Check if computeDistance() was run properly."))
  }
}

#' Access distance matrix from flat structure
#' @noRd
.getDistMatFlat <- function(flat_distances, cellType1, cellType2, slide = NULL,
                           returnTranspose = FALSE, verbose = TRUE) {
  
  # Create the expected flat name
  flat_name <- .createDistMatrixName(cellType1, cellType2, slide)
  
  # Try direct access
  if (flat_name %in% names(flat_distances)) {
    dist_matrix <- flat_distances[[flat_name]]
    if (returnTranspose) {
      return(t(dist_matrix))
    } else {
      return(dist_matrix)
    }
  }
  
  # Try symmetric access (swap cellType1 and cellType2)
  symmetric_name <- .createDistMatrixName(cellType2, cellType1, slide)
  
  if (symmetric_name %in% names(flat_distances)) {
    dist_matrix <- flat_distances[[symmetric_name]]
    if (verbose) {
      message(paste("Using transpose of distance matrix for", cellType2, "->", cellType1))
    }
    if (returnTranspose) {
      return(dist_matrix)  # Don't transpose again
    } else {
      return(t(dist_matrix))
    }
  }
  
  # If we get here, distance matrix not found
  .throwDistMatNotFoundError(cellType1, cellType2, slide)
} 