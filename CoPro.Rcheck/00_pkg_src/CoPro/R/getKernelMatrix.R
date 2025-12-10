#' Get Kernel Matrix
#' 
#' A unified accessor function for kernel matrices that handles the complex
#' nested structure and provides symmetric access when needed.
#' 
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @param sigma Sigma value for kernel selection
#' @param cellType1 First cell type name
#' @param cellType2 Second cell type name  
#' @param slide Slide ID (required for CoProMulti objects, ignored for CoProSingle)
#' @param returnTranspose If TRUE, forces return of transpose when accessing symmetric matrices
#' @param verbose Whether to print detailed error messages
#'
#' @return Kernel matrix as a numeric matrix
#' @export
#' @rdname getKernelMatrix
#' @aliases getKernelMatrix,CoProSingle-method
#' @aliases getKernelMatrix,CoProMulti-method
setGeneric("getKernelMatrix", 
           function(object, sigma, cellType1, cellType2, slide = NULL, 
                   returnTranspose = FALSE, verbose = TRUE) 
             standardGeneric("getKernelMatrix"))

#' @rdname getKernelMatrix
#' @export
setMethod("getKernelMatrix", "CoProSingle",
          function(object, sigma, cellType1, cellType2, slide = NULL, 
                   returnTranspose = FALSE, verbose = TRUE) {
            .getKernelMatrixSingle(object, sigma, cellType1, cellType2, 
                                 returnTranspose, verbose)
          })

#' @rdname getKernelMatrix  
#' @export
setMethod("getKernelMatrix", "CoProMulti",
          function(object, sigma, cellType1, cellType2, slide = NULL, 
                   returnTranspose = FALSE, verbose = TRUE) {
            .getKernelMatrixMulti(object, sigma, cellType1, cellType2, 
                                slide, returnTranspose, verbose)
          })

#' Internal function for single slide kernel matrix access
#' @noRd
.getKernelMatrixSingle <- function(object, sigma, cellType1, cellType2, 
                                 returnTranspose = FALSE, verbose = TRUE) {
  
  # Input validation
  .validateKernelInputs(object, sigma, cellType1, cellType2, verbose)
  
  # Check if kernelMatrices exist
  if (length(object@kernelMatrices) == 0) {
    stop("No kernel matrices found. Run computeKernelMatrix() first.")
  }
  
  # Always use flat structure access
  return(.getKernelMatrixFlat(object@kernelMatrices, sigma, cellType1, cellType2, 
                             slide = NULL, returnTranspose, verbose))
}

#' Internal function for multi slide kernel matrix access
#' @noRd
.getKernelMatrixMulti <- function(object, sigma, cellType1, cellType2, slide,
                                returnTranspose = FALSE, verbose = TRUE) {
  
  # Input validation
  .validateKernelInputs(object, sigma, cellType1, cellType2, verbose)
  
  if (is.null(slide)) {
    stop("slide parameter is required for CoProMulti objects")
  }
  
  if (!slide %in% getSlideList(object)) {
    stop(paste("Slide", slide, "not found. Available slides:", 
               paste(getSlideList(object), collapse = ", ")))
  }
  
  # Check if kernelMatrices exist
  if (length(object@kernelMatrices) == 0) {
    stop("No kernel matrices found. Run computeKernelMatrix() first.")
  }
  
  # Always use flat structure access
  return(.getKernelMatrixFlat(object@kernelMatrices, sigma, cellType1, cellType2, 
                             slide = slide, returnTranspose, verbose))
}



#' Input validation for kernel matrix access
#' @noRd
.validateKernelInputs <- function(object, sigma, cellType1, cellType2, verbose) {
  
  # Check if sigma is valid
  if (!is.numeric(sigma) || length(sigma) != 1 || is.na(sigma) || sigma <= 0) {
    stop("sigma must be a positive numeric value")
  }
  
  # Check if sigma is available
  if (length(object@sigmaValues) > 0 && !sigma %in% object@sigmaValues) {
    available_sigmas <- paste(object@sigmaValues, collapse = ", ")
    stop(paste("Sigma value", sigma, "not in object sigmaValues.", 
               "Available values:", available_sigmas))
  }
  
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

#' Throw standardized kernel not found error
#' @noRd
.throwKernelNotFoundError <- function(sigma, cellType1, cellType2, slide = NULL) {
  if (is.null(slide)) {
    stop(paste("Kernel matrix not found for", cellType1, "->", cellType2, 
               "with sigma =", sigma, ". Check if computeKernelMatrix() was run properly."))
  } else {
    stop(paste("Kernel matrix not found for", cellType1, "->", cellType2, 
               "with sigma =", sigma, "on slide", slide, 
               ". Check if computeKernelMatrix() was run properly."))
  }
}

#' Access kernel matrix from flat structure
#' @noRd
.getKernelMatrixFlat <- function(flat_kernels, sigma, cellType1, cellType2, slide = NULL,
                                returnTranspose = FALSE, verbose = TRUE) {
  
  # Create the expected flat name
  flat_name <- .createKernelMatrixName(sigma, cellType1, cellType2, slide)
  
  # Try direct access
  if (flat_name %in% names(flat_kernels)) {
    kernel_matrix <- flat_kernels[[flat_name]]
    if (returnTranspose) {
      return(t(kernel_matrix))
    } else {
      return(kernel_matrix)
    }
  }
  
  # Try symmetric access (swap cellType1 and cellType2)
  symmetric_name <- .createKernelMatrixName(sigma, cellType2, cellType1, slide)
  
  if (symmetric_name %in% names(flat_kernels)) {
    kernel_matrix <- flat_kernels[[symmetric_name]]
    if (verbose) {
      message(paste("Using transpose of kernel matrix for", cellType2, "->", cellType1))
    }
    if (returnTranspose) {
      return(kernel_matrix)  # Don't transpose again
    } else {
      return(t(kernel_matrix))
    }
  }
  
  # If we get here, kernel not found
  .throwKernelNotFoundError(sigma, cellType1, cellType2, slide)
} 