#' Show method for CoPro objects
#'
#' This method provides a summary of the CoPro object, including the number of cells,
#' genes, slides, and cell types of interest. It also displays the processing steps
#' completed and the parameters used.
#'
#' @param object A CoPro object
#' @return The CoPro object (invisibly)
#' @importFrom methods show
#' @rdname show
#' @aliases show,CoPro-method
#' @export
setMethod("show", "CoPro",
          function(object) {
            # Detect if this is a multi-slide object
            is_multi <- inherits(object, "CoProMulti") || inherits(object, "CoProm")
            
            # Header
            if (is_multi) {
              cat("'CoProMulti' object for multi-slide spatial coordinated progression detection\n")
            } else {
              cat("'CoProSingle' object for spatial coordinated progression detection\n")
            }
            cat("------------------------\n")

            # Main metrics
            cat(sprintf("Number of cells: %d\n", nrow(object@normalizedData)))
            cat(sprintf("Number of genes: %d\n", ncol(object@normalizedData)))
            
            # Handle slide information appropriately
            if (is_multi) {
              cat(sprintf("Number of slides: %d\n", length(object@slideList)))
              if (length(object@slideList) > 0) {
                cat(sprintf("Slide IDs: %s\n", paste(object@slideList, collapse = ", ")))
              }
            } else {
              cat("Number of slides: 1 (single-slide analysis)\n")
            }
            
            # Cell type information
            if (length(object@cellTypesOfInterest) > 0) {
              cat(sprintf("Cell types of interest: %s\n", 
                         paste(object@cellTypesOfInterest, collapse = ", ")))
            }

            # Processing status
            cat("\nProcessing steps completed:\n")
            if(length(object@pcaGlobal) != 0) cat("- PCA\n")
            if(length(object@distances) != 0) cat("- Distance computation\n")
            if(length(object@kernelMatrices) != 0) cat("- Kernel matrix computation\n")
            if(length(object@skrCCAOut) != 0) cat("- skrCCA\n")
            if(length(object@cellScores) != 0) cat("- Cell and gene scores\n")
            if(length(object@normalizedCorrelation) != 0) cat("- Normalized correlation\n")
            
            # Parameters
            if(length(object@nPCA) > 0) {
              cat(sprintf("\nParameters:\n- nPCA: %d\n", object@nPCA))
            }
            if(length(object@nCC) > 0) {
              cat(sprintf("- nCC: %d\n", object@nCC))
            }
            if(length(object@sigmaValues) > 0) {
              cat(sprintf("- Sigma values: %s\n", 
                         paste(round(object@sigmaValues, 4), collapse = ", ")))
            }
            if(length(object@sigmaValueChoice) > 0) {
              cat(sprintf("- Optimal sigma: %g\n", object@sigmaValueChoice))
            }

            # Additional information
            if(length(object@metaData) > 0 && ncol(object@metaData) > 0) {
              cat("\nAvailable metadata fields:\n")
              cat(paste("-", names(object@metaData), collapse = "\n"))
              cat("\n")
            }

            invisible(x = object)
          }
)
