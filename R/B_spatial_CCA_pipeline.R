




#' computeGeneAndCellScores
#' @importFrom stats setNames
#' @param object A `CoPro` object containing CCA results
#' and kernel matrices.
#'
#' @return A `CoPro` object with gene and cell score computed
#' @export
#'
setGeneric(
  "computeGeneAndCellScores",
  function(object) standardGeneric("computeGeneAndCellScores")
)

#' @rdname computeGeneAndCellScores
#' @aliases computeGeneAndCellScores,CoPro-method
#' @export
setMethod(
  "computeGeneAndCellScores", "CoPro",
  function(object) {
    ## Check for required components
    if (length(object@skrCCAOut) == 0) {
      stop("CCA results are not available. Please run CCA first.")
    }
    if (length(object@kernelMatrices) == 0) {
      stop(paste("Kernel matrices are not available.",
           "Please compute the kernel matrices first."))
    }
    ## choose cell types
    if (length(object@cellTypesOfInterest) != 0) {
      cts <- object@cellTypesOfInterest
    } else {
      warning(paste(
        "no cell type of interest specified,",
        "using all cell types to run the analysis"
      ))
      cts <- unique(object@cellTypesSub)
    }

    ## check sigmaValues
    if (length(object@sigmaValues) == 0) {
      stop("sigmaValues is empty, please specify")
    }

    sigmaValues <- object@sigmaValues
    nCC <- object@nCC

    ## load whether the PCs are being scaled prior to CCA

    if (length(object@scalePCs) == 0) {
      stop("object@scalePCs not specified")
    }
    scalePCs <- object@scalePCs

    PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

    sigma_names <- paste("sigma", sigmaValues, sep = "_")

    ## cell scores and gene scores are both by cell types
    cellScores <- setNames(
      vector(mode = "list", length = length(sigmaValues)),
      sigma_names
    )
    for (t in sigma_names) {
      cellScores[[t]] <- setNames(
        vector(mode = "list", length = length(cts)),
        cts
      )
    }


    for (t in sigma_names) {
      for (i in cts) {
        cellScores[[t]][[i]] <- matrix(
          nrow = sum(object@cellTypesSub == i),
          ncol = nCC
        )
        colnames(cellScores[[t]][[i]]) <- paste0("CC_", 1:nCC)
        rownames(cellScores[[t]][[i]]) <- rownames(object@normalizedDataSub)[
          object@cellTypesSub == i
        ]
      }

    }

    ## gene scores

    geneScores <- setNames(
      vector(mode = "list", length = length(sigmaValues)),
      sigma_names
    )
    for (t in sigma_names) {
      geneScores[[t]] <- setNames(
        vector(mode = "list", length = length(cts)),
        cts
      )
    }

    for (t in sigma_names) {
      for (i in cts) {
        geneScores[[t]][[i]] <- matrix(
          nrow = ncol(object@normalizedDataSub),
          ncol = nCC
        )
        colnames(geneScores[[t]][[i]]) <- paste0("CC_", 1:nCC)
        rownames(geneScores[[t]][[i]]) <- colnames(object@normalizedDataSub)
      }
    }


    ## go over all cell types, then over all sigma values
    for (tt in seq_along(sigmaValues)) {
      t <- sigma_names[tt]
      for (i in cts) {
        for (cc_index in seq_len(nCC)) {
          cc_name <- paste0("CC_", cc_index)
          w_1 <- object@skrCCAOut[[t]][[i]][, cc_index, drop = FALSE]
          if (scalePCs) {
            geneScores[[t]][[i]][, cc_name] <- as.vector(
              matrix(w_1 * object@pcaGlobal[[i]]$sdev, nrow = 1) %*%
                t(object@pcaGlobal[[i]]$rotation)
            )
          } else {
            geneScores[[t]][[i]][, cc_name] <- as.vector(
              matrix(w_1, nrow = 1) %*%
                t(object@pcaGlobal[[i]]$rotation)
            )
          }
          cellScores[[t]][[i]][, cc_name] <- as.vector(PCmats[[i]] %*% w_1)
        }

      }
    }

    ## save cellscores and gene scores
    object@cellScores <- cellScores
    object@geneScores <- geneScores

    ## add cell score information to the cell metadata
    meta_t <- stats::setNames(
      vector(mode = "list", length = length(cts)),
      cts
    )
    for (i in cts) {
      meta_t[[i]] <- object@metaDataSub[object@cellTypesSub == i, ]
      for (t in sigma_names){
        for (cc_index in seq_len(nCC)) {
          cc_name <- paste0("CC_", cc_index)
          cellScoreColName <- paste0("cellScore_", t, "_cc_index_", cc_index)
          meta_t[[i]][, cellScoreColName] <-
            cellScores[[t]][[i]][rownames(meta_t[[i]]), cc_name]
        }

      }

    }

    ## combine each cell type meta.data back
    names(meta_t) <- NULL
    meta_all <- do.call(rbind, meta_t)[rownames(object@metaDataSub), ]
    object@metaDataSub <- meta_all

    return(object)
  }
)



#' Assign distance matrix manually
#'
#' @param object A `CoPro` object
#' @param distanceList A list object that contains all pairwise distances
#' between any two pairs of cells.
#'
#' @return A `CoPro` object with specified
#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
#'
setGeneric("assignDistanceManually",
           function(object,
                    distanceList) standardGeneric("assignDistanceManually")
)


#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
setMethod(
  "assignDistanceManually", "CoPro",
  function(object, distanceList) {
    if (!is.list(distanceList)) {
      stop(paste(
        "distanceList must be a nested list object with names",
        "specified by cell types"
      ))
    }

    ## choose cell types
    if (length(object@cellTypesOfInterest) != 0) {
      cts <- object@cellTypesOfInterest
    } else {
      warning(paste(
        "no cell type of interest specified,",
        "using all cell types to run the analysis"
      ))
      cts <- unique(object@cellTypesSub)
    }

    if (names(distanceList) != cts) {
      stop(paste(
        "The names of distanceList do not match cell types",
        "of interest"
      ))
    }

    for (i in cts) {
      if (names(distanceList[[i]]) != cts) {
        stop(paste("The names of distanceList[[", i,
          "]] do not match cell types ",
          "of interest",
          sep = ""
        ))
      }
    }

    object@distances <- distanceList
    return(object)
  }
)



setMethod("show", "CoPro",
          function(object) {
            # Header
            cat("'CoPro' object for spatial coordinated progression detection\n")
            cat("------------------------\n")

            # Main metrics
            cat(sprintf("Number of cells: %d\n", nrow(object@normalizedData)))
            cat(sprintf("Number of genes: %d\n", ncol(object@normalizedData)))

            # Processing status
            cat("\nProcessing steps completed:\n")
            if(length(object@pcaGlobal) != 0) cat("- PCA\n")
            if(length(object@skrCCAOut) != 0) cat("- skrCCA\n")

            # Additional information
            if(length(object@metaData) > 0) {
              cat("\nAvailable metadata fields:\n")
              cat(paste("-", names(object@metaData), collapse = "\n"))
              cat("\n")
            }

      invisible(x = object)
          }
)
