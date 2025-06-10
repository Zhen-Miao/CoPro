


#' computePCA with irlba package
#' @importFrom stats setNames
#' @importFrom irlba prcomp_irlba
#' @param object A `CoPro` object
#' @param nPCA Number of Pcs
#' @param center Whether to center data before PCA
#' @param scale. Whether to scale data by sd before PCA
#'
#' @return A `CoPro` object
#'
#' @rdname computePCA
#' @aliases computePCA,CoPro-method
#'
#' @export
#'
setGeneric("computePCA",
           function(object, nPCA = 40, center = TRUE,
                    scale. = TRUE) standardGeneric("computePCA")
)

#' @rdname computePCA
#' @importFrom stats setNames
#' @aliases computePCA,CoPro-method
#' @export
setMethod(
  "computePCA", "CoPro",
  function(object, nPCA = 40, center = TRUE, scale. = TRUE) {
    ## choose cell types
    if (length(object@cellTypesOfInterest) != 0) {
      cts <- object@cellTypesOfInterest
    } else {
      warning("no cell type of interest specified,
                      using all cell types to run the analysis")
      cts <- unique(object@cellTypesSub)
    }

    ## PCA results will be saved under the name of cell types
    object@pcaGlobal <- setNames(
      rep(list(), length = length(cts)),
      cts
    )

    ## iterate over cell types
    for (i in cts) {
      ## cell type specific subset
      subD <- as.matrix(object@normalizedDataSub[object@cellTypesSub == i, ])

      ## center and scale the data using our own function, because
      ## for genes with mostly zero expression, this will not scale
      ## them up too much

      if (center & scale.) {
        scaledData <- center_scale_matrix_opt(subD)
      } else if (center) {
        scaledData <- t(t(subD) - colMeans(subD))
      } else {
        warning(paste("It is not recommended to skip both centering,",
          "and scaling of the data, unless the data has been centered and",
          "scaled when creating the CoPro object.",
          sep = " "
        ))
        scaledData <- subD
      }

      ## PCA, on the matrix that is already centered and scaled
      pca <- prcomp_irlba(scaledData, center = FALSE, scale. = FALSE, n = nPCA)
      object@pcaGlobal[[i]] <- pca
    }

    ## return
    return(object)
  }
)



#' Compute Normalized Correlation (approximation)
#'
#' This method calculates the normalized correlation between pairs of cell types
#' based on CCA weights and the respective kernel matrix. It uses
#' the spectral norm of the kernel matrix for normalization.
#'
#' @param object A `CoPro` object containing CCA results and kernel matrices.
#' @param tol tolerance for approximate SVD calculation
#' @return The `CoPro` object with the normalized correlation value
#' between any pair of cell types
#' added as a new slot, `normalizedCorrelation`.
#' @export
#'
setGeneric(
  "computeNormalizedCorrelation",
  function(object, tol = 1e-4) standardGeneric("computeNormalizedCorrelation")
)


#' @rdname computeNormalizedCorrelation
#' @aliases computeNormalizedCorrelation,CoPro-method
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @export
setMethod(
  "computeNormalizedCorrelation", "CoPro",
  function(object, tol = 1e-4) {
    ## Check for required components
    if (length(object@skrCCAOut) == 0) {
      stop("CCA results are not available. Please run CCA first.")
    }
    if (length(object@kernelMatrices) == 0) {
      stop(paste(
        "Kernel matrices are not available.",
        "Please compute the kernel matrices first."
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

    ## load whether the PCs are being scaled prior to CCA

    if (length(object@scalePCs) == 0) {
      stop("object@scalePCs not specified")
    }
    scalePCs <- object@scalePCs

    ## check sigmaValues
    if (length(object@sigmaValues) == 0) {
      stop("`sigmaValues` is empty, please specify")
    }

    sigmaValues <- object@sigmaValues

    PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

    pair_cell_types <- combn(cts, 2)

    correlation_value <- vector("list", length = length(sigmaValues))
    sigma_names <- paste("sigma", sigmaValues, sep = "_")
    names(correlation_value) <- sigma_names

    nCC <- object@nCC

    ## calculate all spectral norms
    cat("Calculating spectral norms, ",
        "depending on the data size, this may take a while. \n")
    norm_K12 <- setNames(vector(mode = "list", length = length(sigma_names)),
                         sigma_names)

    for (t in sigma_names) {
      norm_K12[[t]] <- setNames(vector(mode = "list", length = length(cts)),
                                cts)
      for (i in cts) {
        norm_K12[[t]][[i]] <- setNames(vector(mode = "list",
                                              length = length(cts)), cts)
      }
    }

    for (t in sigma_names) {
      for (pp in seq_len(ncol(pair_cell_types))) {
        cellType1 <- pair_cell_types[1, pp]
        cellType2 <- pair_cell_types[2, pp]
        K <- object@kernelMatrices[[t]][[cellType1]][[cellType2]]
        ## Calculate the spectral norm of the kernel matrix
        svd_result <- irlba::irlba(K, nv = 1, tol = tol)
        norm_K12[[t]][[cellType1]][[cellType2]] <- svd_result$d[1]
      }
    }

    cat("Finished calculating spectral norms \n")


    for (tt in seq_along(sigmaValues)) {
      t <- sigma_names[tt]
      correlation_value[[t]] <- data.frame(
        sigmaValues = sigmaValues[tt],
        cellType1 = rep(pair_cell_types[1, ], times = nCC),
        cellType2 = rep(pair_cell_types[2, ], times = nCC),
        CC_index = rep(x = 1:nCC, each = ncol(pair_cell_types)),
        normalizedCorrelation = numeric(length = ncol(pair_cell_types) * nCC),
        stringsAsFactors = FALSE
      )
      for (pp in seq_len(ncol(pair_cell_types))) {
        for (cc_index in seq_len(nCC)) {
          cellType1 <- pair_cell_types[1, pp]
          cellType2 <- pair_cell_types[2, pp]

          w_1 <- object@skrCCAOut[[t]][[cellType1]][, cc_index, drop = FALSE]
          w_2 <- object@skrCCAOut[[t]][[cellType2]][, cc_index, drop = FALSE]

          A <- PCmats[[cellType1]]
          B <- PCmats[[cellType2]]

          A_w1 <- A %*% w_1
          B_w2 <- B %*% w_2

          ## get pre-calculated spectral norm
          K <- object@kernelMatrices[[t]][[cellType1]][[cellType2]]
          norm_K12_sel <- norm_K12[[t]][[cellType1]][[cellType2]]

          ## Calculate normalized correlation
          correlation_value[[t]]$"normalizedCorrelation"[
            pp + (cc_index - 1) * ncol(pair_cell_types)] <-
            (t(A_w1) %*% K %*% B_w2) /
            (sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K12_sel)
        }
      }
    }

    ## Store the result in the object
    object@normalizedCorrelation <- correlation_value

    ## obtain the sigma value with the highest
    ## normalized correlation
    ncorr <- do.call(rbind, correlation_value)
    ncorr$ct12 <- paste(ncorr$cellType1, ncorr$cellType2, sep = "-")

    # Calculate the mean of column 2 for each unique value in column 1
    ## only for cc_index == 1
    meanCorr <- tapply(
      ncorr$"normalizedCorrelation"[ncorr$"CC_index" == 1],
      ncorr$"sigmaValues"[ncorr$"CC_index" == 1], mean
    )

    # Find the value of column 1 with the highest mean in column 2
    sigmaValueChoice <- as.numeric(names(which.max(meanCorr)))
    object@sigmaValueChoice <- sigmaValueChoice

    ## Return the modified object
    return(object)
  }
)


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
