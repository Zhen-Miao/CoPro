
# Define a virtual class that is a union of 'matrix' and 'sparseMatrix'
library(Matrix)
setClassUnion("matrixOrSparseMatrix", c("matrix", "dgCMatrix",
                                        "dgTMatrix"))



#' S4 object
#'
#' @slot normalizedData A `matrix` object to store normalized data.
#' @slot normalizedDataSub A `matrix` object to store the subset of
#' the normalized data. The subset only contain relevant cell types
#' of interest specified by the user.
#' @slot locationData A `data.frame` object to store the location. It should
#' either contain two columns named by "x" and "y", or three columns named by
#' "x", "y", and "z". No other names allowed
#' @slot locationDataSub A `data.frame` object to store the subset of
#' the location data.
#' @slot metaData A `data.frame` object to store metadata for each cell.
#' @slot metaDataSub A `data.frame` object to store the subset of
#' the meta data.
#' @slot cellTypes A `vector` object with elements being character. It should
#' match the number of cells in the data matrix and each represents a cell type
#' label of a cell.
#' @slot cellTypesSub A `vector` object with elements being character. It stores
#' the subset of  the cell type labels.
#' @slot cellTypesOfInterest A `vector` object with elements being character.
#' Specifies the cell types of interest. It will be used to subset the dataset
#' @slot pcaResults A `list` object to store the PCA output for each cell type.
#' @slot distances A `list` object to store the pairwise distances between any
#' two cell types of interest.
#' @slot geneList A `vector` object with elements being character. To store the
#' gene names.
#' @slot kernelMatrices A `list` object. To store the kernel matrix generated
#' from the distance matrices.
#' @slot sigmaSquares A `vector` object with elements being numeric. To store
#' a set of sigma values used for generating the kernel matrix.
#' @slot nPCA A single numeric value. Number of PCs in to retain for downstream
#' analyses.
#' @slot skrCCAOut A `list` object. Output from the skrCCA.
#' @slot cellScores A `matrix` object. Cell scores for each cell type.
#' @slot geneScores A `matrix` object. Gene scores for each cell type.
#' @slot scalePCs A `logical` value. Whether to scale each PC before computing
#' @slot normalizedCorrelation A `list` object. Normalized correlation values
#' for each sigma value.
#' skrCCA
#'
#' @export
#'
setClass("CoPro",
  slots = list(

    ## cell by gene data matrix
    normalizedData = "matrixOrSparseMatrix",
    normalizedDataSub = "matrixOrSparseMatrix",

    ## location data
    locationData = "data.frame",
    locationDataSub = "data.frame",

    ## cell by meta.data matrix
    metaData = "data.frame",
    metaDataSub = "data.frame",

    ## cell type labels for each cell
    cellTypes = "character",
    cellTypesSub = "character",

    ## cell types of interest
    cellTypesOfInterest = "character",

    ## pca results of normalized data matrix
    pcaResults = "list",

    ## cell by cell distances
    distances = "list",

    ## geneList
    geneList = "character",
    kernelMatrices = "list",
    sigmaSquares = "numeric",
    nPCA = "numeric",
    scalePCs = "logical",

    ## skr CCA output
    skrCCAOut = "list",
    cellScores = "list",
    geneScores = "list",
    normalizedCorrelation = "list",
    sigmaSquaredChoice = "numeric"
  )
)

#' Function to create a new object
#' @importFrom methods new
#' @param normalizedData A `matrix` object to store normalized data.
#' @param locationData A `data.frame` object to store the location. It should
#' either contain two columns named by "x" and "y", or three columns named by
#' "x", "y", and "z". No other names allowed
#' @param metaData A `data.frame` object to store metadata for each cell.
#' @param cellTypes A `vector` object with elements being character. It should
#' match the number of cells in the data matrix and each represents a cell type
#' label of a cell.
#' @rdname newCoPro
#' @aliases newCoPro,CoPro-method
#' @return A `CoPro` object
#' @export
#'
setGeneric(
  "newCoPro",
  function(normalizedData, locationData, metaData, cellTypes)
    standardGeneric("newCoPro")
)


#' @rdname newCoPro
#' @aliases newCoPro,CoPro-method
#' @export
setMethod(
  "newCoPro", signature(
    "matrixOrSparseMatrix", "data.frame",
    "data.frame", "character"
  ),
  function(normalizedData, locationData, metaData, cellTypes) {
    ## check dimension of input
    if (length(cellTypes) != nrow(normalizedData) |
      nrow(normalizedData) != nrow(metaData) |
      nrow(normalizedData) != nrow(locationData)) {
      stop("input data do not match dimensionality, please check")
    }

    ## check the format of location data
    if (!all(tolower(colnames(locationData)) %in% c("x", "y", "z"))) {
      stop(paste("locationData should only contain x, y, (or z)",
        "axis info and colnames should be named accordingly",
        sep = " "
      ))
    }

    ## check cellTypes are characters
    if (!is.character(cellTypes)) {
      stop("Cell types must be character strings.")
    }


    colnames(locationData) <- tolower(colnames(locationData))

    ## check cell id and gene names
    if (is.null(rownames(metaData)) | is.null(rownames(normalizedData)) |
      is.null(rownames(locationData))) {
      stop(paste("please make sure the rownames of data,",
        "metaData, and locationData are cell barcodes",
        sep = " "
      ))
    } else if (any((rownames(metaData) != rownames(normalizedData)) |
      any(rownames(locationData) != rownames(normalizedData)))) {
      stop(paste("please make sure the cell barcodes match,",
        "between data, metaData,and locationData",
        sep = " "
      ))
    }

    ## check gene names
    if (is.null(colnames(normalizedData))) {
      stop("please make sure colnames of data are gene names")
    }

    geneList <- colnames(normalizedData)


    ## create new object
    new("CoPro",
      normalizedData = normalizedData,
      metaData = metaData,locationData = locationData,
      cellTypes = cellTypes, geneList = geneList
    )
  }
)



#' subsetData
#'
#' Take a subset of the original matrix based on cell types of interest. The
#' original data are stored without being thrown away
#'
#' @param object A `CoPro` object
#' @param cellTypesOfInterest Input cell types of interest as a vector of
#' characters for subsetting the data
#'
#' @rdname subsetData
#' @aliases subsetData,CoPro-method
#' @return A `CoPro` object with subset slots
#' @export
#'
setGeneric("subsetData", function(object, cellTypesOfInterest)
  standardGeneric("subsetData")
)

#' @rdname subsetData
#' @aliases subsetData,CoPro-method
#' @export
setMethod("subsetData", "CoPro", function(object, cellTypesOfInterest) {

  if (length(cellTypesOfInterest) < 2) {
    stop("at least two cell types are needed for this analysis")
  }

  if (!all(cellTypesOfInterest %in% object@cellTypes)) {
    stop("some cellTypesOfInterest are not in cellTypes, please check")
  }

  subsetIndices <- object@cellTypes %in% cellTypesOfInterest

  if (sum(subsetIndices) < 10) {
    stop("Fewer than 10 cells from the subset, please check the
         cellTypesOfInterest")
  }

  object@cellTypesOfInterest <- cellTypesOfInterest

  ## subset the data
  object@normalizedDataSub <- object@normalizedData[subsetIndices, ]
  object@metaDataSub <- object@metaData[subsetIndices, ]
  object@locationDataSub <- object@locationData[subsetIndices, ]
  object@cellTypesSub <- object@cellTypes[subsetIndices]

  object
})

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
setGeneric("computePCA", function(object, nPCA = 40,
                                  center = TRUE, scale. = TRUE)
  standardGeneric("computePCA")
)

#' @rdname computePCA
#' @importFrom stats setNames
#' @aliases computePCA,CoPro-method
#' @export
setMethod(
  "computePCA", "CoPro",
  function(object, nPCA = 40, center = TRUE, scale. = TRUE) {
    ## choose cell types
    if (lenght(object@cellTypesOfInterest) != 0) {
      cts <- object@cellTypesOfInterest
    } else {
      warning("no cell type of interest specified,
                      using all cell types to run the analysis")
      cts <- unique(object@cellTypesSub)
    }

    ## PCA results will be saved under the name of cell types
    object@pcaResults <- setNames(
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
          "and scaling of the data",
          sep = " "
        ))
        scaledData <- subD
      }

      ## PCA, on the matrix that is already centered and scaled
      pca <- prcomp_irlba(scaledData, center = FALSE, scale. = FALSE, n = nPCA)
      object@pcaResults[[i]] <- pca
    }

    ## return
    return(object)
  }
)

#' computeDistance between pairs of cell types
#'
#' @importFrom fields rdist
#' @importFrom utils combn
#' @importFrom stats setNames
#' @param object A `CoPro` object
#' @param distType Type of distance to compute: "Euclidean2D",
#'  "Euclidean3D", or "Morphology-Aware"
#' @param xDistScale Scale for x distance
#' @param yDistScale Scale for y distance
#' @param zDistScale Scale for z distance
#'
#' @return `CoPro` object with distance matrix computed
#' @export
#' @note To-do: add morphology-aware kernel
setGeneric("computeDistance",
           function(object, distType =
                      c("Euclidean2D", "Euclidean3D","Morphology-Aware"),
                    xDistScale = 1, yDistScale = 1,
                    zDistScale = 1) standardGeneric("computeDistance")
)

#' @rdname computeDistance
#' @aliases computeDistance,CoPro-method
#' @importFrom utils combn
#' @importFrom stats setNames
#' @importFrom fields rdist
#' @export
setMethod(
  "computeDistance", "CoPro",
  function(object, distType = c(
             "Euclidean2D", "Euclidean3D",
             "Morphology-Aware"
           ),
           xDistScale = 1,
           yDistScale = 1, zDistScale = 1) {
    ## match arg
    distType <- match.arg(distType)

    ## choose cell types
    if (length(object@cellTypesOfInterest) != 0) {
      cts <- object@cellTypesOfInterest
    } else {
      warning(paste("no cell type of interest specified,",
        "using all cell types to run the analysis",
        sep = " "
      ))
      cts <- unique(object@cellTypesSub)
    }

    ## check dist
    if (distType == "Euclidean3D") {
      if (!all(c("x", "y", "z") %in% object@locationDataSub)) {
        stop(paste("please make sure x, y, z are all available to run",
          "3D Euclidean distance calcuation",
          sep = " "
        ))
      }
    }

    ## distances are between any two cell types

    distances <- setNames(rep(list(), length = length(cts)), cts)
    for (i in cts) {
      distances[[i]] <- setNames(rep(list(), length = length(cts)), cts)
    }

    # n_mat <- length(cts)

    pair_cell_types <- combn(cts, 2)
    ct_ind_sub <- object@cellTypesSub

    for (pp in seq_len(ncol(pair_cell_types))) {
      i <- pair_cell_types[1, pp]
      j <- pair_cell_types[2, pp]
      if (distType == "Euclidean2D") {
        mat1 <- cbind(
          object@locationDataSub$x[ct_ind_sub == i] * xDistScale,
          object@locationDataSub$y[ct_ind_sub == i] * yDistScale
        )
        mat2 <- cbind(
          object@locationDataSub$x[ct_ind_sub == j] * xDistScale,
          object@locationDataSub$y[ct_ind_sub == j] * yDistScale
        )
      } else if (distType == "Euclidean3D") {
        mat1 <- cbind(
          object@locationDataSub$x[ct_ind_sub == i] * xDistScale,
          object@locationDataSub$y[ct_ind_sub == i] * yDistScale,
          object@locationDataSub$z[ct_ind_sub == i] * zDistScale
        )
        mat2 <- cbind(
          object@locationDataSub$x[ct_ind_sub == j] * xDistScale,
          object@locationDataSub$y[ct_ind_sub == j] * yDistScale,
          object@locationDataSub$z[ct_ind_sub == j] * zDistScale
        )
      } else if (distType == "Morphology-Aware") {
        stop("morphology-aware kernel is not availabe at this moment")
      }

      ## compute distance
      distances_ij <- fields::rdist(mat1, mat2)
      if (any(distances_ij == 0)) {
        warning(paste("Zero distances detected, replacing with",
          "the smallest non-zero distances, please",
          "consider checking the location of cells",
          sep = " "
        ))
        distances_ij[distances_ij == 0] <-
          min(distances_ij[distances_ij != 0])
      }

      ## save the distances
      distances[[i]][[j]] <- distances_ij
    }


    object@distances <- distances
    return(object)
  }
)

#' Compute Kernel Matrix for CoPro
#'
#' This method calculates the kernel matrices for pairs of cell types based on
#' their distances and a range of sigma square values.
#' The matrices are adjusted by clipping the upper quantile of
#'  the values to reduce the effect of outliers.
#' The results are stored within the object.
#'
#' @importFrom utils combn
#' @param object A `CoPro` object.
#' @param sigmaSquares A vector of sigma square values used for kernel calculation.
#' @param lowerLimit The lower limit for the kernel function, default is 0.05.
#' @param upperQuantile The quantile used for clipping the kernel values, default is 0.8.
#' @return The `CoPro` object with computed kernel matrices added. The kernel
#' matrices are organized into a three-layer nested list object. The first layer
#' is indexed by the sigma value, and the second and the third layers are cell
#' types
#' @export
setGeneric("computeKernelMatrix",
           function(object, sigmaSquares,lowerLimit = 0.05, upperQuantile = 0.8,
                    verbose = TRUE) standardGeneric("computeKernelMatrix")
)


#' @rdname computeKernelMatrix
#' @aliases computeKernelMatrix,CoPro-method
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @export
setMethod(
  "computeKernelMatrix", "CoPro",
  function(object, sigmaSquares,
           lowerLimit = 0.05, upperQuantile = 0.8, verbose = TRUE) {
    ## make sure distance matrix exist
    if (length(object@distances) == 0) {
      stop("Please run computeDistance before computing kernel")
    }


    cts <- object@cellTypesOfInterest
    n_mat <- length(cts)
    if (n_mat < 2) {
      stop("At least two cell types are needed to compute kernel matrices.")
    }

    if (length(sigmaSquares) == 0) {
      warning("No Sigma specified, setting to the 5% quantile of cell distance")
      dist12 <- object@distances[[1]][[2]]
      sigmaSquares <- quantile(dist12[dist12 > 0], 0.05)
    }

    if (!is.numeric(sigmaSquares)) {
      stop("SigmaSquares must be numeric values or a vector")
    }

    ## save the sigmaSquares
    object@sigmaSquares <- sigmaSquares

    ## Initialize the list of kernel matrices
    kernel_mat <- vector("list", length(sigmaSquares))
    sigma_names <- paste("sigma", sigmaSquares, sep = "_")
    names(kernel_mat) <- sigma_names
    for (t in sigma_names) {
      kernel_mat[[t]] <- setNames(vector("list", n_mat), cts)
      for (i in cts) {
        kernel_mat[[t]][[i]] <- setNames(vector("list", n_mat), cts)
      }
    }

    pair_cell_types <- combn(cts, 2)

    for (tt in seq_along(sigmaSquares)) {
      t <- sigma_names[tt]
      sigma_square_choose <- sigmaSquares[tt]

      for (pp in seq_len(ncol(pair_cell_types))) {
        i <- pair_cell_types[1, pp]
        j <- pair_cell_types[2, pp]

        kernel_current <- kernel_from_distance(
          sigma_square = sigma_square_choose,
          dist_mat = object@distances[[i]][[j]],
          lower_limit = lowerLimit
        )

        ## print info
        if(verbose){
          cat(paste("Current Sigma value is ",sigma_square_choose ))
          cat(paste("Quantiles of N_neighbors for cell type", i, "\n"))
          print(quantile(rowSums(kernel_current != 0)))
          cat(paste("Quantiles of N_neighbors for cell type", j, "\n"))
          print(quantile(colSums(kernel_current != 0)))
          cat("\n")
        }


        ## Clipping large values
        upper_clip <- quantile(kernel_current[kernel_current != 0], upperQuantile)
        kernel_current[kernel_current >= upper_clip] <- upper_clip

        kernel_mat[[t]][[i]][[j]] <- kernel_current
      }
    }

    object@kernelMatrices <- kernel_mat
    return(object)
  }
)


#' runSkrCCA
#' @importFrom stats setNames
#' @param object A CoPro object
#' @param scalePCs Whether to scale each PCs to a uniform variance before
#' running the program
#' @param maxIter Maximum iterations
#'
#' @return CoPro object with distnace matrix computed
#' @export
#'
setGeneric("runSkrCCA",
           function(object,scalePCs = TRUE,maxIter = 200)
             standardGeneric("runSkrCCA")
)

#' @rdname runSkrCCA
#' @aliases runSkrCCA,CoPro-method
#' @importFrom stats setNames
#' @export
setMethod(
  "runSkrCCA", "CoPro",
  function(object,
           scalePCs = TRUE,
           maxIter = 200) {
    ## check whether the kernel matrix is available
    if (length(object@kernelMatrices) == 0) {
      stop("Kernel matrix is empty, please run computeKernelMatrix first")
    }

    ## record whether PCs have been scaled
    if (length(object@scalePCs) == 0) {
      object@scalePCs <- scalePCs
    } else if (object@scalePCs != scalePCs) {
      stop("Previously set scalePCs was different from the function input")
    }

    ## check sigmaSquares
    if (length(object@sigmaSquares) == 0) {
      stop("sigmaSquares is empty, please specify")
    } else {
      sigmaSquares <- object@sigmaSquares
    }

    ## choose cell types
    if (length(object@cellTypesOfInterest) != 0) {
      cts <- object@cellTypesOfInterest
    } else {
      warning("no cell type of interest specified,
                      using all cell types to run the analysis")
      cts <- unique(object@cellTypesSub)
    }

    ## get PC matrices
    allPCs <- object@pcaResults
    PCmats <- setNames(
      vector("list", length = length(allPCs)),
      names(allPCs)
    )

    ## optionally, scale the PCs before running CCA
    if (scalePCs) {
      for (i in names(allPCs)) {
        pca_A_sd <- allPCs[[i]]$sdev
        PCmats[[i]] <- scale(allPCs[[i]]$x,
          center = FALSE,
          scale = pca_A_sd
        )
      }
    } else {
      for (i in names(allPCs)) {
        PCmats[[i]] <- allPCs[[i]]$x
      }
    }

    ## run across different sigma values
    cca_out <- vector("list", length = length(sigmaSquares))
    sigma_names <- paste("sigma", sigmaSquares, sep = "_")
    names(cca_out) <- sigma_names

    ## for loop to run the analysis
    for (tt in seq_along(sigmaSquares)) {
      t <- sigma_names[tt]
      cca_result <- optimize_bilinear_multi(
        X_list = PCmats,
        K_list = object@kernelMatrices[[t]],
        max_iter = maxIter
      )
      names(cca_result) <- cts
      cca_out[[t]] <- cca_result
    }

    object@skrCCAOut <- cca_out
    return(object)
  }
)


#' Compute Normalized Correlation for CoPro
#'
#' This method calculates the normalized correlation between pairs of cell types
#' based on CCA weights and the respective kernel matrix. It uses the spectral norm
#' of the kernel matrix for normalization.
#'
#' @param object A `CoPro` object containing CCA results and kernel matrices.
#' @return The `CoPro` object with the normalized correlation value
#' between any pair of cell types
#' added as a new slot, `normalizedCorrelation`.
#' @export
#'
setGeneric("computeNormalizedCorrelation",
           function(object) standardGeneric("computeNormalizedCorrelation")
)


#' @rdname computeNormalizedCorrelation
#' @aliases computeNormalizedCorrelation,CoPro-method
#' @importFrom utils combn
#' @export
setMethod(
  "computeNormalizedCorrelation", "CoPro",
  function(object) {
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

    ## check sigmaSquares
    if (length(object@sigmaSquares) == 0) {
      stop("sigmaSquares is empty, please specify")
    }

    sigmaSquares <- object@sigmaSquares


    ## get PC matrices
    allPCs <- object@pcaResults
    PCmats <- setNames(
      vector("list", length = length(allPCs)),
      names(allPCs)
    )

    ## optionally, scale the PCs before running CCA
    if (scalePCs) {
      for (i in names(allPCs)) {
        pca_A_sd <- allPCs[[i]]$sdev
        PCmats[[i]] <- scale(allPCs[[i]]$x,
          center = FALSE,
          scale = pca_A_sd
        )
      }
    } else {
      for (i in names(allPCs)) {
        PCmats[[i]] <- allPCs[[i]]$x
      }
    }


    pair_cell_types <- combn(cts, 2)

    correlation_value <- vector("list", length = length(sigmaSquares))
    sigma_names <- paste("sigma", sigmaSquares, sep = "_")
    names(correlation_value) <- sigma_names


    for (tt in seq_along(sigmaSquares)) {
      t <- sigma_names[tt]
      correlation_value[[t]] <- data.frame(
        sigmaSquares = sigmaSquares[tt],
        cellType1 = pair_cell_types[1, ],
        cellType2 = pair_cell_types[2, ],
        normalizedCorrelation = numeric(length = ncol(pair_cell_types)),
        stringsAsFactors = FALSE
      )
      for (pp in seq_len(ncol(pair_cell_types))) {
        cellType1 <- pair_cell_types[1, pp]
        cellType2 <- pair_cell_types[2, pp]

        w_1 <- object@skrCCAOut[[t]][[cellType1]]
        w_2 <- object@skrCCAOut[[t]][[cellType2]]

        A <- PCmats[[cellType1]]
        B <- PCmats[[cellType2]]
        K <- object@kernelMatrices[[t]][[cellType1]][[cellType2]]

        A_w1 <- A %*% w_1
        B_w2 <- B %*% w_2

        ## Calculate the spectral norm of the kernel matrix
        norm_K12 <- norm(K, type = "2")

        ## Calculate normalized correlation
        correlation_value[[t]]$normalizedCorrelation[pp] <-
          (t(A_w1) %*% K %*% B_w2) /
            (sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K12)
      }
    }

    ## Store the result in the object
    object@normalizedCorrelation <- correlation_value

    ## obtain the sigmaSqured value with the highest
    ## normalized correlation
    ncorr <- do.call(rbind, correlation_value)
    ncorr$ct12 <- paste(ncorr$cellType1, ncorr$cellType2, sep = "-")

    # Calculate the mean of column 2 for each unique value in column 1
    meanCorr <- tapply(ncorr$normalizedCorrelation,
                       ncorr$sigmaSquares, mean)

    # Find the value of column 1 with the highest mean in column 2
    sigmaSquaredChoice <- as.numeric(names(which.max(meanCorr)))
    object@sigmaSquaredChoice <- sigmaSquaredChoice

    ## Return the modified object
    return(object)
  }
)



#' computeGeneAndCellScores
#'
#' @param object A `CoPro` object containing CCA results
#' and kernel matrices.
#'
#' @return A `CoPro` object with gene and cell score computed
#' @export
#'
setGeneric("computeGeneAndCellScores",
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
      stop("Kernel matrices are not available. Please compute the kernel matrices first.")
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

    ## check sigmaSquares
    if (length(object@sigmaSquares) == 0) {
      stop("sigmaSquares is empty, please specify")
    }

    sigmaSquares <- object@sigmaSquares

    ## load whether the PCs are being scaled prior to CCA

    if (length(object@scalePCs) == 0) {
      stop("object@scalePCs not specified")
    }
    scalePCs <- object@scalePCs

    ## get PC matrices
    allPCs <- object@pcaResults
    PCmats <- setNames(
      vector("list", length = length(allPCs)),
      names(allPCs)
    )

    ## optionally, scale the PCs before running CCA
    if (scalePCs) {
      for (i in names(allPCs)) {
        pca_A_sd <- allPCs[[i]]$sdev
        PCmats[[i]] <- scale(allPCs[[i]]$x,
          center = FALSE,
          scale = pca_A_sd
        )
      }
    } else {
      for (i in names(allPCs)) {
        PCmats[[i]] <- allPCs[[i]]$x
      }
    }

    # pair_cell_types <- combn(cts, 2)

    sigma_names <- paste("sigma", sigmaSquares, sep = "_")

    ## cell scores and gene scores are both by cell types
    cellScores <- setNames(
      vector(mode = "list", length = length(cts)),
      cts
    )
    for (i in cts) {
      cellScores[[i]] <- matrix(
        nrow = sum(object@cellTypesSub == i),
        ncol = length(sigmaSquares)
      )
      colnames(cellScores[[i]]) <- sigma_names
      rownames(cellScores[[i]]) <- rownames(object@normalizedDataSub)[
        object@cellTypesSub == i
      ]
    }

    ## gene scores
    geneScores <- setNames(
      vector(mode = "list", length = length(cts)),
      cts
    )
    for (i in cts) {
      geneScores[[i]] <- matrix(
        nrow = ncol(object@normalizedDataSub),
        ncol = length(sigmaSquares)
      )
      colnames(geneScores[[i]]) <- sigma_names
      rownames(geneScores[[i]]) <- colnames(object@normalizedDataSub)
    }

    ## go over all cell types, then over all sigma values
    for (tt in seq_along(sigmaSquares)) {
      t <- sigma_names[tt]
      for (i in cts) {
        w_1 <- object@skrCCAOut[[t]][[i]]
        if (scalePCs) {
          geneScores[[i]][, t] <- as.vector(
            matrix(w_1 * allPCs[[i]]$sdev, nrow = 1) %*%
              t(allPCs[[i]]$rotation)
          )
        } else {
          geneScores[[i]][, t] <- as.vector(
            matrix(w_1, nrow = 1) %*%
              t(allPCs[[i]]$rotation)
          )
        }
        cellScores[[i]][, t] <- as.vector(PCmats[[i]] %*% w_1)
      }
    }

    ## change the column names to make it more informative
    for (i in cts) {
      colnames(cellScores[[i]]) <-
        paste0("cellScore_", colnames(cellScores[[i]]))
      colnames(geneScores[[i]]) <-
        paste0("geneScore_", colnames(geneScores[[i]]))
    }

    ## save cellscores and gene scores
    object@cellScores <- cellScores
    object@geneScores <- geneScores

    ## add cell score information to the cell metadata
    meta_t <- stats::setNames(vector(mode = "list", length = length(cts)),
                             cts)
    for (t in cts) {
      meta_t[[t]] <- object@metaDataSub[object@cellTypesSub == t,]
      meta_t[[t]] <- cbind(meta_t[[t]], cellScores[[t]][rownames(loc_t[[t]]),])
    }

    ## combine each cell type meta.data back
    names(meta_t) <- NULL
    meta_all <- do.call(rbind, meta_t)[rownames(object@metaDataSub),]
    object@metaDataSub <- meta_all

    return(object)
  }
)
