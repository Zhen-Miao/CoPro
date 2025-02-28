# Define a virtual class that is a union of 'matrix' and 'sparseMatrix'

setClassUnion("matrixOrSparseMatrix", c(
  "matrix", "dgCMatrix", "dgTMatrix"
))

setClassUnion("factorOrCharacter", c("factor", "character"))
setClassUnion("matrixOrDataFrame", c("matrix", "data.frame"))


#' CoPro object of spatial transcriptomics data
#' @import Matrix
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
#' @slot sigmaValues A `vector` object with elements being numeric. To store
#' a set of sigma values used for generating the kernel matrix.
#' @slot nPCA A single numeric value. Number of PCs to retain for downstream
#' analyses.
#' @slot nCC A single numeric value. Number of canonical components to retain
#' for downstream analyses.
#' @slot scalePCs A `logical` value. Whether to scale each PC before computing
#'  skrCCA
#' @slot skrCCAOut A `list` object. Output from the skrCCA.
#' @slot skrCCAPermuOut A `list` object. Output from the skrCCA after
#' permutation. This helps establish the null distribution
#' @slot cellPermu A `list` object that stores the cell permutation labels
#' @slot nPermu A `numeric` value specifying the number of permutations
#'  conducted.
#' @slot cellScores A `matrix` object. Cell scores for each cell type.
#' @slot geneScores A `matrix` object. Gene scores for each cell type.
#' @slot geneScoresTest A `list` object. Tested gene scores
#' @slot normalizedCorrelation A `list` object. Normalized correlation values
#' for each sigma value.
#' @slot normalizedCorrelationPermu A `list` object.
#'  Normalized correlation values
#' for each sigma value after permutation
#' @slot sigmaValueChoice A `numeric` value. The optimal sigma squared based
#' on the median normalized correlation value.
#' @name CoPro-class
#' @export
setClass("CoPro",
  slots = list(

    ## cell by gene data matrix
    normalizedData = "matrixOrSparseMatrix",
    normalizedDataSub = "matrixOrSparseMatrix",

    ## location data
    locationData = "matrixOrDataFrame",
    locationDataSub = "matrixOrDataFrame",

    ## cell by meta.data matrix
    metaData = "data.frame",
    metaDataSub = "data.frame",

    ## cell type labels for each cell
    cellTypes = "factorOrCharacter",
    cellTypesSub = "factorOrCharacter",

    ## cell types of interest
    cellTypesOfInterest = "character",

    ## pca results of normalized data matrix
    pcaResults = "list",

    ## cell by cell distances
    distances = "list",

    ## geneList
    geneList = "character",
    kernelMatrices = "list",
    sigmaValues = "numeric",
    nPCA = "numeric",
    scalePCs = "logical",
    nCC = "numeric",

    ## skr CCA output
    skrCCAOut = "list",
    cellScores = "list",
    geneScores = "list",
    geneScoresTest = "list",
    normalizedCorrelation = "list",
    sigmaValueChoice = "numeric",

    ## permutation output
    nPermu = "numeric",
    skrCCAPermuOut = "list",
    cellPermu = "list",
    normalizedCorrelationPermu = "list"
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
  function(normalizedData, locationData, metaData,
           cellTypes)  standardGeneric("newCoPro")
)


#' @rdname newCoPro
#' @aliases newCoPro,CoPro-method
#' @export
setMethod(
  "newCoPro", signature(
    "matrixOrSparseMatrix", "matrixOrDataFrame",
    "data.frame", "factorOrCharacter"
  ),
  function(normalizedData, locationData, metaData, cellTypes) {
    ## check dimension of input
    if (length(cellTypes) != nrow(normalizedData) |
      nrow(normalizedData) != nrow(metaData) |
      nrow(normalizedData) != nrow(locationData)) {
      stop("input data do not match dimensionality, please check")
    }

    ## check the format of location data
    if (!all(c("x", "y") %in% tolower(colnames(locationData)))) {
      stop(paste("locationData should contain x, y, (or z)",
                 "axis info and colnames should be named accordingly"
      ))
    }
    colnames(locationData) <- tolower(colnames(locationData))

    ## convert cellTypes to characters
    if (!is.character(cellTypes)) {
      cellTypes <- as.character(cellTypes)
    }

    ## convert locationData to data.frame
    if (is.matrix(locationData)) {
      locationData <- as.data.frame(locationData)
    }

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
    methods::new("CoPro",
      normalizedData = normalizedData,
      metaData = metaData, locationData = locationData,
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
setGeneric("subsetData",
           function(object,
                    cellTypesOfInterest) standardGeneric("subsetData")
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

  return(object)
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
          "and scaling of the data, unless the data has been centered and",
          "scaled when creating the CoPro object.",
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
#' @importFrom stats setNames quantile
#' @param object A `CoPro` object
#' @param distType Type of distance to compute: "Euclidean2D",
#'  "Euclidean3D", or "Morphology-Aware"
#' @param xDistScale Scale for x distance
#' @param yDistScale Scale for y distance
#' @param zDistScale Scale for z distance
#' @param verbose Whether to print info about the quantile of the distance
#' @param normalizeDistance Whether to normalize distance? The normalization
#'  will make sure that the 0.01% cell-cell distance will become 0.01, thus
#'  ensuring no matter which input scale is used for the distance matrix,
#'  the output will roughly be in mm^3. This ensures that the kernel sizes
#'  from 0.001 to 0.1 will make sense. Default = TRUE
#' @param truncateLowDist Whether to truncate small distances so that the cells
#'  that are nearly overlapping with each other do not have a super small
#'  distance. Default = TRUE.
#' @return `CoPro` object with distance matrix computed
#' @export
#' @note To-do: add morphology-aware kernel
setGeneric(
  "computeDistance",
  function(object, distType =
             c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
           xDistScale = 1, yDistScale = 1,
           zDistScale = 1, normalizeDistance = TRUE, truncateLowDist = TRUE,
           verbose = TRUE) standardGeneric("computeDistance")
)

#' @rdname computeDistance
#' @aliases computeDistance,CoPro-method
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @importFrom fields rdist
#' @export
setMethod(
  "computeDistance", "CoPro",
  function(object, distType = c(
             "Euclidean2D", "Euclidean3D",
             "Morphology-Aware"
           ),
           xDistScale = 1, yDistScale = 1, zDistScale = 1,
           normalizeDistance = TRUE, truncateLowDist = TRUE,
           verbose = TRUE) {
    ## match arg
    distType <- match.arg(distType)

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

    ## check dist
    if (distType == "Euclidean3D") {
      if (!all(c("x", "y", "z") %in% object@locationDataSub)) {
        stop(paste(
          "please make sure x, y, z are all available to run",
          "3D Euclidean distance calcuation"
        ))
      }
    }

    ## distances are between any two cell types
    distances <- setNames(rep(list(), length = length(cts)), cts)
    for (i in cts) {
      distances[[i]] <- setNames(rep(list(), length = length(cts)), cts)
    }

    pair_cell_types <- combn(cts, 2)
    ct_ind_sub <- object@cellTypesSub

    ## notify users if normalizeDistance = TRUE
    if (normalizeDistance) {
      cat("normalizeDistance is set to TRUE, so distance will be",
          "normalized, so that 0.01 percentile distance will be scaled",
          "to 0.01\n")
    }
    dist_1percentile <- vector(mode = "numeric",
                               length = ncol(pair_cell_types))

    ## calculate the distances
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
          "for potential errors"
        ))
        distances_ij[distances_ij == 0] <-
          min(distances_ij[distances_ij != 0])
      }

      dist_1percentile[pp] <- quantile(distances_ij[distances_ij != 0], 1e-4)

      if (truncateLowDist) {
        distances_ij[distances_ij < dist_1percentile[pp]] <-
          dist_1percentile[pp]
      }

      ## save the distances
      distances[[i]][[j]] <- distances_ij
      if (verbose) {
        cat("quantile of the distances between", i, "and", j, "is: \n")
        print(quantile(distances_ij))
      }
    }

    min_1percentile <- min(dist_1percentile)

    if (normalizeDistance) {
      cat("The scaling factor for normalizing distance is",
          0.01 / min_1percentile, "\n")
      for (pp in seq_len(ncol(pair_cell_types))) {
        i <- pair_cell_types[1, pp]
        j <- pair_cell_types[2, pp]

        distances_ij <- distances[[i]][[j]]
        distances_ij <- distances_ij / min_1percentile * 0.01
        distances[[i]][[j]] <- distances_ij
      }
    }


    object@distances <- distances
    return(object)
  }
)


.CheckSigmaValuesToRemove <- function(kernel_current, lowerLimit,
                                      sigma_choose, sigmaValues, i, j,
                                      minAveCellNeighor) {
  n_cell1 <- nrow(kernel_current)
  n_cell2 <- ncol(kernel_current)

  minPropZero <- minAveCellNeighor * min(n_cell1, n_cell2) / (n_cell1 * n_cell2)
  if (mean(kernel_current <= lowerLimit) < minPropZero) {
    warning(paste("Kernel matrix for cell types", i, "and", j,
                  "with sigma =", sigma_choose,
                  "contains almost all zeros. Specifically, more than",
                  minPropZero*100,"% total counts are zero" ))
    if (length(sigmaValues) == 1) {
      stop(paste("Only one sigma value is specified,",
                 "which resulted in all Gaussian kernel being small.",
                 "Please provide a larger sigma value"))
    }else {
      warning(paste("Dropping sigma value of ",
                    sigma_choose,
                    "because all Gaussian kernel values are too small,",
                    "which will not produce meaningful results."))
      return(TRUE)
    }
  }else if (all(is.na(kernel_current))) {
    warning(paste("Kernel matrix for cell types", i, "and", j,
                  "with sigma =", sigma_choose,
                  "contains all NA."))
    if (length(sigmaValues) == 1) {
      stop(paste("Only one sigma value is specified,",
                 "which resulted in all Gaussian kernel being NA."))
    }else {
      warning(paste("Dropping sigma value of ",
                    sigma_choose,
                    "because all Gaussian kernel values are NA."))
      return(TRUE)
    }
  }

  return(FALSE)
}


#' Compute Kernel Matrix for CoPro
#'
#' This method calculates the kernel matrices for pairs of cell types based on
#' their distances and a range of sigma values.
#' The formula of calculating kernel matrix is:
#' \deqn{K(x, y) = \exp\left(-\frac{\|x-y\|^2}{2 \sigma^2}\right)}
#' The matrices are adjusted by clipping the upper quantile of
#'  the values to reduce the effect of outliers. The results are stored
#'  within the object.
#'
#' @importFrom utils combn
#' @param object A `CoPro` object.
#' @param sigmaValues A vector of sigma values used for kernel calculation.
#' @param lowerLimit The lower limit for the kernel function, default is 1e-7.
#' @param upperQuantile The quantile used for clipping the kernel values,
#' default is 0.85.
#' @param verbose Whether to output the progress and related information
#' @param normalizeKernel Whether to normalize the kernel matrix?
#' Default = FALSE. Note that normalization will not affect any downstream
#' analyses, it is for numerical stability and easier interpretation only.
#' @param minAveCellNeighor What is the minimum average number of cell in the
#'  neighbor? This step is to help set up the expected sparsity of the
#'  kernel matrix. If a kernel sigma value is too small, this result in too
#'  few neighbors for most cells, resulting in an overly-sparse matrix that
#'  makes the parameter estimation hard. Thus, the sigma values that results in
#'  an overly-sparse matrix will be removed for later analysis.
#' @param rowNormalizeKernel Whether the kernel matrix will be row-wise
#' normalized? Note that row or column wise normalization will result in an
#' asymmetric result in skrCCA inference.
#' @param colNormalizeKernel Whether the kernel matrix will be column-wise
#' normalized? Note that row or column wise normalization will result in an
#' asymmetric result in skrCCA inference.
#' @return The `CoPro` object with computed kernel matrices added. The kernel
#' matrices are organized into a three-layer nested list object. The first layer
#' is indexed by the sigma value, and the second and the third layers are cell
#' types
#' @export
#' @note To-do: Shall we include row or column normalization of the kernel?
setGeneric(
  "computeKernelMatrix",
  function(object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
           normalizeKernel = FALSE, minAveCellNeighor = 2,
           rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
           verbose = TRUE) standardGeneric("computeKernelMatrix"))

#' @rdname computeKernelMatrix
#' @aliases computeKernelMatrix,CoPro-method
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @export
setMethod(
  "computeKernelMatrix", "CoPro",
  function(object, sigmaValues,
           lowerLimit = 1e-7, upperQuantile = 0.85, normalizeKernel = FALSE,
           minAveCellNeighor = 2,
           rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
           verbose = TRUE) {
    ## make sure distance matrix exist
    if (length(object@distances) == 0) {
      stop("Please run computeDistance before computing kernel")
    }

    if (rowNormalizeKernel && colNormalizeKernel){
      stop("Cannot do both row-wise and column-wise normalization.")
    }

    cts <- object@cellTypesOfInterest
    n_mat <- length(cts)
    if (n_mat < 2) {
      stop("At least two cell types are needed to compute kernel matrices.")
    }

    if (length(sigmaValues) == 0) {
      warning("No Sigma specified, setting to the 5% quantile of cell distance")
      dist12 <- object@distances[[1]][[2]]
      sigmaValues <- quantile(dist12[dist12 > 0], 0.05)
    }

    if (!is.numeric(sigmaValues)) {
      stop("sigmaValues must be numeric values or a vector")
    }

    ## notify users if normalizeDistance = TRUE
    if (normalizeKernel) {
      cat("normalizeKernel is set to TRUE. Kernel matrix will be",
          "normalized so that median row sums of kernel will be",
          "1 \n")
    }

    ## save the sigmaValues
    object@sigmaValues <- sigmaValues

    ## Initialize the list of kernel matrices
    kernel_mat <- vector("list", length(sigmaValues))
    sigma_names <- paste("sigma", sigmaValues, sep = "_")
    names(kernel_mat) <- sigma_names
    for (t in sigma_names) {
      kernel_mat[[t]] <- setNames(vector("list", n_mat), cts)
      for (i in cts) {
        kernel_mat[[t]][[i]] <- setNames(vector("list", n_mat), cts)
      }
    }

    pair_cell_types <- combn(cts, 2)

    ## sigma values to leave out
    sigmaValuesToRemove <- vector(mode = "logical",
                                  length = length(sigmaValues))
    names(sigmaValuesToRemove) <- sigma_names

    for (tt in seq_along(sigmaValues)) {
      t <- sigma_names[tt]
      sigma_choose <- sigmaValues[tt]
      cat("current sigma value is\n", sigma_choose)

      for (pp in seq_len(ncol(pair_cell_types))) {
        i <- pair_cell_types[1, pp]
        j <- pair_cell_types[2, pp]

        kernel_current <- kernel_from_distance(
          sigma = sigma_choose,
          dist_mat = object@distances[[i]][[j]],
          lower_limit = lowerLimit
        )

        sigmaValuesToRemove[t] <- .CheckSigmaValuesToRemove(
          kernel_current = kernel_current, lowerLimit = lowerLimit,
          minAveCellNeighor = minAveCellNeighor, sigma_choose = sigma_choose,
          sigmaValues = sigmaValues, i = i, j = j)

        if (sigmaValuesToRemove[t]) {
          kernel_mat[[t]][[i]][[j]] <- kernel_current
          next
        }

        ## Clipping large values
        upper_clip <- quantile(kernel_current[kernel_current >= lowerLimit],
                               upperQuantile)
        kernel_current[kernel_current >= upper_clip] <- upper_clip


        if ((normalizeKernel && !rowNormalizeKernel) && !colNormalizeKernel) {
          ## calculate row sum
          rs_kernel <- rowSums(kernel_current)
          kernel_current <- kernel_current / median(rs_kernel[rs_kernel > 1e-5])

        }else if (rowNormalizeKernel) {
          ## calculate row sum
          rs_kernel <- rowSums(kernel_current)
          cat("quantile of kernel matrix rowSums \n")
          cat(quantile(rs_kernel))
          cat("\n")
          nz_ind <- rs_kernel > 1e-4
          kernel_current[nz_ind,] <- kernel_current[nz_ind,] / rs_kernel[nz_ind]
        }else if (colNormalizeKernel) {
          kernel_current <- t(kernel_current)
          ## calculate col sum
          rs_kernel <- rowSums(kernel_current)
          cat("quantile of kernel matrix colSums \n")
          cat(quantile(rs_kernel))
          cat("\n")
          nz_ind <- rs_kernel > 1e-4
          kernel_current[nz_ind,] <- kernel_current[nz_ind,] / rs_kernel[nz_ind]
          kernel_current <- t(kernel_current)
        }

        ## print info
        if (verbose) {
          cat("Current Sigma value is", sigma_choose)
          cat("\n")
          cat("Quantiles of N_neighbors for cell type", i, "\n")
          cat(quantile(rowSums(kernel_current >= lowerLimit)))
          cat("\n")
          cat("Quantiles of N_neighbors for cell type", j, "\n")
          cat(quantile(colSums(kernel_current >= lowerLimit)))
          cat("\n")
        }

        # ## remove small values
        kernel_current[kernel_current < lowerLimit] <- 0
        kernel_mat[[t]][[i]][[j]] <- kernel_current
      }
    }

    ## Remove kernel matrices and sigmaValues that were marked for removal
    if (any(sigmaValuesToRemove)) {
      cat("removing", sum(sigmaValuesToRemove), "sigmaValues values", "\n")
      kernel_mat <- kernel_mat[!sigmaValuesToRemove]
      object@sigmaValues <- object@sigmaValues[!sigmaValuesToRemove]
    }

    object@kernelMatrices <- kernel_mat
    return(object)
  }
)

## get PC matrices
.getAllPCMats <- function(allPCs, scalePCs) {

  if (length(allPCs) == 0) {
    stop("PCA results do not exist, run computePCA() first.")
  }

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
  return(PCmats)
}

#' runSkrCCA
#' @importFrom stats setNames
#' @param object A CoPro object
#' @param scalePCs Whether to scale each PCs to a uniform variance before
#' running the program
#' @param nCC Number of canonical vectors to compute, default = 2
#' @param tol Tolerance for termination, default = 1e-5
#' @param maxIter Maximum iterations
#'
#' @return CoPro object with distnace matrix computed
#' @export
#'
setGeneric(
  "runSkrCCA",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           maxIter = 200) standardGeneric("runSkrCCA"))


#' @rdname runSkrCCA
#' @aliases runSkrCCA,CoPro-method
#' @importFrom stats setNames
#' @export
setMethod(
  "runSkrCCA", "CoPro",
  function(object,
           scalePCs = TRUE, nCC = 2, tol = 1e-5,
           maxIter = 200) {
    ## check whether the kernel matrix is available
    if (length(object@kernelMatrices) == 0) {
      stop("Kernel matrix is empty, please run computeKernelMatrix first")
    }

    ## record whether each PC will been rescaled
    if (length(object@scalePCs) == 0) {
      object@scalePCs <- scalePCs
    }

    ## check sigmaValues
    if (length(object@sigmaValues) == 0) {
      stop("sigmaValues is empty, please specify")
    } else {
      sigmaValues <- object@sigmaValues
    }

    ## choose cell types
    if (length(object@cellTypesOfInterest) != 0) {
      cts <- object@cellTypesOfInterest
    } else {
      warning("no cell type of interest specified,
                      using all cell types to run the analysis")
      cts <- unique(object@cellTypesSub)
    }

    PCmats <- .getAllPCMats(allPCs = object@pcaResults, scalePCs = scalePCs)

    ## run across different sigma values
    cca_out <- vector("list", length = length(sigmaValues))
    sigma_names <- paste("sigma", sigmaValues, sep = "_")
    names(cca_out) <- sigma_names

    ## for loop to run the analysis
    for (tt in seq_along(sigmaValues)) {
      t <- sigma_names[tt]
      cca_result <- optimize_bilinear_multi(
        X_list = PCmats,
        K_list = object@kernelMatrices[[t]],
        max_iter = maxIter, tol = tol
      )
      names(cca_result) <- cts

      if (nCC == 1) {
        cca_out[[t]] <- cca_result
      }else {
        cca_result_n <- optimize_bilinear_multi_n(
          X_list = PCmats, K_list = object@kernelMatrices[[t]],
          w_list = cca_result,
          cellTypesOfInterest = cts, nCC = nCC,
          max_iter = maxIter, tol = tol
        )
        cca_out[[t]] <- cca_result_n
      }

    }

    object@skrCCAOut <- cca_out
    object@nCC <- nCC
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

    PCmats <- .getAllPCMats(allPCs = object@pcaResults, scalePCs = scalePCs)

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

    PCmats <- .getAllPCMats(allPCs = object@pcaResults, scalePCs = scalePCs)

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
              matrix(w_1 * object@pcaResults[[i]]$sdev, nrow = 1) %*%
                t(object@pcaResults[[i]]$rotation)
            )
          } else {
            geneScores[[t]][[i]][, cc_name] <- as.vector(
              matrix(w_1, nrow = 1) %*%
                t(object@pcaResults[[i]]$rotation)
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
            if(!is.null(object@pcaResults)) cat("- PCA\n")
            if(!is.null(object@skrCCAOut)) cat("- skrCCA\n")

            # Additional information
            if(length(object@metaData) > 0) {
              cat("\nAvailable metadata fields:\n")
              cat(paste("-", names(object@metaData), collapse = "\n"))
              cat("\n")
            }

            invisible(x = object)
          }
)
