

#' subsetDataOne
#'
#' @param object A `CoPro` object
#' @param cellTypesOfInterest Input cell types of interest as a vector of
#' characters for subsetting the data
#'
#' @returns The CoPro object
#' @export
subsetDataOne <- function(object, cellTypesOfInterest) {

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
}

#' computeDistance for one cell type
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
computeDistanceOne <- function(
    object,
    distType = c("Euclidean2D", "Euclidean3D","Morphology-Aware"),
    xDistScale = 1, yDistScale = 1, zDistScale = 1,
    normalizeDistance = TRUE, truncateLowDist = TRUE,
    verbose = TRUE){
  ## check input
  if (!is(object, "CoPro")) {
    stop("Input object must be a CoPro object")
  }

  ## match arg
  distType <- match.arg(distType)

  ## make sure only one cell type is present
  if (length(object@cellTypesOfInterest) == 1) {
    cts <- object@cellTypesOfInterest
  } else {
    if (length(unique(object@cellTypesSub)) == 1) {
      cts <- unique(object@cellTypesSub)
    } else {
      stop(paste("Only run this function when there is one cell type.",
                 "Please run computeDistance() instead."))
    }
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

  ## distances with itself
  distances <- setNames(rep(list(), length = length(cts)), cts)
  for (i in cts) {
    distances[[i]] <- setNames(rep(list(), length = length(cts)), cts)
  }


  ct_ind_sub <- object@cellTypesSub

  ## notify users if normalizeDistance = TRUE
  if (normalizeDistance) {
    cat("normalizeDistance is set to TRUE, so distance will be",
        "normalized, so that 0.01 percentile distance will be scaled",
        "to 0.01\n")
  }
  dist_1percentile <- vector(mode = "numeric",
                             length = 1)

  if (distType == "Euclidean2D") {
    mat1 <- cbind(
      object@locationDataSub$x[ct_ind_sub == cts] * xDistScale,
      object@locationDataSub$y[ct_ind_sub == cts] * yDistScale
    )

  } else if (distType == "Euclidean3D") {
    mat1 <- cbind(
      object@locationDataSub$x[ct_ind_sub == cts] * xDistScale,
      object@locationDataSub$y[ct_ind_sub == cts] * yDistScale,
      object@locationDataSub$z[ct_ind_sub == cts] * zDistScale
    )
  } else if (distType == "Morphology-Aware") {
    stop("morphology-aware kernel is not availabe at this moment")
  }

  ## compute distance
  distances_ij <- fields::rdist(mat1)
  diag(distances_ij) <- Inf

  if (any(distances_ij == 0)) {
    warning(paste("Zero distances detected, replacing with",
                  "the smallest non-zero distances, please",
                  "consider checking the location of cells",
                  "for potential errors"
    ))
    distances_ij[distances_ij == 0] <-
      min(distances_ij[distances_ij != 0])
  }

  dist_1percentile[1] <- quantile(distances_ij[distances_ij != 0], 1e-4)

  if (truncateLowDist) {
    distances_ij[distances_ij < dist_1percentile[1]] <-
      dist_1percentile[1]
  }

  ## save the distances
  distances[[cts]][[cts]] <- distances_ij
  if (verbose) {
    cat("quantile of the distances between", cts, "and", cts, "is: \n")
    print(quantile(distances_ij))
  }


  min_1percentile <- min(dist_1percentile)

  if (normalizeDistance) {
    cat("The scaling factor for normalizing distance is",
        0.01 / min_1percentile, "\n")

    distances_ij <- distances_ij / min_1percentile * 0.01
    distances[[cts]][[cts]] <- distances_ij
  }

  object@distances <- distances
  return(object)

}

#' Compute Kernel Matrix for one cell type
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
#' @return The `CoPro` object with computed kernel matrices added. The kernel
#' matrices are organized into a three-layer nested list object. The first layer
#' is indexed by the sigma value, and the second and the third layers are cell
#' types
#' @export
#' @note To-do: Shall we include row or column normalization of the kernel?
computeKernelMatrixOne <- function(object, sigmaValues,
         lowerLimit = 1e-7, upperQuantile = 0.85, normalizeKernel = FALSE,
         minAveCellNeighor = 2, verbose = TRUE) {
  ## make sure distance matrix exist
  if (length(object@distances) == 0) {
    stop("Please run computeDistance before computing kernel")
  }

  cts <- object@cellTypesOfInterest
  n_mat <- length(cts)

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

  ## sigma values to leave out
  sigmaValuesToRemove <- vector(mode = "logical",
                                length = length(sigmaValues))
  names(sigmaValuesToRemove) <- sigma_names

  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    sigma_choose <- sigmaValues[tt]


    kernel_current <- kernel_from_distance(
        sigma = sigma_choose,
        dist_mat = object@distances[[cts]][[cts]],
        lower_limit = lowerLimit)


    sigmaValuesToRemove[t] <- .CheckSigmaValuesToRemove(
        kernel_current = kernel_current, lowerLimit = lowerLimit,
        minAveCellNeighor = minAveCellNeighor, sigma_choose = sigma_choose,
        sigmaValues = sigmaValues, i = cts, j = cts)

    if (sigmaValuesToRemove[t]) {
        kernel_mat[[t]][[cts]][[cts]] <- kernel_current
        next
      }

      ## Clipping large values
      upper_clip <- quantile(kernel_current[kernel_current >= lowerLimit],
                             upperQuantile)
      kernel_current[kernel_current >= upper_clip] <- upper_clip

      if (normalizeKernel) {
        rs_kernel <- rowSums(kernel_current)
        kernel_current <- kernel_current / median(rs_kernel[rs_kernel != 0])
      }

      ## print info
      if (verbose) {
        cat("Current Sigma value is", sigma_choose)
        cat("\n")
        cat("Quantiles of N_neighbors for cell type", cts, "\n")
        cat(quantile(rowSums(kernel_current >= lowerLimit)))
        cat("\n")
        cat("Quantiles of N_neighbors for cell type", cts, "\n")
        cat(quantile(colSums(kernel_current >= lowerLimit)))
        cat("\n")
      }

      # ## remove small values
      kernel_current[kernel_current < lowerLimit] <- 0
      kernel_mat[[t]][[cts]][[cts]] <- kernel_current

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


#' SkrCCA optimization function for multiple groups
#'
#' @param X_list List of data matrices (cell by PC matrix)
#' @param K_list Kernel matrices between any two pairs of matrices
#' @param max_iter Maximum number of iterations
#' @param tol tolerance of accuracy
#'
#' @return `w_list`, a list of weight vectors for each cell type
#' @noRd
optimize_bilinear_multi_One <- function(X_list, K_list, max_iter = 1000,
                                    tol = 1e-5) {
  n_mat <- 1
  n_features <- ncol(X_list[[1]])

  # Initialize w1 by its right singular vector
  w_list <- list()
  w_list[[1]] <- svd(X_list[[1]])$v[, 1, drop = FALSE]


  # Iterative refinement
  iter <- 0
  diff_i <- vector(length = n_mat)
  while (iter <= max_iter) {
    ## after going over all j, we will update i
      ## initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = n_mat, nrow = n_features)

      K11 <- K_list[[1]][[1]]
      ## compute Y_ij
      X1 <- X_list[[1]]
      Y <- t(X1) %*% K11 %*% X1
      w2 <- w_list[[1]]

      w_i_left[, 1] <- Y %*% w2

      w_i_new <- normalize_vec(rowSums(w_i_left))
      diff_i[1] <- max(abs(w_i_new - w_list[[1]]))
      w_list[[1]] <- w_i_new

    if (mean(diff_i) <= tol) {
      print(paste("convergence reached at", iter, "iterations", sep = " "))
      break
    }
    iter <- iter + 1
  }

  if (iter > max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  ## convert w_list into matrix
  w_list[[1]] <- matrix(w_list[[1]], ncol = 1)

  ## return results
  return(w_list)
}




bilinear_w_from_Y_resi_One <- function(w_list_new, Y_resi,
                                    n_features, max_iter, tol) {
  # Iterative refinement
  iter <- 0
  while (iter < max_iter) {
    ## after going over all j, we will update i

      ## initialize the w_i_left matrix
      w_i_left <- matrix(data = 0, ncol = 1, nrow = n_features)
      diff_i <- vector(length = 1)

      w2 <- w_list_new[[1]]
      Y <- Y_resi[[1]][[1]] ## make sure all Y_resi exist except for i == j
      w_i_left[, 1] <- Y %*% w2

      w_i_new <- normalize_vec(rowSums(w_i_left))
      diff_i[1] <- max(abs(w_i_new - w_list_new[[1]]))
      w_list_new[[1]] <- w_i_new

    if (mean(diff_i) <= tol) {
      print(paste("convergence reached at", iter, "iterations", sep = " "))
      break
    }
    iter <- iter + 1
  }

  if (iter == max_iter) {
    warning("Maximum number of iterations reached without convergence")
  }

  return(w_list_new)
}



#' Run multi version of skrCCA to detect the third to n_th component
#'
#' @importFrom stats setNames
#'
#' @param X_list List of data matrices (cell by PC matrix)
#' @param K_list Kernel matrices between any two pairs of matrices
#' @param w_list A list of weights for each data matrix, pre-computed and to
#'  be regressed out
#' @param cellTypesOfInterest A vector specifying cell types of interest
#' @param nCC number of canonical vectors to compute, need to be greater than
#'  two, default = 2.
#' @param max_iter Maximum number of iterations
#' @param tol tolerance of accuracy
#'
#' @return A list of weights
#' @noRd
optimize_bilinear_multi_n_One <- function(X_list, K_list, w_list,
                                      cellTypesOfInterest,
                                      nCC = 2,
                                      max_iter = 1000,
                                      tol = 1e-5) {

  n_mat <- 1
  n_features <- ncol(X_list[[1]])
  cts <- cellTypesOfInterest

  ## check if w_list has the same length of x_list
  if (n_mat != length(w_list) || n_mat != length(cellTypesOfInterest) ) {
    stop(paste("the input X_list, w_list, cellTypesOfInterest",
               "are of different length!"))
  }

  ## check the dimension of the w_list
  if (length(dim(w_list[[1]])) == 0) {
    stop("the input w_list must be a matrix!")
  }

  ## get all pairs for X, K, W to obtain Y_resi
  ## let Y_resi be the same structure as K
  Y_resi <- setNames(vector(mode = "list", length = n_mat), cts)
  Y_resi[[cts]] <- setNames(vector(mode = "list", length = n_mat), cts)

  ## Y is between any two cell types, and each time, we update the Y_resi
  ## to regress out all previous canonical vectors

  for (qq in 1:(nCC - 1)) { ## qq: existing w_qq to be regressed out

    ## step 1: obtain or update Y_resi
      K11 <- K_list[[cts]][[cts]]

      ## compute Y_ij
      X1 <- X_list[[cts]]
      w1 <- w_list[[1]][, qq, drop = FALSE]

      if (qq == 1) {
        Y1 <- t(X1) %*% K11 %*% X1
      }else {
        Y1 <- Y_resi[[cts]][[cts]]
      }

      Y_resi[[cts]][[cts]] <- Y1 - ((t(w1) %*% Y1 %*% w1)[1, 1] * (w1 %*% t(w1)))

    ## step 2: initialize w_list_new by SVD
    w_list_new <- rep(list(), length = n_mat)
    w_list_new[[1]] <- svd(t(Y_resi[[cts]][[cts]]))$v[, 1, drop = FALSE]

    ## step 3: Iterative refinement
    w_list_qq <- bilinear_w_from_Y_resi_One(w_list_new = w_list_new,
                                        Y_resi = Y_resi,
                                        n_features = n_features,
                                        max_iter = max_iter, tol = tol)

    ## step 4: add w_list_qq to columns in w_list
    w_list[[1]] <- cbind(w_list[[1]], w_list_qq[[1]])

  }

  ## return results
  return(w_list)
}

#' runSkrCCA for one cell type
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
runSkrCCAOne <- function(object,
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
    cca_result <- optimize_bilinear_multi_One(
      X_list = PCmats,
      K_list = object@kernelMatrices[[t]],
      max_iter = maxIter, tol = tol
    )
    names(cca_result) <- cts

    if (nCC == 1) {
      cca_out[[t]] <- cca_result
    }else {
      cca_result_n <- optimize_bilinear_multi_n_One(
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

#' Compute Normalized Correlation for one type
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
computeNormalizedCorrelationOne <- function(object, tol = 1e-4) {
  ## Check for required components
  if (length(object@skrCCAOut) == 0) {
    stop("CCA results are not available. Please run CCA first.")
  }

  ## choose cell types
  cts <- object@cellTypesOfInterest

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

  # pair_cell_types <- combn(cts, 2)

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
    norm_K12[[t]] <- setNames(vector(mode = "list", length = 1), cts)
    norm_K12[[t]][[cts]] <- setNames(vector(mode = "list", length = 1), cts)
  }

  for (t in sigma_names) {
      K <- object@kernelMatrices[[t]][[cts]][[cts]]
      ## Calculate the spectral norm of the kernel matrix
      svd_result <- irlba::irlba(K, nv = 1, tol = tol)
      norm_K12[[t]][[cts]][[cts]] <- svd_result$d[1]
  }

  cat("Finished calculating spectral norms \n")


  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    correlation_value[[t]] <- data.frame(
      sigmaValues = sigmaValues[tt],
      cellType1 = rep(cts, times = nCC),
      cellType2 = rep(cts, times = nCC),
      CC_index = rep(x = 1:nCC, each = 1),
      normalizedCorrelation = numeric(length = nCC),
      stringsAsFactors = FALSE
    )

  for (cc_index in seq_len(nCC)) {

    A <- PCmats[[cts]]
    w_1 <- object@skrCCAOut[[t]][[cts]][, cc_index, drop = FALSE]

    A_w1 <- A %*% w_1

    ## get pre-calculated spectral norm
    K <- object@kernelMatrices[[t]][[cts]][[cts]]
    norm_K12_sel <- norm_K12[[t]][[cts]][[cts]]

    ## Calculate normalized correlation
    correlation_value[[t]]$"normalizedCorrelation"[cc_index] <-
          (t(A_w1) %*% K %*% A_w1) /
          (sqrt(sum(A_w1^2)) * sqrt(sum(A_w1^2)) * norm_K12_sel)
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

#' computeGeneAndCellScores for one cell type
#' @importFrom stats setNames
#' @param object A `CoPro` object containing CCA results
#' and kernel matrices.
#'
#' @return A `CoPro` object with gene and cell score computed
#' @export
#'
computeGeneAndCellScoresOne <- function(object) {
  ## Check for required components
  if (length(object@skrCCAOut) == 0) {
    stop("CCA results are not available. Please run CCA first.")
  }
  if (length(object@kernelMatrices) == 0) {
    stop(paste("Kernel matrices are not available.",
               "Please compute the kernel matrices first."))
  }

  ## choose cell types
  cts <- object@cellTypesOfInterest

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
      cellScores[[t]][[cts]] <- matrix(
        nrow = sum(object@cellTypesSub == cts),
        ncol = nCC
      )
      colnames(cellScores[[t]][[cts]]) <- paste0("CC_", 1:nCC)
      rownames(cellScores[[t]][[cts]]) <- rownames(object@normalizedDataSub)[
        object@cellTypesSub == cts
      ]
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

      geneScores[[t]][[cts]] <- matrix(
        nrow = ncol(object@normalizedDataSub),
        ncol = nCC
      )
      colnames(geneScores[[t]][[cts]]) <- paste0("CC_", 1:nCC)
      rownames(geneScores[[t]][[cts]]) <- colnames(object@normalizedDataSub)

  }


  ## go over all cell types, then over all sigma values
  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
      for (cc_index in seq_len(nCC)) {
        cc_name <- paste0("CC_", cc_index)
        w_1 <- object@skrCCAOut[[t]][[cts]][, cc_index, drop = FALSE]
        if (scalePCs) {
          geneScores[[t]][[cts]][, cc_name] <- as.vector(
            matrix(w_1 * object@pcaResults[[cts]]$sdev, nrow = 1) %*%
              t(object@pcaResults[[cts]]$rotation)
          )
        } else {
          geneScores[[t]][[cts]][, cc_name] <- as.vector(
            matrix(w_1, nrow = 1) %*%
              t(object@pcaResults[[cts]]$rotation)
          )
        }
        cellScores[[t]][[cts]][, cc_name] <- as.vector(PCmats[[cts]] %*% w_1)
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

  meta_t <- object@metaDataSub[object@cellTypesSub == cts, ]
  for (t in sigma_names){
      for (cc_index in seq_len(nCC)) {
        cc_name <- paste0("CC_", cc_index)
        cellScoreColName <- paste0("cellScore_", t, "_cc_index_", cc_index)

        meta_t[, cellScoreColName] <-
          cellScores[[t]][[cts]][rownames(meta_t[[cts]]), cc_name]
      }
    }

  ## combine each cell type meta.data back
  meta_all <- meta_t
  object@metaDataSub <- meta_all

  return(object)
}


#' Retrieve the Correlation within one cell types
#'
#' @param object A `CoPro` object
#' @param cellTypeA Cell type label for the cell type of interest
#' @param sigmaValueChoice A particular sigma squared value
#' for the correlation
#' @param ccIndex Canonical vector index, default = 1
#'
#' @return A data.frame with two columns, AK and B, where AK represents the
#' cell score of cell type A times the kernel matrix, and B represents the
#' cell score of cell type B.
#' @export
getCorrOneType <- function(object, cellTypeA, ccIndex = 1,
                            sigmaValueChoice) {
  ## check input
  if (length(cellTypeA) != 1) {
    stop("Must give a single cellTypeA for correlation plot")
  }

  ## choose cell types

  cts <- object@cellTypesOfInterest
  if (cellTypeA != cts){
    stop("cellTypeA must be in cellTypesSub.")
  }


  ## set sigmaValueChoice
  if (is.null(sigmaValueChoice)) {
    if (length(object@sigmaValueChoice) == 0) {
      stop(paste(
        "sigmaValueChoice is not given,",
        "and NormalizedCorrelation not computed,",
        "please either specify a particular sigmaValueChoice or",
        "run computeNormalizedCorrelation()"
      ))
    }else {
      warning(paste(
        "sigmaValueChoice is not given",
        "default set to the value with highest",
        "normalized correlation."
      ))
      sigmaValueChoice <- object@sigmaValueChoice
    }
  }

  if (!(sigmaValueChoice %in% object@sigmaValues)) {
    stop("sigmaValueChoice does not exist in the list of sigmaValues")
  }

  ## make sure normalizedCorrelation exists
  if (length(object@cellScores) == 0 ||
      length(object@geneScores) == 0) {
    stop(paste(
      "cellScores slot does not exist,",
      "run computeGeneAndCellScores first"
    ))
  }

  ## load the cellScores and kernel matrix
  sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")

  x1 <- t(object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = FALSE])
  x2 <- object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = TRUE]
  ktemp <- object@kernelMatrices[[sigma_name]][[cellTypeA]][[cellTypeA]]

  df <- data.frame(AK = (x1 %*% ktemp)[1, , drop = TRUE], B = x2)
  return(df)
}





