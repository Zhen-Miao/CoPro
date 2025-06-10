


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

  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)

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
  if(length(cts) != 1){
    stop("Only run this function when there is exactly one cell type")
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
            matrix(w_1 * object@pcaGlobal[[cts]]$sdev, nrow = 1) %*%
              t(object@pcaGlobal[[cts]]$rotation)
          )
        } else {
          geneScores[[t]][[cts]][, cc_name] <- as.vector(
            matrix(w_1, nrow = 1) %*%
              t(object@pcaGlobal[[cts]]$rotation)
          )
        }
        cellScores[[t]][[cts]][, cc_name] <- as.vector(PCmats[[cts]] %*% w_1)
      }


  }

  ## save cellscores and gene scores
  object@cellScores <- cellScores
  object@geneScores <- geneScores

  ## add cell score information to the cell metadata

  meta_t <- object@metaDataSub[object@cellTypesSub == cts, ]
  for (t in sigma_names){
      for (cc_index in seq_len(nCC)) {
        cc_name <- paste0("CC_", cc_index)
        cellScoreColName <- paste0("cellScore_", t, "_cc_index_", cc_index)

        meta_t[, cellScoreColName] <-
          cellScores[[t]][[cts]][rownames(meta_t), cc_name]
      }
    }

  ## combine each cell type meta.data back
  meta_all <- meta_t
  object@metaDataSub <- meta_all

  return(object)
}

