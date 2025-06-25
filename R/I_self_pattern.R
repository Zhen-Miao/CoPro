


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

