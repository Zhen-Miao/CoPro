#' computeGeneAndCellScores
#' @importFrom stats setNames
#' @param object A `CoPro` or `CoProMulti` object containing CCA results
#' and kernel matrices.
#'
#' @return A `CoPro` or `CoProMulti` object with gene and cell scores computed
#' @export
#'
setGeneric(
  "computeGeneAndCellScores",
  function(object) standardGeneric("computeGeneAndCellScores")
)

.checkInputGAC <- function(object) {
    ## determine if the object is a CoPro object or CoProMulti object
    if(is(object, "CoPro")) {
      is_multi <- FALSE
    } else if(is(object, "CoProMulti")) {
      is_multi <- TRUE
    } else {
      stop("object is not a CoPro or CoProMulti object")
    }
    if(is_multi) {
      slides <- object@slideList
    } else {
      slides <- NULL
    }

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

    ## check pcaGlobal
    if(length(object@pcaGlobal) == 0) {
      stop("pcaGlobal is empty, please specify")
    }

    ## check nCC
    if(length(object@nCC) == 0) {
      stop("nCC is empty, please specify")
    }
    

    ## load whether the PCs are being scaled prior to CCA

    if (length(object@scalePCs) == 0) {
      stop("object@scalePCs not specified")
    }
    scalePCs <- object@scalePCs

    return(list(cts = cts, sigmaValues = sigmaValues, scalePCs = scalePCs, 
                slides = slides, is_multi = is_multi))

}

.initializeCSGS <- function(cts, sigma_names, nCC, 
          cellTypesSub, cellNamesSub, geneNamesSub) {
    ## cell scores and gene scores are both by cell types
    cellScores <- setNames(
      vector(mode = "list", length = length(sigma_names)),
      sigma_names
    )
    for (t in sigma_names) {
      cellScores[[t]] <- setNames(
        vector(mode = "list", length = length(cts)),
        cts
      )
    }

    ## Initialize geneScores with same structure as cellScores (different dimensions will be set below)
    geneScores <- cellScores

    for (t in sigma_names) {
      for (i in cts) {
        cellScores[[t]][[i]] <- .createScoreMatrix(
          nrows = sum(cellTypesSub == i),
          ncols = nCC,
          row_names = cellNamesSub[cellTypesSub == i]
        )

        geneScores[[t]][[i]] <- .createScoreMatrix(
          nrows = length(geneNamesSub),
          ncols = nCC,
          row_names = geneNamesSub
        )
      }

    }
    return(list(cellScores = cellScores, geneScores = geneScores))
}

.initializeCSGSMulti <- function(cts, sigma_names, nCC, slides) {
    ## Initialize structure for multi-slide data
    cellScores <- setNames(
      vector(mode = "list", length = length(sigma_names)),
      sigma_names
    )
    
    for (t in sigma_names) {
      cellScores[[t]] <- setNames(
        vector(mode = "list", length = length(slides)),
        slides
      )
      for (sID in slides) {
        cellScores[[t]][[sID]] <- setNames(
          vector(mode = "list", length = length(cts)),
          cts
        )
      }
    }

    geneScores <- setNames(
      vector(mode = "list", length = length(sigma_names)),
      sigma_names
    )
    for (t in sigma_names) {
      geneScores[[t]] <- setNames(
        vector(mode = "list", length = length(cts)),
        cts
      )
    }

    return(list(cellScores = cellScores, geneScores = geneScores))
}

.CSToMeta <- function(object, cellScores, cts, sigma_names, nCC) {
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
    meta_all <- do.call(rbind, meta_t)
    cell_names <- rownames(object@metaDataSub)
    row_names <- rownames(meta_all)

    ## extensive checking on the rownames of meta_all
    ## as this is a very common error
    if(is.null(row_names)) {
      stop("rownames of meta_all is null, code error")
    } else if (!any(row_names %in% cell_names)) {
      warning("rownames of meta_all are not in the rownames of metaDataSub, ",
              "some cell scores will be lost")
    } else {
      meta_all <- meta_all[cell_names, ]
    }

    return(meta_all)
}

.CSToMetaMulti <- function(object, cellScores, cts, sigma_names, nCC, slides) {
    ## Add cell score information to metadata for multi-slide data
    for (t in sigma_names) {
      for (sID in slides) {
        scores_slide <- cellScores[[t]][[sID]]
        
        for (ct in cts) {
          scores_ct <- scores_slide[[ct]]
          if (!is.null(scores_ct) && nrow(scores_ct) > 0) {
            meta_indices_ct_slide <- object@metaDataSub$slideID == sID & object@cellTypesSub == ct
            cells_ct_slide <- rownames(object@metaDataSub)[meta_indices_ct_slide]

            # Ensure scores match metadata rows
            if (all(rownames(scores_ct) %in% cells_ct_slide) && 
                length(rownames(scores_ct)) == length(cells_ct_slide)) {
              scores_ct_aligned <- scores_ct[cells_ct_slide, , drop = FALSE]
              for (cc_index in seq_len(nCC)) {
                cc_colname <- paste0("cellScore_", t, "_cc_index_", cc_index)
                object@metaDataSub[cells_ct_slide, cc_colname] <- 
                  scores_ct_aligned[, paste0("CC_", cc_index)]
              }
            } else {
              warning(paste("Mismatch between cell score rownames and metadata for:", sID, ct))
            }
          }
        }
      }
    }
    return(object@metaDataSub)
}

.computeGACCore <- function(object, cts, sigmaValues, scalePCs) {
    
    # some values that are needed
    nCC <- object@nCC
    PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
    sigma_names <- paste("sigma", sigmaValues, sep = "_")

    csgs <- .initializeCSGS(cts, sigma_names, nCC, 
           cellTypesSub = object@cellTypesSub, 
           cellNamesSub = rownames(object@normalizedDataSub),
           geneNamesSub = colnames(object@normalizedDataSub))

    cellScores <- csgs$cellScores
    geneScores <- csgs$geneScores


    ## go over all cell types, then over all sigma values
    for (tt in seq_along(sigmaValues)) {
      t <- sigma_names[tt]
      for (i in cts) {
        for (cc_index in seq_len(nCC)) {
          cc_name <- paste0("CC_", cc_index)
          w_1 <- object@skrCCAOut[[t]][[i]][, cc_index, drop = FALSE]

          ## get the sdev of the PCs
          if(!scalePCs) {
            sdev_use <- 1L
          } else {
            sdev_use <- object@pcaGlobal[[i]]$sdev
          }

          ## compute the gene scores
          geneScores[[t]][[i]][, cc_name] <- as.vector(
              matrix(w_1 * sdev_use, nrow = 1) %*%
                t(object@pcaGlobal[[i]]$rotation)
            )

          ## compute the cell scores -- independent of scalePCs
          cellScores[[t]][[i]][, cc_name] <- as.vector(PCmats[[i]] %*% w_1)
        }

      }
    }

    ## save cellscores and gene scores
    object@cellScores <- cellScores
    object@geneScores <- geneScores

    meta_all <- .CSToMeta(object, cellScores, cts, sigma_names, nCC)
    object@metaDataSub <- meta_all

    return(object)
}

.checkXW <- function(X, W, sID, ct) {
  if (!is.null(X) && !is.null(W) && nrow(X) > 0 && ncol(X) == nrow(W)) {
    return(TRUE)
  } else if (is.null(X) || is.null(W)) {
    warning("Missing PCA or Weight data for cell score calc: ", sID, ", ", ct)
  } else if (nrow(X) == 0) {
    warning("No cells for cell score calc: ", sID, ", ", ct)
  } else {
    warning("Dimension mismatch for cell score calc: ", sID, ", ", ct, " ", nrow(X), " ", nrow(W))
    }
  return(FALSE)
}

.createScoreMatrix <- function(nrows, ncols, row_names = NULL, col_names = NULL, fill_value = NA) {
  ## Helper function to create matrices with consistent naming
  mat <- matrix(fill_value, nrow = nrows, ncol = ncols)
  if (is.null(col_names)) {
    colnames(mat) <- paste0("CC_", 1:ncols)
  } else {
    colnames(mat) <- col_names
  }
  if (!is.null(row_names)) {
    rownames(mat) <- row_names
  }
  return(mat)
}

.computeGeneScores <- function(W_ct, pca_obj_ct, scalePCs, nCC, gene_names) {
  ## Helper function to compute gene scores
  gene_score_mat <- .createScoreMatrix(length(gene_names), nCC, row_names = gene_names)
  
  if (!is.null(pca_obj_ct) && !is.null(W_ct)) {
    rotation <- pca_obj_ct$rotation
    sdev <- pca_obj_ct$sdev
    
    for (cc in 1:nCC) {
      w_cc <- W_ct[, cc, drop = FALSE]
      
      ## get the sdev of the PCs
      if (!scalePCs) {
        sdev_use <- 1L
      } else {
        sdev_use <- sdev
      }
      
      ## compute the gene scores
      gene_scores_cc <- as.vector(
        matrix(w_cc * sdev_use, nrow = 1) %*% t(rotation)
      )
      gene_score_mat[, paste0("CC_", cc)] <- gene_scores_cc
    }
  }
  
  return(gene_score_mat)
}

.computeGACMultiCore <- function(object, cts, sigmaValues, scalePCs, slides) {
  
  # Initialize variables
  nCC <- object@nCC
  sigma_names <- paste("sigma", sigmaValues, sep = "_")
  
  # Initialize data structures
  csgs <- .initializeCSGSMulti(cts, sigma_names, nCC, slides)
  cellScores <- csgs$cellScores
  geneScores <- csgs$geneScores
  
  ## Iterate over all sigma values and slides
  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    W_list <- object@skrCCAOut[[t]] # Shared weights for this sigma
    
    # Calculate Cell Scores (Slide-Specific)
    for (sID in slides) {
      X_list_slide <- object@pcaResults[[sID]]
      meta_slide <- object@metaDataSub[object@metaDataSub$slideID == sID, ]
      celltype_slide <- object@cellTypesSub[object@metaDataSub$slideID == sID]
      
      for (ct in cts) {
        X_ct <- X_list_slide[[ct]]
        W_ct <- W_list[[ct]] # Shared weight matrix for cell type ct
        check_XW <- .checkXW(X_ct, W_ct, sID, ct)
        cell_ids_ct_slide <- rownames(meta_slide)[celltype_slide == ct]
        
        if (check_XW) {
          scores_mat <- X_ct %*% W_ct
          colnames(scores_mat) <- paste0("CC_", 1:nCC)
          rownames(scores_mat) <- cell_ids_ct_slide
        } else {
          # Create empty matrix if no cells or data mismatch
          scores_mat <- .createScoreMatrix(length(cell_ids_ct_slide), nCC, row_names = cell_ids_ct_slide)
        }
        
        cellScores[[t]][[sID]][[ct]] <- scores_mat
      }
    }
    
    # Calculate Gene Scores (Shared across slides)
    for (ct in cts) {
      pca_obj_ct <- object@pcaGlobal[[ct]]
      W_ct <- W_list[[ct]]
      gene_score_mat <- .computeGeneScores(W_ct, pca_obj_ct, scalePCs, nCC, object@geneList)
      geneScores[[t]][[ct]] <- gene_score_mat
    }
  }
  
  ## Save cell scores and gene scores
  object@cellScores <- cellScores
  object@geneScores <- geneScores
  
  ## Add cell scores to metadata
  meta_all <- .CSToMetaMulti(object, cellScores, cts, sigma_names, nCC, slides)
  object@metaDataSub <- meta_all
  
  return(object)
}

#' @rdname computeGeneAndCellScores
#' @aliases computeGeneAndCellScores,CoPro-method
#' @export
setMethod(
  "computeGeneAndCellScores", "CoPro",
  function(object) {
    input_check <- .checkInputGAC(object)
    cts <- input_check$cts
    sigmaValues <- input_check$sigmaValues
    scalePCs <- input_check$scalePCs

    object <- .computeGACCore(object, cts, sigmaValues, scalePCs)
    return(object)
  }
)

#' @rdname computeGeneAndCellScores
#' @aliases computeGeneAndCellScores,CoProMulti-method
#' @export
setMethod(
  "computeGeneAndCellScores", "CoProMulti",
  function(object) {
    input_check <- .checkInputGAC(object)
    cts <- input_check$cts
    sigmaValues <- input_check$sigmaValues
    scalePCs <- input_check$scalePCs
    slides <- input_check$slides

    object <- .computeGACMultiCore(object, cts, sigmaValues, scalePCs, slides)
    return(object)
  }
)


