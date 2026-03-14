#' computeGeneAndCellScores
#' @importFrom stats setNames
#' @importFrom methods slotNames
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
    if(is(object, "CoProSingle")) {
      is_multi <- FALSE
    } else if(is(object, "CoProMulti")) {
      is_multi <- TRUE
    } else {
      stop("object is not a CoProSingle or CoProMulti object")
    }
    if(is_multi) {
      slides <- getSlideList(object)
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

.initializeCSGS <- function(cts, sigma_names, nCC, object) {
    ## Initialize flat structure for cell scores and gene scores
    cellScores <- list()
    geneScores <- list()
    
    ## Get required data
    cellNamesSub <- rownames(object@metaDataSub)
    cellTypesSub <- object@cellTypesSub
    geneNamesSub <- if("geneList" %in% slotNames(object)) object@geneList else colnames(object@normalizedDataSub)
    
    # Create flat structures with informative names
    for (t in sigma_names) {
      sigma_val <- as.numeric(gsub("sigma_", "", t))
      
      for (ct in cts) {
        # Cell scores with flat names
        cell_flat_name <- .createCellScoresName(sigma_val, ct, slide = NULL)
        cellScores[[cell_flat_name]] <- .createScoreMatrix(
          nrows = sum(cellTypesSub == ct),
          ncols = nCC,
          row_names = cellNamesSub[cellTypesSub == ct]
        )
        
        # Gene scores with flat names (we can create similar naming for gene scores)
        gene_flat_name <- .createGeneScoresName(sigma_val, ct, slide = NULL)
        geneScores[[gene_flat_name]] <- .createScoreMatrix(
          nrows = length(geneNamesSub),
          ncols = nCC,
          row_names = geneNamesSub
        )
      }
    }

    return(list(cellScores = cellScores, geneScores = geneScores))
}

.CSToMeta <- function(object, cellScores, cts, sigma_names, nCC) {
    ## Function to add cell score information to metadata using flat structure
    
    for (t in sigma_names) {
      sigma_val <- as.numeric(gsub("sigma_", "", t))
      
      for (ct in cts) {
        cell_flat_name <- .createCellScoresName(sigma_val, ct, slide = NULL)
        scores_ct <- cellScores[[cell_flat_name]]
        
        if (!is.null(scores_ct) && nrow(scores_ct) > 0) {
          # Get all cells of this cell type
          meta_indices_ct <- object@cellTypesSub == ct
          cells_ct <- rownames(object@metaDataSub)[meta_indices_ct]

          # Ensure scores match metadata rows
          if (all(rownames(scores_ct) %in% cells_ct) && 
              length(rownames(scores_ct)) == length(cells_ct)) {
            scores_ct_aligned <- scores_ct[cells_ct, , drop = FALSE]
            
            for (cc_index in seq_len(nCC)) {
              cc_colname <- paste0("cellScore_", t, "_cc_index_", cc_index)
              object@metaDataSub[cells_ct, cc_colname] <- 
                scores_ct_aligned[, paste0("CC_", cc_index)]
            }
          } else {
            warning(paste("Mismatch between cell score rownames and metadata for cell type:", ct))
          }
        }
      }
    }
    return(object@metaDataSub)
}

.computeCellScores <- function(X_data, W_ct, cell_ids, nCC) {
  ## Helper function to compute cell scores for a given X matrix and weight matrix
  if (is.null(X_data) || is.null(W_ct) || length(cell_ids) == 0) {
    return(.createScoreMatrix(length(cell_ids), nCC, row_names = cell_ids))
  }
  
  scores_mat <- X_data %*% W_ct
  colnames(scores_mat) <- paste0("CC_", 1:nCC)
  rownames(scores_mat) <- cell_ids
  return(scores_mat)
}

.computeGACCore <- function(object, cts, sigmaValues, scalePCs) {
    
    # some values that are needed
    nCC <- object@nCC
    PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
    sigma_names <- paste("sigma", sigmaValues, sep = "_")

    csgs <- .initializeCSGS(cts, sigma_names, nCC, object)

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
          gene_flat_name <- .createGeneScoresName(sigmaValues[tt], i, slide = NULL)
          geneScores[[gene_flat_name]][, cc_name] <- as.vector(
              matrix(w_1 * sdev_use, nrow = 1) %*%
                t(object@pcaGlobal[[i]]$rotation)
            )

          ## compute the cell scores -- independent of scalePCs
          cell_flat_name <- .createCellScoresName(sigmaValues[tt], i, slide = NULL)
          cellScores[[cell_flat_name]][, cc_name] <- as.vector(PCmats[[i]] %*% w_1)
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
  if (ncols == 0) {
    stop("ncols must be > 0 for score matrix creation")
  }
  mat <- matrix(fill_value, nrow = nrows, ncol = ncols)
  if (is.null(col_names)) {
    colnames(mat) <- paste0("CC_", seq_len(ncols))
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
  
  # Initialize data structures - now aggregated across slides
    csgs <- .initializeCSGS(cts, sigma_names, nCC, object)
    cellScores <- csgs$cellScores
    geneScores <- csgs$geneScores

  ## Iterate over all sigma values and slides to compute and aggregate scores
  for (tt in seq_along(sigmaValues)) {
    t <- sigma_names[tt]
    W_list <- object@skrCCAOut[[t]] # Shared weights for this sigma
    
    # Calculate Cell Scores (Slide-Specific, then aggregated)
    for (sID in slides) {
      X_list_slide <- object@pcaResults[[sID]]
      slide_indices <- .getSlideIndices(object, sID)
      meta_slide <- object@metaDataSub[slide_indices, ]
      celltype_slide <- object@cellTypesSub[slide_indices]
      
      for (ct in cts) {
        X_ct <- X_list_slide[[ct]]
        W_ct <- W_list[[ct]] # Shared weight matrix for cell type ct
        check_XW <- .checkXW(X_ct, W_ct, sID, ct)
        cell_ids_ct_slide <- rownames(meta_slide)[celltype_slide == ct]
        
        if (check_XW && length(cell_ids_ct_slide) > 0) {
          # Compute scores for this slide and cell type
          scores_mat_slide <- X_ct %*% W_ct
          colnames(scores_mat_slide) <- paste0("CC_", 1:nCC)
          rownames(scores_mat_slide) <- cell_ids_ct_slide
          
          # Add these scores to the aggregated matrix (flat structure)
          cell_flat_name <- .createCellScoresName(sigmaValues[tt], ct, slide = NULL)
          cellScores[[cell_flat_name]][cell_ids_ct_slide, ] <- scores_mat_slide
        }
      }
    }
    
    # Calculate Gene Scores (Shared across slides)
    for (ct in cts) {
      pca_obj_ct <- object@pcaGlobal[[ct]]
      W_ct <- W_list[[ct]]
      gene_score_mat <- .computeGeneScores(W_ct, pca_obj_ct, scalePCs, nCC, object@geneList)
      gene_flat_name <- .createGeneScoresName(sigmaValues[tt], ct, slide = NULL)
      geneScores[[gene_flat_name]] <- gene_score_mat
    }
  }
  
  ## Save cell scores and gene scores
  object@cellScores <- cellScores
  object@geneScores <- geneScores
  
  ## Add cell scores to metadata
  meta_all <- .CSToMeta(object, cellScores, cts, sigma_names, nCC)
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


