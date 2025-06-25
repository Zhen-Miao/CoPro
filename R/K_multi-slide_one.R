


#' Compute Gene and Cell Scores for Multi-Slide Data
#'
#' Calculates slide-specific cell scores and potentially shared gene scores.
#'
#' @importFrom stats setNames
#' @param object A `CoProMulti` object with skrCCA results.
#' @param sigmaChoice A specific sigma value to use. If NULL (default), uses `object@sigmaValueChoice`. If that's also NULL, uses the first available sigma.
#'
#' @return `CoProMulti` object with `cellScores` and `geneScores` slots populated.
#' @export
#' @rdname computeGeneAndCellScoresMultiOne
#' @aliases computeGeneAndCellScoresMultiOne,CoProm-method
setGeneric("computeGeneAndCellScoresMultiOne", function(
    object, sigmaChoice = NULL) standardGeneric("computeGeneAndCellScoresMultiOne"))

#' @rdname computeGeneAndCellScoresMultiOne
#' @aliases computeGeneAndCellScoresMultiOne,CoProm-method
#' @importFrom stats setNames
#' @export
setMethod("computeGeneAndCellScoresMultiOne", "CoProMulti", function(
    object, sigmaChoice = NULL) {

  # --- Input Checks ---
  if (length(object@skrCCAOut) == 0) stop("skrCCA results missing.")
  if (length(object@pcaResults) == 0) stop("PCA results missing.")
  cts <- object@cellTypesOfInterest
  slides <- object@slideList
  nCC <- object@nCC
  if(is.null(nCC) || length(nCC)==0) nCC <- ncol(object@skrCCAOut[[1]][[1]]) # Infer nCC

  # Determine sigma to use
  if (is.null(sigmaChoice)) {
    sigmaChoice <- object@sigmaValueChoice
    if (is.null(sigmaChoice) || length(sigmaChoice) == 0) {
      if(length(object@skrCCAOut) > 0) {
        sigmaChoice <- as.numeric(gsub("sigma_","", names(object@skrCCAOut)[1]))
        warning("sigmaValueChoice not set or found. Using first available sigma: ", sigmaChoice)
      } else {
        stop("Cannot compute scores: No skrCCA results available to choose sigma.")
      }
    }
  }
  sig_name <- paste0("sigma_", sigmaChoice)
  if (!sig_name %in% names(object@skrCCAOut)) {
    stop(paste("Chosen sigma", sigmaChoice, "not found in skrCCA results."))
  }

  W_list <- object@skrCCAOut[[sig_name]] # Shared weights

  # --- Calculate Cell Scores (Slide-Specific) ---
  cell_scores_all <- setNames(vector("list", length = length(slides)), slides)
  for(sID in slides) {
    cell_scores_slide <- setNames(vector("list", length = length(cts)), cts)
    X_list_slide <- object@pcaResults[[sID]]
    meta_slide <- object@metaDataSub[object@metaDataSub$slideID == sID, ]
    celltype_slide <- object@cellTypesSub[object@metaDataSub$slideID == sID]

    ct <- cts
    X_ct <- X_list_slide[[ct]]
    W_ct <- W_list[[ct]] # Shared weight matrix for cell type ct

      if(!is.null(X_ct) && !is.null(W_ct) && nrow(X_ct)>0 && ncol(X_ct) == nrow(W_ct)){
        scores_mat <- X_ct %*% W_ct
        colnames(scores_mat) <- paste0("CC_", 1:nCC)
        rownames(scores_mat) <- rownames(meta_slide)[celltype_slide == ct]

        cell_scores_slide[[ct]] <- scores_mat
      } else {
        # Create empty matrix if no cells or data mismatch
        cell_ids_ct_slide <- rownames(meta_slide)[celltype_slide == ct]
        cell_scores_slide[[ct]] <- matrix(NA, nrow=length(cell_ids_ct_slide), ncol=nCC,
                                          dimnames=list(cell_ids_ct_slide, paste0("CC_", 1:nCC)))
        if(is.null(X_ct) || is.null(W_ct)) warning("Missing PCA or Weight data for cell score calc: ", sID, ", ", ct)
        else if(nrow(X_ct) == 0) warning("No cells for cell score calc: ", sID, ", ", ct)
        else warning("Dimension mismatch for cell score calc: ", sID, ", ", ct)
      }

    cell_scores_all[[sID]] <- cell_scores_slide
  }
  # Store under the chosen sigma
  object@cellScores[[sig_name]] <- cell_scores_all


  # --- Calculate Gene Scores (Shared - requires PCA loadings) ---
  # This requires access to the PCA rotation matrices, which were not explicitly
  # stored in the proposed pcaResults structure. We need to adapt computePCAMulti
  # to store these, or recompute/access them from the integrated object if possible.
  # Assuming PCA objects per cell type (run on integrated data) are accessible:
  # For simplicity, let's assume a function getPCALoadings(object, cellType) exists.

  gene_scores_sigma <- setNames(vector("list", length=length(cts)), cts)
  scalePCs <- object@scalePCs # Should be FALSE if PCA done on integrated/scaled data

       ct <- cts
         pca_obj_ct <- object@pcaGlobal[[ct]]
         W_ct <- W_list[[ct]]
           gene_score_mat <- matrix(NA, nrow=length(object@geneList), ncol=nCC,
                                      dimnames=list(object@geneList, paste0("CC_", 1:nCC)))

           if(!is.null(pca_obj_ct) && !is.null(W_ct)){
              rotation <- pca_obj_ct$rotation # Assuming standard prcomp structure
              sdev <- pca_obj_ct$sdev

              for(cc in 1:nCC){
                   w_cc <- W_ct[, cc, drop=FALSE]
                   # Project weights back to gene space
                   # Note: scaling depends on whether input to PCA was scaled in runSkrCCA AND scalePCs setting
                   # If PCA input was scaled (e.g. scale.=TRUE in prcomp) OR scalePCs=TRUE, need sdev.
                   if (scalePCs) { # This logic might need refinement based on exact PCA steps
                       gene_scores_cc <- rotation %*% (w_cc * sdev[seq_len(nrow(w_cc))]) # Check dimensions match
                   } else {
                       gene_scores_cc <- rotation %*% w_cc
                   }
                    gene_score_mat[, paste0("CC_", cc)] <- gene_scores_cc[,1]
               }
           }
            gene_scores_sigma[[ct]] <- gene_score_mat


  # Store under the chosen sigma
  object@geneScores[[sig_name]] <- gene_scores_sigma

  # --- Add Cell Scores to metaDataSub ---
  # Create columns for the chosen sigma
  if(!is.null(object@cellScores[[sig_name]])){
    for(sID in slides){
      scores_slide <- object@cellScores[[sig_name]][[sID]]
      # meta_indices_slide <- which(object@metaDataSub$slideID == sID)

      ct <- cts
        scores_ct <- scores_slide[[ct]]
        if(!is.null(scores_ct) && nrow(scores_ct)>0){
          meta_indices_ct_slide <- object@metaDataSub$slideID == sID & object@cellTypesSub == ct # Use cellType from metaDataSub
          cells_ct_slide <- rownames(object@metaDataSub)[meta_indices_ct_slide]

          # Ensure scores match metadata rows
          if(all(rownames(scores_ct) %in% cells_ct_slide) && length(rownames(scores_ct)) == length(cells_ct_slide)){
            scores_ct_aligned <- scores_ct[cells_ct_slide, , drop=FALSE] # Align order
            for(cc in 1:nCC){
              cc_colname <- paste0("CC", cc, "_Score_", sig_name)
              object@metaDataSub[cells_ct_slide, cc_colname] <- scores_ct_aligned[, paste0("CC_", cc)]
            }
          } else {
            warning(paste("Mismatch between cell score rownames and metadata for:", sID, ct))
          }
        }

    }
  }
  return(object)
})





