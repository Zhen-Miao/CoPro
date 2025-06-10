#--- Pipeline Functions ---


#' Compute PCA on Integrated Multi-Slide Data
#'
#' Performs PCA on the integrated data stored within the `CoProm` object.
#' Assumes integration has created a common space across slides.
#'
#' @importFrom stats setNames prcomp
#' @importFrom irlba prcomp_irlba
#' @param object A `CoProm` object with the `integratedData` slot populated.
#' @param nPCA Number of principal components to compute for each cell type.
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'  Default is "raw"
#' @param center_per_slide After the global PCA, do we do center per slide
#'  again? By default this is set to FALSE
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#'
#' @return A `CoProm` object with the `pcaResults` slot populated.
#'         `pcaResults` structure: `list(slideID = list(cellType = pc_matrix))`.
#' @export
#' @rdname computePCAMulti
#' @aliases computePCAMulti,CoProm-method
setGeneric("computePCAMulti",
           function(object, nPCA = 40,
                    center = TRUE, scale. = TRUE,
                    dataUse = "raw",
                    center_per_slide = FALSE) standardGeneric("computePCAMulti"))

#' @rdname computePCAMulti
#' @aliases computePCAMulti,CoProm-method
#' @importFrom stats setNames prcomp predict
#' @importFrom irlba prcomp_irlba
#' @export
setMethod("computePCAMulti",
          "CoProMulti",
          function(object, nPCA = 40, center = TRUE, scale. = TRUE,
                   dataUse = "raw",
                   center_per_slide = FALSE) {
  ## arg match
  if(!(dataUse %in% c("raw", "integrated"))) stop("dataUse must be raw or integrated")
  ## Check if integrated data exists
  if (dataUse == "integrated" && length(object@integratedData) == 0) {
    stop("integratedData slot is empty. Run integrateSlidesMulti first.")
  }
  cts <- object@cellTypesOfInterest
  if(length(cts) == 0) {
    stop("cellTypesOfInterest not set. Run subsetDataMulti first.")
  }
  slides <- object@slideList

  # Initialize pcaResults structure
  pca_results_all <- setNames(vector("list", length = length(slides)), slides)
  for(sID in slides) {
    pca_results_all[[sID]] <- setNames(vector("list", length = length(cts)), cts)
  }

  pca_global <- setNames(vector("list", length = length(cts)), cts)

  # Perform PCA per cell type on the integrated data
  for (ct in cts) {
    message(paste("Performing PCA for cell type:", ct))

    if (dataUse == "integrated" && !ct %in% names(object@integratedData)) {
      warning(paste("No integrated data found for cell type:", ct, "- Skipping PCA."))
      next
    }

    if(dataUse == "integrated"){
      integrated_mat_ct <- object@integratedData[[ct]]
      # Needs verification based on integrateSlidesMulti output
    }else{ ## raw
      integrated_mat_ct <- object@normalizedDataSub[object@cellTypesSub == ct,]
    }

    # Ensure it's a matrix
    if (!is.matrix(integrated_mat_ct) && !inherits(integrated_mat_ct, "Matrix")) {
      message("converting data matrix into dense matrix")
      integrated_mat_ct <- as.matrix(integrated_mat_ct)
    }

    # Check dimensions
    if(nrow(integrated_mat_ct) != sum(object@cellTypesSub == ct)) {
      stop(paste("Integrated data dimensions mismatch for cell type:", ct))
    }

    if (center & scale.) {
      scaledData <- center_scale_matrix_opt(integrated_mat_ct)
      message("data centered and scaled")
    } else if (center) {
      scaledData <- t(t(integrated_mat_ct) - colMeans(integrated_mat_ct))
    } else {
      warning(paste("It is not recommended to skip both centering,",
                    "and scaling of the data, unless the data has been centered and",
                    "scaled when creating the CoPro object.",
                    sep = " "
      ))
      scaledData <- integrated_mat_ct
    }

    # Perform PCA on the combined integrated data for this cell type
    pca_ct <- prcomp_irlba(scaledData, n = nPCA,
                             center = FALSE, scale. = FALSE)
    message("PCA computed")
    pca_global[[ct]] <- pca_ct

    # Project each slide's data onto the shared PCs
    for (sID in slides) {
      slide_ID_ct <- object@metaDataSub$slideID[object@cellTypesSub == ct]
      row_names_ct <- rownames(object@metaDataSub)[object@cellTypesSub == ct]
      pca_sub <- pca_ct$x[slide_ID_ct == sID, , drop = FALSE]
      rownames(pca_sub) <- row_names_ct[slide_ID_ct == sID]

      if(center_per_slide){
        cat("Centering per slide","/n")
        pca_sub <- scale(pca_sub,center = TRUE, scale = FALSE)
      }

      pca_results_all[[sID]][[ct]] <- pca_sub

    }
  } # End loop over cell types

  object@pcaGlobal <- pca_global
  object@pcaResults <- pca_results_all
  object@nPCA <- nPCA
  # Default scalePCs set to True
  object@scalePCs <- TRUE

  return(object)
})








#' Compute Normalized Correlation for Multi-Slide Data
#'
#' Calculates normalized correlation, potentially per slide or aggregated.
#' Needs careful definition of the desired output.
#'
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @importFrom stats aggregate
#' @param object A `CoProm` object with skrCCA results.
#' @param tol Tolerance for SVD approximation of spectral norm.
#' @param calculationMode "aggregate" or "perSlide". Determines how correlation is computed/reported.
#'
#' @return `CoProm` object with `normalizedCorrelation` slot populated.
#' @export
#' @rdname computeNormalizedCorrelationMulti
#' @aliases computeNormalizedCorrelationMulti,CoProm-method
setGeneric("computeNormalizedCorrelationMulti", function(
    object, tol = 1e-4, calculationMode = c("aggregate", "perSlide")
) standardGeneric("computeNormalizedCorrelationMulti"))

#' @rdname computeNormalizedCorrelationMulti
#' @aliases computeNormalizedCorrelationMulti,CoProm-method
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @export
setMethod("computeNormalizedCorrelationMulti", "CoProMulti", function(
    object, tol = 1e-4, calculationMode = c("aggregate", "perSlide")) {

  calculationMode <- match.arg(calculationMode)

  # --- Input Checks ---
  if (length(object@skrCCAOut) == 0) stop("skrCCA results missing. Run runSkrCCAMulti.")
  if (length(object@pcaResults) == 0) stop("PCA results missing.")
  if (length(object@kernelMatrices) == 0) stop("Kernel matrices missing.")
  cts <- object@cellTypesOfInterest
  if (length(cts) < 2) stop("Need at least two cell types.")
  slides <- object@slideList
  sigmas_run <- names(object@skrCCAOut)
  if (length(sigmas_run) == 0) stop("No skrCCA results found.")
  nCC <- object@nCC
  if(is.null(nCC) || length(nCC)==0) nCC <- ncol(object@skrCCAOut[[1]][[1]]) # Infer nCC

  # --- Precompute Spectral Norms (Per Slide, Per Sigma) ---
  message("Calculating spectral norms (can take time)...")
  norm_K_all <- setNames(vector("list", length = length(sigmas_run)), sigmas_run)
  for (sig_name in sigmas_run) {
    norm_K_sigma <- setNames(vector("list", length = length(slides)), slides)
    K_list_sigma <- object@kernelMatrices[[sig_name]]
    pair_cell_types <- combn(cts, 2)

    for (sID in slides) {
      norm_K_slide <- setNames(vector("list", length = length(cts)), cts)
      for (ct_i in cts) norm_K_slide[[ct_i]] <- setNames(vector("list", length=length(cts)), cts)

      K_list_slide <- K_list_sigma[[sID]]

      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]
        K <- K_list_slide[[ct_i]][[ct_j]]
        if (!is.null(K) && nrow(K) > 0 && ncol(K) > 0) {
          # Handle potential sparse matrices if irlba supports them well
          svd_d <- tryCatch({
            irlba::irlba(K, nv = 1, nu = 0, tol = tol)$d[1] # Only need singular value
          }, error = function(e) {
            warning(paste("SVD failed for K[",sID,",",ct_i,",",ct_j,"]:", e$message))
            NA # Return NA if SVD fails
          })
          norm_K_slide[[ct_i]][[ct_j]] <- svd_d
          norm_K_slide[[ct_j]][[ct_i]] <- svd_d # Store transpose too
        } else {
          norm_K_slide[[ct_i]][[ct_j]] <- NA
          norm_K_slide[[ct_j]][[ct_i]] <- NA
        }
      }
      norm_K_sigma[[sID]] <- norm_K_slide
    }
    norm_K_all[[sig_name]] <- norm_K_sigma
  }
  message("Finished calculating spectral norms.")

  # --- Calculate Correlation ---
  correlation_results <- setNames(vector("list", length = length(sigmas_run)), sigmas_run)

  for (sig_name in sigmas_run) {
    W_list_sigma <- object@skrCCAOut[[sig_name]] # Shared weights

    if (calculationMode == "perSlide") {
      correlation_per_slide <- setNames(vector("list", length=length(slides)), slides)
      for(sID in slides) {
        df_slide <- data.frame(
          sigmaValue = character(),
          slideID = character(),
          cellType1 = character(), cellType2 = character(),
          CC_index = integer(), normalizedCorrelation = numeric(),
          stringsAsFactors = FALSE
        )
        X_list_slide <- object@pcaResults[[sID]]
        K_list_slide <- object@kernelMatrices[[sig_name]][[sID]]
        norm_K_slide <- norm_K_all[[sig_name]][[sID]]

        if(is.null(X_list_slide) || is.null(K_list_slide)) next # Skip if data missing

        for (pp in seq_len(ncol(pair_cell_types))) {
          ct_i <- pair_cell_types[1, pp]
          ct_j <- pair_cell_types[2, pp]

          X_i <- X_list_slide[[ct_i]]
          X_j <- X_list_slide[[ct_j]]
          K_ij <- K_list_slide[[ct_i]][[ct_j]]
          norm_K_ij <- norm_K_slide[[ct_i]][[ct_j]]

          if(is.null(X_i) || is.null(X_j) || is.null(K_ij) ||
             is.na(norm_K_ij) || norm_K_ij < 1e-9 ||
             nrow(X_i)==0 || nrow(X_j)==0) next # Skip if data missing/invalid

          for (cc in 1:nCC) {
            w_i <- W_list_sigma[[ct_i]][, cc, drop = FALSE]
            w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]

            Xiw <- X_i %*% w_i
            Xjw <- X_j %*% w_j

            numerator <- (t(Xiw) %*% K_ij %*% Xjw)[1, 1]
            denom_norm <- sqrt(sum(Xiw^2)) * sqrt(sum(Xjw^2)) * norm_K_ij

            norm_corr_val <- ifelse(abs(denom_norm) < 1e-9, 0, numerator / denom_norm)

            df_slide <- rbind(df_slide, data.frame(
              sigmaValue = as.numeric(gsub("sigma_", "", sig_name)),
              slideID = sID,
              cellType1 = ct_i, cellType2 = ct_j,
              CC_index = cc, normalizedCorrelation = norm_corr_val,
              stringsAsFactors = FALSE
            ))
          } # end CC loop
        } # end pair loop
        correlation_per_slide[[sID]] <- df_slide
      } # end slide loop
      correlation_results[[sig_name]] <- do.call(rbind, correlation_per_slide)

    } else { # calculationMode == "aggregate"
      df_agg <- data.frame(
        sigmaValue = numeric(),
        cellType1 = character(), cellType2 = character(),
        CC_index = integer(), aggregateCorrelation = numeric(),
        stringsAsFactors = FALSE
      )
      # Calculate sum(w_i^T X_iq^T K_ijq X_jq w_j) / ( sum(sqrt(sum(Xiw_q^2))) * sum(sqrt(sum(Xjw_q^2))) * mean(normK_ijq) ) ?
      # Or sum(numerator_q) / sum(denominator_q) ?
      # Simpler: calculate objective value achieved by weights and normalize somehow?
      # Let's calculate sum(num_q) / ( sqrt(sum(||Xiw_q||^2)) * sqrt(sum(||Xjw_q||^2)) * mean(normK_ijq) )
      # This needs careful definition based on desired interpretation.

      message("Aggregate correlation calculation needs precise definition. Placeholder implementation.")
      # Placeholder: Calculate mean per-slide correlation
      temp_res <- computeNormalizedCorrelationMulti(object, tol=tol, calculationMode="perSlide")
      df_per_slide <- temp_res@normalizedCorrelation[[sig_name]]
      if(!is.null(df_per_slide) && nrow(df_per_slide) > 0){
        df_agg <- aggregate(normalizedCorrelation ~ sigmaValue + cellType1 + cellType2 + CC_index,
                            data=df_per_slide, FUN = mean)
        colnames(df_agg)[colnames(df_agg) == "normalizedCorrelation"] <- "aggregateCorrelation"
      } else {
        df_agg <- data.frame( # Empty df if no per-slide results
          sigmaValue = numeric(), cellType1 = character(), cellType2 = character(),
          CC_index = integer(), aggregateCorrelation = numeric(), stringsAsFactors = FALSE
        )
      }

      correlation_results[[sig_name]] <- df_agg
    }
  } # End sigma loop

  object@normalizedCorrelation <- correlation_results

  # --- Choose Optimal Sigma (based on aggregate or mean per-slide CC1 correlation) ---
  if(length(correlation_results) > 0) {
    first_sigma_res <- correlation_results[[1]]
    if(!is.null(first_sigma_res) && nrow(first_sigma_res) > 0) {
      corr_col_name <- ifelse("aggregateCorrelation" %in% names(first_sigma_res), "aggregateCorrelation", "normalizedCorrelation")
      all_corrs <- do.call(rbind, correlation_results)
      cc1_corrs <- all_corrs[all_corrs$CC_index == 1, ]

      if(nrow(cc1_corrs) > 0) {
        # Average correlation across pairs and potentially slides for CC1 for each sigma
        mean_corr_per_sigma <- tapply(cc1_corrs[[corr_col_name]], cc1_corrs$sigmaValue, mean, na.rm=TRUE)
        if(any(!is.na(mean_corr_per_sigma))){
          object@sigmaValueChoice <- as.numeric(names(which.max(mean_corr_per_sigma)))
        } else {
          warning("Could not determine optimal sigma: All correlations were NA.")
          object@sigmaValueChoice <- numeric() # Empty
        }

      } else {
        warning("Could not determine optimal sigma: No CC1 correlation values found.")
        object@sigmaValueChoice <- numeric() # Empty
      }
    } else {
      warning("Could not determine optimal sigma: Correlation results are empty.")
      object@sigmaValueChoice <- numeric() # Empty
    }
  } else {
    object@sigmaValueChoice <- numeric() # Empty if no results
  }

  return(object)
})


#' Compute Gene and Cell Scores for Multi-Slide Data
#'
#' Calculates slide-specific cell scores and potentially shared gene scores.
#'
#' @importFrom stats setNames
#' @param object A `CoProm` object with skrCCA results.
#' @param sigmaChoice A specific sigma value to use. If NULL (default), uses `object@sigmaValueChoice`. If that's also NULL, uses the first available sigma.
#'
#' @return `CoProm` object with `cellScores` and `geneScores` slots populated.
#' @export
#' @rdname computeGeneAndCellScoresMulti
#' @aliases computeGeneAndCellScoresMulti,CoProm-method
setGeneric("computeGeneAndCellScoresMulti", function(
    object, sigmaChoice = NULL) standardGeneric("computeGeneAndCellScoresMulti"))

#' @rdname computeGeneAndCellScoresMulti
#' @aliases computeGeneAndCellScoresMulti,CoProm-method
#' @importFrom stats setNames
#' @export
setMethod("computeGeneAndCellScoresMulti", "CoProMulti", function(
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

    for (ct in cts) {
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

       for(ct in cts){
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
                       gene_scores_cc <- rotation %*% (w_cc * sdev[1:nrow(w_cc)]) # Check dimensions match
                   } else {
                       gene_scores_cc <- rotation %*% w_cc
                   }
                    gene_score_mat[, paste0("CC_", cc)] <- gene_scores_cc[,1]
               }
           }
            gene_scores_sigma[[ct]] <- gene_score_mat
       }

  # Store under the chosen sigma
  object@geneScores[[sig_name]] <- gene_scores_sigma

  # --- Add Cell Scores to metaDataSub ---
  # Create columns for the chosen sigma
  if(!is.null(object@cellScores[[sig_name]])){
    for(sID in slides){
      scores_slide <- object@cellScores[[sig_name]][[sID]]
      # meta_indices_slide <- which(object@metaDataSub$slideID == sID)

      for(ct in cts){
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
  }
  return(object)
})
