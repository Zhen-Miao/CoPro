

#' Compute Distances within each Slide -- one cell type
#'
#' Calculates pairwise cell distances for specified cell types within each slide.
#'
#' @importFrom fields rdist
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @param object A `CoProm` object.
#' @param distType Type of distance ("Euclidean2D", "Euclidean3D").
#' @param xDistScale Scaling factors for coordinates.
#' @param yDistScale Scaling factors for coordinates.
#' @param zDistScale Scaling factors for coordinates.
#' @param normalizeDistance Normalize distances based on 0.01 percentile?
#' @param truncateLowDist Truncate very small distances?
#' @param verbose Print progress messages?
#'
#' @return `CoProm` object with the `distances` slot populated:
#'         `list(slideID = list(cellType1 = list(cellType2 = dist_matrix)))`.
#' @export
#' @rdname computeDistanceMultiOne
#' @aliases computeDistanceMultiOne,CoProm-method
setGeneric("computeDistanceMultiOne", function(
    object, distType = c("Euclidean2D", "Euclidean3D"),
    xDistScale = 1, yDistScale = 1, zDistScale = 1,
    normalizeDistance = TRUE, truncateLowDist = TRUE,
    verbose = TRUE) standardGeneric("computeDistanceMultiOne"))

#' @rdname computeDistanceMultiOne
#' @aliases computeDistanceMultiOne,CoProm-method
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @importFrom fields rdist
#' @export
setMethod("computeDistanceMultiOne", "CoProm", function(
    object, distType = c("Euclidean2D", "Euclidean3D"),
    xDistScale = 1, yDistScale = 1, zDistScale = 1,
    normalizeDistance = TRUE, truncateLowDist = TRUE,
    verbose = TRUE) {

  distType <- match.arg(distType)
  cts <- object@cellTypesOfInterest
  if (length(cts) > 1) {
    stop(paste("Only run this function when there is one cell type.",
               "Please run computeDistanceMulti() instead."))
  }
  slides <- object@slideList

  distances_all <- setNames(vector("list", length = length(slides)), slides)
  global_min_1percentile <- Inf # To normalize across slides

  for (sID in slides) {
    if (verbose) message(paste("Computing distances for slide:", sID))

    # Initialize structure for this slide
    distances_slide <- setNames(vector("list", length = length(cts)), cts)

    distances_slide[[cts]] <- setNames(vector("list", length = length(cts)), cts)


    # Get subset indices for the current slide
    slide_indices <- which(object@metaDataSub$slideID == sID)
    locationData_slide <- object@locationDataSub[slide_indices, , drop = FALSE]
    # cellTypes_slide <- object@cellTypesSub[slide_indices]
    # cellIDs_slide <- rownames(locationData_slide)

    # pair_cell_types <- combn(cts, 2)
    dist_1percentile_slide <- vector(mode = "numeric", length = 1)
    # all_zero_pairs <- c() # Track pairs with all zero distances

    ct_i <- ct_j <- cts

    loc_i <- locationData_slide

      if (distType == "Euclidean2D") {
        mat1 <- cbind(loc_i$x * xDistScale, loc_i$y * yDistScale)
      } else { # Euclidean3D
        if (!"z" %in% colnames(loc_i)) stop("z coordinate missing for 3D distance.")
        mat1 <- cbind(loc_i$x * xDistScale, loc_i$y * yDistScale, loc_i$z * zDistScale)
      }

      distances_ij <- fields::rdist(mat1)
      rownames(distances_ij) <- colnames(distances_ij) <- rownames(loc_i)
      diag(distances_ij) <- Inf

      min_dist_nonzero <- min(distances_ij[distances_ij > 0], Inf)

      if (any(distances_ij == 0)) {
        warning(paste("Zero distances detected between", ct_i, "and", ct_j, "in slide", sID,
                      ". Replacing with smallest non-zero distance:", min_dist_nonzero))
        distances_ij[distances_ij == 0] <- min_dist_nonzero
      }

      if(is.infinite(min_dist_nonzero)){
        stop(paste("All distances are zero or missing between",
                      ct_i, "and", ct_j, "in slide", sID))
        # dist_1percentile_slide[pp] <- NA
        # all_zero_pairs <- c(all_zero_pairs, pp)
      } else {
        percentile_choice <- min(1e-3, 2/(max(nrow(distances_ij), ncol(distances_ij))))
        dist_1percentile_slide <- quantile(distances_ij[distances_ij > 0], percentile_choice)
        global_min_1percentile <- min(global_min_1percentile, dist_1percentile_slide, na.rm = TRUE)

        if (truncateLowDist) {
          distances_ij[distances_ij < dist_1percentile_slide] <- dist_1percentile_slide
        }
      }
      distances_slide[[ct_i]][[ct_j]] <- distances_ij
      if (verbose && !is.infinite(min_dist_nonzero)) {
        cat("Slide:", sID, ", Pair:", ct_i, "-", ct_j, "\n")
        print(quantile(distances_ij))
      }

    distances_all[[sID]] <- distances_slide
  } # End slide loop (sID)

  # Normalize across slides if requested
  if (normalizeDistance) {
    if(is.infinite(global_min_1percentile)){
      warning("Cannot normalize distances - no valid non-zero distances found across slides.")
    } else {
      scaling_factor <- 0.01 / global_min_1percentile
      if (verbose) cat("Global distance scaling factor:", scaling_factor, "\n")

      for (sID in slides) {
          ct_i <- ct_j <- cts
          # Avoid scaling if matrix is NA/empty or had all zeros
          if(!is.null(distances_all[[sID]][[ct_i]][[ct_j]])){
            distances_all[[sID]][[ct_i]][[ct_j]] <- distances_all[[sID]][[ct_i]][[ct_j]] * scaling_factor
          }
      }
    }
  }

  object@distances <- distances_all
  return(object)
})


#' Compute Kernel Matrices for Multi-Slide Data
#'
#' Calculates kernel matrices for each slide based on distances and sigma values.
#'
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @param object A `CoProm` object with the `distances` slot populated.
#' @param sigmaValues A numeric vector of sigma values.
#' @param lowerLimit Lower threshold for kernel values.
#' @param upperQuantile Quantile for clipping upper kernel values.
#' @param normalizeKernel Normalize kernel matrix values? (Not recommended as it breaks CCA properties if done row/col-wise).
#' @param minAveCellNeighor Minimum average neighbors check (applied per slide).
#' @param verbose Print progress?
#'
#' @return `CoProm` object with `kernelMatrices` slot populated:
#'         `list(sigma = list(slideID = list(cellType1 = list(cellType2 = kernel_matrix))))`.
#' @export
#' @rdname computeKernelMatrixMultiOne
#' @aliases computeKernelMatrixMultiOne,CoProm-method
setGeneric("computeKernelMatrixMultiOne", function(
    object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
    normalizeKernel = FALSE, minAveCellNeighor = 2,
    verbose = TRUE) standardGeneric("computeKernelMatrixMultiOne"))

#' @rdname computeKernelMatrixMultiOne
#' @aliases computeKernelMatrixMultiOne,CoProm-method
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @export
setMethod("computeKernelMatrixMultiOne", "CoProm", function(
    object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
    normalizeKernel = FALSE, minAveCellNeighor = 2,
    verbose = TRUE) {

  if (length(object@distances) == 0) stop("Distances not computed. Run computeDistanceMulti first.")
  if (normalizeKernel) warning("normalizeKernel=TRUE is generally not recommended for CCA-based methods.")

  cts <- object@cellTypesOfInterest
  if (length(cts) != 1) stop("Only one cellTypesOfInterest is allowed.")
  slides <- object@slideList
  if (length(sigmaValues) == 0) stop("sigmaValues must be provided.")
  object@sigmaValues <- sigmaValues # Store provided sigma values

  kernel_matrices_all <- setNames(vector("list", length(sigmaValues)),
                                  paste0("sigma_", sigmaValues))

  for (tt in seq_along(sigmaValues)) {
    sigma_val <- sigmaValues[tt]
    sigma_name <- names(kernel_matrices_all)[tt]
    if (verbose) message(paste("Computing kernels for sigma =", sigma_val))

    kernels_for_sigma <- setNames(vector("list", length = length(slides)), slides)
    sigma_valid_across_slides <- TRUE # Track if this sigma works for all slides

    for (sID in slides) {
      # Initialize structure for this slide
      kernels_slide <- setNames(vector("list", length = length(cts)), cts)
      for (i in cts) {
        kernels_slide[[i]] <- setNames(vector("list", length = length(cts)), cts)
      }

      dist_slide <- object@distances[[sID]]
      if(length(dist_slide) == 0) {
        warning(paste("No distance matrix found for slide", sID))
        sigma_valid_across_slides <- FALSE
        next # Skip slide if no distances
      }

      sigma_valid_for_slide <- TRUE
      ct_i <- ct_j <- cts
      dist_mat_ij <- dist_slide[[ct_i]][[ct_j]]

        if (length(dist_mat_ij) == 0 || all(is.na(dist_mat_ij)) ||
            nrow(dist_mat_ij) == 0 || ncol(dist_mat_ij) == 0 ) {
          if(verbose) message(paste("  Skipping kernel for",
                                    ct_i, "-", ct_j, "in slide", sID,
                                    "(missing/empty distance matrix)"))
          kernels_slide[[ct_i]][[ct_j]] <- matrix(nrow=nrow(dist_mat_ij), ncol=ncol(dist_mat_ij)) # Store empty/NA
          next
        }

        kernel_current <- kernel_from_distance( # Assumes kernel_from_distance is available
          sigma = sigma_val,
          dist_mat = dist_mat_ij,
          lower_limit = lowerLimit
        )

        # Check sparsity (using helper .CheckSigmaValuesToRemove - needs access or re-implementation)
        # Simplified check: if matrix is nearly all zero, flag sigma for this slide
        n_cells1 <- nrow(kernel_current)
        n_cells2 <- ncol(kernel_current)
        min_count_nonzero <- minAveCellNeighor * min(n_cells1, n_cells2)
        if (sum(kernel_current > lowerLimit) < min_count_nonzero) {
          warning(paste("Kernel matrix for", ct_i,"-", ct_j, "in slide", sID,
                        "with sigma =", sigma_val,
                        "is too sparse. Flagging sigma."))
          sigma_valid_for_slide <- FALSE
          # Store the sparse kernel anyway? Or NULL? Let's store it for now.
          # kernels_slide[[ct_i]][[ct_j]] <- NULL # Or keep the sparse matrix
        }

        if (anyNA(kernel_current)) {
          warning(paste("NA values in kernel for", ct_i, "-", ct_j, "in slide",
                        sID, "sigma =", sigma_val))
          # Decide how to handle NAs, maybe imputation or flagging
        }

        # Clipping large values
        upper_clip <- quantile(kernel_current[kernel_current >= lowerLimit & !is.na(kernel_current)], upperQuantile, na.rm = TRUE)
        if(!is.na(upper_clip)){
          kernel_current[kernel_current >= upper_clip & !is.na(kernel_current)] <- upper_clip
        }


        # Normalization (apply cautiously)
        if (normalizeKernel) {
          # Example: Scale by median row sum (use with caution)
          rs_kernel <- rowSums(kernel_current, na.rm=TRUE)
          median_rs <- median(rs_kernel[rs_kernel > 1e-5], na.rm=TRUE)
          if(!is.na(median_rs) && median_rs > 1e-6){
            kernel_current <- kernel_current / median_rs
          }
        }

      # Ensure values below lowerLimit are zero
      kernel_current[kernel_current < lowerLimit & !is.na(kernel_current)] <- 0
      kernels_slide[[ct_i]][[ct_j]] <- kernel_current


      if(!sigma_valid_for_slide){
        sigma_valid_across_slides <- FALSE
        warning(paste("Sigma value", sigma_val, "invalid for slide", sID, "due to sparse kernels."))
      }
      kernels_for_sigma[[sID]] <- kernels_slide

    } # End slide loop (sID)

    if(sigma_valid_across_slides){
      kernel_matrices_all[[sigma_name]] <- kernels_for_sigma
    } else {
      kernel_matrices_all[[sigma_name]] <- NULL # Remove sigma if invalid for any slide
      warning(paste("Removing sigma value", sigma_val, "as it was invalid for one or more slides."))
    }

  } # End sigma loop (tt)

  # Filter out NULL entries (invalid sigmas)
  kernel_matrices_all <- kernel_matrices_all[!sapply(kernel_matrices_all, is.null)]
  object@sigmaValues <- as.numeric(gsub("sigma_", "", names(kernel_matrices_all))) # Update valid sigmas

  object@kernelMatrices <- kernel_matrices_all
  return(object)
})


# X_list_all: list(slideID = list(cellType = pc_matrix))
.scalePCMatsMulti <- function(X_list_all, pca_object, slides, cts){

  # Initialize pcaResults structure
  X_list_scaled <- setNames(vector("list", length = length(slides)), slides)
  for(sID in slides) {
    X_list_scaled[[sID]] <- setNames(vector("list", length = length(cts)), cts)
  }

  ct <- cts
    pca_A_sd <- pca_object[[ct]]$sdev
    for(sID in slides){
      X_list_scaled[[sID]][[ct]] <- scale(X_list_all[[sID]][[ct]],
                                          center = FALSE, scale = pca_A_sd)
    }
  return(X_list_scaled)
}

#' Run Multi-Slide Spatially Kernel Restricted CCA (skrCCA)
#'
#' Performs skrCCA on the `CoProm` object using data integrated across slides.
#' Calls the multi-slide optimization functions.
#'
#' @importFrom stats setNames
#' @param object A `CoProm` object with `pcaResults` and `kernelMatrices` populated.
#' @param nCC Number of canonical components to compute.
#' @param sigmaChoice A specific sigma value to use. If NULL (default), uses all valid sigma values from `object@sigmaValues`.
#' @param tol Tolerance for optimization convergence.
#' @param maxIter Maximum iterations for optimization.
#' @param n_cores Number of cores for parallel computation within optimization.
#' @param scalePCs Whether to scale each PC before computing skrCCA
#'
#' @return `CoProm` object with `skrCCAOut` slot populated:
#'         `list(sigma = list(cellType = weight_matrix))` (shared weights).
#' @export
#' @rdname runSkrCCAMultiOne
#' @aliases runSkrCCAMultiOne,CoProm-method
setGeneric("runSkrCCAMultiOne", function(
    object, nCC = 2, sigmaChoice = NULL, tol = 1e-5, scalePCs,
    maxIter = 200, n_cores = 1) standardGeneric("runSkrCCAMultiOne"))

#' @rdname runSkrCCAMultiOne
#' @aliases runSkrCCAMultiOne,CoProm-method
#' @importFrom stats setNames
#' @export
setMethod("runSkrCCAMultiOne", "CoProm", function(
    object, nCC = 2, sigmaChoice = NULL, tol = 1e-5, scalePCs,
    maxIter = 200, n_cores = 1) {

  # --- Input Checks ---
  if (length(object@pcaResults) == 0) stop("PCA results missing. Run computePCAMulti.")
  if (length(object@kernelMatrices) == 0) stop("Kernel matrices missing. Run computeKernelMatrixMulti.")
  cts <- object@cellTypesOfInterest
  if (length(cts) != 1) stop("Only one cellTypesOfInterest needed for skrCCA.")
  slides <- object@slideList

  # Determine which sigma values to use
  if (is.null(sigmaChoice)) {
    sigmas_to_run <- object@sigmaValues
    if (length(sigmas_to_run) == 0) stop("No valid sigma values found in object.")
  } else {
    if (!paste0("sigma_", sigmaChoice) %in% names(object@kernelMatrices)) {
      stop(paste("Chosen sigma value", sigmaChoice, "not found or was invalid."))
    }
    sigmas_to_run <- sigmaChoice
  }
  sigma_names_run <- paste0("sigma_", sigmas_to_run)

  # --- Prepare Input Lists for Optimization ---
  # X_list_all: list(slideID = list(cellType = pc_matrix))
  X_list_all <- object@pcaResults # Assumes correct structure from computePCAMulti

  if(scalePCs){
    X_list_all <- .scalePCMatsMulti(X_list_all = X_list_all,
                                    pca_object = object@pcaGlobal,
                                    slides = slides,
                                    cts = cts )
  }

  # Check consistency
  if(!all(slides %in% names(X_list_all))) stop("Slide mismatch in pcaResults")
  for(sID in slides){
    if(!all(cts %in% names(X_list_all[[sID]]))) stop(paste("Cell type mismatch in pcaResults for slide", sID))
  }


  # --- Run Optimization per Sigma ---
  cca_out_all_sigmas <- setNames(vector("list", length = length(sigmas_to_run)), sigma_names_run)

  for (idx in seq_along(sigmas_to_run)) {
    sig_val <- sigmas_to_run[idx]
    sig_name <- sigma_names_run[idx]
    message(paste("Running skrCCA for sigma =", sig_val))

    # K_list_all: list(slideID = list(cellType1 = list(cellType2 = kernel_matrix)))
    K_list_all_sigma <- object@kernelMatrices[[sig_name]]
    # Check consistency
    if(!all(slides %in% names(K_list_all_sigma))) stop(paste("Slide mismatch in kernelMatrices for sigma", sig_val))

    # Run optimization for the first component
    # Ensure optimize_bilinear_multi_slides is loaded/available
    cca_result_1 <- optimize_bilinear_multi_slides_one(
      X_list_all = X_list_all,
      K_list_all = K_list_all_sigma,
      max_iter = maxIter, tol = tol, n_cores = n_cores
    )
    # Name the result list using cell types
    names(cca_result_1) <- cts # Uses implicit order, assumes optimization returns unnamed list


    if (nCC == 1) {
      cca_out_all_sigmas[[sig_name]] <- cca_result_1
    } else {
      # Run optimization for subsequent components
      # Ensure optimize_bilinear_multi_n_slides_one is loaded/available
      cca_result_n <- optimize_bilinear_multi_n_slides_one(
        X_list_all = X_list_all,
        K_list_all = K_list_all_sigma,
        w_list = cca_result_1, # Pass named list from previous step
        nCC = nCC,
        max_iter = maxIter, tol = tol # n_cores not directly used here but in helper
      )
      # optimize_bilinear_multi_n_slides should return named list
      if(is.null(names(cca_result_n))) names(cca_result_n) <- cts # Ensure names

      cca_out_all_sigmas[[sig_name]] <- cca_result_n
    }
  } # End sigma loop

  object@skrCCAOut <- cca_out_all_sigmas
  object@nCC <- nCC
  return(object)
})


#' Compute Normalized Correlation for Multi-Slide Data one cell type
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
#' @rdname computeNormalizedCorrelationMultiOne
#' @aliases computeNormalizedCorrelationMultiOne,CoProm-method
setGeneric("computeNormalizedCorrelationMultiOne", function(
    object, tol = 1e-4, calculationMode = c("aggregate", "perSlide")
) standardGeneric("computeNormalizedCorrelationMultiOne"))

#' @rdname computeNormalizedCorrelationMultiOne
#' @aliases computeNormalizedCorrelationMultiOne,CoProm-method
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @export
setMethod("computeNormalizedCorrelationMultiOne", "CoProm", function(
    object, tol = 1e-4, calculationMode = c("aggregate", "perSlide")) {

  calculationMode <- match.arg(calculationMode)

  # --- Input Checks ---
  if (length(object@skrCCAOut) == 0) stop("skrCCA results missing. Run runSkrCCAMulti.")
  if (length(object@pcaResults) == 0) stop("PCA results missing.")
  if (length(object@kernelMatrices) == 0) stop("Kernel matrices missing.")
  cts <- object@cellTypesOfInterest
  if (length(cts) != 1) stop("Need one cell type.")
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
    # pair_cell_types <- combn(cts, 2)

    for (sID in slides) {
      norm_K_slide <- setNames(vector("list", length = length(cts)), cts)
      for (ct_i in cts) norm_K_slide[[ct_i]] <- setNames(vector("list", length=length(cts)), cts)

      K_list_slide <- K_list_sigma[[sID]]
      ct_i <- ct_j <- cts

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

        if(length(X_list_slide) == 0 || length(K_list_slide) == 0) next # Skip if data missing

        ct_i <- ct_j <- cts

          X_i <- X_j <- X_list_slide[[ct_j]]
          K_ij <- K_list_slide[[ct_i]][[ct_j]]
          norm_K_ij <- norm_K_slide[[ct_i]][[ct_j]]

          if(is.null(X_i) || is.null(K_ij) ||
             is.na(norm_K_ij) || norm_K_ij < 1e-9 ||
             nrow(X_i)==0 || nrow(X_j)==0) next # Skip if data missing/invalid

          for (cc in 1:nCC) {
            w_i <- w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]
            Xiw <- Xjw <- X_j %*% w_j

            numerator <- (t(Xiw) %*% K_ij %*% Xjw)[1, 1]
            denom_norm <- sum(Xiw^2) * norm_K_ij

            norm_corr_val <- ifelse(abs(denom_norm) < 1e-9, 0, numerator / denom_norm)

            df_slide <- rbind(df_slide, data.frame(
              sigmaValue = as.numeric(gsub("sigma_", "", sig_name)),
              slideID = sID,
              cellType1 = ct_i, cellType2 = ct_j,
              CC_index = cc, normalizedCorrelation = norm_corr_val,
              stringsAsFactors = FALSE
            ))
          } # end CC loop
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
#' @rdname computeGeneAndCellScoresMultiOne
#' @aliases computeGeneAndCellScoresMultiOne,CoProm-method
setGeneric("computeGeneAndCellScoresMultiOne", function(
    object, sigmaChoice = NULL) standardGeneric("computeGeneAndCellScoresMultiOne"))

#' @rdname computeGeneAndCellScoresMultiOne
#' @aliases computeGeneAndCellScoresMultiOne,CoProm-method
#' @importFrom stats setNames
#' @export
setMethod("computeGeneAndCellScoresMultiOne", "CoProm", function(
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
