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

#' Prepare PC matrices with optional scaling
#' @param pc_data Either pcaGlobal (single slide) or pcaResults (multi-slide)
#' @param pca_global Global PCA objects for scaling (only needed if scalePCs = TRUE)
#' @param scalePCs Whether to scale PCs
#' @param slides Slide names (NULL for single slide)
#' @param cts Cell type names
#' @return Prepared PC matrices
#' @noRd
.preparePCMatrices <- function(pc_data, pca_global = NULL, scalePCs = TRUE, 
                              slides = NULL, cts = NULL) {
  
  if (is.null(slides)) {
    # Single slide case
    if (length(pc_data) == 0) {
      stop("PCA results do not exist, run computePCA() first.")
    }
    
    PCmats <- setNames(vector("list", length = length(pc_data)), names(pc_data))
    
    if (scalePCs) {
      for (i in names(pc_data)) {
        pca_A_sd <- pc_data[[i]]$sdev
        PCmats[[i]] <- scale(pc_data[[i]]$x, center = FALSE, scale = pca_A_sd)
      }
    } else {
      for (i in names(pc_data)) {
        PCmats[[i]] <- pc_data[[i]]$x
      }
    }
    return(PCmats)
    
  } else {
    # Multi-slide case
    if (length(pc_data) == 0) {
      stop("PCA results missing. Run computePCAMulti.")
    }
    
    if (scalePCs) {
      if (is.null(pca_global) || length(pca_global) == 0) {
        stop("Cannot scale PCs: pcaGlobal slot is empty. Ensure computePCAMulti was run successfully.")
      }
      
      missing_cts <- cts[!cts %in% names(pca_global)]
      if (length(missing_cts) > 0) {
        stop(paste("Cannot scale PCs: PCA objects missing for cell types:", 
                  paste(missing_cts, collapse = ", ")))
      }
      
      # Scale matrices
      X_list_scaled <- setNames(vector("list", length = length(slides)), slides)
      for (sID in slides) {
        X_list_scaled[[sID]] <- setNames(vector("list", length = length(cts)), cts)
      }
      
      for (ct in cts) {
        pca_ct_obj <- pca_global[[ct]]
        if (is.null(pca_ct_obj) || !"sdev" %in% names(pca_ct_obj)) {
          stop(paste("Invalid PCA object for cell type:", ct))
        }
        
        pca_A_sd <- pca_ct_obj$sdev
        if (is.null(pca_A_sd) || length(pca_A_sd) == 0) {
          stop(paste("Invalid sdev for cell type:", ct))
        }
        
        for (sID in slides) {
          if (!sID %in% names(pc_data) || !ct %in% names(pc_data[[sID]]) || 
              is.null(pc_data[[sID]][[ct]])) {
            stop(paste("Missing data for slide:", sID, "cell type:", ct))
          }
          
          X_list_scaled[[sID]][[ct]] <- scale(pc_data[[sID]][[ct]], 
                                             center = FALSE, scale = pca_A_sd)
        }
      }
      return(X_list_scaled)
    } else {
      # No scaling needed, return as-is but validate structure
      if (!all(slides %in% names(pc_data))) {
        stop("Slide mismatch in pcaResults")
      }
      for (sID in slides) {
        if (!all(cts %in% names(pc_data[[sID]]))) {
          stop(paste("Cell type mismatch in pcaResults for slide", sID))
        }
      }
      return(pc_data)
    }
  }
}

#' Unified input validation for skrCCA
#' @param object CoPro or CoProMulti object
#' @param scalePCs Whether to scale PCs
#' @param nCC Number of components
#' @param tol Tolerance
#' @param maxIter Maximum iterations
#' @param sigmaChoice Specific sigma choice (multi-slide only)
#' @param n_cores Number of cores (multi-slide only)
#' @return List with validated parameters
#' @noRd
.validateSkrCCAInputs <- function(object, scalePCs, nCC, tol, maxIter, 
                                 sigmaChoice = NULL, n_cores = 1) {
  
  # Check kernel matrices
  if (length(object@kernelMatrices) == 0) {
    stop("Kernel matrix is empty, please run computeKernelMatrix first")
  }
  
  # Record PC scaling preference
  if (length(object@scalePCs) == 0) {
    object@scalePCs <- scalePCs
  }
  
  # Check sigma values
  if (length(object@sigmaValues) == 0) {
    stop("sigmaValues is empty, please specify")
  }
  
  # Determine cell types of interest
  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    warning("no cell type of interest specified, using all cell types to run the analysis")
    cts <- unique(object@cellTypesSub)
  }
  
  # Handle multi-slide specific validation
  is_multi <- inherits(object, "CoProMulti")
  if (is_multi) {
    # Determine sigma values to run
    if (is.null(sigmaChoice)) {
      sigmas_to_run <- object@sigmaValues
      if (length(sigmas_to_run) == 0) {
        stop("No valid sigma values found in object.")
      }
    } else {
      if (!paste0("sigma_", sigmaChoice) %in% names(object@kernelMatrices)) {
        stop(paste("Chosen sigma value", sigmaChoice, "not found or was invalid."))
      }
      sigmas_to_run <- sigmaChoice
    }
    
    return(list(
      cts = cts, 
      sigmaValues = object@sigmaValues,
      sigmas_to_run = sigmas_to_run,
      n_cores = n_cores,
      is_multi = TRUE
    ))
  } else {
    return(list(
      cts = cts, 
      sigmaValues = object@sigmaValues,
      sigmas_to_run = object@sigmaValues,
      n_cores = 1,
      is_multi = FALSE
    ))
  }
}

#' Core skrCCA execution logic (unified for single and multi-slide)
#' @param object CoPro or CoProMulti object
#' @param validation_result Result from .validateSkrCCAInputs
#' @param scalePCs Whether to scale PCs
#' @param nCC Number of components
#' @param tol Tolerance
#' @param maxIter Maximum iterations
#' @return Updated object with skrCCA results
#' @noRd
.runSkrCCAUnified <- function(object, validation_result, scalePCs, nCC, tol, maxIter) {
  
  cts <- validation_result$cts
  sigmas_to_run <- validation_result$sigmas_to_run
  n_cores <- validation_result$n_cores
  is_multi <- validation_result$is_multi
  
  # Determine sigma names to run
  if (is_multi) {
    sigma_names_run <- paste0("sigma_", sigmas_to_run)
  } else {
    sigma_names_run <- paste("sigma", sigmas_to_run, sep = "_")
  }
  
  # Prepare PC matrices
  if (is_multi) {
    slides <- getSlideList(object)
    X_list_all <- .preparePCMatrices(
      pc_data = object@pcaResults,
      pca_global = object@pcaGlobal,
      scalePCs = scalePCs,
      slides = slides,
      cts = cts
    )
  } else {
    PCmats <- .preparePCMatrices(
      pc_data = object@pcaGlobal,
      scalePCs = scalePCs
    )
  }
  
  # Initialize output structure
  cca_out <- setNames(vector("list", length = length(sigmas_to_run)), sigma_names_run)
  
  # Run optimization for each sigma value
  for (idx in seq_along(sigmas_to_run)) {
    sig_val <- sigmas_to_run[idx]
    sig_name <- sigma_names_run[idx]
    
    if (is_multi) {
      message(paste("Running skrCCA for sigma =", sig_val))
      
      # Validate kernel matrix structure for multi-slide
      K_list_all_sigma <- object@kernelMatrices[[sig_name]]
      if (!all(slides %in% names(K_list_all_sigma))) {
        stop(paste("Slide mismatch in kernelMatrices for sigma", sig_val))
      }
      
      # Run multi-slide optimization
      cca_result_1 <- optimize_bilinear_multi_slides(
        X_list_all = X_list_all,
        K_list_all = K_list_all_sigma,
        max_iter = maxIter, 
        tol = tol, 
        n_cores = n_cores
      )
      
    } else {
      # Run single-slide optimization
      cca_result_1 <- optimize_bilinear(
        X_list = PCmats,
        K_list = object@kernelMatrices[[sig_name]],
        max_iter = maxIter, 
        tol = tol
      )
    }
    
    # Ensure proper naming
    names(cca_result_1) <- cts
    
    # Handle multiple components
    if (nCC == 1) {
      cca_out[[sig_name]] <- cca_result_1
    } else {
      if (is_multi) {
        cca_result_n <- optimize_bilinear_n_multi_slides(
          X_list_all = X_list_all,
          K_list_all = K_list_all_sigma,
          w_list = cca_result_1,
          cellTypesOfInterest = cts,
          nCC = nCC,
          max_iter = maxIter, 
          tol = tol,
          n_cores = n_cores
        )
      } else {
        cca_result_n <- optimize_bilinear_n(
          X_list = PCmats, 
          K_list = object@kernelMatrices[[sig_name]],
          w_list = cca_result_1,
          cellTypesOfInterest = cts, 
          nCC = nCC,
          max_iter = maxIter, 
          tol = tol
        )
      }
      
      # Ensure proper naming
      if (is.null(names(cca_result_n))) names(cca_result_n) <- cts
      cca_out[[sig_name]] <- cca_result_n
    }
  }
  
  # Update object and return
  object@skrCCAOut <- cca_out
  object@nCC <- nCC
  return(object)
}

#' runSkrCCA
#' @importFrom stats setNames
#' @param object A CoPro object
#' @param scalePCs Whether to scale each PCs to a uniform variance before
#' running the program
#' @param nCC Number of canonical vectors to compute, default = 2
#' @param tol Tolerance for termination, default = 1e-5
#' @param maxIter Maximum iterations
#' @param sigmaChoice Specific sigma value to use (CoProMulti only)
#' @param n_cores Number of cores for parallel processing (CoProMulti only)
#'
#' @return CoPro object with distance matrix computed
#' @export
#'
setGeneric(
  "runSkrCCA",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1) standardGeneric("runSkrCCA"))

#' @rdname runSkrCCA
#' @aliases runSkrCCA,CoPro-method
#' @export
setMethod(
  "runSkrCCA", "CoPro",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1) {
    
    validation_result <- .validateSkrCCAInputs(
      object, scalePCs, nCC, tol, maxIter, sigmaChoice, n_cores
    )
    
    .runSkrCCAUnified(object, validation_result, scalePCs, nCC, tol, maxIter)
  }
)

#' @rdname runSkrCCA
#' @aliases runSkrCCA,CoProMulti-method
#' @export
setMethod(
  "runSkrCCA", "CoProMulti",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1) {
    
    validation_result <- .validateSkrCCAInputs(
      object, scalePCs, nCC, tol, maxIter, sigmaChoice, n_cores
    )
    
    .runSkrCCAUnified(object, validation_result, scalePCs, nCC, tol, maxIter)
  }
)
