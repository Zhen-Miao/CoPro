#' Prepare PC matrices with optional scaling
#' @note The structure of pc_data and pca_global has been messy,
#'  but for now, let us make them consistent. pca_global is the direct output 
#'  of irlba::prcomp_irlba, and pcaResults is only used for multi-slide, where
#'  the pca matrix in each slide is saved here for convenience.
#' @param pc_data only used for mulit-slide, which is saved in pcaResults
#' @param pca_global Global PCA objects for scaling, which is the direct output 
#' of irlba::prcomp_irlba, saved in pcaGlobal slot.
#' @param scalePCs Whether to scale PCs
#' @param slides Slide names (NULL for single slide)
#' @param cts Cell type names
#' @return Prepared PC matrices
#' @noRd
.preparePCMatrices <- function(pc_data = NULL, pca_global, scalePCs = TRUE, 
                              slides = NULL, cts) {
  
  if (is.null(slides)) {
    # Single slide case
    if (length(pca_global) == 0) {
      stop("PCA results do not exist, run computePCA() first.")
    }
    # check if all cell types are in pcaGlobal
    if(!all(cts %in% names(pca_global))) {
      stop(paste("Cell types not found in pcaGlobal:",
       paste(cts[!cts %in% names(pca_global)], collapse = ", ")))
    }
    
    PCmats <- setNames(vector("list", length = length(cts)), cts)
    
    if (scalePCs) {
      for (ct in cts) {
        pca_A_sd <- pca_global[[ct]]$sdev
        PCmats[[ct]] <- scale(pca_global[[ct]]$x, center = FALSE, scale = pca_A_sd)
      }
    } else {
      for (ct in cts) {
        PCmats[[ct]] <- pca_global[[ct]]$x
      }
    }
    return(PCmats)
    
  } else {
    # Multi-slide case
    # pc_data is pcaResults
    if (length(pc_data) == 0) {
      stop("PCA results missing. Run computePCAMulti.")
    }

    
    if (scalePCs) {
      if (is.null(pca_global) || length(pca_global) == 0) {
        stop(paste("Cannot scale PCs: pcaGlobal slot is empty.",
                   "Ensure computePCAMulti was run successfully."))
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
                                 sigmaChoice = NULL, n_cores = 1,
                                 step_size = 1) {
  
  # check if scalePCs is a logical value
  if (!is.logical(scalePCs) || length(scalePCs) != 1) {
    stop("scalePCs must be a logical value")
  }
  
  # Validate numeric parameters
  if (!is.numeric(nCC) || length(nCC) != 1 || nCC <= 0 || nCC != as.integer(nCC)) {
    stop("nCC must be a positive integer")
  }
  
  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("tol must be a positive number")
  }
  
  if (!is.numeric(maxIter) || length(maxIter) != 1 || maxIter <= 0 || maxIter != as.integer(maxIter)) {
    stop("maxIter must be a positive integer")
  }
  
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0 || n_cores != as.integer(n_cores)) {
    stop("n_cores must be a positive integer")
  }

  if (!is.numeric(step_size) || length(step_size) != 1 || step_size <= 0 || step_size > 1) {
    stop("step_size must be a single numeric value in (0, 1]")
  }

  # Check kernel matrices
  if (length(object@kernelMatrices) == 0) {
    stop("Kernel matrix is empty, please run computeKernelMatrix first")
  }
  
  # Check sigma values
  if (length(object@sigmaValues) == 0) {
    stop("sigmaValues is empty, please specify")
  }
  
  # Determine cell types of interest
  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    warning(paste("no cell type of interest specified,",
                  "using all cell types to run the analysis"))
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
      # Check if any kernel matrix name contains this sigma value
      sigma_prefix <- paste0("sigma", sigmaChoice)
      has_sigma <- any(grepl(paste0("(^|\\|)", sigma_prefix, "(\\||$)"),
                             names(object@kernelMatrices)))
      if (!has_sigma) {
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
    # For single-slide, warn if multi-slide parameters are used
    if (!is.null(sigmaChoice)) {
      warning("sigmaChoice parameter is ignored for single-slide CoPro objects")
    }
    if (n_cores != 1) {
      warning("n_cores parameter is ignored for single-slide CoPro objects")
    }
    
    return(list(
      cts = cts, 
      sigmaValues = object@sigmaValues,
      sigmas_to_run = object@sigmaValues,
      n_cores = 1,
      is_multi = FALSE
    ))
  }
}



.check_input_transferred_weight_1 <- function(transferred_weight_1, cts, nPCA){
  ## make sure it is a list object
  if(any(names(transferred_weight_1) != cts)){
    stop("transferred_weight_1 must be a list object with cell type names")
  }
  ## check within each element
  for(ct in cts){
    tw <- transferred_weight_1[[ct]]
    if(!is.matrix(tw)){
      stop(paste0("transferred_weight_1[[", ct, "]] must be a matrix"))
    }
    if(ncol(tw) != 1){
      stop(paste0("transferred_weight_1[[", ct, "]] must have one column"))
    }
    if(nrow(tw) != nPCA){
      stop(paste0("transferred_weight_1[[", ct, "]] must have nrow=nPCA"))
    }
    tws <- sum(as.vector(tw)^2)
    if(abs(tws-1) > 1e-3){
      stop(paste0("transferred_weight_1[[", ct, "]][,1] must be a unit vector"))
    }

  }

  return(0)
}

#' Prepare data matrices for skrCCA optimization
#' @param object CoPro or CoProMulti object
#' @param is_multi Whether this is multi-slide
#' @param scalePCs Whether to scale PCs
#' @param cts Cell types of interest
#' @return List containing prepared matrices (X_list_all for multi-slide, PCmats for single-slide)
#' @noRd
.prepareDataMatrices <- function(object, is_multi, scalePCs, cts) {
  # When scalePCs = FALSE, build sdev2_list so the optimizer can enforce
  # the correct CCA constraint w'(X'X)w = 1 via weighted normalization.
  # When scalePCs = TRUE, X is already whitened so sdev2_list = NULL
  # and the standard ||w|| = 1 constraint is equivalent.
  if (!scalePCs) {
    sdev2_list <- setNames(
      lapply(cts, function(ct) object@pcaGlobal[[ct]]$sdev^2),
      cts
    )
  } else {
    sdev2_list <- NULL
  }

  if (is_multi) {
    slides <- getSlideList(object)
    X_list_all <- .preparePCMatrices(
      pc_data = object@pcaResults,
      pca_global = object@pcaGlobal,
      scalePCs = scalePCs,
      slides = slides,
      cts = cts
    )
    return(list(X_list_all = X_list_all, slides = slides, sdev2_list = sdev2_list))
  } else {
    PCmats <- .preparePCMatrices(
      pca_global = object@pcaGlobal,
      scalePCs = scalePCs,
      cts = cts
    )
    return(list(PCmats = PCmats, sdev2_list = sdev2_list))
  }
}

#' Run optimization for a single sigma value
#' @param object CoPro or CoProMulti object
#' @param sig_val Sigma value
#' @param sig_name Sigma name
#' @param data_matrices Prepared data matrices
#' @param transferred_weight_1 Transferred weights or NULL
#' @param is_multi Whether multi-slide
#' @param cts Cell types of interest
#' @param nCC Number of components
#' @param maxIter Maximum iterations
#' @param tol Tolerance
#' @param n_cores Number of cores
#' @param step_size Step size for damped power iteration
#' @return Optimization result or NULL if failed
#' @noRd
.runSingleSigmaOptimization <- function(object, sig_val, sig_name, data_matrices,
                                       transferred_weight_1, is_multi, cts, nCC,
                                       maxIter, tol, n_cores, step_size = 1) {
  
  tryCatch({
    # Get first component
    cca_result_1 <- .getFirstComponent(
      object, sig_name, data_matrices, transferred_weight_1,
      is_multi, cts, maxIter, tol, n_cores, step_size
    )

    # Handle multiple components if needed
    if (nCC == 1) {
      return(cca_result_1)
    } else {
      return(.getSubsequentComponents(
        object, sig_name, data_matrices, cca_result_1,
        is_multi, cts, nCC, maxIter, tol, n_cores, step_size
      ))
    }
    
  }, error = function(e) {
    # Log the error with full context
    full_error_msg <- paste(conditionMessage(e), 
                   if (!is.null(e$call)) paste("in", deparse(e$call)[1]))
    warning(paste("Optimization failed for sigma =",
                 sig_val, ":", full_error_msg), 
            immediate. = TRUE)
    return(NULL)
  })
}



#' Get first component for optimization
#' @param object CoPro or CoProMulti object
#' @param sig_name Sigma name
#' @param data_matrices Prepared data matrices
#' @param transferred_weight_1 Transferred weights or NULL
#' @param is_multi Whether multi-slide
#' @param cts Cell types of interest
#' @param maxIter Maximum iterations
#' @param tol Tolerance
#' @param n_cores Number of cores
#' @param step_size Step size for damped power iteration
#' @return First component result
#' @noRd
.getFirstComponent <- function(object, sig_name, data_matrices, transferred_weight_1,
                              is_multi, cts, maxIter, tol, n_cores,
                              step_size = 1) {
  
  # Use transferred weights if provided
  if (!is.null(transferred_weight_1)) {
    return(transferred_weight_1)
  }
  
  # Extract sigma value from sigma name
  sigma_val <- as.numeric(gsub("sigma_", "", sig_name))
  
  if (is_multi) {
    # Run multi-slide optimization with flat kernels
    slides <- data_matrices$slides
    return(optimize_bilinear_multi_slides(
      X_list_all = data_matrices$X_list_all,
      flat_kernels = object@kernelMatrices,
      sigma = sigma_val,
      slides = slides,
      max_iter = maxIter,
      tol = tol,
      n_cores = n_cores,
      step_size = step_size,
      sdev2_list = data_matrices$sdev2_list
    ))

  } else {
    # Run single-slide optimization with flat kernels
    return(optimize_bilinear(
      X_list = data_matrices$PCmats,
      flat_kernels = object@kernelMatrices,
      sigma = sigma_val,
      max_iter = maxIter,
      tol = tol,
      step_size = step_size,
      sdev2_list = data_matrices$sdev2_list
    ))
  }
}

#' Get subsequent components for optimization
#' @param object CoPro or CoProMulti object
#' @param sig_name Sigma name
#' @param data_matrices Prepared data matrices
#' @param cca_result_1 First component result
#' @param is_multi Whether multi-slide
#' @param cts Cell types of interest
#' @param nCC Number of components
#' @param maxIter Maximum iterations
#' @param tol Tolerance
#' @param n_cores Number of cores
#' @param step_size Step size for damped power iteration
#' @return All components result
#' @noRd
.getSubsequentComponents <- function(object, sig_name, data_matrices, cca_result_1,
                                    is_multi, cts, nCC, maxIter, tol, n_cores,
                                    step_size = 1) {
  
  # Extract sigma value from sigma name
  sigma_val <- as.numeric(gsub("sigma_", "", sig_name))
  
  if (is_multi) {
    # Run multi-slide optimization with flat kernels
    slides <- data_matrices$slides
    cca_result_n <- optimize_bilinear_n_multi_slides(
      X_list_all = data_matrices$X_list_all,
      flat_kernels = object@kernelMatrices,
      sigma = sigma_val,
      slides = slides,
      w_list = cca_result_1,
      cellTypesOfInterest = cts,
      nCC = nCC,
      max_iter = maxIter,
      tol = tol,
      n_cores = n_cores,
      step_size = step_size,
      sdev2_list = data_matrices$sdev2_list
    )
  } else {
    # Run single-slide optimization with flat kernels
    cca_result_n <- optimize_bilinear_n(
      X_list = data_matrices$PCmats,
      flat_kernels = object@kernelMatrices,
      sigma = sigma_val,
      w_list = cca_result_1,
      cellTypesOfInterest = cts,
      nCC = nCC,
      max_iter = maxIter,
      tol = tol,
      step_size = step_size,
      sdev2_list = data_matrices$sdev2_list
    )
  }
  
  # Ensure proper naming
  if (is.null(names(cca_result_n))) names(cca_result_n) <- cts
  return(cca_result_n)
}

#' Process optimization results and generate summary
#' @param cca_out Results from all sigma optimizations
#' @param sigma_names_run Names of sigma values that were run
#' @return List with successful and failed sigma information
#' @noRd
.processOptimizationResults <- function(cca_out, sigma_names_run) {
  
  successful_sigmas <- character(0)
  failed_sigmas <- character(0)
  
  for (sig_name in sigma_names_run) {
    if (is.null(cca_out[[sig_name]])) {
      failed_sigmas <- c(failed_sigmas, sig_name)
    } else {
      successful_sigmas <- c(successful_sigmas, sig_name)
    }
  }
  
  # Summary reporting
  if (length(failed_sigmas) > 0) {
    warning(paste("Optimization failed for", length(failed_sigmas), "sigma value(s):", 
                 paste(failed_sigmas, collapse = ", ")), immediate. = TRUE)
  }
  
  if (length(successful_sigmas) > 0) {
    message(paste("Optimization succeeded for", length(successful_sigmas), "sigma value(s):", 
                 paste(successful_sigmas, collapse = ", ")))
  } else {
    stop("All optimization attempts failed. Please check your data and parameters.")
  }
  
  return(list(successful = successful_sigmas, failed = failed_sigmas))
}

#' Core skrCCA execution logic (unified for single and multi-slide)
#' @param object CoPro or CoProMulti object
#' @param validation_result Result from .validateSkrCCAInputs
#' @param scalePCs Whether to scale PCs
#' @param nCC Number of components
#' @param tol Tolerance
#' @param transferred_weight_1 If we use cross-slide weight transfer function,
#'  the transferred weight on each PC. Otherwise, the value should be set to NULL.
#' @param maxIter Maximum iterations
#' @return Updated object with skrCCA results
#' @noRd
.runSkrCCAUnified <- function(object, validation_result,
 scalePCs, nCC, tol, transferred_weight_1, maxIter, step_size = 1) {
  
  # Extract parameters from validation result
  cts <- validation_result$cts
  sigmas_to_run <- validation_result$sigmas_to_run
  n_cores <- validation_result$n_cores
  is_multi <- validation_result$is_multi
  
  # Determine sigma names to run
  sigma_names_run <- paste("sigma", sigmas_to_run, sep = "_")
  
  # Prepare data matrices
  data_matrices <- .prepareDataMatrices(object, is_multi, scalePCs, cts)
  
  # Initialize output structure
  cca_out <- setNames(vector("list", length = length(sigmas_to_run)), sigma_names_run)
  
  # Run optimization for each sigma value
  for (idx in seq_along(sigmas_to_run)) {
    sig_val <- sigmas_to_run[idx]
    sig_name <- sigma_names_run[idx]
    message(paste("Running skrCCA for sigma =", sig_val))
    
    # Run optimization for this sigma
    cca_out[[sig_name]] <- .runSingleSigmaOptimization(
      object, sig_val, sig_name, data_matrices, transferred_weight_1,
      is_multi, cts, nCC, maxIter, tol, n_cores, step_size
    )
  }
  
  # Process results and generate summary
  .processOptimizationResults(cca_out, sigma_names_run)
  
  # Update object and return
  object@skrCCAOut <- cca_out
  object@nCC <- nCC
  return(object)
}

#' runSkrCCA
#' @importFrom stats setNames
#' @importFrom parallel mclapply
#' @param object A CoPro object
#' @param scalePCs Whether to scale each PCs to a uniform variance before
#' running the program
#' @param nCC Number of canonical vectors to compute, default = 2
#' @param tol Tolerance for termination, default = 1e-5
#' @param transferred_weight_1 If we use cross-slide weight transfer function,
#'  the transferred weight on each PC. Otherwise, the value should be set to NULL.
#' @param maxIter Maximum iterations
#' @param sigmaChoice Specific sigma value to use (CoProMulti only, ignored for CoPro)
#' @param n_cores Number of cores for parallel processing (CoProMulti only, ignored for CoPro)
#' @param step_size Step size for damped power iteration. Default 1 (standard
#'   power iteration). Values in (0,1) blend old and new weights for smoother
#'   convergence, which can help with many cells or many CCs.
#'
#' @return CoPro object with skrCCA results computed
#' @export
#'
setGeneric(
  "runSkrCCA",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           transferred_weight_1 = NULL,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1,
           step_size = 1) standardGeneric("runSkrCCA"))

#' @rdname runSkrCCA
#' @aliases runSkrCCA,CoPro-method
#' @export
setMethod(
  "runSkrCCA", "CoPro",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           transferred_weight_1 = NULL,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1,
           step_size = 1) {

    # validate inputs
    validation_result <- .validateSkrCCAInputs(
      object, scalePCs, nCC, tol, maxIter, sigmaChoice, n_cores, step_size
    )
    # validate transferred_weight_1
    if(!is.null(transferred_weight_1)){
      .check_input_transferred_weight_1(transferred_weight_1,
      cts = validation_result$cts, nPCA = object@nPCA)
    }
    # run skrCCA
    object <- .runSkrCCAUnified(object, validation_result, scalePCs,
     nCC, tol, transferred_weight_1, maxIter, step_size)
    object@scalePCs = scalePCs
    return(object)
  }
)

#' @rdname runSkrCCA
#' @aliases runSkrCCA,CoProMulti-method
#' @export
setMethod(
  "runSkrCCA", "CoProMulti",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
          transferred_weight_1 = NULL,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1,
           step_size = 1) {

    # validate inputs
    validation_result <- .validateSkrCCAInputs(
      object, scalePCs, nCC, tol, maxIter, sigmaChoice, n_cores, step_size
    )
    # validate transferred_weight_1
    if(!is.null(transferred_weight_1)){
      .check_input_transferred_weight_1(transferred_weight_1,
      cts = validation_result$cts, nPCA = object@nPCA)
    }
    # run skrCCA
    object <- .runSkrCCAUnified(object, validation_result, scalePCs,
     nCC, tol, transferred_weight_1, maxIter, step_size)
    object@scalePCs = scalePCs
    return(object)
  }
)
