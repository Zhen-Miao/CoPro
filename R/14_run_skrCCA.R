#' Record a per-slide view of the global PC scores
#'
#' `pcaResults[[slide]][[ct]]` is, except when `center_per_slide = TRUE`,
#' exactly a set of rows of `pcaGlobal[[ct]]$x`. Storing the values again
#' doubles the memory the PC scores occupy for no new information. Store the
#' row indices and rebuild the slice on demand in `.resolvePCSlice()`.
#' @noRd
.newPCSlice <- function(rows) {
  list(type = "pcaSlice", rows = as.integer(rows))
}

#' Is this `pcaResults` entry a stored view rather than a materialized matrix?
#' @noRd
.isPCSlice <- function(entry) {
  is.list(entry) && identical(entry$type, "pcaSlice")
}

#' Materialize one `pcaResults` entry
#'
#' Accepts either representation. A matrix is returned untouched, which covers
#' both the `center_per_slide = TRUE` case -- where the slice is re-centered
#' and so is genuinely different data -- and objects saved before slices became
#' views.
#' @param entry One `pcaResults[[slide]][[ct]]` element.
#' @param scores The cell-by-PC matrix the indices refer to, normally
#'   `pca_global[[ct]]$x` or a scaled copy of it.
#' @return A cell-by-PC matrix.
#' @noRd
.resolvePCSlice <- function(entry, scores) {
  if (is.null(entry)) return(NULL)
  if (!is.list(entry) || !identical(entry$type, "pcaSlice")) return(entry)
  if (is.null(scores)) {
    stop("Cannot resolve per-slide PC scores: the global PCA is missing.")
  }
  scores[entry$rows, , drop = FALSE]
}

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
      stop("PCA results missing. Run computePCA.")
    }

    
    if (scalePCs) {
      if (is.null(pca_global) || length(pca_global) == 0) {
        stop(paste("Cannot scale PCs: pcaGlobal slot is empty.",
                   "Ensure computePCA was run successfully."))
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

        # Scale the global scores once per cell type, then take slide views of
        # the result, instead of calling scale() once per (slide, cell type).
        # Slices that are stored as matrices -- center_per_slide, or a legacy
        # object -- are scaled individually as before.
        scaled_global <- scale(pca_ct_obj$x, center = FALSE, scale = pca_A_sd)

        for (sID in slides) {
          if (!sID %in% names(pc_data) || !ct %in% names(pc_data[[sID]]) ||
              is.null(pc_data[[sID]][[ct]])) {
            stop(paste("Missing data for slide:", sID, "cell type:", ct))
          }

          entry <- pc_data[[sID]][[ct]]
          X_list_scaled[[sID]][[ct]] <- if (.isPCSlice(entry)) {
            .resolvePCSlice(entry, scaled_global)
          } else {
            scale(entry, center = FALSE, scale = pca_A_sd)
          }
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
      resolved <- pc_data
      for (sID in slides) {
        for (ct in cts) {
          resolved[[sID]][[ct]] <- .resolvePCSlice(
            pc_data[[sID]][[ct]], pca_global[[ct]]$x
          )
        }
      }
      return(resolved)
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
                                 step_size = 1, objective = "sumcov",
                                 slideWeight = NULL,
                                 minCellsPerSlide = .min_cells_per_slide) {

  # check if scalePCs is a logical value
  if (!is.logical(scalePCs) || length(scalePCs) != 1) {
    stop("scalePCs must be a logical value")
  }

  objective <- match.arg(objective, c("sumcov", "sumcor"))

  # slideWeight only means something when a per-slide denominator exists.
  # Silently ignoring it under sumcov would let a user believe they had asked
  # for a weighting that never took effect.
  if (identical(objective, "sumcov") && !is.null(slideWeight)) {
    stop("slideWeight applies only to objective = \"sumcor\". Under ",
         "\"sumcov\" each slide is already weighted by its own cell count and ",
         "score scale, which is what sumcor exists to change.")
  }
  if (identical(objective, "sumcor")) {
    slideWeight <- match.arg(
      if (is.null(slideWeight)) "size" else slideWeight,
      c("size", "equal")
    )
  }

  if (!is.numeric(minCellsPerSlide) || length(minCellsPerSlide) != 1 ||
      minCellsPerSlide < 0) {
    stop("minCellsPerSlide must be a single non-negative number")
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
    
    # Across slides the criteria genuinely differ. With a single slide they
    # usually do not -- see `.sumcorSingleSlideIsExact()` for exactly when --
    # so one-slide objects go down the exact SUMCOV solvers rather than
    # iterating to the same answer.
    n_slides <- length(getSlideList(object))
    use_sumcor <- identical(objective, "sumcor") && n_slides > 1
    if (identical(objective, "sumcor") && n_slides <= 1) {
      .reportSingleSlideSumcor(object, cts, slideWeight, "Only one slide present")
    }

    return(list(
      cts = cts,
      sigmaValues = object@sigmaValues,
      sigmas_to_run = sigmas_to_run,
      n_cores = n_cores,
      is_multi = TRUE,
      objective = objective,
      slideWeight = slideWeight,
      minCellsPerSlide = minCellsPerSlide,
      use_sumcor = use_sumcor
    ))
  } else {
    # For single-slide, warn if multi-slide parameters are used
    if (!is.null(sigmaChoice)) {
      warning("sigmaChoice parameter is ignored for single-slide CoPro objects")
    }
    if (n_cores != 1) {
      warning("n_cores parameter is ignored for single-slide CoPro objects")
    }
    
    # Single-slide object: see the multi-slide branch for when sumcor reduces
    # to sumcov here, and when it only approximately does.
    if (identical(objective, "sumcor")) {
      .reportSingleSlideSumcor(object, cts, slideWeight, "Single-slide object")
    }

    return(list(
      cts = cts,
      sigmaValues = object@sigmaValues,
      sigmas_to_run = object@sigmaValues,
      n_cores = 1,
      is_multi = FALSE,
      objective = objective,
      slideWeight = slideWeight,
      minCellsPerSlide = minCellsPerSlide,
      use_sumcor = FALSE
    ))
  }
}

#' Does one-slide `sumcor` reduce to `sumcov` exactly?
#'
#' `.sumcorSigma()` is the norm \eqn{\|X_i w_i\|}, not a root-mean-square, so
#' with whitened PCs and \eqn{\|w_i\| = 1} the denominators are
#' \eqn{\sqrt{n_i - 1}} rather than 1. A one-slide SUMCOR objective is
#' therefore SUMCOV reweighted by the per-pair constant
#' \eqn{m_{ij} / \sqrt{(n_i-1)(n_j-1)}}, and a per-pair constant leaves the
#' argmax alone only when every pair gets the same one -- so with at most one
#' pair (one or two cell types), or with equal cell counts.
#'
#' Only ever called with a single slide, so the flat `cellTypesSub` vector is
#' that slide's labels for both object classes.
#' @noRd
.sumcorSingleSlideIsExact <- function(object, cts, slideWeight) {
  if (identical(slideWeight, "covariance")) return(TRUE)
  if (length(cts) <= 2L) return(TRUE)
  counts <- as.numeric(table(factor(object@cellTypesSub, levels = cts)))
  isTRUE(all.equal(min(counts), max(counts)))
}

#' Say what a one-slide `sumcor` request is actually going to compute
#'
#' Exact in the common cases; an approximation with three or more cell types at
#' unequal counts, which is nearly every real dataset. How loudly to say so
#' depends on the weighting, because that is what sets the size of the gap. The
#' per-pair constant is \eqn{m_{ij} / \sqrt{(n_i-1)(n_j-1)}}: under `"size"`,
#' \eqn{m_{ij} = \sqrt{n_i n_j}} nearly cancels it and the mismatch is
#' \eqn{1 + O(1/n)} -- immaterial, so a message. Under `"equal"` nothing cancels
#' and the gap can be large, so that one warns.
#' @noRd
.reportSingleSlideSumcor <- function(object, cts, slideWeight, prefix) {
  if (.sumcorSingleSlideIsExact(object, cts, slideWeight)) {
    message(prefix, ": sumcor and sumcov are the same optimization problem ",
            "here, so the exact sumcov solvers are used.")
    return(invisible(NULL))
  }

  detail <- paste0(
    "with ", length(cts), " cell types at unequal cell counts, sumcor and ",
    "sumcov do not share a maximizer on one slide -- the denominators are ",
    "sqrt(n_i - 1), not 1, so the criteria differ by a per-pair constant. ",
    "The exact sumcov solvers are used anyway."
  )

  if (identical(slideWeight, "equal")) {
    warning(prefix, ": ", detail, " Under slideWeight = \"equal\" nothing ",
            "cancels that constant and the difference can be material; call ",
            "optimize_sumcor_pca() directly to optimize the sumcor criterion ",
            "itself, or use slideWeight = \"size\".", call. = FALSE)
  } else {
    message(prefix, ": ", detail, " Under slideWeight = \"size\" the ",
            "cell-count factor nearly cancels it, so the mismatch is ",
            "1 + O(1/n) and immaterial in practice.")
  }
  invisible(NULL)
}

#' Apply the per-slide adequacy rule, or report it
#'
#' Only `sumcor` divides by a per-slide quantity, so only `sumcor` needs thin
#' slides removed: a slide whose weight lands near the null space of its Gram
#' matrix drives `sigma` to the floor and `1/sigma` explodes. Under `sumcov` a
#' thin slide merely contributes little to the summed operator, so it is
#' reported and kept -- dropping it would silently change results computed
#' before this rule existed.
#'
#' @param data_matrices Structure from `.prepareDataMatrices()`.
#' @param cts Cell types of interest.
#' @param validation Structure from `.validateSkrCCAInputs()`.
#' @param nPCA Feature count, used for the rank-deficiency warning.
#' @return `data_matrices`, with slides filtered when `use_sumcor` is `TRUE`,
#'   and `use_sumcor` possibly downgraded if fewer than two slides survive.
#' @noRd
.applySlideAdequacy <- function(data_matrices, cts, validation, nPCA) {
  if (!isTRUE(validation$is_multi)) return(data_matrices)

  counts <- lapply(data_matrices$slides, function(s) {
    vapply(cts, function(ct) {
      X <- data_matrices$X_list_all[[s]][[ct]]
      if (is.null(X)) 0L else nrow(X)
    }, integer(1))
  })
  names(counts) <- data_matrices$slides
  thin <- data_matrices$slides[vapply(
    data_matrices$slides,
    function(s) any(counts[[s]] < validation$minCellsPerSlide),
    logical(1)
  )]

  if (!isTRUE(validation$use_sumcor)) {
    if (length(thin) > 0) {
      message(sprintf(
        paste0("  Note: slide(s) %s have fewer than %d cells for some cell ",
               "type. Under objective = \"sumcov\" they are kept (they simply ",
               "contribute little); objective = \"sumcor\" would drop them."),
        paste(thin, collapse = ", "), validation$minCellsPerSlide
      ))
    }
    return(data_matrices)
  }

  kept <- .dropDegenerateSlides(
    data_matrices$X_list_all, cts,
    minCells = validation$minCellsPerSlide,
    nFeatures = nPCA, what = "runSkrCCA (sumcor)"
  )
  data_matrices$X_list_all <- kept$X_list_all
  data_matrices$slides <- kept$slides
  data_matrices$droppedSlides <- kept$dropped

  if (length(kept$slides) < 2L) {
    warning("runSkrCCA (sumcor): only ", length(kept$slides),
            " slide(s) survived the per-slide cell threshold, and with one ",
            "slide sumcor and sumcov are the same problem. Falling back to the ",
            "exact sumcov solvers.", call. = FALSE)
    data_matrices$use_sumcor <- FALSE
  }
  data_matrices
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
                                       maxIter, tol, n_cores, step_size = 1,
                                       use_sumcor = FALSE,
                                       slideWeight = "size") {

  tryCatch({
    # SUMCOR route. The per-slide operators are built once per sigma and reused
    # for every axis; each is nPC x nPC, so the whole iteration afterwards is
    # small dense algebra with no kernel products in the loop.
    #
    # The exact one- and two-cell-type solvers below do NOT apply here: SUMCOR's
    # operator depends on w through sigma, so there is no fixed matrix to
    # decompose. That is also why the algebraic identity "stacking slides equals
    # summing their operators" stops holding under SUMCOR.
    if (isTRUE(use_sumcor)) {
      ops <- .computeSlideOperators(
        data_matrices$X_list_all, object@kernelMatrices, sig_val,
        data_matrices$slides, cts, n_cores
      )

      if (is.null(transferred_weight_1)) {
        w_first <- optimize_sumcor_pca(
          X_list_all = data_matrices$X_list_all,
          flat_kernels = object@kernelMatrices, sigma = sig_val,
          slides = data_matrices$slides, cell_types = cts,
          slideWeight = slideWeight, sdev2_list = data_matrices$sdev2_list,
          max_iter = maxIter, tol = tol, n_cores = n_cores, ops = ops
        )
      } else {
        w_first <- transferred_weight_1
      }

      # Strip the solver's reporting attributes before they become the seed of
      # the multi-axis matrices.
      w_first <- setNames(
        lapply(cts, function(ct) w_first[[ct]][, 1L, drop = FALSE]), cts
      )
      if (nCC == 1L) return(w_first)

      return(optimize_sumcor_pca_n(
        X_list_all = data_matrices$X_list_all,
        flat_kernels = object@kernelMatrices, sigma = sig_val,
        slides = data_matrices$slides, cell_types = cts,
        w_list = w_first, nCC = nCC, slideWeight = slideWeight,
        sdev2_list = data_matrices$sdev2_list,
        max_iter = maxIter, tol = tol, n_cores = n_cores, ops = ops
      ))
    }

    # The ordinary one- and two-cell-type problems have exact all-axis direct
    # solutions. Form the small PC-space operator once per sigma. Keep the
    # sequential route when the first axis was externally transferred, because
    # later axes must then be conditioned on that supplied direction.
    if (length(cts) <= 2L && is.null(transferred_weight_1)) {
      if (is_multi) {
        Y_resi <- compute_Y_multi_slide(
          data_matrices$X_list_all,
          object@kernelMatrices,
          sig_val,
          data_matrices$slides,
          cts,
          n_cores
        )
      } else {
        Y_resi <- compute_Y_resi(
          data_matrices$PCmats,
          object@kernelMatrices,
          sig_val,
          cts,
          slide = NULL
        )
      }
      if (length(cts) == 1L) {
        return(solve_one_type_eigen(
          Y_resi, cts, nCC = nCC,
          sdev2_list = data_matrices$sdev2_list
        ))
      }
      return(solve_two_type_svd(
        Y_resi, cts, nCC = nCC,
        sdev2_list = data_matrices$sdev2_list
      ))
    }

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

  # Thin slides are dropped under sumcor and only reported under sumcov.
  data_matrices <- .applySlideAdequacy(
    data_matrices, cts, validation_result, object@nPCA
  )
  use_sumcor <- if (!is.null(data_matrices$use_sumcor)) {
    isTRUE(data_matrices$use_sumcor)
  } else {
    isTRUE(validation_result$use_sumcor)
  }
  slideWeight <- validation_result$slideWeight
  if (use_sumcor) {
    message(sprintf(
      "Objective: sumcor (per-slide self-normalized), slideWeight = \"%s\", %d slides.",
      slideWeight, length(data_matrices$slides)
    ))
  }

  # Initialize output structure
  cca_out <- setNames(vector("list", length = length(sigmas_to_run)), sigma_names_run)
  
  # Run optimization for each sigma value
  n_sigmas <- length(sigmas_to_run)
  t_start <- Sys.time()
  for (idx in seq_along(sigmas_to_run)) {
    sig_val <- sigmas_to_run[idx]
    sig_name <- sigma_names_run[idx]
    message(sprintf(
      "Running skrCCA [%d/%d] for sigma = %g ...",
      idx, n_sigmas, sig_val
    ))

    # Run optimization for this sigma
    cca_out[[sig_name]] <- .runSingleSigmaOptimization(
      object, sig_val, sig_name, data_matrices, transferred_weight_1,
      is_multi, cts, nCC, maxIter, tol, n_cores, step_size,
      use_sumcor = use_sumcor, slideWeight = slideWeight
    )
  }
  message(sprintf(
    "skrCCA finished %d sigma value(s) in %.1f s.",
    n_sigmas, as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  ))
  
  # Process results and generate summary
  .processOptimizationResults(cca_out, sigma_names_run)

  # Record what was optimized alongside the weights. An attribute rather than a
  # new S4 slot so objects saved before this existed still load; the reader
  # defaults to "sumcov", which is what those objects were computed under. This
  # mirrors how computeNormalizedCorrelation() records its resolved normalizer.
  attr(cca_out, "ccaObjective") <- list(
    space = "pca",
    objective = if (use_sumcor) "sumcor" else "sumcov",
    requested = validation_result$objective,
    slideWeight = if (use_sumcor) slideWeight else NA_character_,
    sweep = NA_character_,
    slides = if (is_multi) data_matrices$slides else NA_character_,
    droppedSlides = if (is.null(data_matrices$droppedSlides)) {
      character(0)
    } else {
      data_matrices$droppedSlides
    }
  )

  # Update object and return
  object@skrCCAOut <- cca_out
  object@nCC <- nCC
  return(object)
}

#' What objective produced an object's CCA weights
#'
#' Reads the provenance record `runSkrCCA()` attaches to `@skrCCAOut`. Objects
#' computed before the `objective` argument existed carry no record and are
#' reported as `"sumcov"`, which is what they were computed under.
#'
#' @param object A `CoPro` or `CoProMulti` object.
#' @return A list with `space`, `objective`, `requested`, `slideWeight`,
#'   `slides` and `droppedSlides`.
#' @family scores-and-correlation
#' @seealso [runSkrCCA()]
#' @examples
#' \donttest{
#' toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
#' obj <- newCoProSingle(
#'   normalizedData = toy$normalizedData,
#'   locationData   = toy$locationData,
#'   metaData       = toy$metaData,
#'   cellTypes      = toy$cellTypes
#' )
#' obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
#' obj <- computePCA(obj, nPCA = 10)
#' obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE)
#' obj <- runSkrCCA(obj, nCC = 2)
#' getCCAObjective(obj)
#' }
#' @export
getCCAObjective <- function(object) {
  record <- attr(object@skrCCAOut, "ccaObjective")
  if (is.null(record)) {
    return(list(
      space = if (length(grep("^gscca_", names(object@skrCCAOut))) > 0) {
        "gene"
      } else {
        "pca"
      },
      objective = "sumcov",
      requested = NA_character_,
      slideWeight = NA_character_,
      sweep = NA_character_,
      slides = NA_character_,
      droppedSlides = character(0)
    ))
  }
  record
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
#' @param space Which feature space to optimize in. `"pca"` (default) runs the
#'   PC-space optimizer described here. `"gene"` forwards to
#'   [runGeneSpaceCCA()], which needs a single `sigma` -- supply it through
#'   `sigmaChoice` -- and accepts its own arguments (`clip`, `min_prevalence`,
#'   `streaming`, ...) through `...`.
#'
#'   Under `space = "gene"`, `objective`, `tol`, `maxIter` and
#'   `minCellsPerSlide` are forwarded **only when you supply them**, because the
#'   two entry points have different defaults (`objective` `"sumcov"` vs
#'   `"sumcor"`, `tol` `1e-5` vs `1e-6`, `maxIter` `200` vs `3000`,
#'   `minCellsPerSlide` `10` vs `20`). Forwarding a `runSkrCCA()` default would
#'   silently make this call differ from the [runGeneSpaceCCA()] call it claims
#'   to be, so an unset argument keeps gene space's own default. Arguments with
#'   no gene-space meaning (`scalePCs`, `transferred_weight_1`, `n_cores`,
#'   `step_size`, `slideWeight`) are an error rather than silently dropped.
#' @param objective Which canonical criterion to maximize.
#'
#'   `"sumcov"` (default) maximizes the sum of kernel-smoothed cross-covariances
#'   \eqn{\sum_{i<j} w_i' (\sum_s X_i^{(s)'} K_{ij}^{(s)} X_j^{(s)}) w_j} under
#'   \eqn{\|w_i\| = 1}.
#'
#'   `"sumcor"` maximizes the per-slide self-normalized sum, dividing each
#'   slide's cross term by that slide's own score scales.
#'
#'   **With one slide these are usually the same problem**, and `"sumcor"` is
#'   routed to the exact `"sumcov"` solvers. Whitened PCs give
#'   \eqn{X_i'X_i = (n_i-1) I}, so on \eqn{\|w_i\| = 1} the denominators are
#'   \eqn{\sigma_i = \sqrt{n_i - 1}} -- note \eqn{\sigma} is a norm, not a
#'   root-mean-square -- and the objective is SUMCOV reweighted by the
#'   *per-pair* constant \eqn{m_{ij} / \sqrt{(n_i-1)(n_j-1)}}. A per-pair
#'   constant leaves the maximizer alone only when every pair gets the same
#'   one, so the reduction is exact for **one or two cell types at any cell
#'   counts**, and for **three or more only when the counts are equal**.
#'   Outside that, one-slide `"sumcor"` still uses the `"sumcov"` solvers but
#'   warns: the mismatch is \eqn{1 + O(1/n)} under the default
#'   `slideWeight = "size"` and can be material under `"equal"`. Call
#'   [optimize_sumcor_pca()] directly to optimize the criterion itself.
#'
#'   Across slides they always differ. SUMCOV factors exactly as
#'   \eqn{\sum_s \sigma_i^{(s)} \sigma_j^{(s)} \rho_{ij}^{(s)}} -- with no
#'   \eqn{\sqrt{n_i n_j}} factor, since \eqn{\sigma_i \sigma_j \rho_{ij}} is
#'   the SUMCOV term already -- so it already sums per-slide correlations,
#'   weighted by per-slide score scale. That scale factor is what lets a slide
#'   with inflated variance along the canonical direction dominate; `"sumcor"`
#'   removes it, and `slideWeight = "size"` reintroduces the cell-count factor
#'   \eqn{\sqrt{n_i n_j}} on its own.
#' @param slideWeight Per-slide weighting, only valid with
#'   `objective = "sumcor"` (an error otherwise, since under `"sumcov"` the
#'   weighting is fixed by the objective). `"size"` (default) weights slide `s`
#'   by \eqn{\sqrt{n_i^{(s)} n_j^{(s)}}}, so larger slides count for more
#'   without per-slide variance re-entering. `"equal"` weights every slide the
#'   same -- strict Kettenring SUMCOR, matching [runGeneSpaceCCA()].
#' @param minCellsPerSlide Minimum cells per (slide, cell type). Slides below
#'   this are **dropped** under `objective = "sumcor"`, which divides by a
#'   per-slide scale that a near-empty slide drives to its floor. Under
#'   `"sumcov"` they are only reported, not dropped: a thin slide simply
#'   contributes little to the summed operator, and dropping it would change
#'   results computed before this rule existed.
#' @param ... Passed to [runGeneSpaceCCA()] when `space = "gene"`; ignored
#'   otherwise.
#'
#' @section Batch effects:
#' `objective = "sumcor"` removes the per-slide *scale* sensitivity but not a
#' per-slide *mean shift*. PC scores are centered globally
#' (`computePCA(center_per_slide = FALSE)`), so a shared technical shift leaves
#' \eqn{u_i^{(s)} \approx M_i^{(s)} \mathbf{1} + \epsilon} and the numerator
#' picks up \eqn{M_i M_j \mathbf{1}'K\mathbf{1}}, positive whenever both cell
#' types shift the same way. Dividing by `sigma` does not rescue this: for a
#' nonnegative kernel the leading singular pair is close to the Perron vector,
#' so a constant score reaches \eqn{\rho \approx \sigma_{max}(K)} on every
#' slide. Pair multi-slide `"sumcor"` with
#' `computePCA(..., center_per_slide = TRUE)`, which is the half of the fix that
#' addresses the mean shift. This applies to `"sumcov"` too and is independent
#' of the objective choice.
#'
#' @return CoPro object with skrCCA results computed. The objective actually
#'   used is recorded on `@skrCCAOut` and readable with [getCCAObjective()].
#' @family spatial-pipeline
#' @seealso [computePCA()], [computeKernelMatrix()],
#'   [computeNormalizedCorrelation()], [computeGeneAndCellScores()],
#'   [runGeneSpaceCCA()], [getCCAObjective()]
#' @examples
#' \donttest{
#' toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
#' obj <- newCoProSingle(
#'   normalizedData = toy$normalizedData,
#'   locationData   = toy$locationData,
#'   metaData       = toy$metaData,
#'   cellTypes      = toy$cellTypes
#' )
#' obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
#' obj <- computePCA(obj, nPCA = 10)
#' obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
#' obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
#' obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
#' }
#' @export
#'
setGeneric(
  "runSkrCCA",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           transferred_weight_1 = NULL,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1,
           step_size = 1, space = c("pca", "gene"),
           objective = c("sumcov", "sumcor"), slideWeight = NULL,
           minCellsPerSlide = 10, ...) standardGeneric("runSkrCCA"))

#' Forward a `space = "gene"` request to the gene-space implementation
#'
#' Kept in one place so both methods dispatch identically. Gene space needs a
#' single bandwidth rather than the sigma grid the PC-space route sweeps, so
#' `sigmaChoice` supplies it; a one-element `@sigmaValues` is accepted as an
#' unambiguous default.
#' @noRd
.dispatchGeneSpace <- function(object, sigmaChoice, nCC, objective, given,
                               tol, maxIter, minCellsPerSlide, ...) {
  sigma <- sigmaChoice
  if (is.null(sigma)) {
    if (length(object@sigmaValues) == 1L) {
      sigma <- object@sigmaValues
    } else {
      stop("space = \"gene\" analyses one bandwidth at a time. Pass it as ",
           "sigmaChoice; the object currently holds ",
           length(object@sigmaValues), " sigma values (",
           paste(object@sigmaValues, collapse = ", "), ").")
    }
  }

  # These are named formals of runSkrCCA(), so they never reach `...`. Passing
  # one and having it quietly vanish is the dangerous case -- a dropped
  # `transferred_weight_1` looks exactly like a transfer that ran.
  inapplicable <- intersect(
    given,
    c("scalePCs", "transferred_weight_1", "n_cores", "step_size", "slideWeight")
  )
  if (length(inapplicable) > 0L) {
    stop("space = \"gene\" does not use ",
         paste(sprintf("`%s`", inapplicable), collapse = ", "),
         ". Gene space z-standardizes per (slide, cell type) and weights slides ",
         "equally, so there is no PC parametrization, transfer, or slide ",
         "weighting to set here. Drop ",
         if (length(inapplicable) > 1L) "them" else "it",
         " or use space = \"pca\".", call. = FALSE)
  }

  # Forward the ones that do have a gene-space analogue, under its names -- but
  # only when actually supplied, since the two entry points have different
  # defaults (tol 1e-5 vs 1e-6, maxIter 200 vs 3000, minCells 10 vs 20) and
  # forwarding a runSkrCCA default would silently change the gene-space result.
  fwd <- list()
  if ("tol" %in% given) fwd$tol <- tol
  if ("maxIter" %in% given) fwd$max_iter <- maxIter
  if ("minCellsPerSlide" %in% given) fwd$min_cells <- minCellsPerSlide

  # `objective` is the same trap: runSkrCCA() defaults to "sumcov" and
  # runGeneSpaceCCA() to "sumcor", so forwarding the resolved default would make
  # runSkrCCA(space = "gene") optimize a different criterion than the
  # runGeneSpaceCCA() call it claims to forward to. Only an explicit choice
  # travels; otherwise gene space keeps its own default.
  if ("objective" %in% given) fwd$objective <- objective

  do.call(
    runGeneSpaceCCA,
    c(list(object, sigma = sigma, nCC = nCC), fwd, list(...))
  )
}

#' @rdname runSkrCCA
#' @aliases runSkrCCA,CoPro-method
#' @export
setMethod(
  "runSkrCCA", "CoPro",
  function(object, scalePCs = TRUE, nCC = 2, tol = 1e-5,
           transferred_weight_1 = NULL,
           maxIter = 200, sigmaChoice = NULL, n_cores = 1,
           step_size = 1, space = c("pca", "gene"),
           objective = c("sumcov", "sumcor"), slideWeight = NULL,
           minCellsPerSlide = 10, ...) {

    space <- match.arg(space)
    objective <- match.arg(objective)
    if (identical(space, "gene")) {
      return(.dispatchGeneSpace(
        object, sigmaChoice, nCC, objective,
        given = names(match.call())[-1],
        tol = tol, maxIter = maxIter,
        minCellsPerSlide = minCellsPerSlide, ...
      ))
    }

    # validate inputs
    validation_result <- .validateSkrCCAInputs(
      object, scalePCs, nCC, tol, maxIter, sigmaChoice, n_cores, step_size,
      objective = objective, slideWeight = slideWeight,
      minCellsPerSlide = minCellsPerSlide
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
           step_size = 1, space = c("pca", "gene"),
           objective = c("sumcov", "sumcor"), slideWeight = NULL,
           minCellsPerSlide = 10, ...) {

    space <- match.arg(space)
    objective <- match.arg(objective)
    if (identical(space, "gene")) {
      return(.dispatchGeneSpace(
        object, sigmaChoice, nCC, objective,
        given = names(match.call())[-1],
        tol = tol, maxIter = maxIter,
        minCellsPerSlide = minCellsPerSlide, ...
      ))
    }

    # validate inputs
    validation_result <- .validateSkrCCAInputs(
      object, scalePCs, nCC, tol, maxIter, sigmaChoice, n_cores, step_size,
      objective = objective, slideWeight = slideWeight,
      minCellsPerSlide = minCellsPerSlide
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
