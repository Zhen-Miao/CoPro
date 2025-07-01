

.getCellPermu <- function(object, permu_method, nPermu, cts,
                          num_bins_x = 10, num_bins_y = 10) {

  cell_permu <- setNames(vector("list", length = length(cts)), cts)

  if(permu_method == "global")  {
    for(i in cts){
      n_cell <- sum(object@cellTypesSub == i)
      if (i == cts[1]) {
        cell_permu[[i]] <- replicate(nPermu, 1:n_cell)
      }else {
        cell_permu[[i]] <- replicate(nPermu,
                                     sample.int(n = n_cell, replace = FALSE))
      }
    }
  }else if (permu_method == "bin") {

    location_full <- object@locationDataSub
    location_full$"cell_ID" <- rownames(location_full)
    location_full$x_bin <- cut(location_full$x, breaks = num_bins_x, labels = FALSE)
    location_full$y_bin <- cut(location_full$y, breaks = num_bins_y, labels = FALSE)

    for (i in cts){
      n_cell <- sum(object@cellTypesSub == i)
      if (i == cts[1]) {
        cell_permu[[i]] <- replicate(nPermu, 1:n_cell)
      }else {
        cell_loc <- location_full[object@cellTypesSub == i, ]
        cell_permu[[i]] <- matrix(ncol = nPermu, nrow = nrow(cell_loc))

        for (j in seq_len(nPermu)){
          cell_loc_resample <- resample_spatial(location_data = cell_loc,
                                                num_bins_x = num_bins_x,
                                                num_bins_y = num_bins_y)
          cell_permu[[i]][, j] <- match(cell_loc_resample$"cell_ID",
                                       cell_loc$"cell_ID")
        }

      }
    }
  }else{
    stop(paste("the function currently only support global resampling",
         "or bin-wise resampling"))
  }


  return(cell_permu)
}

#' runSkrCCAPermu
#' @importFrom stats setNames
#' @param object A `CoPro` object
#' @param tol Tolerance for termination, default = 1e-5
#' @param maxIter Maximum iterations
#' @param nPermu Number of permutation to run, default = 20
#' @param permu_method Method of permutation. Must be "global" or "bin".
#'  Default = "bin", which will shuffle the cells by bin.
#'  This will conserve local bin structures, so that local covariance
#'  will be maintained
#' @param num_bins_x Number of bins in x for bin-wise permutation, default = 10
#' @param num_bins_y Number of bins in y for bin-wise permutation, default = 10
#' @param verbose Whether to output the progress and related information
#'
#' @return CoPro object with distnace matrix computed
#' @export
#'
runSkrCCAPermu <- function(object, tol = 1e-5, nPermu = 20,
                           maxIter = 200, permu_method = "bin",
                           num_bins_x = 10, num_bins_y = 10, verbose = TRUE) {

  ## check input
  if (!is(object, "CoPro")) {
    stop("Input object must be a CoPro object")
  }

  if (!(permu_method %in% c("bin", "global"))) {
    stop("permu_method must be 'bin' or 'global'. ")
  }

  ## check nPermu input
  object@nPermu <- as.integer(nPermu)

  ## match the arguments from the runSkrCCA() function
  if (length(object@skrCCAOut) == 0) {
    stop("Please run runSkrCCA() before runSkrCCAPermu()")
  }
  nCC <- object@nCC
  scalePCs <- object@scalePCs

  ## choose cell types
  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    warning("no cell type of interest specified,
                      using all cell types to run the analysis")
    cts <- unique(object@cellTypesSub)
  }

  sigmaValueChoice <- object@sigmaValueChoice
  ## set sigmaValueChoice
  if (is.null(sigmaValueChoice)) {
    if (length(object@sigmaValueChoice) == 0) {
      stop(paste(
        "sigmaValueChoice is not given,",
        "and NormalizedCorrelation not computed,",
        "please either specify a particular sigmaValueChoice or",
        "run computeNormalizedCorrelation()"
      ))
    }else {
      warning(paste(
        "sigmaValueChoice is not given",
        "default set to the value with highest",
        "normalized correlation."
      ))
      sigmaValueChoice <- object@sigmaValueChoice
    }
  }

  if (!(sigmaValueChoice %in% object@sigmaValues)) {
    stop("sigmaValueChoice does not exist in the list of sigmaValues")
  }

  s_name <- paste("sigma", sigmaValueChoice, sep = "_")

  ## save output to a list
  cca_permu_out <- vector("list", length = nPermu)
  permu_names <- paste("permu", 1:nPermu, sep = "_")
  names(cca_permu_out) <- permu_names

  ## step 1. generate cell permutations
  cell_permu <- .getCellPermu(object = object, permu_method = permu_method,
                              nPermu = nPermu, cts = cts,
                              num_bins_x = num_bins_x, num_bins_y = num_bins_y)
  if (verbose) {
    cat("Cell permutation finished", "\n")
  }
  object@cellPermu <- cell_permu

  ## get PCA matrices and permute
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
  PCmats2 <- PCmats


  ## for loop to run the analysis
  for (tt in seq_len(nPermu)) {
    t <- permu_names[tt]

    for(i in names(PCmats)){
      PCmats2[[i]] <- PCmats[[i]][cell_permu[[i]][,tt],]
    }

    cca_result <- optimize_bilinear(
      X_list = PCmats2,
      K_list = object@kernelMatrices[[s_name]],
      max_iter = maxIter, tol = tol
    )
    names(cca_result) <- cts

    if (nCC == 1) {
      cca_permu_out[[t]] <- cca_result
    }else {
      cca_result_n <- optimize_bilinear_n(
        X_list = PCmats2,
        K_list = object@kernelMatrices[[s_name]],
        w_list = cca_result,
        cellTypesOfInterest = cts, nCC = nCC,
        max_iter = maxIter, tol = tol
      )
      cca_permu_out[[t]] <- cca_result_n
    }

  }

  object@skrCCAPermuOut <- cca_permu_out

  return(object)
}



#' Compute Normalized Correlation under permutation (approximation)
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
computeNormalizedCorrelationPermu <- function(object, tol = 1e-4) {

  ## check input
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  ## Check for required components
  if (length(object@skrCCAPermuOut) == 0) {
    stop(paste("skrCCAPermuOut is not available.",
               "Please run runSkrCCAPermu() first."))
  }

  ## choose cell types
  cts <- object@cellTypesOfInterest
  nPermu <- object@nPermu


  ## load whether the PCs are being scaled prior to CCA
  if (length(object@scalePCs) == 0) {
    stop("object@scalePCs not specified")
  }
  scalePCs <- object@scalePCs

  ## check sigmaValues
  sigmaValueChoice <- object@sigmaValueChoice
  PCmats <- .getAllPCMats(allPCs = object@pcaGlobal, scalePCs = scalePCs)
  nCC <- object@nCC

  pair_cell_types <- combn(cts, 2)

  ## save output to a list
  correlation_value <- vector("list", length = nPermu)
  permu_names <- paste("permu", 1:nPermu, sep = "_")
  names(correlation_value) <- permu_names
  s_name <- paste("sigma", sigmaValueChoice, sep = "_")

  ## calculate all spectral norms
  cat("Calculating spectral norms, ",
      "depending on the data size, this may take a while. \n")
  norm_K12 <- setNames(vector(mode = "list", length = 1),
                       s_name)

  norm_K12[[s_name]] <-
    setNames(vector(mode = "list", length = length(cts)), cts)
  for (i in cts) {
    norm_K12[[s_name]][[i]] <- setNames(
      vector(mode = "list", length = length(cts)), cts)
  }


  for (pp in seq_len(ncol(pair_cell_types))) {
    cellType1 <- pair_cell_types[1, pp]
    cellType2 <- pair_cell_types[2, pp]
    K <- getKernelMatrix(object, sigma = sigmaValueChoice, 
                         cellType1 = cellType1, cellType2 = cellType2, 
                         verbose = FALSE)
    ## Calculate the spectral norm of the kernel matrix
    svd_result <- irlba::irlba(K, nv = 1, tol = tol)
    norm_K12[[s_name]][[cellType1]][[cellType2]] <- svd_result$d[1]
  }

  cat("Finished calculating spectral norms \n")


  for (tt in seq_len(nPermu)) {
    t <- permu_names[tt]
    correlation_value[[t]] <- data.frame(
      sigmaValues = sigmaValueChoice,
      cellType1 = rep(pair_cell_types[1, ], times = nCC),
      cellType2 = rep(pair_cell_types[2, ], times = nCC),
      CC_index = rep(x = 1:nCC, each = ncol(pair_cell_types)),
      normalizedCorrelation = numeric(length = ncol(pair_cell_types) * nCC),
      stringsAsFactors = FALSE
    )
    for (pp in seq_len(ncol(pair_cell_types))) {
      for (cc_index in seq_len(nCC)) {
        cellType1 <- pair_cell_types[1, pp]
        cellType2 <- pair_cell_types[2, pp]

        w_1 <- object@skrCCAPermuOut[[t]][[cellType1]][, cc_index, drop = FALSE]
        w_2 <- object@skrCCAPermuOut[[t]][[cellType2]][, cc_index, drop = FALSE]

        A <- PCmats[[cellType1]][object@cellPermu[[cellType1]][, tt], ]
        B <- PCmats[[cellType2]][object@cellPermu[[cellType2]][, tt], ]

        A_w1 <- A %*% w_1
        B_w2 <- B %*% w_2

        ## get pre-calculated spectral norm
        K <- getKernelMatrix(object, sigma = sigmaValueChoice, 
                             cellType1 = cellType1, cellType2 = cellType2, 
                             verbose = FALSE)
        norm_K12_sel <- norm_K12[[s_name]][[cellType1]][[cellType2]]

        ## Calculate normalized correlation
        correlation_value[[t]]$"normalizedCorrelation"[
          pp + (cc_index - 1) * ncol(pair_cell_types)] <-
          (t(A_w1) %*% K %*% B_w2) /
          (sqrt(sum(A_w1^2)) * sqrt(sum(B_w2^2)) * norm_K12_sel)
      }
    }
  }

  ## Store the result in the object
  object@normalizedCorrelationPermu <- correlation_value

  ## Return the modified object
  return(object)
}
