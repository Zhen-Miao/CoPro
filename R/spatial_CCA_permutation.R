#' runSkrCCAPermu
#' @importFrom stats setNames
#' @param object A CoPro object
#' @param tol Tolerance for termination, default = 1e-5
#' @param maxIter Maximum iterations
#' @param nPermu Number of permutation to run
#'
#' @return CoPro object with distnace matrix computed
#' @export
#'
setGeneric(
  "runSkrCCAPermu",
  function(object, tol = 1e-5, nPermu,
           maxIter = 200) standardGeneric("runSkrCCAPermu"))


#' @rdname runSkrCCAPermu
#' @aliases runSkrCCAPermu,CoPro-method
#' @importFrom stats setNames
#' @export
setMethod(
  "runSkrCCAPermu", "CoPro",
  function(object, tol = 1e-5, nPermu,
           maxIter = 200) {

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


    ## set sigmaValueChoice
    if (is.null(sigmaValueChoice)) {
      if (length(object@sigmaValueChoice) == 0) {
        stop(paste(
          "sigmaValueChoice is not given,",
          "and NormalizedCorrelation not computed,",
          "please either specify a particular sigmaValueChoice or",
          "run computeNormalizedCorrelation()"
        ))
      }else{
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

    sigmaValueChoice_name <- paste("sigma", sigmaValueChoice, sep = "_")

    ## save output to a list
    cca_permu_out <- vector("list", length = nPermu)
    permu_names <- paste("permu", 1:nPermu, sep = "_")
    names(cca_permu_out) <- permu_names

    ## step 1. generate cell permutations
    cell_permu <- setNames(vector("list", length = length(cts)), cts)
    for(i in cts){
      n_cell <- sum(object@cellTypesSub == i)
      cell_permu[[i]] <- replicate(nPermu,
                                   sample.int(n = n_cell, replace = FALSE))
    }

    object@cellPermu <- cell_permu

    ## get PCA matrices and permute
    PCmats <- .getAllPCMats(allPCs = object@pcaResults, scalePCs = scalePCs)
    PCmats2 <- PCmats


    ## for loop to run the analysis
    for (tt in seq_len(nPermu)) {
      t <- permu_names[tt]

      for(i in names(PCmats)){
        PCmats2[[i]] <- PCmats[[i]][cell_permu[[i]][,tt],]
      }

      cca_result <- optimize_bilinear_multi(
        X_list = PCmats2,
        K_list = object@kernelMatrices[[sigmaValueChoice_name]],
        max_iter = maxIter, tol = tol
      )
      names(cca_result) <- cts

      if (nCC == 1) {
        cca_permu_out[[t]] <- cca_result
      }else {
        cca_result_n <- optimize_bilinear_multi_n(
          X_list = PCmats2,
          K_list = object@kernelMatrices[[sigmaValueChoice_name]],
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
)

