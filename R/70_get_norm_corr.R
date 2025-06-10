## get nomralized correlation (all cells and any pair of cell types) for ploting

#' Get normalized correlation vs Sigma squared values
#'
#' Get a data.frame with normalized correlation vs Sigma squared values.
#' This helps to evaluate which sigma squared value to choose for downstream
#' analyses.
#' @importFrom methods is
#' @importFrom methods slot
#'
#' @param object A `CoPro` object
#'
#' @return A `data.frame` with correlation information
#' @rdname getNormalizedCorrelation
#' @aliases getNormCorr,CoProSingle-method
#' @aliases getNormCorr,CoProMulti-method
#' @export
setGeneric("getNormCorr",
           function(object) standardGeneric("getNormCorr")
)

.checkInputNormCorr <- function(object) {
  if ( !( is(object, "CoPro") || is(object, "CoProm") ) ) {
    stop("Input must be a CoPro or CoProm object")
  }
    ## make sure normalizedCorrelation exists
  if (length(object@normalizedCorrelation) == 0) {
    stop(paste(
      "normalizedCorrelation slot does not exist,",
      "run `computeNormalizedCorrelation()` first"
    ))
  }
  return(TRUE)
}

.getNormCorrCore <- function(object) {
  .checkInputNormCorr(object)
    ## organize into a data.frame
  ncorr <- do.call(rbind, slot(object, "normalizedCorrelation"))
  ncorr$"ct12" <- paste(ncorr$"cellType1", ncorr$"cellType2", sep = "-")
  ncorr$"sigmaValues" <- factor(ncorr$"sigmaValues",
    levels = sort(unique(ncorr$"sigmaValues"),
      decreasing = FALSE
    )
  )

  return(ncorr)
  
}


#' @rdname getNormalizedCorrelation
setMethod("getNormCorr", "CoProSingle", function(object) {
  .checkInputNormCorr(object)
  .getNormCorrCore(object = object)
})


#' @rdname getNormalizedCorrelation
setMethod("getNormCorr", "CoProMulti", function(object) {
  .checkInputNormCorr(object)
  .getNormCorrCore(object = object)
})




#' Retrieve the Correlation between two cell types
#'
#' @param object A `CoPro` object
#' @param cellTypeA Cell type label for one cell type
#' @param cellTypeB Cell type label for another cell type
#' @param sigmaValueChoice A particular sigma squared value
#' for the correlation
#' @param ccIndex Canonical vector index, default = 1
#'
#' @return A data.frame with two columns, AK and B, where AK represents the
#' cell score of cell type A times the kernel matrix, and B represents the
#' cell score of cell type B. If the object is a `CoProMulti` object, the
#' data.frame will have an additional column, slideID, indicating the slide ID.
#' @rdname getCorrTwoTypes
#' @aliases getCorrTwoTypes,CoProSingle-method
#' @aliases getCorrTwoTypes,CoProMulti-method
#' @export
setGeneric("getCorrTwoTypes", 
           function(object, cellTypeA, cellTypeB, ccIndex = 1, 
                    sigmaValueChoice = NULL) standardGeneric("getCorrTwoTypes"))



.getCorrTwoTypesSigma <- function(object, cellTypeA, cellTypeB, ccIndex = 1, sigmaValueChoice) {
      ## check input
  if (length(cellTypeA) != 1) {
    stop("Must give a single cellTypeA for correlation plot")
  }
  if (length(cellTypeB) != 1) {
    stop(paste("Must give a single cellTypeB for correlation plot.",
    "if there is only one cell type, please use getCorrOneType() instead."))
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

  ## check cell type A and B in cts
  if (!(cellTypeA %in% cts) || !(cellTypeB %in% cts)) {
    stop(paste(
      "cellTypeA or cellTypeB not in cellTypesSub,",
      "so the correlation plot cannot",
      "be generated."
    ))
  }
  ## make sure cellScores exists
  if (length(object@cellScores) == 0 ||
        length(object@geneScores) == 0) {
    stop(paste(
      "cellScores slot does not exist,",
      "run computeGeneAndCellScores first"
    ))
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


  return(sigmaValueChoice)

}


.getCorrTwoTypesCoreSingle <- function(object, cellTypeA, cellTypeB, ccIndex = 1,
                            sigmaValueChoice) {

  ## load the cellScores and kernel matrix
  sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")

  x1 <- t(object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = FALSE])
  x2 <- object@cellScores[[sigma_name]][[cellTypeB]][, ccIndex, drop = TRUE]
  ktemp <- object@kernelMatrices[[sigma_name]][[cellTypeA]][[cellTypeB]]
  if (length(ktemp) != 0) {
    k <- ktemp
  }else {
    k <- t(object@kernelMatrices[[sigma_name]][[cellTypeB]][[cellTypeA]])
  }

  df <- data.frame(AK = (x1 %*% k)[1, , drop = TRUE], B = x2)
  return(df)
}


.getCorrTwoTypesCoreMulti <- function(object, cellTypeA, cellTypeB, ccIndex = 1,
                            sigmaValueChoice) {

  ## load the cellScores and kernel matrix
  sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")
  df_q <- rep(list(), length = length(object@slideList))
  names(df_q) <- object@slideList

  for (q in object@slideList) {
    x1 <- t(object@cellScores[[sigma_name]][[q]][[cellTypeA]][, ccIndex, drop = FALSE])
    x2 <- object@cellScores[[sigma_name]][[q]][[cellTypeB]][, ccIndex, drop = TRUE]
    ktemp <- object@kernelMatrices[[sigma_name]][[q]][[cellTypeA]][[cellTypeB]]
    if (length(ktemp) != 0) {
      k <- ktemp
    }else {
      k <- t(object@kernelMatrices[[sigma_name]][[q]][[cellTypeB]][[cellTypeA]])
    }
    df_q[[q]] <- data.frame(AK = (x1 %*% k)[1, , drop = TRUE], B = x2, slideID = q)
  }
  df_q <- do.call(rbind, df_q)
  return(df_q)
}

#' @noRd
setMethod("getCorrTwoTypes", "CoProSingle", 
          function(object, cellTypeA, cellTypeB, ccIndex = 1, 
                   sigmaValueChoice = NULL) {
    sigmaValueChoice <- .getCorrTwoTypesSigma(object, cellTypeA, cellTypeB, 
                                               ccIndex = ccIndex, sigmaValueChoice)
    .getCorrTwoTypesCoreSingle(object, cellTypeA, cellTypeB, ccIndex, sigmaValueChoice)
  })

#' @noRd
setMethod("getCorrTwoTypes", "CoProMulti", 
          function(object, cellTypeA, cellTypeB, ccIndex = 1, 
                   sigmaValueChoice = NULL) {
    sigmaValueChoice <- .getCorrTwoTypesSigma(object, cellTypeA, cellTypeB, 
                                              ccIndex = ccIndex, sigmaValueChoice)
    .getCorrTwoTypesCoreMulti(object, cellTypeA, cellTypeB, ccIndex, sigmaValueChoice)
  })


#' Retrieve the Correlation within one cell types
#'
#' @param object A `CoPro` object
#' @param cellTypeA Cell type label for the cell type of interest
#' @param sigmaValueChoice A particular sigma squared value
#' for the correlation
#' @param ccIndex Canonical vector index, default = 1
#'
#' @return A data.frame with two columns, AK and B, where AK represents the
#' cell score of cell type A times the kernel matrix, and B represents the
#' cell score of cell type B.
#' @export
#' @rdname getCorrOneType
#' @aliases getCorrOneType,CoProSingle-method
#' @aliases getCorrOneType,CoProMulti-method
setGeneric("getCorrOneType", 
           function(object, cellTypeA, ccIndex = 1, 
                    sigmaValueChoice = NULL) standardGeneric("getCorrOneType"))


.getCorrOneTypeSigma <- function(object, cellTypeA, ccIndex = 1,
                            sigmaValueChoice) {
 ## check input
  if (length(cellTypeA) != 1) {
    stop("Must give a single cellTypeA for correlation plot")
  }

  ## choose cell types

  cts <- object@cellTypesOfInterest
  if (cellTypeA != cts){
    stop("cellTypeA must be in cellTypesSub.")
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

  ## make sure normalizedCorrelation exists
  if (length(object@cellScores) == 0 ||
      length(object@geneScores) == 0) {
    stop(paste(
      "cellScores slot does not exist,",
      "run computeGeneAndCellScores first"
    ))
  }

  return(sigmaValueChoice)
}






.getCorrOneTypeCoreSingle <- function(object, cellTypeA, ccIndex = 1,
                            sigmaValueChoice) {
  ## load the cellScores and kernel matrix
  sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")

  x1 <- t(object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = FALSE])
  x2 <- object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = TRUE]
  ktemp <- object@kernelMatrices[[sigma_name]][[cellTypeA]][[cellTypeA]]

  df <- data.frame(AK = (x1 %*% ktemp)[1, , drop = TRUE], B = x2)
  return(df)
}

.getCorrOneTypeCoreMulti <- function(object, cellTypeA, ccIndex = 1,
                            sigmaValueChoice) {
  df_q <- rep(list(), length = length(object@slideList))
  names(df_q) <- object@slideList
  sigma_name <- paste("sigma", sigmaValueChoice, sep = "_")
  for (q in object@slideList) {
    x1 <- t(object@cellScores[[sigma_name]][[q]][[cellTypeA]][, ccIndex, drop = FALSE])
    x2 <- object@cellScores[[sigma_name]][[q]][[cellTypeA]][, ccIndex, drop = TRUE]
    ktemp <- object@kernelMatrices[[sigma_name]][[q]][[cellTypeA]][[cellTypeA]]
    df_q[[q]] <- data.frame(AK = (x1 %*% ktemp)[1, , drop = TRUE], B = x2, slideID = q)
  }
  df_q <- do.call(rbind, df_q)
  return(df_q)
}


#' @noRd
setMethod("getCorrOneType", "CoProSingle", 
          function(object, cellTypeA, ccIndex = 1, 
                   sigmaValueChoice = NULL) {
    sigmaValueChoice <- .getCorrOneTypeSigma(object, cellTypeA, ccIndex, sigmaValueChoice)
    .getCorrOneTypeCoreSingle(object, cellTypeA, ccIndex, sigmaValueChoice)
  })

#' @noRd
setMethod("getCorrOneType", "CoProMulti", 
          function(object, cellTypeA, ccIndex = 1, 
                   sigmaValueChoice = NULL) {
    sigmaValueChoice <- .getCorrOneTypeSigma(object, cellTypeA, ccIndex, sigmaValueChoice)
    .getCorrOneTypeCoreMulti(object, cellTypeA, ccIndex, sigmaValueChoice)
  })