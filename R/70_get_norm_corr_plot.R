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
  
  # Check if normalizedCorrelation slot is properly populated
  normCorr <- slot(object, "normalizedCorrelation")
  if (length(normCorr) == 0) {
    stop("normalizedCorrelation slot is empty. Please run computeNormalizedCorrelation() first.")
  }
  
  # Check if the list contains valid data frames
  if (!all(sapply(normCorr, is.data.frame))) {
    stop("normalizedCorrelation slot contains invalid data. Please run computeNormalizedCorrelation() again.")
  }
  
  # Check if any data frames are empty
  if (any(sapply(normCorr, nrow) == 0)) {
    stop("normalizedCorrelation slot contains empty data frames. Please run computeNormalizedCorrelation() again.")
  }
  
  ## organize into a data.frame
  ncorr <- do.call(rbind, normCorr)
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

  # Get the cell scores data for cellTypeA
  cell_score_data_a <- object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = FALSE]
  
  # Ensure it's a matrix before transposing
  if (!is.matrix(cell_score_data_a)) {
    cell_score_data_a <- as.matrix(cell_score_data_a)
  }
  
  x1 <- t(cell_score_data_a)
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
    # Get the cell scores data for cellTypeA - now using aggregated structure
    cell_score_data_a <- object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = FALSE]
    
    # Ensure it's a matrix before transposing
    if (!is.matrix(cell_score_data_a)) {
      cell_score_data_a <- as.matrix(cell_score_data_a)
    }
    
    # Filter for cells from this specific slide
    slide_cells_a <- rownames(object@metaDataSub)[object@metaDataSub$slideID == q & object@cellTypesSub == cellTypeA]
    slide_cells_b <- rownames(object@metaDataSub)[object@metaDataSub$slideID == q & object@cellTypesSub == cellTypeB]
    
    if (length(slide_cells_a) > 0 && length(slide_cells_b) > 0) {
      # Extract data for this slide only
      cell_score_data_a_slide <- cell_score_data_a[slide_cells_a, , drop = FALSE]
      cell_score_data_b_slide <- object@cellScores[[sigma_name]][[cellTypeB]][slide_cells_b, ccIndex, drop = TRUE]
      
      x1 <- t(cell_score_data_a_slide)
      x2 <- cell_score_data_b_slide
      ktemp <- object@kernelMatrices[[sigma_name]][[q]][[cellTypeA]][[cellTypeB]]
      if (length(ktemp) != 0) {
        k <- ktemp
      } else {
        k <- t(object@kernelMatrices[[sigma_name]][[q]][[cellTypeB]][[cellTypeA]])
      }
      df_q[[q]] <- data.frame(AK = (x1 %*% k)[1, , drop = TRUE], B = x2, slideID = q)
    } else {
      # Create empty data frame if no cells in this slide
      df_q[[q]] <- data.frame(AK = numeric(0), B = numeric(0), slideID = character(0))
    }
  }
  df_q <- do.call(rbind, df_q)
  
  # Reorder the data frame to match the original order in metaDataSub
  # Get the cell IDs for the specified cell types in the original order
  cell_ids_meta <- rownames(object@metaDataSub)[object@cellTypesSub %in% c(cellTypeA, cellTypeB)]
  
  # Create a mapping from the combined data to the original order
  # We need to match the rownames of the cell scores to the cell IDs in metaDataSub
  all_cell_ids <- c()
  for (q in object@slideList) {
    # Get cell IDs for cellTypeA and cellTypeB from this slide
    slide_cell_ids_a <- rownames(object@metaDataSub)[object@metaDataSub$slideID == q & object@cellTypesSub == cellTypeA]
    slide_cell_ids_b <- rownames(object@metaDataSub)[object@metaDataSub$slideID == q & object@cellTypesSub == cellTypeB]
    all_cell_ids <- c(all_cell_ids, slide_cell_ids_a, slide_cell_ids_b)
  }
  
  # Create a data frame with the correct order
  # For multi-slide data, we need to match the cell IDs from metaDataSub
  # to the order they appear in the combined data frame
  if (length(all_cell_ids) == nrow(df_q)) {
    # Add rownames to the data frame for ordering
    rownames(df_q) <- all_cell_ids
    
    # Reorder based on the original order in metaDataSub
    common_cells <- intersect(cell_ids_meta, all_cell_ids)
    if (length(common_cells) > 0) {
      df_q <- df_q[common_cells, , drop = FALSE]
      rownames(df_q) <- NULL  # Remove rownames from final output
    }
  }
  
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
  if (length(cts) == 0) {
    warning(paste(
      "no cell type of interest specified,",
      "using all cell types to run the analysis"
    ))
    cts <- unique(object@cellTypesSub)
  }
  
  if (!(cellTypeA %in% cts)) {
    stop(paste(
      "cellTypeA not in cellTypesOfInterest,",
      "so the correlation plot cannot be generated."
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

  ## make sure cellScores exists
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

  # Get the cell scores data
  cell_score_data <- object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = FALSE]
  
  # Ensure it's a matrix before transposing
  if (!is.matrix(cell_score_data)) {
    cell_score_data <- as.matrix(cell_score_data)
  }
  
  x1 <- t(cell_score_data)
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
    # Get the cell scores data for this slide - now using aggregated structure
    slide_cells <- rownames(object@metaDataSub)[object@metaDataSub$slideID == q & object@cellTypesSub == cellTypeA]
    
    if (length(slide_cells) > 0) {
      # Extract data for this slide only
      cell_score_data <- object@cellScores[[sigma_name]][[cellTypeA]][slide_cells, ccIndex, drop = FALSE]
      
      # Ensure it's a matrix before transposing
      if (!is.matrix(cell_score_data)) {
        cell_score_data <- as.matrix(cell_score_data)
      }
      
      x1 <- t(cell_score_data)
      x2 <- cell_score_data[,1, drop = TRUE]
      ktemp <- object@kernelMatrices[[sigma_name]][[q]][[cellTypeA]][[cellTypeA]]
      df_q[[q]] <- data.frame(AK = (x1 %*% ktemp)[1, , drop = TRUE], B = x2, slideID = q)
    } else {
      # Create empty data frame if no cells in this slide
      df_q[[q]] <- data.frame(AK = numeric(0), B = numeric(0), slideID = character(0))
    }
  }
  names(df_q) <- NULL
  df_q <- do.call(rbind, df_q)
  df_q$slideID <- as.factor(df_q$slideID)
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