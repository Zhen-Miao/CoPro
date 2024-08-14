## correlation plot

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
#' @export
getNormCorr <- function(object) {
  ## check input
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  ## make sure normalizedCorrelation exists
  if (length(object@normalizedCorrelation) == 0) {
    stop(paste(
      "normalizedCorrelation slot does not exist,",
      "run computeNormalizedCorrelation first"
    ))
  }

  ## organize into a data.frame
  ncorr <- do.call(rbind, slot(object, "normalizedCorrelation"))
  ncorr$"ct12" <- paste(ncorr$"cellType1", ncorr$"cellType2", sep = "-")
  ncorr$"sigmaSquares" <- factor(ncorr$"sigmaSquares",
    levels = sort(unique(ncorr$"sigmaSquares"),
      decreasing = FALSE
    )
  )

  return(ncorr)
}


#' Get cell score and location information as a data.frame
#'
#' @importFrom stats median
#' @param object A `CoPro` object
#' @param sigmaSquaredChoice A value to specify the sigma squared to
#' use for selecting the particular cell score information
#' @param scoreColorType Should the color be in binary scale or
#' continuous scale?
#'
#' @return A data.frame object with cell scores and their locations
#' @export
getCellScoresInSitu <- function(object, sigmaSquaredChoice,
                                scoreColorType = c("binary", "continuous")) {
  ## check input
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  ## match arg
  scoreColorType <- match.arg(scoreColorType)

  ## make sure normalizedCorrelation exists
  if (length(object@cellScores) == 0 |
    length(object@geneScores) == 0) {
    stop(paste(
      "cellScores slot does not exist,",
      "run computeNormalizedCorrelation first"
    ))
  }

  if (is.null(sigmaSquaredChoice)) {
    stop(paste(
      "sigmaSquaredChoice is not given",
      "default set to the value with highest",
      "normalized correlation."
    ))
    sigmaSquaredChoice <- object@sigmaSquaredChoice
  }

  if (!(sigmaSquaredChoice %in% object@sigmaSquares)) {
    stop("sigmaSquaredChoice does not exist in the list of sigmaSquares")
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

  sigma_name_choice <- paste("cellScore_sigma", sigmaSquaredChoice, sep = "_")

  loc_t <- stats::setNames(
    vector(mode = "list", length = length(cts)),
    cts
  )
  median_score_t <- vector("numeric", length = length(cts))
  names(median_score_t) <- cts

  for (t in cts) {
    loc_t[[t]] <- object@locationDataSub[object@cellTypesSub == t, ]
    loc_t[[t]]$cellScores <- object@cellScores[[t]][
      rownames(loc_t[[t]]),
      sigma_name_choice
    ]
    loc_t[[t]]$cellTypesSub <- t
    median_score_t[t] <- median(object@cellScores[[t]][, sigma_name_choice])
    loc_t[[t]]$cellScores_b <- ifelse(loc_t[[t]]$cellScores > median_score_t[t],
      paste0("high_", t), paste0("low_", t)
    )
  }

  combinations <- expand.grid(c("high", "low"), cts)
  all_binary_scores <- apply(combinations, 1, function(x) paste(x[1], x[2], sep = "_"))

  names(loc_t) <- NULL
  loc_all <- do.call(rbind, loc_t)[rownames(object@locationDataSub), ]

  loc_all$cellScores_b <- factor(loc_all$cellScores_b,
    levels = all_binary_scores
  )

  return(loc_all)
}



#' Assign distance matrix manually
#'
#' @param object A `CoPro` object
#'
#' @return A `CoPro` object with specified
#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
#'
#' @examples
setGeneric("assignDistanceManually",
           function(object,
                    distanceList) standardGeneric("assignDistanceManually")
)


#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
setMethod("assignDistanceManually", "CoPro",
          function(object, distanceList) {

            if(!is.list(distanceList)){
              stop(paste("distanceList must be a nested list object with names",
                         "specified by cell types"))
            }

            ## choose cell types
            if (length(object@cellTypesOfInterest) != 0) {
              cts <- object@cellTypesOfInterest
            } else {
              warning(paste("no cell type of interest specified,",
                            "using all cell types to run the analysis"))
              cts <- unique(object@cellTypesSub)
            }

            if(names(distanceList) != cts){
              stop(paste("The names of distanceList do not match cell types",
                   "of interest"))
            }

            for (i in cts) {
              if(names(distanceList[[i]]) != cts){
                stop(paste("The names of distanceList[[", i,
                           "]] do not match cell types ",
                           "of interest", sep = ""))
              }
            }

            object@distances <- distanceList
            return(object)
}
)
