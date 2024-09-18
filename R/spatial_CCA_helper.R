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
      "run `computeNormalizedCorrelation()` first"
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
#' continuous scale? Need to be either "binary" or "continuous"
#' @param ccIndex Canonical vector index, default = 1
#'
#' @return A data.frame object with cell scores and their locations
#' @export
getCellScoresInSitu <- function(object, sigmaSquaredChoice, ccIndex = 1,
                                scoreColorType = c("binary", "continuous")) {
  ## check input
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  ## match arg
  scoreColorType <- match.arg(scoreColorType)

  ## make sure normalizedCorrelation exists
  if (length(object@cellScores) == 0 ||
    length(object@geneScores) == 0) {
    stop(paste(
      "cellScores slot does not exist,",
      "run `computeGeneAndCellScores()` first"
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

  sigma_name_choice <- paste("sigma", sigmaSquaredChoice, sep = "_")

  loc_t <- stats::setNames(
    vector(mode = "list", length = length(cts)),
    cts
  )
  median_score_t <- vector("numeric", length = length(cts))
  names(median_score_t) <- cts

  for (t in cts) {
    loc_t[[t]] <- object@locationDataSub[object@cellTypesSub == t, ]
    loc_t[[t]]$"cellScores" <- object@cellScores[[sigma_name_choice]][[t]][
      rownames(loc_t[[t]]),
      ccIndex
    ]
    loc_t[[t]]$"cellTypesSub" <- t
    median_score_t[t] <- median(loc_t[[t]]$"cellScores")
    loc_t[[t]]$"cellScores_b" <- ifelse(loc_t[[t]]$"cellScores" > median_score_t[t],
      paste0("high_", t), paste0("low_", t)
    )
  }

  combinations <- expand.grid(c("high", "low"), cts)
  all_binary_scores <- apply(combinations, 1, function(x) paste(x[1], x[2],
                                                                sep = "_"))

  names(loc_t) <- NULL
  loc_all <- do.call(rbind, loc_t)[rownames(object@locationDataSub), ]

  loc_all$cellScores_b <- factor(loc_all$cellScores_b,
    levels = all_binary_scores
  )

  return(loc_all)
}





#' Retrieve the Correlation between two cell types
#'
#' @param object A `CoPro` object
#' @param cellTypeA Cell type label for one cell type
#' @param cellTypeB Cell type label for another cell type
#' @param sigmaSquaredChoice A particular sigma squared value
#' for the correlation
#' @param ccIndex Canonical vector index, default = 1
#'
#' @return A data.frame with two columns, AK and B, where AK represents the
#' cell score of cell type A times the kernel matrix, and B represents the
#' cell score of cell type B.
#' @export
getCorrTwoTypes <- function(object, cellTypeA, cellTypeB, ccIndex = 1,
                            sigmaSquaredChoice) {
  ## check input
  if (length(cellTypeA) != 1) {
    stop("Must give a single cellTypeA for correlation plot")
  }
  if (length(cellTypeB) != 1) {
    stop("Must give a single cellTypeB for correlation plot")
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

  ## set sigmaSquaredChoice
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

  ## make sure normalizedCorrelation exists
  if (length(object@cellScores) == 0 ||
    length(object@geneScores) == 0) {
    stop(paste(
      "cellScores slot does not exist,",
      "run computeGeneAndCellScores first"
    ))
  }

  ## load the cellScores and kernel matrix
  sigma_name <- paste("sigma", sigmaSquaredChoice, sep = "_")

  x1 <- t(object@cellScores[[sigma_name]][[cellTypeA]][, ccIndex, drop = FALSE])
  x2 <- object@cellScores[[sigma_name]][[cellTypeB]][, ccIndex, drop = TRUE]
  ktemp <- object@kernelMatrices[[sigma_name]][[cellTypeA]][[cellTypeB]]
  if (length(ktemp) != 0) {
    k <- ktemp
  }else{
    k <- t(object@kernelMatrices[[sigma_name]][[cellTypeB]][[cellTypeA]])
  }

  df <- data.frame(AK = (x1 %*% k)[1,,drop = TRUE], B = x2)
  return(df)
}
