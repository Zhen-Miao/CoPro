
## correlation plot

plotNormCorr <- function(object){

  ## check input
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  ## make sure normalizedCorrelation exists
  if (length(object@normalizedCorrelation) == 0){
    stop(paste("normalizedCorrelation slot does not exist,",
               "run computeNormalizedCorrelation first"))
  }

  ## organize into a data.frame
  ncorr <- do.call(rbind, slot(object, "normalizedCorrelation"))
  ncorr$ct12 <- paste(ncorr$cellType1, ncorr$cellType2, sep = "-")
  ncorr$sigmaSquares <- factor(ncorr$sigmaSquares,
                               levels = sort(unique(ncorr$sigmaSquares),
                                             decreasing = FALSE))

  p1 <- ggplot(data = ncorr)+
    geom_line(mapping = aes(x = sigmaSquares, y = normalizedCorrelation))+
    facet_wrap(vars(ct12))+
    xlab("Sigma squared")+
    ylab("Norm. Corr.")+
    ggtitle(label = "Norm. Corr. across sigma squared values")+
    theme_minimal()

  return(p1)

}


plotCellScoresInSitu <- function(object, sigmaSquaredChoice,
                                 scoreColorType = c("binary", "continuous") ){
  ## check input
  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }

  ## match arg
  scoreColorType <- match.arg(scoreColorType)

  ## make sure normalizedCorrelation exists
  if (length(object@cellScores) == 0 |
      length(object@geneScores == 0)){
    stop(paste("cellScores slot does not exist,",
               "run computeNormalizedCorrelation first"))
  }

  if(is.null(sigmaSquaredChoice)){
    stop(paste("sigmaSquaredChoice is not given",
               "default set to the value with highest",
               "normalized correlation."))
    sigmaSquaredChoice <- object@sigmaSquaredChoice
  }

  if (!(sigmaSquaredChoice %in% object@sigmaSquares)){
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

  loc_t <- stats::setNames(vector(mode = "list", length = length(cts)),
                           cts)
  median_score_t <- vector("numeric", length = length(cts))
  names(median_score_t) <- cts

  for (t in cts) {
    loc_t[[t]] <- object@locationDataSub[object@cellTypesSub == t,]
    loc_t[[t]]$cellScores <- object@cellScores[[t]][rownames(loc_t[[t]]),
                                               sigma_name_choice]
    loc_t[[t]]$cellTypesSub <- t
    median_score_t[t] <- median(object@cellScores[[t]][,sigma_name_choice])
    loc_t[[t]]$cellScores_b <- ifelse(loc_t[[t]]$cellScores > median_score_t[t],
                                      paste0("high_",t), paste0("low_",t))

  }

  combinations <- expand.grid(c("high", "low"), cts)
  all_binary_scores <- apply(combinations, 1, function(x) paste(x[1], x[2], sep = "-"))

  names(loc_t) <- NULL
  loc_all <- do.call(rbind, loc_t)[rownames(object@locationDataSub),]

  loc_all$cellScores_b = factor(loc_all$cellScores_b,
                                levels = all_binary_scores)

  if (scoreColorType == "binary") {
    p2 <- ggplot(data = loc_all)+
      geom_point(aes(x = x,
                     y = y,
                     color = cellScores_b), size = 0.8)+
      scale_color_discrete(type = c('darkgreen','#d4f8d4','darkred','#ffd8d8'))+
      coord_fixed()+
      theme_minimal()
    return(p2)
  }else{
    p2 <- ggplot(data = loc_all)+
      geom_point(aes(x = x,
                     y = y,
                     color = cellScores), size = 0.8)+
      coord_fixed()+
      theme_minimal()
    return(p2)
  }


}





