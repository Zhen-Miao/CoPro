## remove outliers from a CoPro object



removeOutliers <- function(object, threshold,
                          method = c("zscore", "iqr"),
                          verbose = TRUE) {
  method <- match.arg(method)

  if (verbose) {
    message("Removing outliers using method: ", method)
  }

  ## obtain cell scores
  cs <- object@cellScores

  ## slots to be taken subset
  # @cellTypesSub
  # @cellScores
  # @normalizedDataSub
  # @metaDataSub
  # @slideID if exist
  # @locationDataSub


  if (method == "zscore") {
    z_scores <- abs(scale(object@normalizedDataSub))
    outliers <- which(z_scores > threshold, arr.ind = TRUE)
  } else if (method == "iqr") {
    Q1 <- apply(object@normalizedDataSub, 2, quantile, probs = 0.25)
    Q3 <- apply(object@normalizedDataSub, 2, quantile, probs = 0.75)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - threshold * IQR
    upper_bound <- Q3 + threshold * IQR
    outliers <- which(object@normalizedDataSub < lower_bound | object@normalizedDataSub > upper_bound, arr.ind = TRUE)
  }

  outlier_cells <- unique(outliers[, 1])
  if (length(outlier_cells) > 0) {
    object@normalizedDataSub[outliers] <- NA
    if (verbose) {
      message("Removed ", length(outlier_cells), " outlier cells.")
    }
  } else {
    if (verbose) {
      message("No outliers found.")
    }
  }

  return(object)
}
