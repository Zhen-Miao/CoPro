## remove outliers from a CoPro object



removeOutliers <- function(object, threshold,
                          method = c("zscore", "iqr"),
                          verbose = TRUE) {
  method <- match.arg(method)

  if (verbose) {
    message("Removing outliers using method: ", method)
  }

  ## obtain cell scores
  cs <- object$cellScores

  ## slots to be taken subset
  # @cellTypesSub
  # @cellScores
  # @normalizedDataSub
  # @metaDataSub
  # @slideID if exist
  # @locationDataSub


  if (method == "zscore") {
    z_scores <- abs(scale(object$data))
    outliers <- which(z_scores > threshold, arr.ind = TRUE)
  } else if (method == "iqr") {
    Q1 <- apply(object$data, 2, quantile, probs = 0.25)
    Q3 <- apply(object$data, 2, quantile, probs = 0.75)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - threshold * IQR
    upper_bound <- Q3 + threshold * IQR
    outliers <- which(object$data < lower_bound | object$data > upper_bound, arr.ind = TRUE)
  }

  if (length(outliers) > 0) {
    object$data[outliers] <- NA
    if (verbose) {
      message("Removed ", length(outliers), " outliers.")
    }
  } else {
    if (verbose) {
      message("No outliers found.")
    }
  }

  return(object)
}
