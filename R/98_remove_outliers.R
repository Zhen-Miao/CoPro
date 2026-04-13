## remove outliers from a CoPro object



removeOutliers <- function(object, threshold,
                          method = c("zscore", "iqr"),
                          verbose = TRUE) {
  method <- match.arg(method)

  if (verbose) {
    message("Removing outliers using method: ", method)
  }

  expr <- object@normalizedDataSub
  if (inherits(expr, "Matrix")) {
    if (verbose) {
      message("Converting sparse matrix normalizedDataSub to dense matrix for outlier detection.")
    }
    expr <- as.matrix(expr)
  } else if (is.data.frame(expr)) {
    expr <- as.matrix(expr)
  } else if (!is.matrix(expr)) {
    stop("normalizedDataSub must be a matrix-like object.")
  }
  
  if (method == "zscore") {
    z_scores <- abs(scale(expr))
    outliers <- which(z_scores > threshold, arr.ind = TRUE)
  } else if (method == "iqr") {
    Q1 <- apply(expr, 2, quantile, probs = 0.25)
    Q3 <- apply(expr, 2, quantile, probs = 0.75)
    IQR <- Q3 - Q1
    lower_bound <- Q1 - threshold * IQR
    upper_bound <- Q3 + threshold * IQR
    outliers <- which(expr < lower_bound | expr > upper_bound, arr.ind = TRUE)
  }

  outlier_cells <- unique(outliers[, 1])
  if (length(outlier_cells) > 0) {
    expr[outliers] <- NA
    object@normalizedDataSub <- expr
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
