#' Fit a frozen CoPro score-transfer reference
#'
#' Builds a self-contained reference for transferring fitted CoPro axes to new
#' slides. For each cell type, the function freezes the training-expression
#' mean and standard deviation together with the exact PCA back-projected CoPro
#' weights. [predict()] then applies
#' \deqn{((X_{target} - \mu_{train}) / \sigma_{train}) W}
#' without estimating anything from the target cohort.
#'
#' This frozen log-expression workflow is the recommended default for
#' cross-slide score transfer after an internal benchmark. It preserves a
#' fixed out-of-sample map: adding or removing other target slides does not
#' change an existing cell's score.
#'
#' @param object For `fit_score_reference()`, a fitted `CoProSingle` or
#'   `CoProMulti` object. Run [computeGeneAndCellScores()] first; the exact PCA
#'   back-projected weights in `object@geneScores` are used for scoring. For the
#'   `predict()` method, a `CoProScoreReference` returned by
#'   `fit_score_reference()`.
#' @param sigma Numeric spatial bandwidth identifying the fitted CoPro weights.
#'   By default, uses the single value in `object@sigmaValueChoice`.
#' @param reference_weight How training moments are combined. `"cell_pooled"`
#'   (default) computes ordinary cell-pooled means and sample standard
#'   deviations. `"equal_slide"` averages each training slide's first and
#'   second moments with equal weight.
#'
#' @return `fit_score_reference()` returns a `CoProScoreReference` object
#'   containing the frozen transform, weights, and training provenance.
#'
#' @details
#' The reference and target must use the same normalized-expression
#' representation and modeled gene panel. This function deliberately does not
#' quantile-normalize the target or recompute target-specific moments. For
#' cross-platform quantile normalization, use [getTransferCellScores()] instead.
#'
#' Regression gene scores are useful for interpretation, but they are not the
#' canonical scoring map. This function therefore always uses the exact
#' back-projected weights from [computeGeneAndCellScores()].
#'
#' @seealso [computeGeneAndCellScores()], [getTransferCellScores()]
#' @export
#' @examples
#' genes <- c("g1", "g2")
#' training_x <- cbind(g1 = 1:20, g2 = c(2:11, 1:10))
#' rownames(training_x) <- paste0("train", 1:20)
#' training_type <- rep(c("A", "B"), each = 10)
#' training <- newCoProSingle(
#'   normalizedData = training_x,
#'   locationData = data.frame(
#'     x = 1:20, y = 1:20, row.names = rownames(training_x)
#'   ),
#'   metaData = data.frame(row.names = rownames(training_x)),
#'   cellTypes = training_type
#' )
#' training <- subsetData(training, c("A", "B"))
#'
#' # In a real analysis these weights come from computeGeneAndCellScores().
#' training@geneScores <- list(
#'   "geneScores|sigma0.1|A" = matrix(
#'     c(1, -1), ncol = 1, dimnames = list(genes, "CC_1")
#'   ),
#'   "geneScores|sigma0.1|B" = matrix(
#'     c(0.5, 0.5), ncol = 1, dimnames = list(genes, "CC_1")
#'   )
#' )
#' training@sigmaValueChoice <- 0.1
#' reference <- fit_score_reference(training)
#'
#' target_x <- training_x[1:10, ] + 1
#' rownames(target_x) <- paste0("target", 1:10)
#' target_type <- rep(c("A", "B"), each = 5)
#' target <- newCoProSingle(
#'   normalizedData = target_x,
#'   locationData = data.frame(
#'     x = 1:10, y = 1:10, row.names = rownames(target_x)
#'   ),
#'   metaData = data.frame(row.names = rownames(target_x)),
#'   cellTypes = target_type
#' )
#' target <- subsetData(target, c("A", "B"))
#' transferred <- predict(reference, target)
fit_score_reference <- function(
    object,
    sigma = NULL,
    reference_weight = c("cell_pooled", "equal_slide")) {
  if (!(methods::is(object, "CoProSingle") ||
          methods::is(object, "CoProMulti"))) {
    stop("object must be a CoProSingle or CoProMulti object", call. = FALSE)
  }
  reference_weight <- match.arg(reference_weight)

  cell_types <- as.character(object@cellTypesOfInterest)
  if (!length(cell_types)) {
    stop("Run subsetData() before fitting a score reference.", call. = FALSE)
  }
  if (!nrow(object@normalizedDataSub)) {
    stop("The fitted object has no subsetted expression data.", call. = FALSE)
  }
  if (!length(object@geneScores)) {
    stop(
      "Gene scores are unavailable. Run computeGeneAndCellScores() first.",
      call. = FALSE
    )
  }

  sigma <- .resolve_frozen_score_sigma(object, sigma)
  slide_id <- if (methods::is(object, "CoProMulti")) {
    as.character(getSlideID(object))
  } else {
    rep("single_slide", nrow(object@normalizedDataSub))
  }
  if (length(slide_id) != nrow(object@normalizedDataSub)) {
    stop(
      "Slide IDs are not aligned with the fitted expression rows.",
      call. = FALSE
    )
  }

  references <- stats::setNames(lapply(cell_types, function(cell_type) {
    rows <- which(as.character(object@cellTypesSub) == cell_type)
    if (!length(rows)) {
      stop("No fitted cells are available for cell type ", cell_type, ".",
           call. = FALSE)
    }
    key <- .createGeneScoresName(sigma, cell_type, slide = NULL)
    weights <- object@geneScores[[key]]
    if (is.null(weights)) {
      stop(
        "Missing exact gene scores for ", cell_type, " at sigma ", sigma,
        ". Run computeGeneAndCellScores() for that bandwidth.", call. = FALSE
      )
    }
    weights <- as.matrix(weights)
    if (!is.numeric(weights) || !nrow(weights) || !ncol(weights) ||
          is.null(rownames(weights)) || anyDuplicated(rownames(weights)) ||
          any(!is.finite(weights))) {
      stop("Gene scores for ", cell_type, " are not a finite named matrix.",
           call. = FALSE)
    }
    missing_genes <- setdiff(
      rownames(weights), colnames(object@normalizedDataSub)
    )
    if (length(missing_genes)) {
      stop("Gene-score rows are absent from the fitted expression matrix: ",
           paste(missing_genes, collapse = ", "), call. = FALSE)
    }
    expression <- object@normalizedDataSub[
      rows, rownames(weights), drop = FALSE
    ]
    moments <- .frozen_score_moments(
      expression, slide_id[rows], reference_weight
    )
    list(
      genes = rownames(weights),
      center = moments$center,
      scale = moments$scale,
      weights = weights,
      n_training_cells = length(rows),
      training_slides = unique(slide_id[rows])
    )
  }), cell_types)
  component_names <- lapply(references, function(reference) {
    colnames(reference$weights)
  })
  components_missing <- any(vapply(component_names, is.null, logical(1)))
  components_differ <- !all(vapply(
    component_names, identical, logical(1), component_names[[1L]]
  ))
  if (components_missing || components_differ) {
    stop(
      paste(
        "Exact gene scores must use the same named components for every",
        "cell type."
      ),
      call. = FALSE
    )
  }

  structure(
    list(
      version = "1.0.0",
      transform = "frozen_log_center_scale",
      sigma = sigma,
      reference_weight = reference_weight,
      cell_types = cell_types,
      references = references,
      training_cell_ids = rownames(object@normalizedDataSub),
      training_slides = unique(slide_id),
      package_version = as.character(utils::packageVersion("CoPro")),
      fitted_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    class = "CoProScoreReference"
  )
}

.resolve_frozen_score_sigma <- function(object, sigma) {
  if (is.null(sigma)) {
    chosen <- as.numeric(object@sigmaValueChoice)
    if (length(chosen) != 1L || !is.finite(chosen) || chosen <= 0) {
      stop(
        paste(
          "sigma must be supplied when object@sigmaValueChoice is not one",
          "positive finite value."
        ),
        call. = FALSE
      )
    }
    return(chosen)
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) ||
        sigma <= 0) {
    stop("sigma must be one positive finite number.", call. = FALSE)
  }
  as.numeric(sigma)
}

.frozen_score_moments <- function(expression, slide_id, reference_weight) {
  if (identical(reference_weight, "cell_pooled")) {
    center <- as.numeric(colMeans(expression))
    scale <- if (.is_bpcells(expression)) {
      sqrt(as.numeric(BPCells::colVars(expression)))
    } else {
      as.numeric(.columnSds(expression))
    }
  } else {
    slides <- unique(slide_id)
    first <- vapply(slides, function(slide) {
      as.numeric(colMeans(expression[slide_id == slide, , drop = FALSE]))
    }, numeric(ncol(expression)))
    second <- vapply(slides, function(slide) {
      block <- expression[slide_id == slide, , drop = FALSE]
      as.numeric(colMeans(block ^ 2))
    }, numeric(ncol(expression)))
    if (is.null(dim(first))) {
      first <- matrix(first, nrow = ncol(expression))
      second <- matrix(second, nrow = ncol(expression))
    }
    center <- rowMeans(first)
    scale <- sqrt(pmax(rowMeans(second) - center ^ 2, 0))
  }
  scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- 1
  names(center) <- names(scale) <- colnames(expression)
  list(center = center, scale = scale)
}

#' Predict scores from a frozen CoPro reference
#'
#' @param newdata A subsetted `CoProSingle` or `CoProMulti` target object using
#'   the same normalized-expression representation and genes as the reference.
#'   No PCA, kernel, or CCA fit is required on the target.
#' @param aggregate If `FALSE` (default), return a named list of score matrices
#'   by cell type. If `TRUE`, return one matrix in target-cell order.
#' @param chunk_size Positive integer number of target cells scored at once.
#' @param ... Unused.
#'
#' @return `predict.CoProScoreReference()` returns target cell scores as a named
#'   list or an aggregated matrix.
#'
#' @rdname fit_score_reference
#' @method predict CoProScoreReference
#' @export
predict.CoProScoreReference <- function(
    object, newdata, aggregate = FALSE, chunk_size = 20000L, ...) {
  if (!inherits(object, "CoProScoreReference")) {
    stop("object must be a CoProScoreReference", call. = FALSE)
  }
  if (!(methods::is(newdata, "CoProSingle") ||
          methods::is(newdata, "CoProMulti"))) {
    stop("newdata must be a CoProSingle or CoProMulti object", call. = FALSE)
  }
  if (!is.logical(aggregate) || length(aggregate) != 1L || is.na(aggregate)) {
    stop("aggregate must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
        !is.finite(chunk_size) || chunk_size < 1 ||
        chunk_size > .Machine$integer.max || chunk_size != floor(chunk_size)) {
    stop("chunk_size must be one positive integer.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)

  target_types <- as.character(newdata@cellTypesOfInterest)
  if (!setequal(target_types, object$cell_types)) {
    stop(
      paste(
        "Reference and target cellTypesOfInterest must contain the same",
        "cell types."
      ),
      call. = FALSE
    )
  }
  expression <- newdata@normalizedDataSub
  if (!nrow(expression)) {
    stop("The target has no subsetted expression data.", call. = FALSE)
  }
  target_labels <- as.character(newdata@cellTypesSub)
  if (length(target_labels) != nrow(expression)) {
    stop(
      "Target cell types are not aligned with expression rows.",
      call. = FALSE
    )
  }

  scores <- stats::setNames(lapply(object$cell_types, function(cell_type) {
    rows <- which(target_labels == cell_type)
    reference <- object$references[[cell_type]]
    missing_genes <- setdiff(reference$genes, colnames(expression))
    if (length(missing_genes)) {
      stop(
        "Target is missing frozen-reference genes for ", cell_type, ": ",
        paste(missing_genes, collapse = ", "), call. = FALSE
      )
    }
    result <- matrix(
      NA_real_, nrow = length(rows), ncol = ncol(reference$weights),
      dimnames = list(
        rownames(expression)[rows], colnames(reference$weights)
      )
    )
    if (!length(rows)) return(result)

    starts <- seq.int(1L, length(rows), by = chunk_size)
    for (start in starts) {
      local <- start:min(start + chunk_size - 1L, length(rows))
      block_rows <- rows[local]
      block <- as.matrix(expression[
        block_rows, reference$genes, drop = FALSE
      ])
      block <- sweep(block, 2L, reference$center, "-")
      block <- sweep(block, 2L, reference$scale, "/")
      result[local, ] <- block %*% reference$weights
    }
    result
  }), object$cell_types)

  if (aggregate) {
    scores <- do.call(rbind, unname(scores))
    scores <- scores[rownames(expression), , drop = FALSE]
  }
  attr(scores, "reference_sigma") <- object$sigma
  attr(scores, "reference_weight") <- object$reference_weight
  scores
}

#' @param x A `CoProScoreReference` object.
#' @param ... Unused.
#' @rdname fit_score_reference
#' @method print CoProScoreReference
#' @export
print.CoProScoreReference <- function(x, ...) {
  cat(
    "CoPro frozen score reference\n",
    "  sigma: ", format(x$sigma), "\n",
    "  cell types: ", paste(x$cell_types, collapse = ", "), "\n",
    "  reference weighting: ", x$reference_weight, "\n",
    "  training cells: ", length(x$training_cell_ids), "\n",
    sep = ""
  )
  invisible(x)
}
