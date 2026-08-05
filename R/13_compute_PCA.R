#' Compute PCA on Single- or Multi-Slide Data
#'
#' Performs PCA on the normalized data stored within a `CoProSingle` or
#' `CoProMulti` object. Multi-slide data must share the same gene columns; the
#' default preprocessing removes slide-level location and scale internally.
#'
#' @importFrom stats setNames prcomp
#' @importFrom irlba irlba prcomp_irlba
#' @param object A `CoProSingle` or `CoProMulti` object.
#' @param nPCA Number of principal components to compute for each cell type.
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'   Default is "raw". For single slide, this argument is ignored.
#' @param center_per_slide For multi-slide data, apply `center` and `scale.`
#'   within every (slide, cell type) block *before* fitting one joint PCA per
#'   cell type. This removes between-slide location and scale effects from the
#'   covariance that chooses the shared loading matrix. Default `TRUE` for
#'   `CoProMulti`. Set `FALSE` to recover the legacy pooled preprocessing.
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#' @param scalePCs Whether to scale (whiten) PCs by their standard deviation
#'   before downstream CCA optimization. Default `TRUE`. With `FALSE`, CoPro
#'   carries the diagonal PC-variance metric through observed, permutation,
#'   and deflation calculations, so this is a supported reparameterization.
#'
#' @details For cell type \eqn{i} and slide \eqn{s}, the recommended multi-slide
#' path forms
#' \deqn{Z_i^{(s)} = (X_i^{(s)} - 1\mu_i^{(s)'})D_i^{(s)-1},}
#' stacks the \eqn{Z_i^{(s)}} by rows, and computes one truncated SVD
#' \eqn{Z_i = U_i\Delta_i V_i'}. Every slide therefore uses the same loading
#' matrix \eqn{V_i}, with score block \eqn{Z_i^{(s)}V_i}. Because each
#' standardized gene column has zero mean within its block, those PC scores are
#' also centered within slide automatically; no post-PCA recentering is needed.
#'
#' The shared loading is in **within-slide standardized gene coordinates**. If
#' the stored slide scales differ, there is intentionally no single equivalent
#' coefficient vector in raw expression units.
#'
#' @return The input object with `pcaGlobal` populated. For `CoProMulti`,
#'   `pcaResults` has structure
#'   `list(slideID = list(cellType = score-row view))`.
#' @family spatial-pipeline
#' @seealso [computeDistance()], [computeKernelMatrix()], [runSkrCCA()]
#' @examples
#' toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
#' obj <- newCoProSingle(
#'   normalizedData = toy$normalizedData,
#'   locationData   = toy$locationData,
#'   metaData       = toy$metaData,
#'   cellTypes      = toy$cellTypes
#' )
#' obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
#' obj <- computePCA(obj, nPCA = 10)
#' @export
#' @rdname computePCA
setGeneric("computePCA",
           function(object, nPCA = 40,
                    center = TRUE, scale. = TRUE,
                    scalePCs = TRUE,
                    dataUse = "raw",
                    center_per_slide = TRUE) standardGeneric("computePCA"))

# Common input validation function
.validate_pca_params <- function(nPCA, center, scale., scalePCs) {
  if (!is.numeric(nPCA) || nPCA <= 0 || nPCA != as.integer(nPCA)) {
    stop("nPCA must be a positive integer")
  }
  if (!is.logical(center) || length(center) != 1) {
    stop("center must be a single logical value")
  }
  if (!is.logical(scale.) || length(scale.) != 1) {
    stop("scale. must be a single logical value")
  }
  if (!is.logical(scalePCs) || length(scalePCs) != 1) {
    stop("scalePCs must be a single logical value")
  }
}

# Common function to apply centering and scaling
.apply_centering_scaling <- function(mat, center, scale.) {
  if (center && scale.) {
    return(center_scale_matrix_opt(mat))
  } else if (center) {
    # sweep() subtracts in place over one pass. The previous
    # t(t(mat) - colMeans(mat)) spelling made two full transposed copies of a
    # cells-by-genes matrix to achieve the same thing.
    return(sweep(mat, 2L, colMeans(mat), "-"))
  } else if (!scale.) {
    warning(paste(
      "It is not recommended to skip both centering and scaling of the data,",
      "unless the data has been centered and scaled when creating the CoPro object."
    ))
    return(mat)
  } else {
    # Only scaling without centering
    return(scale(mat, center = FALSE, scale = TRUE))
  }
}

.check_pca_input <- function(object, nPCA, center, scale., scalePCs) {
  .validate_pca_params(nPCA, center, scale., scalePCs)

  # Choose cell types
  if (length(object@cellTypesOfInterest) != 0) {
    cts <- object@cellTypesOfInterest
  } else {
    warning(paste(
      "No cell type of interest specified,",
      "using all cell types to run the analysis"
    ))
    cts <- unique(object@cellTypesSub)
  }

  return(cts)
}

.is_bpcells <- function(x) {
  inherits(x, "IterableMatrix")
}

.max_pca_rank <- function(mat) {
  max(0L, min(nrow(mat) - 1L, ncol(mat) - 1L))
}

.resolve_common_pca_rank <- function(matrices, nPCA, cts,
                                     max_ranks = NULL) {
  missing <- cts[vapply(matrices, is.null, logical(1))]
  if (length(missing) > 0L) {
    stop("PCA input is missing for cell type(s): ",
         paste(missing, collapse = ", "))
  }
  if (is.null(max_ranks)) {
    max_ranks <- stats::setNames(
      vapply(matrices, .max_pca_rank, integer(1)), cts
    )
  }
  if (any(max_ranks < 1L)) {
    bad <- names(max_ranks)[max_ranks < 1L]
    stop("PCA requires at least two cells and two genes for every cell type. ",
         "Insufficient data for: ", paste(bad, collapse = ", "))
  }

  common_rank <- min(as.integer(nPCA), min(max_ranks))
  if (common_rank < nPCA) {
    limiting <- names(max_ranks)[max_ranks == min(max_ranks)]
    warning(
      "nPCA (", nPCA, ") exceeds the common feasible rank across cell types. ",
      "Using ", common_rank, " PCs for every cell type; limiting type(s): ",
      paste(limiting, collapse = ", "), "."
    )
  }
  common_rank
}

.max_within_slide_pca_rank <- function(mat, slide_ids, center) {
  irlba_rank <- .max_pca_rank(mat)
  if (!center) return(irlba_rank)
  # Centering each nonempty slide block removes one independent row-space
  # direction per block. This is stricter than pooled centering's N - 1 bound.
  n_blocks <- length(unique(slide_ids))
  max(0L, min(irlba_rank, nrow(mat) - n_blocks))
}

# Compute centering/scaling vectors without materializing a centered sparse
# matrix. irlba applies these vectors inside its matrix products.
.sparse_pca_parameters <- function(mat, center, scale.,
                                   zero_sd_threshold = 1e-3,
                                   nz_propion_threshold = 0.01) {
  n <- nrow(mat)
  means <- as.numeric(Matrix::colMeans(mat))
  sumsq <- as.numeric(Matrix::colSums(mat ^ 2))

  center_arg <- if (center) means else FALSE
  scale_arg <- FALSE
  if (scale.) {
    variance <- if (center) {
      pmax((sumsq - n * means^2) / max(1, n - 1L), 0)
    } else {
      sumsq / max(1, n - 1L)
    }
    scale_values <- sqrt(variance)
    nz_prop <- as.numeric(Matrix::colSums(mat != 0)) / n
    unsafe <- !is.finite(scale_values) |
      scale_values < zero_sd_threshold |
      nz_prop < nz_propion_threshold
    scale_values[unsafe] <- 1
    scale_arg <- scale_values
  }

  list(center = center_arg, scale = scale_arg)
}

.run_pca_irlba <- function(mat, nPCA, center, scale.) {
  if (inherits(mat, "sparseMatrix")) {
    params <- .sparse_pca_parameters(mat, center, scale.)
    message("Input is sparse (", class(mat)[1],
            "), performing implicitly centered/scaled irlba PCA...")
    return(prcomp_irlba(
      mat, n = nPCA, center = params$center, scale. = params$scale
    ))
  }

  scaled_data <- .apply_centering_scaling(mat, center, scale.)
  message("Input is dense (", class(scaled_data)[1],
          "), performing irlba PCA...")
  prcomp_irlba(scaled_data, center = FALSE, scale. = FALSE, n = nPCA)
}

# -------------------------------------------------------------------------
# Within-slide preprocessing for a shared PCA
# -------------------------------------------------------------------------

#' Matrix-free within-slide standardized expression matrix
#'
#' Internal representation of
#' \deqn{Z_s = (X_s - 1\mu_s')\operatorname{diag}(d_s)^{-1}}
#' stacked in the original row order. `irlba()` only needs right and left
#' matrix products, so sparse inputs can stay sparse even though explicit
#' centering would make every zero nonzero.
#'
#' @slot blocks Per-slide expression blocks in their original storage class.
#' @slot rows Original row indices for each block.
#' @slot centers Matrix of per-block gene means.
#' @slot scales Matrix of per-block gene scales.
#' @slot dims Integer vector containing the stacked matrix dimensions.
#' @keywords internal
methods::setClass(
  "CoProWithinSlideMatrix",
  slots = c(
    blocks = "list",
    rows = "list",
    centers = "matrix",
    scales = "matrix",
    dims = "integer"
  )
)

#' @param x A `CoProWithinSlideMatrix` object or conformable numeric operand.
#' @return Matrix dimensions for `dim()`; a dense low-dimensional product for
#'   `%*%`.
#' @rdname CoProWithinSlideMatrix-class
methods::setMethod("dim", "CoProWithinSlideMatrix", function(x) x@dims)

.withinSlideRightMultiply <- function(x, y) {
  y <- as.matrix(y)
  out <- matrix(0, nrow = x@dims[[1L]], ncol = ncol(y))
  for (k in seq_along(x@blocks)) {
    rows <- x@rows[[k]]
    if (length(rows) == 0L) next
    scaled_y <- sweep(y, 1L, x@scales[k, ], "/")
    block_out <- as.matrix(x@blocks[[k]] %*% scaled_y)
    shift <- as.numeric(x@centers[k, ]) %*% scaled_y
    out[rows, ] <- sweep(block_out, 2L, as.numeric(shift), "-")
  }
  out
}

.withinSlideLeftMultiply <- function(x, y) {
  x <- as.matrix(x)
  out <- matrix(0, nrow = nrow(x), ncol = y@dims[[2L]])
  for (k in seq_along(y@blocks)) {
    rows <- y@rows[[k]]
    if (length(rows) == 0L) next
    x_block <- x[, rows, drop = FALSE]
    contribution <- as.matrix(x_block %*% y@blocks[[k]])
    contribution <- contribution -
      rowSums(x_block) %o% as.numeric(y@centers[k, ])
    out <- out + sweep(contribution, 2L, y@scales[k, ], "/")
  }
  out
}

#' @param y A `CoProWithinSlideMatrix` object or conformable numeric operand.
#' @rdname CoProWithinSlideMatrix-class
methods::setMethod(
  "%*%", signature(x = "CoProWithinSlideMatrix", y = "numeric"),
  function(x, y) .withinSlideRightMultiply(x, y)
)
#' @rdname CoProWithinSlideMatrix-class
methods::setMethod(
  "%*%", signature(x = "CoProWithinSlideMatrix", y = "matrix"),
  function(x, y) .withinSlideRightMultiply(x, y)
)
#' @rdname CoProWithinSlideMatrix-class
methods::setMethod(
  "%*%", signature(x = "numeric", y = "CoProWithinSlideMatrix"),
  function(x, y) .withinSlideLeftMultiply(matrix(x, nrow = 1L), y)
)
#' @rdname CoProWithinSlideMatrix-class
methods::setMethod(
  "%*%", signature(x = "matrix", y = "CoProWithinSlideMatrix"),
  function(x, y) .withinSlideLeftMultiply(x, y)
)

.withinSlidePCAParameters <- function(mat, slide_ids, slides, center, scale.,
                                      zero_sd_threshold = 1e-3,
                                      nz_propion_threshold = 0.01) {
  p <- ncol(mat)
  rows <- stats::setNames(lapply(slides, function(s) which(slide_ids == s)), slides)
  centers <- matrix(0, nrow = length(slides), ncol = p,
                    dimnames = list(slides, colnames(mat)))
  scales <- matrix(1, nrow = length(slides), ncol = p,
                   dimnames = list(slides, colnames(mat)))

  for (k in seq_along(slides)) {
    idx <- rows[[k]]
    if (length(idx) == 0L) next
    block <- mat[idx, , drop = FALSE]
    means <- as.numeric(colMeans(block))
    if (center) centers[k, ] <- means

    if (scale.) {
      block_scale <- if (center) {
        if (.is_bpcells(block)) {
          sqrt(as.numeric(BPCells::colVars(block)))
        } else {
          as.numeric(.columnSds(block))
        }
      } else {
        sqrt(as.numeric(colSums(block ^ 2)) / max(1, nrow(block) - 1L))
      }
      nz_prop <- .columnNonzeroFraction(block)
      unsafe <- !is.finite(block_scale) |
        block_scale < zero_sd_threshold |
        nz_prop < nz_propion_threshold
      block_scale[unsafe] <- 1
      scales[k, ] <- block_scale
    }
  }

  list(rows = rows, centers = centers, scales = scales)
}

.materializeWithinSlideMatrix <- function(mat, params) {
  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat),
                dimnames = dimnames(mat))
  for (k in seq_along(params$rows)) {
    rows <- params$rows[[k]]
    if (length(rows) == 0L) next
    block <- as.matrix(mat[rows, , drop = FALSE])
    block <- sweep(block, 2L, params$centers[k, ], "-")
    out[rows, ] <- sweep(block, 2L, params$scales[k, ], "/")
  }
  out
}

.run_within_slide_pca <- function(mat, slide_ids, slides, nPCA,
                                  center, scale.) {
  params <- .withinSlidePCAParameters(
    mat, slide_ids, slides, center = center, scale. = scale.
  )

  if (inherits(mat, "sparseMatrix") || .is_bpcells(mat)) {
    blocks <- lapply(params$rows, function(rows) mat[rows, , drop = FALSE])
    op <- methods::new(
      "CoProWithinSlideMatrix",
      blocks = blocks,
      rows = unname(params$rows),
      centers = params$centers,
      scales = params$scales,
      dims = as.integer(dim(mat))
    )
    message("Input is ", class(mat)[1],
            ", performing implicitly within-slide standardized irlba PCA...")
    sv <- irlba::irlba(op, nv = nPCA, nu = nPCA, fastpath = FALSE)
    x_scores <- sweep(sv$u, 2L, sv$d, "*")
    pca <- list(
      sdev = sv$d / sqrt(max(1, nrow(mat) - 1L)),
      rotation = sv$v,
      x = x_scores,
      center = FALSE,
      scale = FALSE
    )
    class(pca) <- "prcomp"
  } else {
    standardized <- .materializeWithinSlideMatrix(mat, params)
    message("Input is dense (", class(standardized)[1],
            "), performing within-slide standardized irlba PCA...")
    pca <- prcomp_irlba(
      standardized, center = FALSE, scale. = FALSE, n = nPCA
    )
  }

  # The custom matrix-free operator has no dimnames method, so irlba cannot
  # propagate names on that path. Keep its prcomp-compatible output identical
  # to the dense path for downstream gene back-projection and score lookup.
  pc_names <- colnames(pca$x)
  if (is.null(pc_names)) pc_names <- paste0("PC", seq_len(ncol(pca$x)))
  rownames(pca$x) <- rownames(mat)
  colnames(pca$x) <- pc_names
  rownames(pca$rotation) <- colnames(mat)
  colnames(pca$rotation) <- pc_names

  # Extra fields preserve the affine map that defines the shared feature
  # coordinates.  A single raw-unit back-projection does not exist when the
  # per-slide scales differ; the shared loading lives in standardized-gene
  # coordinates, exactly as in gene-space CCA.
  pca$preprocessing <- "within_slide"
  pca$slideCenter <- params$centers
  pca$slideScale <- params$scales
  pca
}

.compute_pca_single <- function(object, nPCA = 40, center = TRUE, scale. = TRUE, scalePCs = TRUE, cts) {
  # PCA results will be saved under the name of cell types
  object@pcaGlobal <- setNames(
    vector("list", length = length(cts)),
    cts
  )

  matrices <- stats::setNames(lapply(cts, function(ct) {
    object@normalizedDataSub[object@cellTypesSub == ct, , drop = FALSE]
  }), cts)
  nPCA_use <- .resolve_common_pca_rank(matrices, nPCA, cts)

  # Iterate over cell types
  for (ct in cts) {
    # Cell type specific subset
    sub_data <- matrices[[ct]]

    # PCA on the matrix that is already centered and scaled
    if (.is_bpcells(sub_data)) {
      scaled_data <- .apply_centering_scaling(sub_data, center, scale.)
      message("Input is BPCell (", class(scaled_data), "), performing BPCell svd...")
      sv <- BPCells::svds(scaled_data, k = nPCA_use, nu = nPCA_use, nv = nPCA_use, threads = 0L)
      x_scores <- sweep(sv$u, 2, sv$d, `*`)

      # add cell names to pca matrix
      cell_ids <- rownames(sub_data)
      if (is.null(cell_ids)) {
        # fallback: use metaDataSub rownames aligned to cellTypesSub
        cell_ids <- rownames(object@metaDataSub)[object@cellTypesSub == ct]
      }
      rownames(x_scores) <- cell_ids
      colnames(x_scores) <- paste0("PC_", seq_len(ncol(x_scores)))
      colnames(sv$v) <- paste0("PC_", seq_len(ncol(sv$v)))

      # prcomp-like object
      pca <- list(
        sdev     = sv$d / sqrt(max(1, nrow(scaled_data) - 1)),
        rotation = sv$v,
        x        = x_scores,
        center   = NULL,
        scale    = NULL
      )
      class(pca) <- "prcomp"
    } else {
      pca <- .run_pca_irlba(sub_data, nPCA_use, center, scale.)
    }
    object@pcaGlobal[[ct]] <- pca
  }

  object@nPCA <- nPCA_use
  object@scalePCs <- scalePCs

  return(object)
}

.check_pca_input_multi <- function(object, nPCA, center, scale., scalePCs,
                                   dataUse, center_per_slide) {
  .validate_pca_params(nPCA, center, scale., scalePCs)

  if (!is.logical(center_per_slide) || length(center_per_slide) != 1L ||
      is.na(center_per_slide)) {
    stop("center_per_slide must be a single non-missing logical value")
  }

  # Validate dataUse argument
  if (!dataUse %in% c("raw", "integrated")) {
    stop("dataUse must be 'raw' or 'integrated'")
  }

  # Check if integrated data exists when needed
  if (dataUse == "integrated" && length(object@integratedData) == 0) {
    stop("integratedData slot is empty. Run integration first.")
  }

  # Check cell types of interest
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0) {
    stop("cellTypesOfInterest not set. Run subsetData first.")
  }

  return(cts)
}

.compute_pca_multi <- function(object, nPCA = 40, center = TRUE, scale. = TRUE,
                               scalePCs = TRUE, dataUse = "raw",
                               center_per_slide = TRUE, cts) {
  slides <- getSlideList(object)

  # Initialize pcaResults structure
  pca_results_all <- setNames(vector("list", length = length(slides)), slides)
  for (slide_id in slides) {
    pca_results_all[[slide_id]] <- setNames(vector("list", length = length(cts)), cts)
  }

  pca_global <- setNames(vector("list", length = length(cts)), cts)

  matrices <- stats::setNames(lapply(cts, function(ct) {
    if (dataUse == "integrated") {
      object@integratedData[[ct]]
    } else {
      object@normalizedDataSub[object@cellTypesSub == ct, , drop = FALSE]
    }
  }), cts)
  within_ranks <- if (center_per_slide) {
    stats::setNames(vapply(cts, function(ct) {
      rows <- object@cellTypesSub == ct
      .max_within_slide_pca_rank(
        matrices[[ct]], getSlideID(object)[rows], center = center
      )
    }, integer(1)), cts)
  } else {
    NULL
  }
  nPCA_use <- .resolve_common_pca_rank(
    matrices, nPCA, cts, max_ranks = within_ranks
  )

  # Perform PCA per cell type on the integrated data
  for (ct in cts) {
    message("Performing PCA for cell type: ", ct)

    if (dataUse == "integrated" && !ct %in% names(object@integratedData)) {
      warning("No integrated data found for cell type: ", ct, " - Skipping PCA.")
      next
    }

    # Get the appropriate data matrix
    if (dataUse == "integrated") {
      mat_ct <- matrices[[ct]]
    } else { # raw
      mat_ct <- matrices[[ct]]
    }

    # Ensure it's a matrix
    if (!.is_bpcells(mat_ct) && !is.matrix(mat_ct) && !inherits(mat_ct, "Matrix")) {
      message("Converting data matrix into dense matrix")
      mat_ct <- as.matrix(mat_ct)
    }

    # Check dimensions
    expected_rows <- sum(object@cellTypesSub == ct)
    if (nrow(mat_ct) != expected_rows) {
      stop("Data dimensions mismatch for cell type: ", ct,
           ". Expected ", expected_rows, " rows, got ", nrow(mat_ct))
    }

    slide_id_ct <- getSlideID(object)[object@cellTypesSub == ct]

    # Recommended path: remove within-slide gene means/scales before the
    # covariance that chooses the single shared loading matrix is formed.
    if (center_per_slide) {
      pca_ct <- .run_within_slide_pca(
        mat_ct, slide_id_ct, slides, nPCA_use, center = center, scale. = scale.
      )
    } else if (.is_bpcells(mat_ct)) {
      scaled_data <- .apply_centering_scaling(mat_ct, center, scale.)
      message("Input is BPCell (", paste(class(scaled_data), collapse = ", "),
              "), performing BPCell svd...")
      sv <- BPCells::svds(scaled_data, k = nPCA_use, nu = nPCA_use, nv = nPCA_use, threads = 0L)
      x_scores <- sweep(sv$u, 2, sv$d, `*`)

      # add cell id to pca matrix
      rownames(x_scores) <- rownames(object@metaDataSub)[object@cellTypesSub == ct]
      colnames(x_scores) <- paste0("PC_", seq_len(ncol(x_scores)))
      colnames(sv$v) <- paste0("PC_", seq_len(ncol(sv$v)))

      # prcomp-like object
      pca_ct <- list(
        sdev     = sv$d / sqrt(max(1, nrow(scaled_data) - 1)),
        rotation = sv$v,
        x        = x_scores,
        center   = NULL,
        scale    = NULL
      )
      class(pca_ct) <- "prcomp"
    } else {
      pca_ct <- .run_pca_irlba(mat_ct, nPCA_use, center, scale.)
    }
    if (!center_per_slide) pca_ct$preprocessing <- "pooled"
    message("PCA computed for cell type: ", ct)

    # Project each slide's data onto the shared PCs
    row_names_ct <- rownames(object@metaDataSub)[object@cellTypesSub == ct]

    # Label the global scores once, so the per-slide views inherit the cell IDs
    # instead of each slice having to be relabelled. prcomp_irlba() does not
    # carry rownames through, which is why the per-slice assignment below used
    # to be necessary.
    rownames(pca_ct$x) <- row_names_ct
    pca_global[[ct]] <- pca_ct

    for (slide_id in slides) {
      rows <- which(slide_id_ct == slide_id)

      # Both preprocessing modes now store rows of the one global score matrix.
      # Under within-slide preprocessing those rows already have zero mean by
      # construction, so no post-PCA materialization or repair is needed.
      pca_results_all[[slide_id]][[ct]] <- .newPCSlice(rows)
    }
  } # End loop over cell types

  object@pcaGlobal <- pca_global
  object@pcaResults <- pca_results_all
  object@nPCA <- nPCA_use
  object@scalePCs <- scalePCs

  return(object)
}

#' @param nPCA Number of principal components to compute for each cell type.
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#' @param scalePCs Whether to scale (whiten) PCs by their standard deviation
#'   before downstream CCA optimization. Default \code{TRUE} (recommended).
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'   Default is "raw". For single slide, this argument is ignored.
#' @param center_per_slide For multi-slide data, apply `center` and `scale.`
#'   within slide before fitting the shared PCA. Default `TRUE`.
#' @rdname computePCA
#' @aliases computePCA,CoProSingle-method
#' @export
setMethod("computePCA", "CoProSingle",
          function(object, nPCA = 40, center = TRUE, scale. = TRUE,
                   scalePCs = TRUE) {
            cts <- .check_pca_input(object, nPCA, center, scale., scalePCs)
            object <- .compute_pca_single(object, nPCA, center, scale., scalePCs, cts)
            return(object)
          })

#' @param nPCA Number of principal components to compute for each cell type.
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'   Default is "raw". For single slide, this argument is ignored.
#' @param center_per_slide Apply `center` and `scale.` within each
#'   (slide, cell type) block before fitting one shared PCA. Default `TRUE`.
#' @rdname computePCA
#' @aliases computePCA,CoProMulti-method
#' @export
setMethod("computePCA", "CoProMulti",
          function(object, nPCA = 40, center = TRUE, scale. = TRUE,
                   scalePCs = TRUE, dataUse = "raw", center_per_slide = TRUE) {
            cts <- .check_pca_input_multi(
              object, nPCA, center, scale., scalePCs, dataUse,
              center_per_slide
            )
            object <- .compute_pca_multi(object, nPCA, center, scale., scalePCs,
                                         dataUse, center_per_slide, cts)
            return(object)
          })
