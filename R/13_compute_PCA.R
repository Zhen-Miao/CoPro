#' Compute PCA on Single- or Multi-Slide Data
#'
#' Performs PCA on the normalized data stored within a `CoProSingle` or
#' `CoProMulti` object. For multi-slide objects, the data is assumed to have
#' already been integrated into a common space across slides.
#'
#' @importFrom stats setNames prcomp
#' @importFrom irlba prcomp_irlba
#' @param object A `CoProMulti` object with the `integratedData` slot populated.
#' @param nPCA Number of principal components to compute for each cell type.
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'   Default is "raw". For single slide, this argument is ignored.
#' @param center_per_slide After the global PCA, do we do center per slide
#'   again? By default this is set to FALSE
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#' @param scalePCs Whether to scale (whiten) PCs by their standard deviation
#'   before downstream CCA optimization. Default `TRUE`. With `FALSE`, CoPro
#'   carries the diagonal PC-variance metric through observed, permutation,
#'   and deflation calculations, so this is a supported reparameterization.
#'
#' @return A `CoProMulti` object with the `pcaResults` slot populated.
#'         `pcaResults` structure: `list(slideID = list(cellType = pc_matrix))`.
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
                    center_per_slide = FALSE) standardGeneric("computePCA"))

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
    return(t(t(mat) - colMeans(mat)))
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

.resolve_common_pca_rank <- function(matrices, nPCA, cts) {
  missing <- cts[vapply(matrices, is.null, logical(1))]
  if (length(missing) > 0L) {
    stop("PCA input is missing for cell type(s): ",
         paste(missing, collapse = ", "))
  }
  max_ranks <- stats::setNames(
    vapply(matrices, .max_pca_rank, integer(1)), cts
  )
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

.check_pca_input_multi <- function(object, nPCA, center, scale., scalePCs, dataUse) {
  .validate_pca_params(nPCA, center, scale., scalePCs)

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
                               scalePCs = TRUE, dataUse = "raw", center_per_slide = FALSE, cts) {
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
  nPCA_use <- .resolve_common_pca_rank(matrices, nPCA, cts)

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

    # Perform PCA on the combined integrated data for this cell type
    if (.is_bpcells(mat_ct)) {
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
    message("PCA computed for cell type: ", ct)
    pca_global[[ct]] <- pca_ct

    # Project each slide's data onto the shared PCs
    slide_id_ct <- getSlideID(object)[object@cellTypesSub == ct]
    row_names_ct <- rownames(object@metaDataSub)[object@cellTypesSub == ct]

    for (slide_id in slides) {
      pca_sub <- pca_ct$x[slide_id_ct == slide_id, , drop = FALSE]
      rownames(pca_sub) <- row_names_ct[slide_id_ct == slide_id]

      if (center_per_slide && nrow(pca_sub) > 0) {
        message("Centering per slide for slide: ", slide_id)
        pca_sub <- scale(pca_sub, center = TRUE, scale = FALSE)
      }

      pca_results_all[[slide_id]][[ct]] <- pca_sub
    }
  } # End loop over cell types

  object@pcaGlobal <- pca_global
  object@pcaResults <- pca_results_all
  object@nPCA <- nPCA_use
  object@scalePCs <- scalePCs

  return(object)
}

#' @param object A `CoProSingle` object with the `normalizedData` slot populated.
#' @param nPCA Number of principal components to compute for each cell type.
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#' @param scalePCs Whether to scale (whiten) PCs by their standard deviation
#'   before downstream CCA optimization. Default \code{TRUE} (recommended).
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'   Default is "raw". For single slide, this argument is ignored.
#' @param center_per_slide After the global PCA, do we do center per slide
#'   again? By default this is set to FALSE
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

#' @param object A `CoProMulti` object with the `normalizedData` slot populated.
#' @param nPCA Number of principal components to compute for each cell type.
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'   Default is "raw". For single slide, this argument is ignored.
#' @param center_per_slide After the global PCA, do we do center per slide
#'   again? By default this is set to FALSE
#' @rdname computePCA
#' @aliases computePCA,CoProMulti-method
#' @export
setMethod("computePCA", "CoProMulti",
          function(object, nPCA = 40, center = TRUE, scale. = TRUE,
                   scalePCs = TRUE, dataUse = "raw", center_per_slide = FALSE) {
            cts <- .check_pca_input_multi(object, nPCA, center, scale., scalePCs, dataUse)
            object <- .compute_pca_multi(object, nPCA, center, scale., scalePCs,
                                         dataUse, center_per_slide, cts)
            return(object)
          })
