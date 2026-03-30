#' Compute PCA on Integrated Multi-Slide Data
#'
#' Performs PCA on the integrated data stored within the `CoProMulti` object.
#' Assumes integration has created a common space across slides.
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
#'   before downstream CCA optimization. Default `TRUE` (recommended). When
#'   PCs are whitened, the unit-norm constraint `||w|| = 1` is equivalent to
#'   the standard CCA constraint `||w'X'Xw|| = 1`. Setting this to `FALSE`
#'   distorts the constraint space and may produce unreliable results.
#'
#' @return A `CoProMulti` object with the `pcaResults` slot populated.
#'         `pcaResults` structure: `list(slideID = list(cellType = pc_matrix))`.
#' @export
#' @rdname computePCAMulti
#' @aliases computePCAMulti,CoProMulti-method
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
  if (!scalePCs) {
    warning(paste(
      "scalePCs = FALSE is not recommended. When PCs are whitened (scaled),",
      "the unit-norm constraint ||w|| = 1 is equivalent to the standard CCA",
      "constraint ||w'X'Xw|| = 1. Without whitening, the optimization may",
      "produce unreliable results because the constraint space is distorted."
    ))
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

.compute_pca_single <- function(object, nPCA = 40, center = TRUE, scale. = TRUE, scalePCs = TRUE, cts) {
  # PCA results will be saved under the name of cell types
  object@pcaGlobal <- setNames(
    vector("list", length = length(cts)),
    cts
  )

  # Iterate over cell types
  for (ct in cts) {
    # Cell type specific subset
    sub_data <- object@normalizedDataSub[object@cellTypesSub == ct, ]

    # Apply centering and scaling
    scaled_data <- .apply_centering_scaling(sub_data, center, scale.)

    # Guard against nPCA exceeding data dimensions
    max_pca <- min(nrow(scaled_data) - 1, ncol(scaled_data))
    if (nPCA >= max_pca) {
      warning(paste0("nPCA (", nPCA, ") exceeds max allowed (", max_pca,
                     ") for cell type '", ct, "'. Reducing to ", max(1, max_pca - 1), "."))
      nPCA_use <- max(1, max_pca - 1)
    } else {
      nPCA_use <- nPCA
    }

    # PCA on the matrix that is already centered and scaled
    if (.is_bpcells(scaled_data)) {
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
      message("Input is dense (", class(scaled_data), "), performing irlba pca...")
      pca <- prcomp_irlba(scaled_data, center = FALSE, scale. = FALSE, n = nPCA_use)
    }
    object@pcaGlobal[[ct]] <- pca
  }

  object@nPCA <- nPCA
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
    stop("cellTypesOfInterest not set. Run subsetDataMulti first.")
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

  # Perform PCA per cell type on the integrated data
  for (ct in cts) {
    message("Performing PCA for cell type: ", ct)

    if (dataUse == "integrated" && !ct %in% names(object@integratedData)) {
      warning("No integrated data found for cell type: ", ct, " - Skipping PCA.")
      next
    }

    # Get the appropriate data matrix
    if (dataUse == "integrated") {
      mat_ct <- object@integratedData[[ct]]
    } else { # raw
      mat_ct <- object@normalizedDataSub[object@cellTypesSub == ct, ]
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

    # Apply centering and scaling
    scaled_data <- .apply_centering_scaling(mat_ct, center, scale.)
    if (center || scale.) {
      message("Data centered and/or scaled")
    }

    # Guard against nPCA exceeding data dimensions
    max_pca <- min(nrow(scaled_data) - 1, ncol(scaled_data))
    if (nPCA >= max_pca) {
      warning(paste0("nPCA (", nPCA, ") exceeds max allowed (", max_pca,
                     ") for cell type '", ct, "'. Reducing to ", max(1, max_pca - 1), "."))
      nPCA_use <- max(1, max_pca - 1)
    } else {
      nPCA_use <- nPCA
    }

    # Perform PCA on the combined integrated data for this cell type
    if (.is_bpcells(scaled_data)) {
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
      pca_ct <- prcomp_irlba(scaled_data, n = nPCA_use,
                             center = FALSE, scale. = FALSE)
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
  object@nPCA <- nPCA
  object@scalePCs <- scalePCs

  return(object)
}

#' Compute PCA on Single-Slide Data
#'
#' This method performs PCA on the normalized data stored within the `CoProSingle` object.
#' It assumes that the data has already been integrated across slides.
#'
#' @param object A `CoProSingle` object with the `normalizedData` slot populated.
#' @param nPCA Number of principal components to compute for each cell type.
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
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

#' Compute PCA on Multi-Slide Data
#'
#' This method performs PCA on the normalized data stored within the `CoProMulti` object.
#' It assumes that the data has already been integrated across slides.
#'
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
