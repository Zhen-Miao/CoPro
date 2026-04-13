## get PC matrices
.getAllPCMats <- function(allPCs, scalePCs) {

  if (length(allPCs) == 0) {
    stop("PCA results do not exist, run computePCA() first.")
  }

  PCmats <- setNames(
    vector("list", length = length(allPCs)),
    names(allPCs)
  )

  ## optionally, scale the PCs before running CCA
  if (scalePCs) {
    for (i in names(allPCs)) {
      pca_A_sd <- allPCs[[i]]$sdev
      PCmats[[i]] <- scale(allPCs[[i]]$x,
                           center = FALSE,
                           scale = pca_A_sd
      )
    }
  } else {
    for (i in names(allPCs)) {
      PCmats[[i]] <- allPCs[[i]]$x
    }
  }
  return(PCmats)
}


#' centering and scaling the matrix
#' @importFrom stats sd
#' @importFrom matrixStats colSds
#' @param matrix Input matrix to be column-centered
#'
#' @return centered and scaled matrix
#' @noRd
center_scale_matrix_opt <- function(input_matrix,
                                    zero_sd_threshold = 1e-3,
                                    nz_propion_threshold = 0.01) {

  if (!.is_bpcells(input_matrix)) {
    # Original behavior for base matrix / Matrix::dgCMatrix
    col_means <- colMeans(input_matrix)
    col_sds   <- apply(input_matrix, 2, sd)
    col_nz    <- colSums(input_matrix != 0) / nrow(input_matrix)

    zero_sd_cols <- which(col_sds < zero_sd_threshold | col_nz < nz_propion_threshold)
    col_sds_safe <- col_sds
    if (length(zero_sd_cols) > 0) col_sds_safe[zero_sd_cols] <- 1.0

    return(scale(input_matrix, center = col_means, scale = col_sds_safe))
  }

  # ---- BPCells path ----

  col_means <- colMeans(input_matrix)
  col_sds <- sqrt(BPCells::colVars(input_matrix))
  col_nz <- colSums(BPCells::binarize(input_matrix)) / nrow(input_matrix)
  zero_sd_cols <- which(col_sds < zero_sd_threshold | col_nz < nz_propion_threshold)
  col_sds_safe <- col_sds
  if (length(zero_sd_cols) > 0) col_sds_safe[zero_sd_cols] <- 1.0

  # Center then scale using BPCells broadcasting (no base::scale())
  centered <- BPCells::add_cols(input_matrix, -col_means)
  scaled   <- BPCells::multiply_cols(centered, 1 / col_sds_safe)

  scaled
}



#' Normalize a vector to unit length
#' @param v Input vector or matrix (if matrix, treated as column vector)
#' @return Normalized vector as column matrix
#' @noRd
normalize_vec <- function(v) {
  v_norm <- sqrt(sum(v^2))

  if (v_norm < 1e-12) {
    warning("Vector has very small norm, may cause numerical issues")
    return(matrix(0, nrow = length(v), ncol = 1))
  }

  normalized <- v / v_norm
  if (is.matrix(v)) {
    return(normalized)
  } else {
    return(matrix(normalized, ncol = 1))
  }
}

#' Normalize w under the CCA constraint w'(X'X)w = 1
#'
#' When \code{sdev2} is \code{NULL} (whitened PCs), this reduces to
#' \code{||w|| = 1}. Otherwise the weighted norm
#' \code{sqrt(sum(w^2 * sdev2))} is used so that the constraint becomes
#' \code{w' diag(sdev^2) w = 1}, which is equivalent to \code{||Xw|| = 1}
#' in the original (unwhitened) PC space.
#'
#' @param v Numeric vector or single-column matrix
#' @param sdev2 Squared sdev vector (same length as v), or NULL for standard norm
#' @return Normalized vector as column matrix
#' @noRd
normalize_vec_weighted <- function(v, sdev2 = NULL) {
  if (is.null(sdev2)) {
    return(normalize_vec(v))
  }
  v_norm <- sqrt(sum(as.numeric(v)^2 * sdev2))
  if (v_norm < 1e-12) {
    warning("Vector has very small weighted norm, may cause numerical issues")
    return(matrix(0, nrow = length(v), ncol = 1))
  }
  normalized <- v / v_norm
  if (is.matrix(v)) {
    return(normalized)
  } else {
    return(matrix(normalized, ncol = 1))
  }
}

#' Apply D-inverse preconditioner to gradient, then normalize under D-norm
#'
#' For the generalized eigenvalue problem \code{Yw = λDw}, the correct power
#' iteration is: (1) compute gradient \code{v = Yw}, (2) apply
#' \code{D^{-1}}: \code{u = v / sdev2}, (3) normalize so that
#' \code{u'Du = 1}. This function performs steps 2--3 in one call.
#' When \code{sdev2} is \code{NULL}, falls back to standard normalization.
#'
#' @param v Gradient vector (single-column matrix)
#' @param sdev2 Squared sdev vector, or NULL for standard norm
#' @return Preconditioned and normalized vector as column matrix
#' @noRd
normalize_gradient_weighted <- function(v, sdev2 = NULL) {
  if (is.null(sdev2)) {
    return(normalize_vec(v))
  }
  # D^{-1} v
  u <- as.numeric(v) / sdev2
  # Normalize under D-norm: sqrt(u' D u) = sqrt(sum(u^2 * sdev2))
  normalize_vec_weighted(matrix(u, ncol = 1), sdev2)
}

# Declare globals used across the package to quiet R CMD check NOTES
utils::globalVariables(c(
  ".computeSpatialCrossCorrelation",
  "getKernelMatrix",
  "getSlideList",
  "x",
  "y"
))

#' Validate that cell types and slide IDs don't contain pipe characters
#' @param cellTypes Character vector of cell type names
#' @param slideIDs Character vector of slide IDs (optional)
#' @return TRUE if valid, stops execution with error if invalid
#' @noRd
.validateSeparatorSafety <- function(cellTypes = NULL, slideIDs = NULL) {

  # Check cell types for pipe characters
  if (!is.null(cellTypes)) {
    cellTypes <- as.character(cellTypes)
    pipe_in_cellTypes <- grepl("\\|", cellTypes)
    if (any(pipe_in_cellTypes)) {
      problematic_types <- cellTypes[pipe_in_cellTypes]
      stop(paste("Cell type names cannot contain pipe characters (|).",
                 "Problematic cell types:", paste(problematic_types, collapse = ", "),
                 "\nPlease rename these cell types to avoid conflicts with internal naming."))
    }
  }

  # Check slide IDs for pipe characters
  if (!is.null(slideIDs)) {
    slideIDs <- as.character(slideIDs)
    pipe_in_slideIDs <- grepl("\\|", slideIDs)
    if (any(pipe_in_slideIDs)) {
      problematic_slides <- slideIDs[pipe_in_slideIDs]
      stop(paste("Slide IDs cannot contain pipe characters (|).",
                 "Problematic slide IDs:", paste(problematic_slides, collapse = ", "),
                 "\nPlease rename these slide IDs to avoid conflicts with internal naming."))
    }
  }

  return(TRUE)
}



#' Assign distance matrix manually
#'
#' @param object A `CoPro` object
#' @param distanceList A list object that contains all pairwise distances
#' between any two pairs of cells.
#'
#' @return A `CoPro` object with specified
#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
#'
setGeneric("assignDistanceManually",
           function(object,
                    distanceList) standardGeneric("assignDistanceManually")
)


#' @rdname assignDistanceManually
#' @aliases assignDistanceManually,CoPro-method
#' @export
setMethod(
  "assignDistanceManually", "CoPro",
  function(object, distanceList) {
    if (!is.list(distanceList)) {
      stop(paste(
        "distanceList must be a nested list object with names",
        "specified by cell types"
      ))
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

    if (!identical(names(distanceList), cts)) {
      stop(paste(
        "The names of distanceList do not match cell types",
        "of interest"
      ))
    }

    for (i in cts) {
      if (!identical(names(distanceList[[i]]), cts)) {
        stop(paste("The names of distanceList[[", i,
          "]] do not match cell types ",
          "of interest",
          sep = ""
        ))
      }
    }

    object@distances <- distanceList
    return(object)
  }
)

#' Get slide IDs from CoPro object
#'
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @return For CoProMulti: the slideID vector; For CoProSingle: character(0)
#' @export
#' @rdname getSlideID
setGeneric("getSlideID", function(object) standardGeneric("getSlideID"))

#' @rdname getSlideID
setMethod("getSlideID", "CoProSingle", function(object) {
  character(0)
})

#' @rdname getSlideID
setMethod("getSlideID", "CoProMulti", function(object) {
  if ("slideID" %in% colnames(object@metaDataSub)) {
    return(object@metaDataSub$slideID)
  } else {
    return(character(0))
  }
})

#' Get slide list from CoPro object
#'
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @return For CoProMulti: the slideList vector; For CoProSingle: character(0)
#' @export
#' @rdname getSlideList
setGeneric("getSlideList", function(object) standardGeneric("getSlideList"))

#' @rdname getSlideList
setMethod("getSlideList", "CoProSingle", function(object) {
  character(0)
})

#' @rdname getSlideList
setMethod("getSlideList", "CoProMulti", function(object) {
  object@slideList
})

#' Check if object is multi-slide
#'
#' @param object A CoPro object (CoProSingle or CoProMulti)
#' @return Logical indicating if object contains multiple slides
#' @export
#' @rdname isMultiSlide
setGeneric("isMultiSlide", function(object) standardGeneric("isMultiSlide"))

#' @rdname isMultiSlide
setMethod("isMultiSlide", "CoProSingle", function(object) {
  FALSE
})

#' @rdname isMultiSlide
setMethod("isMultiSlide", "CoProMulti", function(object) {
  TRUE
})
