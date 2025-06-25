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
  # Calculate column standard deviations
  col_means <- colMeans(input_matrix)
  col_sds <- apply(input_matrix, 2, sd)

  # non-zero proportion
  col_nz <- colSums(input_matrix != 0) / nrow(input_matrix)

  # Identify columns that are not full of zeros (to avoid division by zero)

  zero_sd_cols <- which(col_sds < zero_sd_threshold |
                        col_nz < nz_propion_threshold)

  ## do not scale if the sd is too small, or if the proportion of non-zero
  # values is too low
  col_sds_safe <- col_sds
  if (length(zero_sd_cols) > 0) {
    col_sds_safe[zero_sd_cols] <- 1.0
  }

  scaled_matrix <- scale(input_matrix, center = col_means, scale = col_sds_safe)


  return(scaled_matrix)
}



#' Normalize a vector to unit length
#' @param v Input vector or matrix (if matrix, treated as column vector)
#' @return Normalized vector as column matrix
#' @noRd
normalize_vec <- function(v) {
  if (is.matrix(v)) {
    v_norm <- sqrt(sum(v^2))
  } else {
    v_norm <- sqrt(sum(v^2))
  }
  
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

    if (names(distanceList) != cts) {
      stop(paste(
        "The names of distanceList do not match cell types",
        "of interest"
      ))
    }

    for (i in cts) {
      if (names(distanceList[[i]]) != cts) {
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
  object@slideID
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
