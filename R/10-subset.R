#' subsetData
#'
#' Take a subset of the original matrix based on cell types of interest. The
#' original data are stored without being thrown away
#'
#' @param object A `CoPro` object
#' @param cellTypesOfInterest Input cell types of interest as a vector of
#' characters for subsetting the data
#' @param saveOriginal Logical, whether to save the original data in the
#' subsetted object. Default is FALSE.
#'
#' @rdname subsetData
#' @aliases subsetData,CoProSingle-method
#' @aliases subsetData,CoProMulti-method
#' @return A `CoPro` object with subset slots
#' @export
#'
setGeneric("subsetData",
           function(object,
                    cellTypesOfInterest, 
                    saveOriginal = FALSE) standardGeneric("subsetData")
)

.subset_core <- function(object, ctoi, min_n = 10, saveOriginal = FALSE) {
  if (length(ctoi) < 1) {
    stop("Please specify at least one cell type of interest.")
  }

  if (!all(ctoi %in% object@cellTypes)) {
    stop("some cellTypesOfInterest are not in cellTypes, please check")
  }

  idx <- object@cellTypes %in% ctoi
  if (sum(idx) < min_n) stop("Fewer than ", min_n, " cells after subsetting.")

  ## subset the data
  object@cellTypesOfInterest <- ctoi
  object@normalizedDataSub <- object@normalizedData[idx, , drop = FALSE]
  object@metaDataSub <- object@metaData[idx, , drop = FALSE]
  object@locationDataSub <- object@locationData[idx, , drop = FALSE]
  object@cellTypesSub <- object@cellTypes[idx]
  
  # Only subset slideID if this is a multi-slide object
  if (isMultiSlide(object)) {
    object@slideID <- object@slideID[idx]
  }
  
  if (!saveOriginal) {
    object@normalizedData <- matrix(0)
    object@metaData <- data.frame(0)
    object@locationData <- data.frame(0)
  }

  return(object)

}


setMethod("subsetData", "CoProSingle", function(object, cellTypesOfInterest, saveOriginal = FALSE) {
  .subset_core(object = object, ctoi = cellTypesOfInterest, min_n = 10, saveOriginal = saveOriginal)
})

setMethod("subsetData", "CoProMulti", function(object, cellTypesOfInterest, saveOriginal = FALSE) {
  n_slide = length(object@slideList)
  .subset_core(object, cellTypesOfInterest, min_n = 10 * n_slide, saveOriginal = saveOriginal)
})
