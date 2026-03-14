#' @import Matrix
#' @import methods
# Define necessary class unions from the original CoPro object
setClassUnion("matrixOrSparseMatrix", c("matrix", "dgCMatrix", "dgTMatrix"))
setClassUnion("factorOrCharacter", c("factor", "character"))
setClassUnion("matrixOrDataFrame", c("matrix", "data.frame"))

#' Function to create a new object
#' @importFrom methods new
#' @param normalizedData A `matrix` object to store normalized data.
#' @param locationData A `data.frame` object to store the location. It should
#' either contain two columns named by "x" and "y", or three columns named by
#' "x", "y", and "z". No other names allowed
#' @param metaData A `data.frame` object to store metadata for each cell.
#' @param cellTypes A `vector` object with elements being character. It should
#' match the number of cells in the data matrix and each represents a cell type
#' label of a cell.
#' @rdname newCoProSingle
#' @aliases newCoProSingle,CoPro-method
#' @return A `CoPro` object
#' @export
#'
setGeneric(
  "newCoProSingle",
  function(normalizedData, locationData, metaData,
           cellTypes)  standardGeneric("newCoProSingle")
)


#' @rdname newCoProSingle
#' @aliases newCoProSingle,CoPro-method
#' @export
setMethod(
  "newCoProSingle", signature(
    "matrixOrSparseMatrix", "matrixOrDataFrame",
    "data.frame", "factorOrCharacter"
  ),
  function(normalizedData, locationData, metaData, cellTypes) {
    ## check dimension of input
    if (length(cellTypes) != nrow(normalizedData) ||
        nrow(normalizedData) != nrow(metaData) ||
        nrow(normalizedData) != nrow(locationData)) {
      stop("input data do not match dimensionality, please check")
    }

    ## check the format of location data
    if (!all(c("x", "y") %in% tolower(colnames(locationData)))) {
      stop(paste("locationData should contain x, y, (or z)",
                 "axis info and colnames should be named accordingly"
      ))
    }
    colnames(locationData) <- tolower(colnames(locationData))

    ## convert cellTypes to characters
    if (!is.character(cellTypes)) {
      cellTypes <- as.character(cellTypes)
    }
    
    ## validate cell types don't contain pipe characters
    .validateSeparatorSafety(cellTypes = cellTypes)

    ## convert locationData to data.frame
    if (is.matrix(locationData)) {
      locationData <- as.data.frame(locationData)
    }

    ## check cell id and gene names
    if (is.null(rownames(metaData)) || is.null(rownames(normalizedData)) ||
        is.null(rownames(locationData))) {
      stop(paste("please make sure the rownames of data,",
                 "metaData, and locationData are cell barcodes",
                 sep = " "
      ))
    } else if (any(rownames(metaData) != rownames(normalizedData)) ||
                   any(rownames(locationData) != rownames(normalizedData))) {
      stop(paste("please make sure the cell barcodes match,",
                 "between data, metaData,and locationData",
                 sep = " "
      ))
    }

    ## check gene names
    if (is.null(colnames(normalizedData))) {
      stop("please make sure colnames of data are gene names")
    }

    geneList <- colnames(normalizedData)

    ## create new object
    methods::new("CoProSingle",
                 normalizedData = normalizedData,
                 metaData = metaData, locationData = locationData,
                 cellTypes = cellTypes, geneList = geneList
    )
  }
)


#' Create a new CoProMulti object for Multi-Slide Analysis
#'
#' Initializes a `CoProMulti` object with combined data from multiple slides.
#'
#' @param normalizedData Combined normalized expression matrix (cells x genes) for all slides. Rownames should be unique cell identifiers.
#' @param locationData Combined location data frame (cells x coordinates) for all slides. Rownames must match `normalizedData`. Columns 'x', 'y', (and optionally 'z') required.
#' @param metaData Combined metadata data frame (cells x annotations) for all slides. Rownames must match `normalizedData`.
#' @param cellTypes Combined cell type labels vector for all cells. Length must match `nrow(normalizedData)`.
#' @param slideID Combined slide/sample identifier vector for all cells. Length must match `nrow(normalizedData)`.
#'
#' @return A `CoProMulti` object.
#' @export
#' @rdname newCoProMulti
#' @aliases newCoProMulti,CoProMulti-method
setGeneric(
  "newCoProMulti",
  function(normalizedData, locationData, metaData,
           cellTypes, slideID) standardGeneric("newCoProMulti")
)

#' @rdname newCoProMulti
#' @aliases newCoProMulti,CoProMulti-method
#' @export
setMethod(
  "newCoProMulti", signature(
    "matrixOrSparseMatrix", "matrixOrDataFrame", "data.frame",
    "factorOrCharacter", "factorOrCharacter"
  ),
  function(normalizedData, locationData, metaData, cellTypes, slideID) {

    # --- Input Validation ---
    n_cells <- nrow(normalizedData)
    if (length(cellTypes) != n_cells || nrow(metaData) != n_cells ||
        nrow(locationData) != n_cells || length(slideID) != n_cells) {
      stop("Input data dimensions do not match the number of cells.")
    }

    # Check required columns in locationData
    if (!all(c("x", "y") %in% tolower(colnames(locationData)))) {
      stop("locationData requires columns named 'x' and 'y'.")
    }
    colnames(locationData) <- tolower(colnames(locationData))

    # Check rownames consistency
    if (is.null(rownames(normalizedData)) || is.null(rownames(locationData)) || is.null(rownames(metaData))) {
      stop("Rownames are missing from input data (should be unique cell IDs).")
    }
    if (any(rownames(normalizedData) != rownames(locationData)) || any(rownames(normalizedData) != rownames(metaData))) {
      stop("Rownames mismatch between input data matrices/data frames.")
    }
    if(anyDuplicated(rownames(normalizedData))) {
      stop("Cell IDs (rownames) must be unique across all slides.")
    }

    # Check gene names
    if (is.null(colnames(normalizedData))) {
      stop("colnames of normalizedData (gene names) are missing.")
    }
    geneList <- colnames(normalizedData)

    # Convert factors to characters if necessary
    if (!is.character(cellTypes)) cellTypes <- as.character(cellTypes)
    if (!is.character(slideID)) slideID <- as.character(slideID)
    if (is.matrix(locationData)) locationData <- as.data.frame(locationData) # Ensure data frame
    
    # Validate cell types and slide IDs don't contain pipe characters
    .validateSeparatorSafety(cellTypes = cellTypes, slideIDs = slideID)

    # Get unique slide identifiers
    unique_slides <- unique(slideID)
    if(length(unique_slides) < 2) {
      warning("CoProMulti object created with only one unique slide ID. Multi-slide functions may not be appropriate.")
    }

    # add slideID to metadata, if not already in it
    if("slideID" %in% colnames(metaData)) {
      if(any(metaData[, "slideID", drop=TRUE] != slideID)) stop(
        "metaData contains slideID column, but it does not match slideID")
    }else {
      metaData["slideID"] <- slideID
    }

    # --- Create CoProMulti Object ---
    methods::new("CoProMulti",
                 normalizedData = normalizedData,
                 locationData = locationData,
                 metaData = metaData,
                 cellTypes = cellTypes,
                 slideList = unique_slides,
                 geneList = geneList
    )
  }
)



