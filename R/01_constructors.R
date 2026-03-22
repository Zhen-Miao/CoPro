# Note: class unions are defined in 00-classes.R; not duplicated here

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
    "ANY", "matrixOrDataFrame",
    "data.frame", "factorOrCharacter"
  ),
  function(normalizedData, locationData, metaData, cellTypes) {
    ## validate normalizedData type
    if (!is.matrix(normalizedData) && !inherits(normalizedData, "sparseMatrix") &&
        !inherits(normalizedData, "IterableMatrix")) {
      stop("normalizedData must be a matrix, sparse matrix, or BPCells IterableMatrix")
    }

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
    "ANY", "matrixOrDataFrame", "data.frame",
    "factorOrCharacter", "factorOrCharacter"
  ),
  function(normalizedData, locationData, metaData, cellTypes, slideID) {

    # --- Input Validation ---
    # validate normalizedData type
    if (!is.matrix(normalizedData) && !inherits(normalizedData, "sparseMatrix") &&
        !inherits(normalizedData, "IterableMatrix")) {
      stop("normalizedData must be a matrix, sparse matrix, or BPCells IterableMatrix")
    }

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


#' Internal: partition locations into approximately equal-sized spatial blocks
#'
#' Uses k-means on coordinates x, y, (and z) and recursively subdivides any cluster
#' that still exceeds maxCell.
#'
#' @param locationData data.frame with columns x,y,(optional z)
#' @param n target number of blocks (>= 1)
#' @param maxCell maximum cells per block (>= 1)
#' @return character vector of block labels length nrow(locationData)
.partitionByLocation <- function(locationData, n, maxCell) {
  stopifnot(is.data.frame(locationData))
  cn <- tolower(colnames(locationData))
  if (!all(c("x", "y") %in% cn)) {
    stop("locationData must contain columns named 'x' and 'y' (optionally 'z').")
  }
  # Use x/y/(z if present)
  cols <- c("x", "y")
  if ("z" %in% cn) cols <- c(cols, "z")
  # Reorder/standardize columns
  loc <- locationData[, match(cols, cn), drop = FALSE]
  loc <- as.matrix(loc)
  storage.mode(loc) <- "double"

  m <- nrow(loc)
  if (m == 0) return(character(0))
  if (n <= 1L || m <= maxCell) return(rep("1", m))

  # Scale coords so kmeans isn't dominated by one axis
  loc_sc <- scale(loc)

  # seed is set once by the caller (CreateCoPro) before partitioning begins;
  # do not call set.seed() here so that recursive calls consume the RNG
  # naturally and the global seed is only set once.
  km <- stats::kmeans(loc_sc, centers = n, iter.max = 50, nstart = 5)
  lab <- as.character(km$cluster)

  # Recursive refinement: split any cluster that still exceeds maxCell
  tab <- table(lab)
  if (any(tab > maxCell)) {
    # Prepare new labels
    new_lab <- character(m)
    # Maintain stable ordering of clusters
    cl_ids <- names(tab)
    next_id <- 1L

    for (cl in cl_ids) {
      idx <- which(lab == cl)
      if (length(idx) <= maxCell) {
        new_lab[idx] <- as.character(next_id)
        next_id <- next_id + 1L
      } else {
        # Need to split this cluster further
        extra <- ceiling(length(idx) / maxCell)
        sub <- .partitionByLocation(
          locationData = as.data.frame(loc[idx, , drop = FALSE],
                                       stringsAsFactors = FALSE) |>
            stats::setNames(cols),
          n = extra,
          maxCell = maxCell
        )
        # Map sublabels into global ids
        sub_levels <- unique(sub)
        for (sl in sub_levels) {
          sub_idx <- idx[which(sub == sl)]
          new_lab[sub_idx] <- as.character(next_id)
          next_id <- next_id + 1L
        }
      }
    }
    lab <- new_lab
  }

  lab
}


#' Create a CoPro object, automatically choosing Single vs Multi and splitting large slices.
#'
#' Behavior:
#' 1) one slice and ncell < p  -> newCoProSingle()
#' 2) one slice and ncell > p  -> newCoProMulti() with artificial blocks
#' 3) multi slices, all < p    -> newCoProMulti() with original slideID
#' 4) multi slices, some > p   -> newCoProMulti() with refined slideID (origin + block)
#'
#' @param normalizedData matrix (cells x genes)
#' @param locationData data.frame (cells x coords; needs x,y, optionally z)
#' @param metaData data.frame (cells x annotations)
#' @param cellTypes vector length nrow(normalizedData)
#' @param slideID optional vector length nrow(normalizedData); if NULL, treated as one slice
#' @param maxCell integer max cells per (sub)slice (default 50000)
#' @param plot logical; if TRUE, plot the auto-slicing result (default FALSE)
#' @param seed integer seed for reproducible splitting (default 1)
#' @export
CreateCoPro <- function(normalizedData, locationData, metaData, cellTypes,
                        slideID = NULL, maxCell = 50000L, plot = FALSE, seed = 1L) {

  # Basic dimension checks (match your existing constructor expectations)
  n_cells <- nrow(normalizedData)
  if (length(cellTypes) != n_cells ||
      nrow(metaData) != n_cells ||
      nrow(locationData) != n_cells) {
    stop("Input data do not match dimensionality, please check.")
  }

  # Normalize slideID input
  if (is.null(slideID)) {
    slideID <- rep("slice1", n_cells)
  }
  if (!is.character(slideID)) slideID <- as.character(slideID)

  # Safety check consistent with your existing validation
  .validateSeparatorSafety(cellTypes = as.character(cellTypes), slideIDs = slideID)

  u <- unique(slideID)

  # Fast path: truly single and small -> Single
  if (length(u) == 1L && n_cells <= maxCell) {
    return(newCoProSingle(
      normalizedData = normalizedData,
      locationData   = locationData,
      metaData       = metaData,
      cellTypes      = cellTypes
    ))
  }

  # Check which slides actually need splitting
  oversized <- vapply(u, function(sid) sum(slideID == sid) > maxCell, logical(1))

  if (any(oversized)) {
    message(
      "Slide(s) exceeding cell limit (", maxCell, "): ",
      paste(u[oversized], collapse = ", "),
      ". Performing auto slicing..."
    )
  }

  # Save and restore RNG state so we don't alter the caller's random stream
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)

  t0 <- Sys.time()
  newSlideID <- slideID

  for (sid in u) {
    idx <- which(slideID == sid)
    m <- length(idx)
    if (m > maxCell) {
      n_target <- ceiling(m / maxCell)
      block <- .partitionByLocation(
        locationData = locationData[idx, , drop = FALSE],
        n = n_target,
        maxCell = maxCell
      )
      # Combine original id + block id so origin is preserved
      newSlideID[idx] <- paste0(sid, "_blk", block)
    }
  }

  if (any(oversized)) {
    elapsed <- difftime(Sys.time(), t0, units = "secs")
    n_new <- length(unique(newSlideID[slideID %in% u[oversized]]))
    message(
      "Auto slicing done in ", round(as.numeric(elapsed), 2), "s; ",
      sum(oversized), " oversized slide(s) split into ", n_new, " slices; ",
      "total slides: ", length(unique(newSlideID))
    )
  }

  # visual check
  if (plot) {
    .plotAutoSlicing(locationData, slideID, newSlideID)
  }


  # If after refinement there's still exactly one group and it's small,
  # we could return Single, but your desired behavior says:
  # - if original was one slice AND it was > p => Multi with artificial IDs.
  # - multi slices => Multi anyway.
  # So we only return Single in the earlier fast path.
  return(newCoProMulti(
    normalizedData = normalizedData,
    locationData   = locationData,
    metaData       = metaData,
    cellTypes      = cellTypes,
    slideID        = newSlideID
  ))
}


# a quick visual check for k-mean based splitting
.plotAutoSlicing <- function(locationData, slideID, newSlideID) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping plots.")
    return(invisible(NULL))
  }

  df <- data.frame(
    x = locationData$x,
    y = locationData$y,
    slideID = as.character(slideID),
    newSlideID = as.character(newSlideID)
  )

  for (sid in unique(df$slideID)) {
    d <- df[df$slideID == sid, , drop = FALSE]

    p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, color = newSlideID)) +
      ggplot2::geom_point(size = 0.4, alpha = 0.7) +
      ggplot2::labs(
        title = paste0("Auto-slicing for slice: ", sid),
        color = "Block ID"
      ) +
      ggplot2::theme_bw()

    print(p)
  }

  invisible(NULL)
}
