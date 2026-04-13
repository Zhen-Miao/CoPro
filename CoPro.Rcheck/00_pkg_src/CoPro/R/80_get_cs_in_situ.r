#' Get cell score and location information as a data.frame
#'
#' @importFrom stats median
#' @param object A `CoPro` object
#' @param sigmaValueChoice A value to specify the sigma squared to
#' use for selecting the particular cell score information
#' @param ccIndex Canonical vector index, default = 1
#'
#' @return A data.frame object with cell scores and their locations
#' @rdname getCellScoresInSitu
#' @aliases getCellScoresInSitu,CoProSingle-method
#' @aliases getCellScoresInSitu,CoProMulti-method
#' @export
setGeneric("getCellScoresInSitu",
           function(object, sigmaValueChoice, ccIndex = 1
                    ) standardGeneric("getCellScoresInSitu")
)

.checkInputCs <- function(object, sigmaValueChoice, ccIndex = 1){
  ## make sure cellScores and geneScores exist
  if (length(object@cellScores) == 0 ||
      length(object@geneScores) == 0) {
    stop(paste(
      "cellScores slot does not exist,",
      "run `computeGeneAndCellScores()` first"
    ))
  }

  if (is.null(sigmaValueChoice)) {
    if (length(object@sigmaValueChoice) == 0) {
      stop(paste(
        "sigmaValueChoice is not given,",
        "and sigmaValueChoice not computed,",
        "please either specify a particular sigmaValueChoice or",
        "run computeNormalizedCorrelation()"
      ))
    } else {
      warning(paste(
        "sigmaValueChoice is not given,",
        "using default value with highest normalized correlation."
      ))
      sigmaValueChoice <- object@sigmaValueChoice
    }
  }

  if (!(sigmaValueChoice %in% object@sigmaValues)) {
    stop("sigmaValueChoice does not exist in the list of sigmaValues")
  }

  # Validate ccIndex using flat structure
  cell_types <- object@cellTypesOfInterest
  if (length(cell_types) == 0) {
    cell_types <- unique(object@cellTypesSub)
  }
  
  # Try to find any cell scores matrix to validate ccIndex
  found_valid_matrix <- FALSE
  for (ct in cell_types) {
    flat_name <- .createCellScoresName(sigmaValueChoice, ct, slide = NULL)
    if (flat_name %in% names(object@cellScores)) {
      scores_matrix <- object@cellScores[[flat_name]]
      if (!is.null(scores_matrix) && ncol(scores_matrix) > 0) {
        if (ccIndex > ncol(scores_matrix)) {
          stop("ccIndex exceeds number of canonical components")
        }
        found_valid_matrix <- TRUE
        break
      }
    }
  }
  
  if (!found_valid_matrix) {
    stop("No cell scores found for validation")
  }

  return(TRUE)
}

.getCellScoresInSituCore <- function(object, sigmaValueChoice, ccIndex = 1){
  ## Safely detect single vs multi-slide based on object class
  if (inherits(object, "CoProSingle")) {
    return(.csSingleSlide(object, sigmaValueChoice, ccIndex = ccIndex))
  } else if (inherits(object, "CoProMulti")) {
    return(.csMultiSlide(object, sigmaValueChoice, ccIndex = ccIndex))
  } else {
    # Fallback: try to detect based on slideID slot existence
    has_slideList <- tryCatch({
      slideList_data <- object@slideList
      length(slideList_data) > 0 && length(unique(slideList_data)) > 1
    }, error = function(e) FALSE)
    
    if (has_slideList) {
      return(.csMultiSlide(object, sigmaValueChoice, ccIndex = ccIndex))
    } else {
      return(.csSingleSlide(object, sigmaValueChoice, ccIndex = ccIndex))
    }
  }
}

.csSingleSlide <- function(object, sigmaValueChoice, ccIndex = 1){
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

  loc_t <- stats::setNames(
    vector(mode = "list", length = length(cts)),
    cts
  )
  median_score_t <- vector("numeric", length = length(cts))
  names(median_score_t) <- cts

  for (ct in cts) {
    loc_t[[ct]] <- object@locationDataSub[object@cellTypesSub == ct, , drop = FALSE]
    
    # Check if we have cells for this cell type
    if (nrow(loc_t[[ct]]) == 0) {
      warning(paste("No cells found for cell type", ct))
      next
    }
    
    # Use flat structure to access cell scores
    flat_name <- .createCellScoresName(sigmaValueChoice, ct, slide = NULL)
    if (!flat_name %in% names(object@cellScores)) {
      warning(paste("Cell scores not found for cell type", ct))
      next
    }
    
    cell_scores_matrix <- object@cellScores[[flat_name]]
    if (is.null(cell_scores_matrix)) {
      warning(paste("Cell scores matrix is NULL for cell type", ct))
      next
    }
    
    common_cells <- intersect(rownames(loc_t[[ct]]), rownames(cell_scores_matrix))
    
    if (length(common_cells) == 0) {
      warning(paste("No matching cells between location data and cell scores for cell type", ct))
      next
    }
    
    # Subset to common cells
    loc_t[[ct]] <- loc_t[[ct]][common_cells, , drop = FALSE]
    
    loc_t[[ct]]$"cellScores" <- cell_scores_matrix[common_cells, ccIndex]
    loc_t[[ct]]$"cellTypesSub" <- ct
    median_score_t[ct] <- median(loc_t[[ct]]$"cellScores", na.rm = TRUE)
    loc_t[[ct]]$"cellScores_b" <-
      ifelse(loc_t[[ct]]$"cellScores" > median_score_t[ct],
             paste0("high_", ct), paste0("low_", ct)
             )
  }

  # Remove empty list elements
  loc_t <- loc_t[sapply(loc_t, function(x) !is.null(x) && nrow(x) > 0)]
  
  if (length(loc_t) == 0) {
    stop("No valid data found for any cell type")
  }

  combinations <- expand.grid(c("high", "low"), cts)
  all_binary_scores <- apply(combinations, 1,
                             function(x) paste(x[1], x[2], sep = "_"))

  names(loc_t) <- NULL
  loc_all <- do.call(rbind, loc_t)

  loc_all$cellScores_b <- factor(loc_all$cellScores_b,
    levels = all_binary_scores
  )

  return(loc_all)
}

.csMultiSlide <- function(object, sigmaValueChoice, ccIndex = 1){
  all_slides <- object@slideList
  cts <- object@cellTypesOfInterest
  
  if (length(cts) == 0) {
    warning(paste(
      "no cell type of interest specified,",
      "using all cell types to run the analysis"
    ))
    cts <- unique(object@cellTypesSub)
  }
  
  # Initialize list structure for slide and cell type combinations
  df_by_slide_ct <- vector("list", length = length(all_slides))
  names(df_by_slide_ct) <- all_slides
  
  for (q in all_slides) {
    df_by_slide_ct[[q]] <- vector("list", length = length(cts))
    names(df_by_slide_ct[[q]]) <- cts
  }

  # Note: For multi-slide objects, cell scores are aggregated across slides
  # in the flat structure, so we access them directly without slide-specific keys
  
  # Collect data for each slide and cell type
  for (q in all_slides) {
    for (ct in cts) {
      # Get location data for this slide and cell type
      slide_mask <- object@cellTypesSub == ct & object@metaDataSub$slideID == q
      loc_t <- object@locationDataSub[slide_mask, , drop = FALSE]
      
      # Check if there are cells for this slide and cell type
      if (nrow(loc_t) == 0) {
        next  # Skip if no cells
      }
      
      # Use flat structure to access cell scores (aggregated across slides)
      flat_name <- .createCellScoresName(sigmaValueChoice, ct, slide = NULL)
      if (!flat_name %in% names(object@cellScores)) {
        warning(paste("Cell scores not found for cell type", ct))
        next
      }
      
      cell_scores_matrix <- object@cellScores[[flat_name]]
      
      # Check if the cell scores matrix exists and has the right dimensions
      if (is.null(cell_scores_matrix) || nrow(cell_scores_matrix) == 0) {
        warning(paste("Empty cell scores matrix for cell type", ct))
        next
      }
      
      # Check if ccIndex is valid
      if (ccIndex > ncol(cell_scores_matrix)) {
        stop(paste("ccIndex", ccIndex, "exceeds number of canonical components"))
      }
      
      # Match rownames between location data and cell scores
      common_cells <- intersect(rownames(loc_t), rownames(cell_scores_matrix))
      if (length(common_cells) == 0) {
        warning(paste("No matching cells between location data and cell scores for slide", q, "and cell type", ct))
        next
      }
      
      # Subset to common cells
      loc_t <- loc_t[common_cells, , drop = FALSE]
      
      loc_t$"cellScores" <- cell_scores_matrix[common_cells, ccIndex]
      loc_t$"cellTypesSub" <- ct
      loc_t$"slideID" <- q
      df_by_slide_ct[[q]][[ct]] <- loc_t
    }
  }

  # Flatten the nested structure and combine all data
  all_data_list <- list()
  for (q in all_slides) {
    for (ct in cts) {
      if (!is.null(df_by_slide_ct[[q]][[ct]]) && nrow(df_by_slide_ct[[q]][[ct]]) > 0) {
        all_data_list <- append(all_data_list, list(df_by_slide_ct[[q]][[ct]]))
      }
    }
  }
  
  if (length(all_data_list) == 0) {
    stop("No valid data found across all slides and cell types")
  }

  combined_loc_t <- do.call(rbind, all_data_list)
  
  # Check if we have any cell scores
  if (nrow(combined_loc_t) == 0 || all(is.na(combined_loc_t$cellScores))) {
    stop("No valid cell scores found")
  }
  
  # Calculate median scores per cell type (across all slides)
  median_score_t <- vector("numeric", length = length(cts))
  names(median_score_t) <- cts
  
  for (ct in cts) {
    ct_data <- combined_loc_t[combined_loc_t$cellTypesSub == ct, ]
    if (nrow(ct_data) > 0 && !all(is.na(ct_data$cellScores))) {
      median_score_t[ct] <- median(ct_data$cellScores, na.rm = TRUE)
    } else {
      median_score_t[ct] <- NA
    }
  }

  # Create binary scores per cell type
  combined_loc_t$"cellScores_b" <- NA
  for (ct in cts) {
    ct_mask <- combined_loc_t$cellTypesSub == ct
    if (any(ct_mask) && !is.na(median_score_t[ct])) {
      combined_loc_t$cellScores_b[ct_mask] <- ifelse(
        is.na(combined_loc_t$cellScores[ct_mask]), 
        paste0("unknown_", ct),
        ifelse(combined_loc_t$cellScores[ct_mask] > median_score_t[ct],
               paste0("high_", ct), 
               paste0("low_", ct))
      )
    } else {
      # If no valid median could be calculated, mark as unknown
      combined_loc_t$cellScores_b[ct_mask] <- paste0("unknown_", ct)
    }
  }

  # Create factor levels for all cell types
  combinations <- expand.grid(c("high", "low", "unknown"), cts)
  all_binary_scores <- apply(combinations, 1,
                             function(x) paste(x[1], x[2], sep = "_"))

  combined_loc_t$cellScores_b <- factor(combined_loc_t$cellScores_b,
                                        levels = all_binary_scores
  )

  # Reorder the data frame to match the original order in metaDataSub
  if (is.null(rownames(combined_loc_t)) || length(rownames(combined_loc_t)) == 0) {
    stop("combined_loc_t has no rownames. This indicates a problem with the data structure.")
  }
  
  cell_ids_meta <- rownames(object@metaDataSub)
  common_cells <- intersect(cell_ids_meta, rownames(combined_loc_t))
  if (length(common_cells) > 0) {
    # Reorder to match the original order in metaDataSub
    # common_cells is already ordered by cell_ids_meta (from intersect above)
    # Reorder combined_loc_t rows to match the metaDataSub order
    meta_order <- cell_ids_meta[cell_ids_meta %in% rownames(combined_loc_t)]
    combined_loc_t <- combined_loc_t[meta_order, , drop = FALSE]
  } else {
    warning("No common cell IDs found between combined_loc_t and metaDataSub")
  }

  return(combined_loc_t)
}

#' @rdname getCellScoresInSitu
#' @aliases getCellScoresInSitu,CoProSingle-method
#' @export
setMethod("getCellScoresInSitu", "CoProSingle", function(object, sigmaValueChoice, ccIndex = 1) {
  .checkInputCs(object, sigmaValueChoice, ccIndex = ccIndex)
  .getCellScoresInSituCore(object = object, sigmaValueChoice = sigmaValueChoice, ccIndex = ccIndex)
})

#' @rdname getCellScoresInSitu
#' @aliases getCellScoresInSitu,CoProMulti-method
#' @export
setMethod("getCellScoresInSitu", "CoProMulti", function(object, sigmaValueChoice, ccIndex = 1) {
  .checkInputCs(object, sigmaValueChoice, ccIndex = ccIndex)
  .getCellScoresInSituCore(object = object, sigmaValueChoice = sigmaValueChoice, ccIndex = ccIndex)
})




