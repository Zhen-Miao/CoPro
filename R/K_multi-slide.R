#' @import Matrix
#' @import methods
# Define necessary class unions from the original CoPro object
setClassUnion("matrixOrSparseMatrix", c("matrix", "dgCMatrix", "dgTMatrix"))
setClassUnion("factorOrCharacter", c("factor", "character"))
setClassUnion("matrixOrDataFrame", c("matrix", "data.frame"))

#' CoProm object for Multi-Slide Spatial Transcriptomics Analysis
#'
#' An S4 class to store and manage data for multi-slide spatial transcriptomics
#' analysis using the CoPro methodology.
#'
#' @slot slideID A `factorOrCharacter` vector indicating the slide/sample origin for each cell.
#' @slot slideList A `character` vector storing the unique slide identifiers present.
#' @slot normalizedData A `matrixOrSparseMatrix` object storing the combined normalized data for all cells across slides.
#' @slot locationData A `matrixOrDataFrame` object storing the combined location data.
#' @slot metaData A `data.frame` object storing the combined metadata. Must include columns matching cell barcodes in rows.
#' @slot cellTypes A `factorOrCharacter` vector storing cell type labels for all cells.
#' @slot geneList A `character` vector storing gene names (colnames of `normalizedData`).
#'
#' @slot cellTypesOfInterest A `character` vector specifying the cell types to include in the analysis.
#' @slot normalizedDataSub A `matrixOrSparseMatrix` storing the subset of `normalizedData` for cells matching `cellTypesOfInterest`.
#' @slot locationDataSub A `matrixOrDataFrame` storing the subset of `locationData`.
#' @slot metaDataSub A `data.frame` storing the subset of `metaData`, including the `slideID` column.
#' @slot cellTypesSub A `factorOrCharacter` storing the subset of `cellTypes`.
#'
#' @slot integratedData A `list` object to store the output from data integration (e.g., Seurat integrated assay or corrected values), typically structured by cell type. Format depends on integration method output.
#' @slot pcaResults A `list` object storing PCA results after integration.
#'  Recommended structure: `list(slideID = list(cellType = pc_matrix))`.
#' @slot pcaGlobal A `list` object storing PCA results for each cell type.
#' @slot nPCA A `numeric` value specifying the number of PCs used.
#' @slot scalePCs A `logical` value indicating if PCs were scaled before skrCCA (usually FALSE if PCA done on integrated data).
#'
#' @slot distances A `list` object storing pairwise distances, structured by slide: `list(slideID = list(cellType1 = list(cellType2 = dist_matrix)))`.
#' @slot kernelMatrices A `list` object storing kernel matrices: `list(sigma = list(slideID = list(cellType1 = list(cellType2 = kernel_matrix))))`.
#' @slot sigmaValues A `numeric` vector storing the sigma values used for kernel generation.
#' @slot sigmaValueChoice A `numeric` value storing the chosen sigma based on correlation analysis.
#'
#' @slot nCC A `numeric` value specifying the number of canonical components computed.
#' @slot skrCCAOut A `list` object storing the shared weight vectors from multi-slide skrCCA: `list(sigma = list(cellType = weight_matrix))`.
#' @slot normalizedCorrelation A `list` object storing normalized correlation results. Structure depends on calculation (e.g., `list(sigma = list(slideID = correlation_df))` or `list(sigma = aggregate_correlation_df)`).
#' @slot cellScores A `list` object storing slide-specific cell scores: `list(sigma = list(slideID = list(cellType = cell_score_matrix)))`.
#' @slot geneScores A `list` object storing gene scores (potentially shared): `list(sigma = list(cellType = gene_score_matrix))`.
#'
#' @slot nPermu A `numeric` value specifying the number of permutations (if performed).
#' @slot skrCCAPermuOut A `list` storing permutation results for weights.
#' @slot normalizedCorrelationPermu A `list` storing permutation results for correlation.
#'
#' @export
setClass("CoProm",
         slots = list(
           # Input Data (Combined)
           slideID = "factorOrCharacter",
           slideList = "character",
           normalizedData = "matrixOrSparseMatrix",
           locationData = "matrixOrDataFrame",
           metaData = "data.frame",
           cellTypes = "factorOrCharacter",
           geneList = "character",

           # Subsets based on cellTypesOfInterest
           cellTypesOfInterest = "character",
           normalizedDataSub = "matrixOrSparseMatrix",
           locationDataSub = "matrixOrDataFrame",
           metaDataSub = "data.frame",
           cellTypesSub = "factorOrCharacter",

           # Integration & PCA
           integratedData = "list", # Stores output of integrateSlidesMulti
           pcaResults = "list",     # Stores list(slideID = list(cellType = pc_matrix))
           pcaGlobal = "list",      # Stores list(cellType = pc_result)
           nPCA = "numeric",        # number of PCs
           scalePCs = "logical",    # Whether PCs were scaled before skrCCA

           # Distance & Kernel (Slide-Specific)
           distances = "list",      # Stores list(slideID = list(ct1 = list(ct2 = dist)))
           kernelMatrices = "list", # Stores list(sigma = list(slideID = list(ct1 = list(ct2 = K))))
           sigmaValues = "numeric",
           sigmaValueChoice = "numeric",

           # skrCCA Output (Shared Weights, Slide-Specific Scores)
           nCC = "numeric",
           skrCCAOut = "list",      # Stores list(sigma = list(cellType = W)) - Shared
           normalizedCorrelation = "list", # Store per-slide or aggregate results
           cellScores = "list",     # Stores list(sigma = list(slideID = list(cellType = Scores)))
           geneScores = "list",     # Stores list(sigma = list(cellType = gene_scores)) - Potentially Shared

           # Permutation Slots (Structure TBD based on implementation)
           nPermu = "numeric",
           skrCCAPermuOut = "list",
           normalizedCorrelationPermu = "list"
           # Potentially add cellPermuMulti if needed
         )
)

#' Create a new CoProm object for Multi-Slide Analysis
#'
#' Initializes a `CoProm` object with combined data from multiple slides.
#'
#' @param normalizedData Combined normalized expression matrix (cells x genes) for all slides. Rownames should be unique cell identifiers.
#' @param locationData Combined location data frame (cells x coordinates) for all slides. Rownames must match `normalizedData`. Columns 'x', 'y', (and optionally 'z') required.
#' @param metaData Combined metadata data frame (cells x annotations) for all slides. Rownames must match `normalizedData`.
#' @param cellTypes Combined cell type labels vector for all cells. Length must match `nrow(normalizedData)`.
#' @param slideID Combined slide/sample identifier vector for all cells. Length must match `nrow(normalizedData)`.
#'
#' @return A `CoProm` object.
#' @export
#' @rdname newCoProm
#' @aliases newCoProm,CoProm-method
setGeneric(
  "newCoProm",
  function(normalizedData, locationData, metaData,
           cellTypes, slideID) standardGeneric("newCoProm")
)

#' @rdname newCoProm
#' @aliases newCoProm,CoProm-method
#' @export
setMethod(
  "newCoProm", signature(
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

    # Get unique slide identifiers
    unique_slides <- unique(slideID)
    if(length(unique_slides) < 2) {
      warning("CoProm object created with only one unique slide ID. Multi-slide functions may not be appropriate.")
    }

    # add slideID to metadata, if not already in it
    if("slideID" %in% metaData) {
      if(any(metaData[,"slideID"] != slideID)) stop(
        "metaData contains slideID column, but it does not match slideID")
    }else {
      metaData["slideID"] <- slideID
    }

    # --- Create CoProm Object ---
    methods::new("CoProm",
                 normalizedData = normalizedData,
                 locationData = locationData,
                 metaData = metaData,
                 cellTypes = cellTypes,
                 slideID = slideID,
                 slideList = unique_slides,
                 geneList = geneList
    )
  }
)

#--- Helper Functions ---

#' Normalize a vector to have unit norm
#' @param v numeric vector
#' @return normalized vector
#' @noRd
normalize_vec <- function(v) {
  norm_v <- sqrt(sum(v^2))
  if (norm_v < .Machine$double.eps) {
    return(rep(0, length(v)))
  } else {
    return(v / norm_v)
  }
}


#--- Pipeline Functions ---

#' Subset Data for Multi-Slide Analysis
#'
#' Subsets the combined data in a `CoProm` object based on cell types of interest.
#'
#' @param object A `CoProm` object.
#' @param cellTypesOfInterest A character vector of cell types to keep.
#'
#' @return A `CoProm` object with subset slots populated.
#' @export
#' @rdname subsetDataMulti
#' @aliases subsetDataMulti,CoProm-method
setGeneric("subsetDataMulti",
           function(object,
                    cellTypesOfInterest) standardGeneric("subsetDataMulti"))

#' @rdname subsetDataMulti
#' @aliases subsetDataMulti,CoProm-method
#' @export
setMethod("subsetDataMulti", "CoProm", function(object, cellTypesOfInterest) {
  if (length(cellTypesOfInterest) < 1) { # Allow 1 for integration/PCA
    stop("Please specify at least one cell type of interest.")
  }
  if (!all(cellTypesOfInterest %in% unique(object@cellTypes))) {
    stop("Some cellTypesOfInterest are not present in the object's cellTypes.")
  }

  subsetIndices <- object@cellTypes %in% cellTypesOfInterest

  if (sum(subsetIndices) < 10 * length(object@slideList)) {
    warning("Fewer than 10 cells per slide on average exist in the subset. Check cellTypesOfInterest.")
  }
  if (sum(subsetIndices) == 0) {
    stop("No cells selected with the given cellTypesOfInterest.")
  }


  object@cellTypesOfInterest <- cellTypesOfInterest
  object@normalizedDataSub <- object@normalizedData[subsetIndices, , drop = FALSE]
  object@locationDataSub <- object@locationData[subsetIndices, , drop = FALSE]
  object@metaDataSub <- object@metaData[subsetIndices, , drop = FALSE]
  object@cellTypesSub <- object@cellTypes[subsetIndices]
  # Ensure slideID is retained in metaDataSub if not already there
  if(!"slideID" %in% colnames(object@metaDataSub)) {
    object@metaDataSub$slideID <- object@slideID[subsetIndices]
  }


  return(object)
})

#' Compute PCA on Integrated Multi-Slide Data
#'
#' Performs PCA on the integrated data stored within the `CoProm` object.
#' Assumes integration has created a common space across slides.
#'
#' @importFrom stats setNames prcomp
#' @importFrom irlba prcomp_irlba
#' @param object A `CoProm` object with the `integratedData` slot populated.
#' @param nPCA Number of principal components to compute for each cell type.
#' @param use_irlba Logical, whether to use irlba for faster computation.
#' @param dataUse What data to use, choices between "raw" and "integrated".
#'  Default is "raw"
#' @param center_per_slide After the global PCA, do we do center per slide
#'  again? By default this is set to FALSE
#' @param center Whether to center the matrix before PCA
#' @param scale. Whether to scale the matrix before PCA
#'
#' @return A `CoProm` object with the `pcaResults` slot populated.
#'         `pcaResults` structure: `list(slideID = list(cellType = pc_matrix))`.
#' @export
#' @rdname computePCAMulti
#' @aliases computePCAMulti,CoProm-method
setGeneric("computePCAMulti",
           function(object, nPCA = 40,
                    center = TRUE, scale. = TRUE,
                    dataUse = "raw",
                    center_per_slide = FALSE,
                    use_irlba = TRUE) standardGeneric("computePCAMulti"))

#' @rdname computePCAMulti
#' @aliases computePCAMulti,CoProm-method
#' @importFrom stats setNames prcomp predict
#' @importFrom irlba prcomp_irlba
#' @export
setMethod("computePCAMulti",
          "CoProm",
          function(object, nPCA = 40, center = TRUE, scale. = TRUE,
                   dataUse = "raw",
                   center_per_slide = FALSE, use_irlba = TRUE) {
  ## arg match
  if(!(dataUse %in% c("raw", "integrated"))) stop("dataUse must be raw or integrated")
  ## Check if integrated data exists
  if (dataUse == "integrated" && length(object@integratedData) == 0) {
    stop("integratedData slot is empty. Run integrateSlidesMulti first.")
  }
  cts <- object@cellTypesOfInterest
  if(length(cts) == 0) {
    stop("cellTypesOfInterest not set. Run subsetDataMulti first.")
  }
  slides <- object@slideList

  # Initialize pcaResults structure
  pca_results_all <- setNames(vector("list", length = length(slides)), slides)
  for(sID in slides) {
    pca_results_all[[sID]] <- setNames(vector("list", length = length(cts)), cts)
  }

  pca_global <- setNames(vector("list", length = length(cts)), cts)

  # Perform PCA per cell type on the integrated data
  for (ct in cts) {
    message(paste("Performing PCA for cell type:", ct))

    if (dataUse == "integrated" && !ct %in% names(object@integratedData)) {
      warning(paste("No integrated data found for cell type:", ct, "- Skipping PCA."))
      next
    }

    if(dataUse == "integrated"){
      integrated_mat_ct <- object@integratedData[[ct]]
      # Needs verification based on integrateSlidesMulti output
    }else{ ## raw
      integrated_mat_ct <- object@normalizedDataSub[object@cellTypesSub == ct,]
    }

    # Ensure it's a matrix
    if (!is.matrix(integrated_mat_ct) && !inherits(integrated_mat_ct, "Matrix")) {
      message("converting data matrix into dense matrix")
      integrated_mat_ct <- as.matrix(integrated_mat_ct)
    }

    # Check dimensions
    if(nrow(integrated_mat_ct) != sum(object@cellTypesSub == ct)) {
      stop(paste("Integrated data dimensions mismatch for cell type:", ct))
    }

    if (center & scale.) {
      scaledData <- center_scale_matrix_opt(integrated_mat_ct)
      message("data centered and scaled")
    } else if (center) {
      scaledData <- t(t(integrated_mat_ct) - colMeans(integrated_mat_ct))
    } else {
      warning(paste("It is not recommended to skip both centering,",
                    "and scaling of the data, unless the data has been centered and",
                    "scaled when creating the CoPro object.",
                    sep = " "
      ))
      scaledData <- integrated_mat_ct
    }

    # Perform PCA on the combined integrated data for this cell type
    pca_ct <- prcomp_irlba(scaledData, n = nPCA,
                             center = FALSE, scale. = FALSE)
    message("PCA computed")
    pca_global[[ct]] <- pca_ct

    # Project each slide's data onto the shared PCs
    for (sID in slides) {
      slide_ID_ct <- object@metaDataSub$slideID[object@cellTypesSub == ct]
      row_names_ct <- rownames(object@metaDataSub)[object@cellTypesSub == ct]
      pca_sub <- pca_ct$x[slide_ID_ct == sID, , drop = FALSE]
      rownames(pca_sub) <- row_names_ct[slide_ID_ct == sID]

      if(center_per_slide){
        cat("Centering per slide","/n")
        pca_sub <- scale(pca_sub,center = TRUE, scale = FALSE)
      }

      pca_results_all[[sID]][[ct]] <- pca_sub

    }
  } # End loop over cell types

  object@pcaGlobal <- pca_global
  object@pcaResults <- pca_results_all
  object@nPCA <- nPCA
  # Default scalePCs set to True
  object@scalePCs <- TRUE

  return(object)
})


#' Compute Distances within each Slide
#'
#' Calculates pairwise cell distances for specified cell types within each slide.
#'
#' @importFrom fields rdist
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @param object A `CoProm` object.
#' @param distType Type of distance ("Euclidean2D", "Euclidean3D").
#' @param xDistScale Scaling factors for coordinates.
#' @param yDistScale Scaling factors for coordinates.
#' @param zDistScale Scaling factors for coordinates.
#' @param normalizeDistance Normalize distances based on 0.01 percentile?
#' @param truncateLowDist Truncate very small distances?
#' @param verbose Print progress messages?
#'
#' @return `CoProm` object with the `distances` slot populated:
#'         `list(slideID = list(cellType1 = list(cellType2 = dist_matrix)))`.
#' @export
#' @rdname computeDistanceMulti
#' @aliases computeDistanceMulti,CoProm-method
setGeneric("computeDistanceMulti", function(
    object, distType = c("Euclidean2D", "Euclidean3D"),
    xDistScale = 1, yDistScale = 1, zDistScale = 1,
    normalizeDistance = TRUE, truncateLowDist = TRUE,
    verbose = TRUE) standardGeneric("computeDistanceMulti"))

#' @rdname computeDistanceMulti
#' @aliases computeDistanceMulti,CoProm-method
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @importFrom fields rdist
#' @export
setMethod("computeDistanceMulti", "CoProm", function(
    object, distType = c("Euclidean2D", "Euclidean3D"),
    xDistScale = 1, yDistScale = 1, zDistScale = 1,
    normalizeDistance = TRUE, truncateLowDist = TRUE,
    verbose = TRUE) {

  distType <- match.arg(distType)
  cts <- object@cellTypesOfInterest
  if (length(cts) < 2) {
    stop("At least two cellTypesOfInterest are needed for distance calculations.")
  }
  slides <- object@slideList

  distances_all <- setNames(vector("list", length = length(slides)), slides)
  global_min_1percentile <- Inf # To normalize across slides

  for (sID in slides) {
    if (verbose) message(paste("Computing distances for slide:", sID))

    # Initialize structure for this slide
    distances_slide <- setNames(vector("list", length = length(cts)), cts)
    for (i in cts) {
      distances_slide[[i]] <- setNames(vector("list", length = length(cts)), cts)
    }

    # Get subset indices for the current slide
    slide_indices <- which(object@metaDataSub$slideID == sID)
    locationData_slide <- object@locationDataSub[slide_indices, , drop = FALSE]
    cellTypes_slide <- object@cellTypesSub[slide_indices]
    cellIDs_slide <- rownames(locationData_slide)

    pair_cell_types <- combn(cts, 2)
    dist_1percentile_slide <- vector(mode = "numeric", length = ncol(pair_cell_types))
    all_zero_pairs <- c() # Track pairs with all zero distances

    for (pp in seq_len(ncol(pair_cell_types))) {
      ct_i <- pair_cell_types[1, pp]
      ct_j <- pair_cell_types[2, pp]

      # Indices for cell types within the current slide
      idx_i <- which(cellTypes_slide == ct_i)
      idx_j <- which(cellTypes_slide == ct_j)

      if (length(idx_i) <= 5 || length(idx_j) <= 5) {
        if (verbose) message(paste("Skipping pair", ct_i, "-", ct_j,
                                   "in slide", sID, "(missing cells)"))
        distances_slide[[ct_i]][[ct_j]] <- matrix(nrow=length(idx_i), ncol=length(idx_j)) # Empty/NA matrix
        dist_1percentile_slide[pp] <- NA
        next
      }

      loc_i <- locationData_slide[idx_i, , drop = FALSE]
      loc_j <- locationData_slide[idx_j, , drop = FALSE]

      if (distType == "Euclidean2D") {
        mat1 <- cbind(loc_i$x * xDistScale, loc_i$y * yDistScale)
        mat2 <- cbind(loc_j$x * xDistScale, loc_j$y * yDistScale)
      } else { # Euclidean3D
        if (!"z" %in% colnames(loc_i)) stop("z coordinate missing for 3D distance.")
        mat1 <- cbind(loc_i$x * xDistScale, loc_i$y * yDistScale, loc_i$z * zDistScale)
        mat2 <- cbind(loc_j$x * xDistScale, loc_j$y * yDistScale, loc_j$z * zDistScale)
      }

      distances_ij <- fields::rdist(mat1, mat2)
      rownames(distances_ij) <- rownames(loc_i)
      colnames(distances_ij) <- rownames(loc_j)

      min_dist_nonzero <- min(distances_ij[distances_ij > 0], Inf)

      if (any(distances_ij == 0)) {
        warning(paste("Zero distances detected between", ct_i, "and", ct_j, "in slide", sID,
                      ". Replacing with smallest non-zero distance:", min_dist_nonzero))
        distances_ij[distances_ij == 0] <- min_dist_nonzero
      }

      if(is.infinite(min_dist_nonzero)){
        warning(paste("All distances are zero or missing between",
                      ct_i, "and", ct_j, "in slide", sID))
        dist_1percentile_slide[pp] <- NA
        all_zero_pairs <- c(all_zero_pairs, pp)
      } else {
        percentile_choice <- min(1e-3, 2/(max(nrow(distances_ij), ncol(distances_ij))))
        dist_1percentile_slide[pp] <- quantile(distances_ij[distances_ij > 0], percentile_choice)
        global_min_1percentile <- min(global_min_1percentile, dist_1percentile_slide[pp], na.rm = TRUE)

        if (truncateLowDist) {
          distances_ij[distances_ij < dist_1percentile_slide[pp]] <- dist_1percentile_slide[pp]
        }
      }
      distances_slide[[ct_i]][[ct_j]] <- distances_ij
      if (verbose && !is.infinite(min_dist_nonzero)) {
        cat("Slide:", sID, ", Pair:", ct_i, "-", ct_j, "\n")
        print(quantile(distances_ij))
      }
    } # End pair loop (pp)
    distances_all[[sID]] <- distances_slide
  } # End slide loop (sID)

  # Normalize across slides if requested
  if (normalizeDistance) {
    if(is.infinite(global_min_1percentile)){
      warning("Cannot normalize distances - no valid non-zero distances found across slides.")
    } else {
      scaling_factor <- 0.01 / global_min_1percentile
      if (verbose) cat("Global distance scaling factor:", scaling_factor, "\n")

      for (sID in slides) {
        for (pp in seq_len(ncol(pair_cell_types))) {
          ct_i <- pair_cell_types[1, pp]
          ct_j <- pair_cell_types[2, pp]
          # Avoid scaling if matrix is NA/empty or had all zeros
          if(!is.null(distances_all[[sID]][[ct_i]][[ct_j]]) && !(pp %in% all_zero_pairs)){
            distances_all[[sID]][[ct_i]][[ct_j]] <- distances_all[[sID]][[ct_i]][[ct_j]] * scaling_factor
          }
        }
      }
    }
  }

  object@distances <- distances_all
  return(object)
})


#' Compute Kernel Matrices for Multi-Slide Data
#'
#' Calculates kernel matrices for each slide based on distances and sigma values.
#'
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @param object A `CoProm` object with the `distances` slot populated.
#' @param sigmaValues A numeric vector of sigma values.
#' @param lowerLimit Lower threshold for kernel values.
#' @param upperQuantile Quantile for clipping upper kernel values.
#' @param normalizeKernel Normalize kernel matrix values? (Not recommended as it breaks CCA properties if done row/col-wise).
#' @param minAveCellNeighor Minimum average neighbors check (applied per slide).
#' @param verbose Print progress?
#'
#' @return `CoProm` object with `kernelMatrices` slot populated:
#'         `list(sigma = list(slideID = list(cellType1 = list(cellType2 = kernel_matrix))))`.
#' @export
#' @rdname computeKernelMatrixMulti
#' @aliases computeKernelMatrixMulti,CoProm-method
setGeneric("computeKernelMatrixMulti", function(
    object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
    normalizeKernel = FALSE, minAveCellNeighor = 2,
    verbose = TRUE) standardGeneric("computeKernelMatrixMulti"))

#' @rdname computeKernelMatrixMulti
#' @aliases computeKernelMatrixMulti,CoProm-method
#' @importFrom utils combn
#' @importFrom stats setNames quantile
#' @export
setMethod("computeKernelMatrixMulti", "CoProm", function(
    object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
    normalizeKernel = FALSE, minAveCellNeighor = 2,
    verbose = TRUE) {

  if (length(object@distances) == 0) stop("Distances not computed. Run computeDistanceMulti first.")
  if (normalizeKernel) warning("normalizeKernel=TRUE is generally not recommended for CCA-based methods.")

  cts <- object@cellTypesOfInterest
  if (length(cts) < 2) stop("At least two cellTypesOfInterest are needed.")
  slides <- object@slideList
  if (length(sigmaValues) == 0) stop("sigmaValues must be provided.")
  object@sigmaValues <- sigmaValues # Store provided sigma values

  kernel_matrices_all <- setNames(vector("list", length(sigmaValues)),
                                  paste0("sigma_", sigmaValues))

  for (tt in seq_along(sigmaValues)) {
    sigma_val <- sigmaValues[tt]
    sigma_name <- names(kernel_matrices_all)[tt]
    if (verbose) message(paste("Computing kernels for sigma =", sigma_val))

    kernels_for_sigma <- setNames(vector("list", length = length(slides)), slides)
    sigma_valid_across_slides <- TRUE # Track if this sigma works for all slides

    for (sID in slides) {
      # Initialize structure for this slide
      kernels_slide <- setNames(vector("list", length = length(cts)), cts)
      for (i in cts) {
        kernels_slide[[i]] <- setNames(vector("list", length = length(cts)), cts)
      }

      dist_slide <- object@distances[[sID]]
      if(length(dist_slide) == 0) {
        warning(paste("No distance matrix found for slide", sID))
        sigma_valid_across_slides <- FALSE
        next # Skip slide if no distances
      }

      pair_cell_types <- combn(cts, 2)
      sigma_valid_for_slide <- TRUE

      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]

        dist_mat_ij <- dist_slide[[ct_i]][[ct_j]]

        if (length(dist_mat_ij) == 0 || all(is.na(dist_mat_ij)) ||
            nrow(dist_mat_ij) == 0 || ncol(dist_mat_ij) == 0 ) {
          if(verbose) message(paste("  Skipping kernel for",
                                    ct_i, "-", ct_j, "in slide", sID,
                                    "(missing/empty distance matrix)"))
          kernels_slide[[ct_i]][[ct_j]] <- matrix(nrow=nrow(dist_mat_ij), ncol=ncol(dist_mat_ij)) # Store empty/NA
          next
        }

        kernel_current <- kernel_from_distance( # Assumes kernel_from_distance is available
          sigma = sigma_val,
          dist_mat = dist_mat_ij,
          lower_limit = lowerLimit
        )

        # Check sparsity (using helper .CheckSigmaValuesToRemove - needs access or re-implementation)
        # Simplified check: if matrix is nearly all zero, flag sigma for this slide
        n_cells1 <- nrow(kernel_current)
        n_cells2 <- ncol(kernel_current)
        min_count_nonzero <- minAveCellNeighor * min(n_cells1, n_cells2)
        if (sum(kernel_current > lowerLimit) < min_count_nonzero) {
          warning(paste("Kernel matrix for", ct_i,"-", ct_j, "in slide", sID,
                        "with sigma =", sigma_val,
                        "is too sparse. Flagging sigma."))
          sigma_valid_for_slide <- FALSE
          # Store the sparse kernel anyway? Or NULL? Let's store it for now.
          # kernels_slide[[ct_i]][[ct_j]] <- NULL # Or keep the sparse matrix
        }

        if (anyNA(kernel_current)) {
          warning(paste("NA values in kernel for", ct_i, "-", ct_j, "in slide",
                        sID, "sigma =", sigma_val))
          # Decide how to handle NAs, maybe imputation or flagging
        }

        # Clipping large values
        upper_clip <- quantile(kernel_current[kernel_current >= lowerLimit & !is.na(kernel_current)], upperQuantile, na.rm = TRUE)
        if(!is.na(upper_clip)){
          kernel_current[kernel_current >= upper_clip & !is.na(kernel_current)] <- upper_clip
        }


        # Normalization (apply cautiously)
        if (normalizeKernel) {
          # Example: Scale by median row sum (use with caution)
          rs_kernel <- rowSums(kernel_current, na.rm=TRUE)
          median_rs <- median(rs_kernel[rs_kernel > 1e-5], na.rm=TRUE)
          if(!is.na(median_rs) && median_rs > 1e-6){
            kernel_current <- kernel_current / median_rs
          }
        }

        # Ensure values below lowerLimit are zero
        kernel_current[kernel_current < lowerLimit & !is.na(kernel_current)] <- 0

        kernels_slide[[ct_i]][[ct_j]] <- kernel_current

      } # End pair loop (pp)

      if(!sigma_valid_for_slide){
        sigma_valid_across_slides <- FALSE
        warning(paste("Sigma value", sigma_val, "invalid for slide", sID, "due to sparse kernels."))
      }
      kernels_for_sigma[[sID]] <- kernels_slide

    } # End slide loop (sID)

    if(sigma_valid_across_slides){
      kernel_matrices_all[[sigma_name]] <- kernels_for_sigma
    } else {
      kernel_matrices_all[[sigma_name]] <- NULL # Remove sigma if invalid for any slide
      warning(paste("Removing sigma value", sigma_val, "as it was invalid for one or more slides."))
    }

  } # End sigma loop (tt)

  # Filter out NULL entries (invalid sigmas)
  kernel_matrices_all <- kernel_matrices_all[!sapply(kernel_matrices_all, is.null)]
  object@sigmaValues <- as.numeric(gsub("sigma_", "", names(kernel_matrices_all))) # Update valid sigmas

  object@kernelMatrices <- kernel_matrices_all
  return(object)
})


# X_list_all: list(slideID = list(cellType = pc_matrix))
.scalePCMatsMulti <- function(X_list_all, pca_object, slides, cts){

  # Initialize pcaResults structure
  X_list_scaled <- setNames(vector("list", length = length(slides)), slides)
  for(sID in slides) {
    X_list_scaled[[sID]] <- setNames(vector("list", length = length(cts)), cts)
  }

  for(ct in cts){
    pca_A_sd <- pca_object[[ct]]$sdev
    for(sID in slides){
      X_list_scaled[[sID]][[ct]] <- scale(X_list_all[[sID]][[ct]],
                                          center = FALSE, scale = pca_A_sd)
    }
  }
  return(X_list_scaled)
}

#' Run Multi-Slide Spatially Kernel Restricted CCA (skrCCA)
#'
#' Performs skrCCA on the `CoProm` object using data integrated across slides.
#' Calls the multi-slide optimization functions.
#'
#' @importFrom stats setNames
#' @param object A `CoProm` object with `pcaResults` and `kernelMatrices` populated.
#' @param nCC Number of canonical components to compute.
#' @param sigmaChoice A specific sigma value to use. If NULL (default), uses all valid sigma values from `object@sigmaValues`.
#' @param tol Tolerance for optimization convergence.
#' @param maxIter Maximum iterations for optimization.
#' @param n_cores Number of cores for parallel computation within optimization.
#' @param scalePCs Whether to scale each PC before computing skrCCA
#'
#' @return `CoProm` object with `skrCCAOut` slot populated:
#'         `list(sigma = list(cellType = weight_matrix))` (shared weights).
#' @export
#' @rdname runSkrCCAMulti
#' @aliases runSkrCCAMulti,CoProm-method
setGeneric("runSkrCCAMulti", function(
    object, nCC = 2, sigmaChoice = NULL, tol = 1e-5, scalePCs,
    maxIter = 200, n_cores = 1) standardGeneric("runSkrCCAMulti"))

#' @rdname runSkrCCAMulti
#' @aliases runSkrCCAMulti,CoProm-method
#' @importFrom stats setNames
#' @export
setMethod("runSkrCCAMulti", "CoProm", function(
    object, nCC = 2, sigmaChoice = NULL, tol = 1e-5, scalePCs,
    maxIter = 200, n_cores = 1) {

  # --- Input Checks ---
  if (length(object@pcaResults) == 0) stop("PCA results missing. Run computePCAMulti.")
  if (length(object@kernelMatrices) == 0) stop("Kernel matrices missing. Run computeKernelMatrixMulti.")
  cts <- object@cellTypesOfInterest
  if (length(cts) < 2) stop("At least two cellTypesOfInterest needed for skrCCA.")
  slides <- object@slideList

  # Determine which sigma values to use
  if (is.null(sigmaChoice)) {
    sigmas_to_run <- object@sigmaValues
    if (length(sigmas_to_run) == 0) stop("No valid sigma values found in object.")
  } else {
    if (!paste0("sigma_", sigmaChoice) %in% names(object@kernelMatrices)) {
      stop(paste("Chosen sigma value", sigmaChoice, "not found or was invalid."))
    }
    sigmas_to_run <- sigmaChoice
  }
  sigma_names_run <- paste0("sigma_", sigmas_to_run)

  # --- Prepare Input Lists for Optimization ---
  # X_list_all: list(slideID = list(cellType = pc_matrix))
  X_list_all <- object@pcaResults # Assumes correct structure from computePCAMulti

  if(scalePCs){
    X_list_all <- .scalePCMatsMulti(X_list_all = X_list_all,
                                    pca_object = object@pcaGlobal,
                                    slides = slides,
                                    cts = cts )
  }

  # Check consistency
  if(!all(slides %in% names(X_list_all))) stop("Slide mismatch in pcaResults")
  for(sID in slides){
    if(!all(cts %in% names(X_list_all[[sID]]))) stop(paste("Cell type mismatch in pcaResults for slide", sID))
  }


  # --- Run Optimization per Sigma ---
  cca_out_all_sigmas <- setNames(vector("list", length = length(sigmas_to_run)), sigma_names_run)

  for (idx in seq_along(sigmas_to_run)) {
    sig_val <- sigmas_to_run[idx]
    sig_name <- sigma_names_run[idx]
    message(paste("Running skrCCA for sigma =", sig_val))

    # K_list_all: list(slideID = list(cellType1 = list(cellType2 = kernel_matrix)))
    K_list_all_sigma <- object@kernelMatrices[[sig_name]]
    # Check consistency
    if(!all(slides %in% names(K_list_all_sigma))) stop(paste("Slide mismatch in kernelMatrices for sigma", sig_val))
    for(sID in slides){
      # Need to check pairs, assuming combn structure was used
      pairs_ok <- TRUE
      # Simplified check: check if first pair exists
      if(length(cts)>=2 && is.null(K_list_all_sigma[[sID]][[cts[1]]][[cts[2]]])){
        pairs_ok <- FALSE
      }
      if(!all(cts %in% names(K_list_all_sigma[[sID]])) || !pairs_ok) {
        stop(paste("Cell type/pair mismatch in kernelMatrices for slide", sID, "sigma", sig_val))
      }
    }


    # Run optimization for the first component
    # Ensure optimize_bilinear_multi_slides is loaded/available
    cca_result_1 <- optimize_bilinear_multi_slides(
      X_list_all = X_list_all,
      K_list_all = K_list_all_sigma,
      max_iter = maxIter, tol = tol, n_cores = n_cores
    )
    # Name the result list using cell types
    names(cca_result_1) <- cts # Uses implicit order, assumes optimization returns unnamed list


    if (nCC == 1) {
      cca_out_all_sigmas[[sig_name]] <- cca_result_1
    } else {
      # Run optimization for subsequent components
      # Ensure optimize_bilinear_multi_n_slides is loaded/available
      cca_result_n <- optimize_bilinear_multi_n_slides(
        X_list_all = X_list_all,
        K_list_all = K_list_all_sigma,
        w_list = cca_result_1, # Pass named list from previous step
        cellTypesOfInterest = cts,
        nCC = nCC,
        max_iter = maxIter, tol = tol # n_cores not directly used here but in helper
      )
      # optimize_bilinear_multi_n_slides should return named list
      if(is.null(names(cca_result_n))) names(cca_result_n) <- cts # Ensure names

      cca_out_all_sigmas[[sig_name]] <- cca_result_n
    }
  } # End sigma loop

  object@skrCCAOut <- cca_out_all_sigmas
  object@nCC <- nCC
  return(object)
})


#' Compute Normalized Correlation for Multi-Slide Data
#'
#' Calculates normalized correlation, potentially per slide or aggregated.
#' Needs careful definition of the desired output.
#'
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @importFrom stats aggregate
#' @param object A `CoProm` object with skrCCA results.
#' @param tol Tolerance for SVD approximation of spectral norm.
#' @param calculationMode "aggregate" or "perSlide". Determines how correlation is computed/reported.
#'
#' @return `CoProm` object with `normalizedCorrelation` slot populated.
#' @export
#' @rdname computeNormalizedCorrelationMulti
#' @aliases computeNormalizedCorrelationMulti,CoProm-method
setGeneric("computeNormalizedCorrelationMulti", function(
    object, tol = 1e-4, calculationMode = c("aggregate", "perSlide")
) standardGeneric("computeNormalizedCorrelationMulti"))

#' @rdname computeNormalizedCorrelationMulti
#' @aliases computeNormalizedCorrelationMulti,CoProm-method
#' @importFrom utils combn
#' @importFrom irlba irlba
#' @export
setMethod("computeNormalizedCorrelationMulti", "CoProm", function(
    object, tol = 1e-4, calculationMode = c("aggregate", "perSlide")) {

  calculationMode <- match.arg(calculationMode)

  # --- Input Checks ---
  if (length(object@skrCCAOut) == 0) stop("skrCCA results missing. Run runSkrCCAMulti.")
  if (length(object@pcaResults) == 0) stop("PCA results missing.")
  if (length(object@kernelMatrices) == 0) stop("Kernel matrices missing.")
  cts <- object@cellTypesOfInterest
  if (length(cts) < 2) stop("Need at least two cell types.")
  slides <- object@slideList
  sigmas_run <- names(object@skrCCAOut)
  if (length(sigmas_run) == 0) stop("No skrCCA results found.")
  nCC <- object@nCC
  if(is.null(nCC) || length(nCC)==0) nCC <- ncol(object@skrCCAOut[[1]][[1]]) # Infer nCC

  # --- Precompute Spectral Norms (Per Slide, Per Sigma) ---
  message("Calculating spectral norms (can take time)...")
  norm_K_all <- setNames(vector("list", length = length(sigmas_run)), sigmas_run)
  for (sig_name in sigmas_run) {
    norm_K_sigma <- setNames(vector("list", length = length(slides)), slides)
    K_list_sigma <- object@kernelMatrices[[sig_name]]
    pair_cell_types <- combn(cts, 2)

    for (sID in slides) {
      norm_K_slide <- setNames(vector("list", length = length(cts)), cts)
      for (ct_i in cts) norm_K_slide[[ct_i]] <- setNames(vector("list", length=length(cts)), cts)

      K_list_slide <- K_list_sigma[[sID]]

      for (pp in seq_len(ncol(pair_cell_types))) {
        ct_i <- pair_cell_types[1, pp]
        ct_j <- pair_cell_types[2, pp]
        K <- K_list_slide[[ct_i]][[ct_j]]
        if (!is.null(K) && nrow(K) > 0 && ncol(K) > 0) {
          # Handle potential sparse matrices if irlba supports them well
          svd_d <- tryCatch({
            irlba::irlba(K, nv = 1, nu = 0, tol = tol)$d[1] # Only need singular value
          }, error = function(e) {
            warning(paste("SVD failed for K[",sID,",",ct_i,",",ct_j,"]:", e$message))
            NA # Return NA if SVD fails
          })
          norm_K_slide[[ct_i]][[ct_j]] <- svd_d
          norm_K_slide[[ct_j]][[ct_i]] <- svd_d # Store transpose too
        } else {
          norm_K_slide[[ct_i]][[ct_j]] <- NA
          norm_K_slide[[ct_j]][[ct_i]] <- NA
        }
      }
      norm_K_sigma[[sID]] <- norm_K_slide
    }
    norm_K_all[[sig_name]] <- norm_K_sigma
  }
  message("Finished calculating spectral norms.")

  # --- Calculate Correlation ---
  correlation_results <- setNames(vector("list", length = length(sigmas_run)), sigmas_run)

  for (sig_name in sigmas_run) {
    W_list_sigma <- object@skrCCAOut[[sig_name]] # Shared weights

    if (calculationMode == "perSlide") {
      correlation_per_slide <- setNames(vector("list", length=length(slides)), slides)
      for(sID in slides) {
        df_slide <- data.frame(
          sigmaValue = character(),
          slideID = character(),
          cellType1 = character(), cellType2 = character(),
          CC_index = integer(), normalizedCorrelation = numeric(),
          stringsAsFactors = FALSE
        )
        X_list_slide <- object@pcaResults[[sID]]
        K_list_slide <- object@kernelMatrices[[sig_name]][[sID]]
        norm_K_slide <- norm_K_all[[sig_name]][[sID]]

        if(is.null(X_list_slide) || is.null(K_list_slide)) next # Skip if data missing

        for (pp in seq_len(ncol(pair_cell_types))) {
          ct_i <- pair_cell_types[1, pp]
          ct_j <- pair_cell_types[2, pp]

          X_i <- X_list_slide[[ct_i]]
          X_j <- X_list_slide[[ct_j]]
          K_ij <- K_list_slide[[ct_i]][[ct_j]]
          norm_K_ij <- norm_K_slide[[ct_i]][[ct_j]]

          if(is.null(X_i) || is.null(X_j) || is.null(K_ij) ||
             is.na(norm_K_ij) || norm_K_ij < 1e-9 ||
             nrow(X_i)==0 || nrow(X_j)==0) next # Skip if data missing/invalid

          for (cc in 1:nCC) {
            w_i <- W_list_sigma[[ct_i]][, cc, drop = FALSE]
            w_j <- W_list_sigma[[ct_j]][, cc, drop = FALSE]

            Xiw <- X_i %*% w_i
            Xjw <- X_j %*% w_j

            numerator <- (t(Xiw) %*% K_ij %*% Xjw)[1, 1]
            denom_norm <- sqrt(sum(Xiw^2)) * sqrt(sum(Xjw^2)) * norm_K_ij

            norm_corr_val <- ifelse(abs(denom_norm) < 1e-9, 0, numerator / denom_norm)

            df_slide <- rbind(df_slide, data.frame(
              sigmaValue = as.numeric(gsub("sigma_", "", sig_name)),
              slideID = sID,
              cellType1 = ct_i, cellType2 = ct_j,
              CC_index = cc, normalizedCorrelation = norm_corr_val,
              stringsAsFactors = FALSE
            ))
          } # end CC loop
        } # end pair loop
        correlation_per_slide[[sID]] <- df_slide
      } # end slide loop
      correlation_results[[sig_name]] <- do.call(rbind, correlation_per_slide)

    } else { # calculationMode == "aggregate"
      df_agg <- data.frame(
        sigmaValue = numeric(),
        cellType1 = character(), cellType2 = character(),
        CC_index = integer(), aggregateCorrelation = numeric(),
        stringsAsFactors = FALSE
      )
      # Calculate sum(w_i^T X_iq^T K_ijq X_jq w_j) / ( sum(sqrt(sum(Xiw_q^2))) * sum(sqrt(sum(Xjw_q^2))) * mean(normK_ijq) ) ?
      # Or sum(numerator_q) / sum(denominator_q) ?
      # Simpler: calculate objective value achieved by weights and normalize somehow?
      # Let's calculate sum(num_q) / ( sqrt(sum(||Xiw_q||^2)) * sqrt(sum(||Xjw_q||^2)) * mean(normK_ijq) )
      # This needs careful definition based on desired interpretation.

      message("Aggregate correlation calculation needs precise definition. Placeholder implementation.")
      # Placeholder: Calculate mean per-slide correlation
      temp_res <- computeNormalizedCorrelationMulti(object, tol=tol, calculationMode="perSlide")
      df_per_slide <- temp_res@normalizedCorrelation[[sig_name]]
      if(!is.null(df_per_slide) && nrow(df_per_slide) > 0){
        df_agg <- aggregate(normalizedCorrelation ~ sigmaValue + cellType1 + cellType2 + CC_index,
                            data=df_per_slide, FUN = mean)
        colnames(df_agg)[colnames(df_agg) == "normalizedCorrelation"] <- "aggregateCorrelation"
      } else {
        df_agg <- data.frame( # Empty df if no per-slide results
          sigmaValue = numeric(), cellType1 = character(), cellType2 = character(),
          CC_index = integer(), aggregateCorrelation = numeric(), stringsAsFactors = FALSE
        )
      }

      correlation_results[[sig_name]] <- df_agg
    }
  } # End sigma loop

  object@normalizedCorrelation <- correlation_results

  # --- Choose Optimal Sigma (based on aggregate or mean per-slide CC1 correlation) ---
  if(length(correlation_results) > 0) {
    first_sigma_res <- correlation_results[[1]]
    if(!is.null(first_sigma_res) && nrow(first_sigma_res) > 0) {
      corr_col_name <- ifelse("aggregateCorrelation" %in% names(first_sigma_res), "aggregateCorrelation", "normalizedCorrelation")
      all_corrs <- do.call(rbind, correlation_results)
      cc1_corrs <- all_corrs[all_corrs$CC_index == 1, ]

      if(nrow(cc1_corrs) > 0) {
        # Average correlation across pairs and potentially slides for CC1 for each sigma
        mean_corr_per_sigma <- tapply(cc1_corrs[[corr_col_name]], cc1_corrs$sigmaValue, mean, na.rm=TRUE)
        if(any(!is.na(mean_corr_per_sigma))){
          object@sigmaValueChoice <- as.numeric(names(which.max(mean_corr_per_sigma)))
        } else {
          warning("Could not determine optimal sigma: All correlations were NA.")
          object@sigmaValueChoice <- numeric() # Empty
        }

      } else {
        warning("Could not determine optimal sigma: No CC1 correlation values found.")
        object@sigmaValueChoice <- numeric() # Empty
      }
    } else {
      warning("Could not determine optimal sigma: Correlation results are empty.")
      object@sigmaValueChoice <- numeric() # Empty
    }
  } else {
    object@sigmaValueChoice <- numeric() # Empty if no results
  }

  return(object)
})


#' Compute Gene and Cell Scores for Multi-Slide Data
#'
#' Calculates slide-specific cell scores and potentially shared gene scores.
#'
#' @importFrom stats setNames
#' @param object A `CoProm` object with skrCCA results.
#' @param sigmaChoice A specific sigma value to use. If NULL (default), uses `object@sigmaValueChoice`. If that's also NULL, uses the first available sigma.
#'
#' @return `CoProm` object with `cellScores` and `geneScores` slots populated.
#' @export
#' @rdname computeGeneAndCellScoresMulti
#' @aliases computeGeneAndCellScoresMulti,CoProm-method
setGeneric("computeGeneAndCellScoresMulti", function(
    object, sigmaChoice = NULL) standardGeneric("computeGeneAndCellScoresMulti"))

#' @rdname computeGeneAndCellScoresMulti
#' @aliases computeGeneAndCellScoresMulti,CoProm-method
#' @importFrom stats setNames
#' @export
setMethod("computeGeneAndCellScoresMulti", "CoProm", function(
    object, sigmaChoice = NULL) {

  # --- Input Checks ---
  if (length(object@skrCCAOut) == 0) stop("skrCCA results missing.")
  if (length(object@pcaResults) == 0) stop("PCA results missing.")
  cts <- object@cellTypesOfInterest
  slides <- object@slideList
  nCC <- object@nCC
  if(is.null(nCC) || length(nCC)==0) nCC <- ncol(object@skrCCAOut[[1]][[1]]) # Infer nCC

  # Determine sigma to use
  if (is.null(sigmaChoice)) {
    sigmaChoice <- object@sigmaValueChoice
    if (is.null(sigmaChoice) || length(sigmaChoice) == 0) {
      if(length(object@skrCCAOut) > 0) {
        sigmaChoice <- as.numeric(gsub("sigma_","", names(object@skrCCAOut)[1]))
        warning("sigmaValueChoice not set or found. Using first available sigma: ", sigmaChoice)
      } else {
        stop("Cannot compute scores: No skrCCA results available to choose sigma.")
      }
    }
  }
  sig_name <- paste0("sigma_", sigmaChoice)
  if (!sig_name %in% names(object@skrCCAOut)) {
    stop(paste("Chosen sigma", sigmaChoice, "not found in skrCCA results."))
  }

  W_list <- object@skrCCAOut[[sig_name]] # Shared weights

  # --- Calculate Cell Scores (Slide-Specific) ---
  cell_scores_all <- setNames(vector("list", length = length(slides)), slides)
  for(sID in slides) {
    cell_scores_slide <- setNames(vector("list", length = length(cts)), cts)
    X_list_slide <- object@pcaResults[[sID]]
    meta_slide <- object@metaDataSub[object@metaDataSub$slideID == sID, ]
    celltype_slide <- object@cellTypesSub[object@metaDataSub$slideID == sID]

    for (ct in cts) {
      X_ct <- X_list_slide[[ct]]
      W_ct <- W_list[[ct]] # Shared weight matrix for cell type ct

      if(!is.null(X_ct) && !is.null(W_ct) && nrow(X_ct)>0 && ncol(X_ct) == nrow(W_ct)){
        scores_mat <- X_ct %*% W_ct
        colnames(scores_mat) <- paste0("CC_", 1:nCC)
        rownames(scores_mat) <- rownames(meta_slide)[celltype_slide == ct]

        cell_scores_slide[[ct]] <- scores_mat
      } else {
        # Create empty matrix if no cells or data mismatch
        cell_ids_ct_slide <- rownames(meta_slide)[celltype_slide == ct]
        cell_scores_slide[[ct]] <- matrix(NA, nrow=length(cell_ids_ct_slide), ncol=nCC,
                                          dimnames=list(cell_ids_ct_slide, paste0("CC_", 1:nCC)))
        if(is.null(X_ct) || is.null(W_ct)) warning("Missing PCA or Weight data for cell score calc: ", sID, ", ", ct)
        else if(nrow(X_ct) == 0) warning("No cells for cell score calc: ", sID, ", ", ct)
        else warning("Dimension mismatch for cell score calc: ", sID, ", ", ct)
      }
    }
    cell_scores_all[[sID]] <- cell_scores_slide
  }
  # Store under the chosen sigma
  object@cellScores[[sig_name]] <- cell_scores_all


  # --- Calculate Gene Scores (Shared - requires PCA loadings) ---
  # This requires access to the PCA rotation matrices, which were not explicitly
  # stored in the proposed pcaResults structure. We need to adapt computePCAMulti
  # to store these, or recompute/access them from the integrated object if possible.
  # Assuming PCA objects per cell type (run on integrated data) are accessible:
  # For simplicity, let's assume a function getPCALoadings(object, cellType) exists.

  gene_scores_sigma <- setNames(vector("list", length=length(cts)), cts)
  scalePCs <- object@scalePCs # Should be FALSE if PCA done on integrated/scaled data

       for(ct in cts){
         pca_obj_ct <- object@pcaGlobal[[ct]]
         W_ct <- W_list[[ct]]
           gene_score_mat <- matrix(NA, nrow=length(object@geneList), ncol=nCC,
                                      dimnames=list(object@geneList, paste0("CC_", 1:nCC)))

           if(!is.null(pca_obj_ct) && !is.null(W_ct)){
              rotation <- pca_obj_ct$rotation # Assuming standard prcomp structure
              sdev <- pca_obj_ct$sdev

              for(cc in 1:nCC){
                   w_cc <- W_ct[, cc, drop=FALSE]
                   # Project weights back to gene space
                   # Note: scaling depends on whether input to PCA was scaled in runSkrCCA AND scalePCs setting
                   # If PCA input was scaled (e.g. scale.=TRUE in prcomp) OR scalePCs=TRUE, need sdev.
                   if (scalePCs) { # This logic might need refinement based on exact PCA steps
                       gene_scores_cc <- rotation %*% (w_cc * sdev[1:nrow(w_cc)]) # Check dimensions match
                   } else {
                       gene_scores_cc <- rotation %*% w_cc
                   }
                    gene_score_mat[, paste0("CC_", cc)] <- gene_scores_cc[,1]
               }
           }
            gene_scores_sigma[[ct]] <- gene_score_mat
       }

  # Store under the chosen sigma
  object@geneScores[[sig_name]] <- gene_scores_sigma

  # --- Add Cell Scores to metaDataSub ---
  # Create columns for the chosen sigma
  if(!is.null(object@cellScores[[sig_name]])){
    for(sID in slides){
      scores_slide <- object@cellScores[[sig_name]][[sID]]
      # meta_indices_slide <- which(object@metaDataSub$slideID == sID)

      for(ct in cts){
        scores_ct <- scores_slide[[ct]]
        if(!is.null(scores_ct) && nrow(scores_ct)>0){
          meta_indices_ct_slide <- object@metaDataSub$slideID == sID & object@cellTypesSub == ct # Use cellType from metaDataSub
          cells_ct_slide <- rownames(object@metaDataSub)[meta_indices_ct_slide]

          # Ensure scores match metadata rows
          if(all(rownames(scores_ct) %in% cells_ct_slide) && length(rownames(scores_ct)) == length(cells_ct_slide)){
            scores_ct_aligned <- scores_ct[cells_ct_slide, , drop=FALSE] # Align order
            for(cc in 1:nCC){
              cc_colname <- paste0("CC", cc, "_Score_", sig_name)
              object@metaDataSub[cells_ct_slide, cc_colname] <- scores_ct_aligned[, paste0("CC_", cc)]
            }
          } else {
            warning(paste("Mismatch between cell score rownames and metadata for:", sID, ct))
          }
        }
      }
    }
  }
  return(object)
})
