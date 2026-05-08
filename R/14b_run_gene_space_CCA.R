#' Prepare standardized gene expression matrices per slide per cell type
#'
#' Extracts expression from the CoPro object, filters genes by prevalence,
#' clips extreme values, and standardizes (center + scale) per slide per cell type.
#'
#' @param object A CoProMulti object
#' @param clip Clipping method: \code{"quantile"} (98th pctile) or a numeric threshold
#' @param min_prevalence Minimum fraction of cells expressing a gene
#' @param min_cells Minimum absolute number of cells expressing a gene
#' @param cts Cell types of interest
#' @param slides Slide IDs
#' @return List with \code{Z_by_slide}, \code{genes}, and \code{slides}
#' @importFrom stats quantile sd
#' @noRd
.prepareGeneSpaceData <- function(object, clip = "quantile",
                                  min_prevalence = 0.008, min_cells = 20,
                                  cts, slides) {
  expr <- object@normalizedDataSub
  cell_types_vec <- object@cellTypesSub
  slide_ids <- getSlideID(object)
  cell_names <- rownames(object@metaDataSub)

  # Gene filtering: prevalence across all cells of interest
  gene_expressed <- colMeans(expr > 0)
  gene_n_expressed <- colSums(expr > 0)
  keep_genes <- (gene_expressed >= min_prevalence) & (gene_n_expressed >= min_cells)
  genes <- colnames(expr)[keep_genes]

  if (length(genes) == 0) {
    stop("No genes pass the prevalence filter. Lower min_prevalence or min_cells.")
  }

  expr <- expr[, genes, drop = FALSE]

  # Clip expression values
  if (is.character(clip) && clip == "quantile") {
    q98 <- quantile(as.numeric(expr[expr > 0]), 0.98)
    expr[expr > q98] <- q98
  } else if (is.numeric(clip)) {
    expr[expr > clip] <- clip
  }

  # Build Z_by_slide: standardize per slide per cell type
  Z_by_slide <- setNames(vector("list", length(slides)), slides)

  for (s in slides) {
    Z_by_slide[[s]] <- setNames(vector("list", length(cts)), cts)
    s_idx <- slide_ids == s

    for (ct in cts) {
      ct_idx <- s_idx & (cell_types_vec == ct)
      if (sum(ct_idx) < 3) next

      Z <- as.matrix(expr[ct_idx, genes, drop = FALSE])
      Z <- scale(Z, center = TRUE, scale = TRUE)
      # Replace NaN columns (zero-variance genes on this slide) with 0
      Z[is.nan(Z)] <- 0
      rownames(Z) <- cell_names[ct_idx]
      Z_by_slide[[s]][[ct]] <- Z
    }
  }

  # Filter slides that have all cell types present
  valid_slides <- slides[vapply(slides, function(s) {
    all(vapply(cts, function(ct) !is.null(Z_by_slide[[s]][[ct]]), logical(1)))
  }, logical(1))]

  if (length(valid_slides) == 0) {
    stop("No slides have all cell types present after filtering.")
  }

  Z_by_slide <- Z_by_slide[valid_slides]

  list(Z_by_slide = Z_by_slide, genes = genes, slides = valid_slides)
}

#' Precompute per-slide G x G covariance matrices
#'
#' For each slide, computes self-covariance (C_self) and kernel-smoothed
#' cross-covariance (C_cross) matrices in gene space.
#'
#' @param Z_by_slide Standardized expression matrices per slide per cell type
#' @param flat_kernels Flat list of kernel matrices from CoPro object
#' @param sigma Sigma value for kernel lookup
#' @param slides Slide IDs
#' @param cell_types Cell type names
#' @return List with \code{C_self} and \code{C_cross}
#' @noRd
.precomputeCovarianceMatrices <- function(Z_by_slide, flat_kernels, sigma,
                                         slides, cell_types) {
  C_self <- setNames(vector("list", length(slides)), slides)
  C_cross <- setNames(vector("list", length(slides)), slides)

  pairs <- combn(cell_types, 2, simplify = FALSE)

  for (s in slides) {
    C_self[[s]] <- setNames(vector("list", length(cell_types)), cell_types)
    C_cross[[s]] <- list()

    for (ct in cell_types) {
      Z <- Z_by_slide[[s]][[ct]]
      C_self[[s]][[ct]] <- crossprod(Z) / nrow(Z)
    }

    for (pair in pairs) {
      ct_i <- pair[1]
      ct_j <- pair[2]
      key <- paste0(ct_i, "-", ct_j)

      Z_i <- Z_by_slide[[s]][[ct_i]]
      Z_j <- Z_by_slide[[s]][[ct_j]]
      K_ij <- get_kernel_matrix_flat(flat_kernels, sigma, ct_i, ct_j, slide = s)

      # Align kernel rows/cols to Z cell order
      cells_i <- rownames(Z_i)
      cells_j <- rownames(Z_j)
      K_sub <- K_ij[cells_i, cells_j, drop = FALSE]

      n_i <- nrow(Z_i)
      n_j <- nrow(Z_j)
      C_cross[[s]][[key]] <- crossprod(Z_i, K_sub %*% Z_j) / sqrt(n_i * n_j)
    }
  }

  list(C_self = C_self, C_cross = C_cross)
}

#' Store gene-space CCA results into CoPro object slots
#'
#' Computes cell scores from gene weights and stores gene weights, cell scores,
#' and raw weight vectors into the standard CoPro object slots.
#'
#' @param object CoPro object
#' @param w_list Weight matrices per cell type (G x nCC)
#' @param Z_by_slide Standardized expression per slide per cell type
#' @param sigma Sigma value used
#' @param cts Cell types
#' @param nCC Number of components
#' @param genes Gene names (matching columns of Z)
#' @param slides Slide IDs
#' @return Updated CoPro object
#' @noRd
.storeGeneSpaceCCAResults <- function(object, w_list, Z_by_slide, sigma,
                                      cts, nCC, genes, slides) {
  sigma_name <- paste("sigma", sigma, sep = "_")
  cell_types_vec <- object@cellTypesSub
  slide_ids <- getSlideID(object)
  cell_names <- rownames(object@metaDataSub)

  # Store raw weight vectors in skrCCAOut
  cca_out <- list()
  cca_out[[sigma_name]] <- w_list
  object@skrCCAOut <- cca_out
  object@nCC <- nCC

  # Compute gene scores: for gene-space CCA, weights ARE gene scores directly
  geneScores <- list()
  for (ct in cts) {
    gene_flat_name <- .createGeneScoresName(sigma, ct, slide = NULL)
    gs_mat <- .createScoreMatrix(length(genes), nCC, row_names = genes)
    for (cc in seq_len(nCC)) {
      gs_mat[, paste0("CC_", cc)] <- w_list[[ct]][, cc]
    }
    geneScores[[gene_flat_name]] <- gs_mat
  }
  object@geneScores <- geneScores

  # Compute cell scores: Z %*% w, per-slide z-normalized
  cellScores <- list()
  for (ct in cts) {
    cell_flat_name <- .createCellScoresName(sigma, ct, slide = NULL)
    ct_idx <- cell_types_vec == ct
    ct_cells <- cell_names[ct_idx]
    cs_mat <- .createScoreMatrix(length(ct_cells), nCC, row_names = ct_cells)

    for (cc in seq_len(nCC)) {
      all_scores <- rep(NA_real_, length(ct_cells))
      names(all_scores) <- ct_cells

      for (s in slides) {
        Z_s <- Z_by_slide[[s]][[ct]]
        if (is.null(Z_s)) next
        s_cells <- rownames(Z_s)
        raw <- as.numeric(Z_s %*% w_list[[ct]][, cc, drop = FALSE])
        # Per-slide z-normalization
        s_mean <- mean(raw)
        s_sd <- sd(raw)
        if (s_sd > 0) {
          raw <- (raw - s_mean) / s_sd
        }
        all_scores[s_cells] <- raw
      }
      cs_mat[, paste0("CC_", cc)] <- all_scores
    }
    cellScores[[cell_flat_name]] <- cs_mat
  }
  object@cellScores <- cellScores

  # Add scores to metadata
  for (ct in cts) {
    cell_flat_name <- .createCellScoresName(sigma, ct, slide = NULL)
    scores_ct <- cellScores[[cell_flat_name]]
    ct_cells <- rownames(scores_ct)

    for (cc in seq_len(nCC)) {
      col_name <- paste0("cellScore_", sigma_name, "_cc_index_", cc)
      object@metaDataSub[ct_cells, col_name] <- scores_ct[ct_cells, paste0("CC_", cc)]
    }
  }

  object
}


#' Run gene-space canonical correlation analysis
#'
#' Batch-robust CCA that operates directly in gene space instead of PCA space.
#' For multi-slide data, uses average per-slide canonical correlation where
#' each slide's contribution is normalized by its own score variance, preventing
#' batch-level mean shifts from inflating the objective.
#'
#' Requires kernel matrices to already be computed via
#' \code{\link{computeKernelMatrix}}. Does NOT require
#' \code{\link{computePCA}}.
#'
#' @param object A CoPro object with kernel matrices computed.
#' @param sigma Sigma value to use (single numeric value). Must match a sigma
#'   for which kernels were computed.
#' @param nCC Number of canonical components (default 2).
#' @param clip Clipping method for gene expression: \code{"quantile"} for
#'   98th percentile (default), or a numeric value for a fixed threshold
#'   (e.g., 1 for count data).
#' @param min_prevalence Minimum fraction of cells expressing a gene
#'   (default 0.008 = 0.8 percent).
#' @param min_cells Minimum number of cells expressing a gene (default 20).
#' @param max_iter Maximum iterations per component (default 3000).
#' @param tol Convergence tolerance (default 1e-6).
#' @param verbose Print progress messages (default TRUE).
#'
#' @return The CoPro object with gene weights in \code{geneScores},
#'   cell scores in \code{cellScores}, and weight vectors in \code{skrCCAOut}.
#'
#' @details
#' The objective maximized is:
#' \deqn{f_{avg}(w) = \frac{1}{S} \sum_{s=1}^{S} \sum_{A<B}
#'   \frac{w_A^\top C_{AB}^{(s)} w_B}{\sigma_A^{(s)} \sigma_B^{(s)}}}
#' where \eqn{\sigma_A^{(s)} = \sqrt{w_A^\top C_{AA}^{(s)} w_A}} is the
#' per-slide score standard deviation. Subsequent components use Gram-Schmidt
#' deflation in weight space.
#'
#' @family spatial-pipeline
#' @seealso [runSkrCCA()], [computeKernelMatrix()]
#' @export
setGeneric(
  "runGeneSpaceCCA",
  function(object, sigma, nCC = 2, clip = "quantile",
           min_prevalence = 0.008, min_cells = 20,
           max_iter = 3000, tol = 1e-6,
           verbose = TRUE) standardGeneric("runGeneSpaceCCA")
)

#' @rdname runGeneSpaceCCA
#' @aliases runGeneSpaceCCA,CoProMulti-method
#' @export
setMethod(
  "runGeneSpaceCCA", "CoProMulti",
  function(object, sigma, nCC = 2, clip = "quantile",
           min_prevalence = 0.008, min_cells = 20,
           max_iter = 3000, tol = 1e-6,
           verbose = TRUE) {

    # Validate inputs
    if (length(object@kernelMatrices) == 0) {
      stop("Kernel matrices not found. Run computeKernelMatrix() first.")
    }
    if (!is.numeric(sigma) || length(sigma) != 1) {
      stop("sigma must be a single numeric value.")
    }
    if (!is.numeric(nCC) || nCC < 1) {
      stop("nCC must be a positive integer.")
    }

    cts <- if (length(object@cellTypesOfInterest) > 0) {
      object@cellTypesOfInterest
    } else {
      unique(object@cellTypesSub)
    }

    slides <- getSlideList(object)

    if (verbose) message("=== Gene-Space CCA (P1b) ===")

    # Step 1: Prepare gene-space data
    if (verbose) message("Step 1: Preparing gene-space data...")
    gsd <- .prepareGeneSpaceData(
      object, clip = clip,
      min_prevalence = min_prevalence, min_cells = min_cells,
      cts = cts, slides = slides
    )

    if (verbose) {
      message(sprintf("  Genes: %d, Slides: %d, Cell types: %s",
                      length(gsd$genes), length(gsd$slides),
                      paste(cts, collapse = ", ")))
    }

    # Step 2: Precompute G x G covariance matrices
    if (verbose) message("Step 2: Precomputing covariance matrices...")
    covmats <- .precomputeCovarianceMatrices(
      gsd$Z_by_slide, object@kernelMatrices, sigma,
      gsd$slides, cts
    )

    # Step 3: Run P1b optimization
    if (verbose) message("Step 3: P1b power iteration...")

    if (verbose) message("  Finding CC 1 ...")
    w_list <- optimize_genespace_avg_corr(
      C_self_slide = covmats$C_self,
      C_cross_slide = covmats$C_cross,
      slides = gsd$slides,
      cell_types = cts,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    )

    if (nCC > 1) {
      w_list <- optimize_genespace_avg_corr_n(
        C_self_slide = covmats$C_self,
        C_cross_slide = covmats$C_cross,
        slides = gsd$slides,
        cell_types = cts,
        w_list = w_list,
        nCC = nCC,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose
      )
    }

    # Step 4: Store results
    if (verbose) message("Step 4: Storing results...")
    object <- .storeGeneSpaceCCAResults(
      object, w_list, gsd$Z_by_slide, sigma,
      cts, nCC, gsd$genes, gsd$slides
    )

    if (verbose) message("Done.")
    object
  }
)
