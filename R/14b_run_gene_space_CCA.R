# Per-slide minimum number of cells per cell type. Below this the
# G x G covariance from the slide is too noisy to be useful (rank <= n-1).
# Slides failing this for any requested cell type are dropped with a warning.
.min_cells_per_slide <- 10

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
  if (is.character(clip) && length(clip) == 1 && clip == "quantile") {
    q98 <- quantile(as.numeric(expr[expr > 0]), 0.98)
    expr[expr > q98] <- q98
  } else if (is.numeric(clip) && length(clip) == 1) {
    expr[expr > clip] <- clip
  } else {
    stop("clip must be \"quantile\" or a single numeric threshold. Got: ",
         deparse(clip))
  }

  # Build Z_by_slide: standardize per slide per cell type
  Z_by_slide <- setNames(vector("list", length(slides)), slides)

  for (s in slides) {
    Z_by_slide[[s]] <- setNames(vector("list", length(cts)), cts)
    s_idx <- slide_ids == s

    for (ct in cts) {
      ct_idx <- s_idx & (cell_types_vec == ct)
      # Require at least .min_cells_per_slide cells per (slide, cell-type).
      # The threshold is a numerical-stability floor: with n cells, the
      # G x G cross-covariance from this slide has rank <= n-1, so for
      # G >> n the G x G estimate is dominated by sampling noise and the
      # per-slide sigma floor (1e-12) starts to bite. 10 keeps the slide-
      # level covariance well-defined while still admitting modest cell
      # populations.
      if (sum(ct_idx) < .min_cells_per_slide) next

      Z <- as.matrix(expr[ct_idx, genes, drop = FALSE])
      Z <- scale(Z, center = TRUE, scale = TRUE)
      # Per-slide zero-variance handling: a gene with no expression variance
      # on this (slide, ct) pair becomes a NaN column after scale(); we set
      # it to 0 so the slide contributes nothing for that gene. The gene is
      # still kept in the global gene set if it passes the prevalence filter,
      # so its weight will be driven by other slides where it has variance.
      Z[is.nan(Z)] <- 0
      rownames(Z) <- cell_names[ct_idx]
      Z_by_slide[[s]][[ct]] <- Z
    }
  }

  # Filter slides that have all cell types present (per-slide min-cell threshold)
  valid_slides <- slides[vapply(slides, function(s) {
    all(vapply(cts, function(ct) !is.null(Z_by_slide[[s]][[ct]]), logical(1)))
  }, logical(1))]

  dropped <- setdiff(slides, valid_slides)
  if (length(dropped) > 0) {
    for (s in dropped) {
      missing_cts <- cts[vapply(cts, function(ct) is.null(Z_by_slide[[s]][[ct]]), logical(1))]
      warning(sprintf("Slide '%s' dropped: cell type(s) %s have fewer than %d cells.",
                       s, paste(missing_cts, collapse = ", "),
                       .min_cells_per_slide))
    }
  }

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
#' @importFrom utils combn
#' @noRd
.precomputeCovarianceMatrices <- function(Z_by_slide, flat_kernels, sigma,
                                         slides, cell_types) {
  C_self <- setNames(vector("list", length(slides)), slides)
  C_cross <- setNames(vector("list", length(slides)), slides)

  pairs <- combn(cell_types, 2, simplify = FALSE)

  # Normalization convention: each side of every covariance is scaled by
  # 1/sqrt(n_side). For self-covariance both sides are the same cell type
  # so 1/sqrt(n)*1/sqrt(n) = 1/n. For cross-covariance the two sides have
  # different counts, giving 1/sqrt(n_i * n_j). The ratio in the objective
  # rho = w'C_cross w / (sigma_i * sigma_j) is therefore dimensionless and
  # invariant to per-slide cell counts.
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

      n_i <- nrow(Z_i)
      n_j <- nrow(Z_j)
      C_cross[[s]][[key]] <- crossprod(Z_i, K_ij %*% Z_j) / sqrt(n_i * n_j)
    }
  }

  list(C_self = C_self, C_cross = C_cross)
}

#' Resolve streaming distance / kernel argument lists with defaults
#' @importFrom utils modifyList
#' @noRd
.resolveStreamingArgs <- function(distanceArgs, kernelArgs) {
  d_def <- list(
    distType = "Euclidean2D",
    xDistScale = 1, yDistScale = 1, zDistScale = 1,
    normalizeDistance = TRUE,
    normalizeTarget = 0.01,
    normalizationScope = "global",  # or "per_slide"
    truncateLowDist = TRUE,
    knn_k = 10,
    geodesic_threshold = 10,
    geodesic_cutoff = 7
  )
  k_def <- list(
    lowerLimit = 1e-7,
    upperQuantile = 0.85,
    normalizeKernel = FALSE,
    minAveCellNeighor = 2,
    rowNormalizeKernel = FALSE,
    colNormalizeKernel = FALSE
  )

  unknown_d <- setdiff(names(distanceArgs), names(d_def))
  unknown_k <- setdiff(names(kernelArgs), names(k_def))
  if (length(unknown_d) > 0) {
    stop("Unknown distanceArgs: ", paste(unknown_d, collapse = ", "),
         ". Allowed: ", paste(names(d_def), collapse = ", "))
  }
  if (length(unknown_k) > 0) {
    stop("Unknown kernelArgs: ", paste(unknown_k, collapse = ", "),
         ". Allowed: ", paste(names(k_def), collapse = ", "))
  }

  d <- modifyList(d_def, distanceArgs)
  k <- modifyList(k_def, kernelArgs)

  d$distType <- match.arg(d$distType,
                          c("Euclidean2D", "Euclidean3D", "Morphology-Aware"))
  d$normalizationScope <- match.arg(d$normalizationScope,
                                    c("per_slide", "global"))
  if (any(c(d$xDistScale, d$yDistScale, d$zDistScale) <= 0)) {
    stop("Distance scales must be positive.")
  }
  if (!is.numeric(d$normalizeTarget) || length(d$normalizeTarget) != 1 ||
      !is.finite(d$normalizeTarget) || d$normalizeTarget <= 0) {
    stop("normalizeTarget must be a positive finite scalar.")
  }
  if (k$rowNormalizeKernel && k$colNormalizeKernel) {
    stop("Cannot do both row-wise and column-wise normalization.")
  }
  if (k$lowerLimit <= 0 || k$lowerLimit >= 1) {
    stop("lowerLimit must be in (0, 1).")
  }
  if (k$upperQuantile <= 0 || k$upperQuantile >= 1) {
    stop("upperQuantile must be in (0, 1).")
  }

  list(dist = d, kern = k)
}

#' Compute one slide's pairwise distance matrix, on demand
#' Used by the streaming covariance precompute. No slot storage.
#' @noRd
.streamingPairDistance <- function(object, slide, ct_i, ct_j, dist_args) {
  d <- dist_args
  mat1 <- .getCoordinateMatrix(object, ct_i, d$distType,
                               d$xDistScale, d$yDistScale, d$zDistScale,
                               slideID = slide)
  mat2 <- .getCoordinateMatrix(object, ct_j, d$distType,
                               d$xDistScale, d$yDistScale, d$zDistScale,
                               slideID = slide)

  distances_ij <- fields::rdist(mat1, mat2)
  dimnames(distances_ij) <- list(rownames(mat1), rownames(mat2))

  if (d$distType == "Morphology-Aware") {
    slide_indices <- which(.getSlideIndices(object, slide))
    all_coords <- cbind(
      object@locationDataSub$x[slide_indices] * d$xDistScale,
      object@locationDataSub$y[slide_indices] * d$yDistScale
    )
    rownames(all_coords) <- rownames(object@locationDataSub)[slide_indices]
    knn_adj <- .computeKnnGraph(all_coords, k = d$knn_k)
    geodesic_all <- .computeGeodesicDistance(knn_adj)
    idx_i <- match(rownames(mat1), rownames(all_coords))
    idx_j <- match(rownames(mat2), rownames(all_coords))
    distances_ij <- .applyMorphologyFilter(
      d_E = distances_ij,
      d_g = geodesic_all[idx_i, idx_j, drop = FALSE],
      geodesic_threshold = d$geodesic_threshold,
      geodesic_cutoff = d$geodesic_cutoff,
      verbose = FALSE
    )
  }

  distances_ij
}

#' Streaming per-slide covariance precompute (no slot storage)
#'
#' Single-pass over slides. For each slide and each cell-type pair, computes
#' the pairwise distance, derives a per-slide scaling factor from the
#' minimum low-percentile across that slide's pairs, builds the kernel,
#' and reduces to a G x G cross-covariance. The n x n distance and kernel
#' matrices are released before the next slide is processed, so the
#' resident n x n footprint stays at one slide's worth of pairs at a time.
#'
#' Differs from \code{.precomputeCovarianceMatrices + computeDistance +
#' computeKernelMatrix}: distance normalization is computed per slide
#' rather than globally across slides. For datasets where all slides
#' share the same coordinate scale (e.g., the same imaging modality
#' across patients), the per-slide factors are nearly identical and the
#' cross-covariances differ only at the level of the per-slide percentile
#' jitter.
#'
#' @param object CoProMulti object (already subsetted to cells of interest)
#' @param Z_by_slide Standardized expression list from \code{.prepareGeneSpaceData}
#' @param sigma Numeric scalar, kernel bandwidth
#' @param slides Slide IDs to process
#' @param cell_types Cell type names
#' @param distanceArgs Named list overriding distance defaults
#' @param kernelArgs Named list overriding kernel defaults
#' @param verbose Logical
#' @return List with \code{C_self} and \code{C_cross}, identical structure
#'   to \code{.precomputeCovarianceMatrices}.
#' @importFrom utils combn
#' @noRd
.streamingCovariancePrecompute <- function(object, Z_by_slide, sigma,
                                           slides, cell_types,
                                           distanceArgs = list(),
                                           kernelArgs = list(),
                                           verbose = TRUE) {
  args <- .resolveStreamingArgs(distanceArgs, kernelArgs)
  d <- args$dist
  k <- args$kern

  .check_dist_type(d$distType, object)

  pairs <- combn(cell_types, 2, simplify = FALSE)
  scope <- d$normalizationScope

  C_self <- setNames(vector("list", length(slides)), slides)
  C_cross <- setNames(vector("list", length(slides)), slides)

  # See .precomputeCovarianceMatrices for the 1/n vs 1/sqrt(n_i*n_j) note.
  for (s in slides) {
    C_self[[s]] <- setNames(vector("list", length(cell_types)), cell_types)
    for (ct in cell_types) {
      Z <- Z_by_slide[[s]][[ct]]
      C_self[[s]][[ct]] <- crossprod(Z) / nrow(Z)
    }
  }

  # ---- Phase 1: low-percentile estimation ----
  # In all cases, each (slide, pair) distance is materialized once, its low
  # percentile is recorded, then the matrix is freed. The only difference
  # between scopes is what we do with the percentiles afterwards.
  per_slide_percentiles <- vector("list", length(slides))
  names(per_slide_percentiles) <- slides

  if (isTRUE(d$normalizeDistance)) {
    if (verbose && scope == "global") {
      message("  Streaming phase 1: per-pair percentiles (global scope)...")
    }
    for (s in slides) {
      pcts <- numeric(length(pairs))
      for (pp in seq_along(pairs)) {
        ct_i <- pairs[[pp]][1]
        ct_j <- pairs[[pp]][2]
        dist_mat <- .streamingPairDistance(object, s, ct_i, ct_j, d)
        proc <- .processDistanceMatrix(dist_mat, d$truncateLowDist)
        pcts[pp] <- proc$percentile
        rm(dist_mat, proc)
        gc(verbose = FALSE, full = TRUE)
      }
      per_slide_percentiles[[s]] <- pcts
    }

    if (scope == "global") {
      all_pcts <- unlist(per_slide_percentiles)
      finite_pcts <- all_pcts[is.finite(all_pcts) & !is.na(all_pcts)]
      if (length(finite_pcts) == 0) {
        stop("Streaming: no valid distance percentile across slides.")
      }
      global_scaling <- d$normalizeTarget / min(finite_pcts)
      if (verbose) {
        message(sprintf("  Streaming GLOBAL scaling factor = %g", global_scaling))
      }
    }
  }

  # ---- Phase 2: per-slide kernel + G x G reduction ----
  for (s in slides) {
    if (verbose) message(sprintf("  Streaming slide: %s", s))

    if (isTRUE(d$normalizeDistance)) {
      if (scope == "per_slide") {
        finite_pcts <- per_slide_percentiles[[s]]
        finite_pcts <- finite_pcts[is.finite(finite_pcts) & !is.na(finite_pcts)]
        if (length(finite_pcts) == 0) {
          stop(sprintf("Streaming: no valid distance percentile for slide '%s'.", s))
        }
        slide_scaling <- d$normalizeTarget / min(finite_pcts)
        if (verbose) {
          message(sprintf("    per-slide scaling factor = %g", slide_scaling))
        }
      } else {  # "global"
        slide_scaling <- global_scaling
      }
    } else {
      slide_scaling <- 1
    }

    # Per-pair distance + kernel + reduction; each pair's n x n matrices
    # live only until that pair's G x G covariance is computed.
    C_cross[[s]] <- list()
    for (pp in seq_along(pairs)) {
      ct_i <- pairs[[pp]][1]
      ct_j <- pairs[[pp]][2]
      key <- paste0(ct_i, "-", ct_j)

      dist_mat <- .streamingPairDistance(object, s, ct_i, ct_j, d)
      proc <- .processDistanceMatrix(dist_mat, d$truncateLowDist)
      dist_mat <- proc$distances * slide_scaling
      rm(proc)

      kernel_mat <- kernel_from_distance(
        sigma = sigma, dist_mat = dist_mat,
        lower_limit = k$lowerLimit
      )
      rm(dist_mat)
      gc(verbose = FALSE, full = TRUE)

      kernel_mat <- .processKernelMatrix(
        kernel_mat, k$lowerLimit, k$upperQuantile,
        k$normalizeKernel, k$rowNormalizeKernel, k$colNormalizeKernel
      )

      Z_i <- Z_by_slide[[s]][[ct_i]]
      Z_j <- Z_by_slide[[s]][[ct_j]]
      n_i <- nrow(Z_i)
      n_j <- nrow(Z_j)

      C_cross[[s]][[key]] <- crossprod(Z_i, kernel_mat %*% Z_j) /
        sqrt(n_i * n_j)

      rm(kernel_mat)
      gc(verbose = FALSE, full = TRUE)
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
  # Gene-space CCA stores its weights in @skrCCAOut under a "gscca_"-prefixed
  # key so they cannot collide with runSkrCCA's "sigma_"-prefixed keys. The
  # two CCA flavors live in different spaces (PCA vs gene), so feeding
  # gene-space weights to computeGeneAndCellScores() (which applies a PCA
  # back-projection via @pcaGlobal$rotation) would silently produce garbage.
  # The prefix prevents that by making the keys non-overlapping; the guard
  # in .checkInputGAC catches the misuse with a clear error.
  gscca_name <- paste0("gscca_", sigma_name)
  cell_types_vec <- object@cellTypesSub
  slide_ids <- getSlideID(object)
  cell_names <- rownames(object@metaDataSub)

  # Store raw weight vectors in skrCCAOut (merge, don't overwrite)
  cca_out <- object@skrCCAOut
  cca_out[[gscca_name]] <- w_list
  object@skrCCAOut <- cca_out
  object@nCC <- nCC

  # Compute gene scores: for gene-space CCA, weights ARE gene scores directly
  geneScores <- object@geneScores
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
  cellScores <- object@cellScores
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
#' @param streaming Logical. If \code{TRUE}, fuse distance + kernel +
#'   covariance reduction into a per-slide loop and free the n x n matrices
#'   between slides. Bypasses \code{\link{computeDistance}} and
#'   \code{\link{computeKernelMatrix}}: \code{object@@distances} and
#'   \code{object@@kernelMatrices} are not populated. Default \code{FALSE}.
#' @param distanceArgs Named list of distance parameters passed through to
#'   the streaming path (e.g., \code{distType}, \code{normalizeDistance},
#'   \code{normalizeTarget}, \code{truncateLowDist}, scaling factors,
#'   morphology-aware parameters). Ignored when \code{streaming = FALSE}.
#' @param kernelArgs Named list of kernel parameters passed through to the
#'   streaming path (e.g., \code{lowerLimit}, \code{upperQuantile},
#'   \code{normalizeKernel}). Ignored when \code{streaming = FALSE}.
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
#' deflation in weight space. The power iteration uses a frozen-sigma
#' surrogate (sigma values held fixed at the previous iterate when computing
#' each weight update), making the algorithm an ALS-style alternating
#' maximization rather than exact coordinate ascent.
#'
#' Memory scales as
#' \eqn{O(G^2 \times S \times (C + C(C-1)/2))}
#' for precomputed covariance matrices: \eqn{S} slides, each storing \eqn{C}
#' self-covariances and \eqn{C(C-1)/2} cross-covariances of size
#' \eqn{G \times G}. For example G=5000, S=10, C=3 gives \eqn{10 \times 6}
#' matrices of \eqn{200} MB each, approximately 12 GB.
#'
#' Per-slide gene handling: each (slide, cell type) expression matrix is
#' independently centered and scaled. A gene with zero variance on a
#' particular (slide, cell type) pair is set to zero on that slide and
#' contributes nothing from it; the gene is still retained globally if it
#' passes the prevalence filter, with its weight driven by other slides
#' where it has variance. Slides where any requested cell type has fewer
#' than 10 cells are dropped with a warning.
#'
#' Storage: gene-space CCA results live in \code{@@skrCCAOut} under the
#' \code{gscca_sigma_<value>} key (distinct from \code{runSkrCCA}'s
#' \code{sigma_<value>} keys) so the two CCA flavors do not collide.
#' \code{computeGeneAndCellScores()} only operates on \code{runSkrCCA}
#' (PCA-space) outputs; gene-space CCA already populates
#' \code{@@geneScores} and \code{@@cellScores} directly, so no further call
#' is needed.
#'
#' Streaming mode (\code{streaming = TRUE}) reduces peak memory to roughly
#' one slide's pairwise n x n footprint at a time. Distance normalization
#' defaults to \code{normalizationScope = "global"} (a single factor across
#' all slides, matching the slot-based pipeline's semantics under
#' \code{normalizeDistance = TRUE}); under a fixed RNG seed the streaming
#' result is bit-identical to the slot-based path. Pass
#' \code{distanceArgs = list(normalizationScope = "per_slide")} to opt into
#' per-slide factors instead -- useful when slides have different coordinate
#' scales, but be aware the per-slide percentile jitter can perturb
#' degenerate canonical components on heterogeneous datasets. Use
#' \code{streaming = FALSE} when downstream code reads
#' \code{object@@distances} or \code{object@@kernelMatrices}.
#'
#' @family spatial-pipeline
#' @seealso [runSkrCCA()], [computeKernelMatrix()]
#' @export
setGeneric(
  "runGeneSpaceCCA",
  function(object, sigma, nCC = 2, clip = "quantile",
           min_prevalence = 0.008, min_cells = 20,
           max_iter = 3000, tol = 1e-6,
           streaming = FALSE,
           distanceArgs = list(),
           kernelArgs = list(),
           verbose = TRUE) standardGeneric("runGeneSpaceCCA")
)

#' @rdname runGeneSpaceCCA
#' @aliases runGeneSpaceCCA,CoPro-method
#' @export
setMethod(
  "runGeneSpaceCCA", "CoPro",
  function(object, sigma, nCC = 2, clip = "quantile",
           min_prevalence = 0.008, min_cells = 20,
           max_iter = 3000, tol = 1e-6,
           streaming = FALSE,
           distanceArgs = list(),
           kernelArgs = list(),
           verbose = TRUE) {
    stop("runGeneSpaceCCA requires a CoProMulti object (multi-slide data). ",
         "Got: ", class(object)[1])
  }
)

#' @rdname runGeneSpaceCCA
#' @aliases runGeneSpaceCCA,CoProMulti-method
#' @export
setMethod(
  "runGeneSpaceCCA", "CoProMulti",
  function(object, sigma, nCC = 2, clip = "quantile",
           min_prevalence = 0.008, min_cells = 20,
           max_iter = 3000, tol = 1e-6,
           streaming = FALSE,
           distanceArgs = list(),
           kernelArgs = list(),
           verbose = TRUE) {

    # Validate inputs
    if (!is.logical(streaming) || length(streaming) != 1 || is.na(streaming)) {
      stop("streaming must be a single logical value.")
    }
    if (!is.list(distanceArgs)) {
      stop("distanceArgs must be a named list.")
    }
    if (!is.list(kernelArgs)) {
      stop("kernelArgs must be a named list.")
    }
    if (!streaming) {
      if (length(object@kernelMatrices) == 0) {
        stop("Kernel matrices not found. Run computeKernelMatrix() first, ",
             "or call with streaming = TRUE.")
      }
      if (!is.numeric(sigma) || length(sigma) != 1) {
        stop("sigma must be a single numeric value.")
      }
      if (!sigma %in% object@sigmaValues) {
        stop(sprintf(
          "sigma = %g not found in object@sigmaValues. Available: %s",
          sigma, paste(object@sigmaValues, collapse = ", ")
        ))
      }
      if (length(distanceArgs) > 0 || length(kernelArgs) > 0) {
        warning("distanceArgs / kernelArgs are ignored when streaming = FALSE.")
      }
    } else {
      if (!is.numeric(sigma) || length(sigma) != 1) {
        stop("sigma must be a single numeric value.")
      }
    }
    if (!is.numeric(nCC) || length(nCC) != 1 || nCC < 1 || nCC != as.integer(nCC)) {
      stop("nCC must be a positive integer.")
    }

    cts <- if (length(object@cellTypesOfInterest) > 0) {
      object@cellTypesOfInterest
    } else {
      unique(object@cellTypesSub)
    }

    if (length(cts) < 2) {
      stop("runGeneSpaceCCA requires at least 2 cell types. Found: ",
           paste(cts, collapse = ", "))
    }

    slides <- getSlideList(object)

    if (verbose) message("=== Gene-Space CCA ===")

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
    if (streaming) {
      if (verbose) message("Step 2: Streaming covariance precompute...")
      covmats <- .streamingCovariancePrecompute(
        object, gsd$Z_by_slide, sigma,
        gsd$slides, cts,
        distanceArgs = distanceArgs,
        kernelArgs = kernelArgs,
        verbose = verbose
      )
      if (!sigma %in% object@sigmaValues) {
        object@sigmaValues <- c(object@sigmaValues, sigma)
      }
    } else {
      if (verbose) message("Step 2: Precomputing covariance matrices...")
      covmats <- .precomputeCovarianceMatrices(
        gsd$Z_by_slide, object@kernelMatrices, sigma,
        gsd$slides, cts
      )
    }

    # Step 3: Power iteration for canonical components
    if (nCC > length(gsd$genes)) {
      stop(sprintf("nCC (%d) exceeds number of genes after filtering (%d).",
                    nCC, length(gsd$genes)))
    }

    if (verbose) message("Step 3: Power iteration for canonical components...")

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
