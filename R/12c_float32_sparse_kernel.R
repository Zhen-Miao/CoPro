# =============================================================================
# Direct float32 sparse kernels
# -----------------------------------------------------------------------------
# CoProFloat32SparseMatrix is a compact row-compressed kernel:
#   p   : zero-based row pointer (integer)
#   j   : zero-based column indices (integer)
#   x   : IEEE-754 float32 values packed into a raw vector
#   Dim : stored matrix dimensions
#   symmetric: whether the strict upper triangle represents a symmetric matrix
#
# Kernels remain encoded through multiplication. The compiled X' K X operator
# accumulates in float32 across disjoint row blocks and returns only the small
# PC-space result as an ordinary R double matrix.
# =============================================================================

#' Test whether an object is a CoPro float32 sparse kernel
#' @param x An R object.
#' @return A single logical value.
#' @noRd
.isFloat32SparseKernel <- function(x) {
  inherits(x, "CoProFloat32SparseMatrix")
}

#' Construct a validated float32 CSR kernel
#' @noRd
.newFloat32SparseKernel <- function(x, dimnames = NULL) {
  required <- c("p", "j", "x", "Dim", "transposed")
  if (!is.list(x) || !all(required %in% names(x))) {
    stop("Invalid float32 sparse-kernel payload.")
  }
  if (is.null(x$symmetric)) x$symmetric <- FALSE
  if (length(x$p) != x$Dim[1] + 1L ||
      length(x$j) * 4 != length(x$x)) {
    stop("Inconsistent float32 sparse-kernel dimensions.")
  }
  if (isTRUE(x$symmetric) && x$Dim[1] != x$Dim[2]) {
    stop("A symmetric float32 sparse kernel must be square.")
  }
  x$Dimnames <- dimnames
  class(x) <- c("CoProFloat32SparseMatrix", "list")
  x
}

#' Dimensions of an encoded CoPro float32 sparse matrix
#' @param x A `CoProFloat32SparseMatrix`.
#' @export
dim.CoProFloat32SparseMatrix <- function(x) {
  if (isTRUE(x$transposed)) rev(x$Dim) else x$Dim
}

#' Transpose an encoded CoPro float32 sparse matrix without copying its values
#' @param x A `CoProFloat32SparseMatrix`.
#' @method t CoProFloat32SparseMatrix
#' @export
t.CoProFloat32SparseMatrix <- function(x) {
  if (isTRUE(x$symmetric)) return(x)
  x$transposed <- !isTRUE(x$transposed)
  if (!is.null(x$Dimnames)) x$Dimnames <- rev(x$Dimnames)
  x
}

#' Dimension names of an encoded CoPro float32 sparse matrix
#' @param x A `CoProFloat32SparseMatrix`.
#' @export
dimnames.CoProFloat32SparseMatrix <- function(x) {
  x$Dimnames
}

#' Number of represented nonzero entries
#' @noRd
.float32KernelNnz <- function(x) {
  length(x$j)
}

#' CPU cores actually available to this process
#'
#' `parallel::detectCores()` reports the whole machine, which on a shared HPC
#' node is not what this job was granted -- spawning that many threads
#' oversubscribes cores allocated to other jobs. This honors the common
#' scheduler allocations first (SLURM, SGE/Grid Engine, PBS/Torque, LSF) and
#' only falls back to the machine core count on an unmanaged workstation.
#' @noRd
.detectAllocatedCores <- function() {
  env_cores <- function(name) {
    value <- suppressWarnings(as.integer(Sys.getenv(name, "")))
    if (length(value) == 1L && !is.na(value) && value >= 1L) value else NA_integer_
  }
  ## Scheduler allocations, most job-specific first.
  for (name in c("SLURM_CPUS_PER_TASK", "NSLOTS", "PBS_NUM_PPN",
                 "LSB_DJOB_NUMPROC", "PBS_NP")) {
    cores <- env_cores(name)
    if (!is.na(cores)) return(cores)
  }
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(cores) || cores < 1L) {
    cores <- suppressWarnings(parallel::detectCores())
  }
  if (is.na(cores) || cores < 1L) 1L else as.integer(cores)
}

#' Number of worker threads for float32 sparse operators
#'
#' Resolution order: an explicit `options(CoPro.float32Threads=)` always wins;
#' otherwise the count is derived automatically from the cores actually
#' allocated to this process (see [.detectAllocatedCores]), then capped by
#' `OMP_NUM_THREADS` when a cluster sets it, by CRAN's core limit during
#' `R CMD check`, and by a default ceiling of eight beyond which the
#' memory-bandwidth-bound operators stop scaling. Raise the ceiling for very
#' large jobs with `options(CoPro.float32Threads = <n>)`.
#' @noRd
.float32KernelThreads <- function() {
  option <- getOption("CoPro.float32Threads", NULL)
  if (!is.null(option)) {
    option <- suppressWarnings(as.integer(option))
    if (length(option) != 1L || is.na(option) || option < 1L) {
      stop("options(CoPro.float32Threads=) must be a positive integer.")
    }
    return(option)
  }

  cores <- .detectAllocatedCores()

  ## A cluster that limits OpenMP threads (often to the SLURM allocation) should
  ## limit ours too.
  omp <- suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", "")))
  if (length(omp) == 1L && !is.na(omp) && omp >= 1L) {
    cores <- min(cores, omp)
  }

  ## Respect CRAN's check-time core limit.
  ceiling_threads <- 8L
  check_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(check_limit) && !identical(tolower(check_limit), "false")) {
    ceiling_threads <- 2L
  }

  max(1L, min(ceiling_threads, cores))
}

#' Apply an encoded float32 sparse kernel to a dense matrix
#' @noRd
.float32KernelMatMult <- function(K, X, n_threads = .float32KernelThreads()) {
  if (!.isFloat32SparseKernel(K)) {
    return(K %*% X)
  }
  if (is.null(dim(X))) X <- matrix(X, ncol = 1L)
  float32_csr_matmul_cpp(
    K$p, K$j, K$x, K$Dim, X,
    n_threads = n_threads,
    transpose = isTRUE(K$transposed),
    symmetric = isTRUE(K$symmetric)
  )
}

#' Form X_left' K X_right without decoding K
#' @noRd
.kernelXKY <- function(
    X_left, K, X_right,
    n_threads = .float32KernelThreads()) {
  if (!.isFloat32SparseKernel(K)) {
    return(as.matrix(crossprod(X_left, K %*% X_right)))
  }
  if (isTRUE(K$transposed) && !isTRUE(K$symmetric)) {
    return(t(float32_csr_xky_cpp(
      K$p, K$j, K$x, K$Dim,
      X_right, X_left, n_threads, symmetric = FALSE
    )))
  }
  float32_csr_xky_cpp(
    K$p, K$j, K$x, K$Dim,
    X_left, X_right, n_threads,
    symmetric = isTRUE(K$symmetric)
  )
}

#' Row and column sums for an encoded float32 sparse kernel
#' @noRd
.float32KernelSums <- function(K) {
  sums <- float32_csr_sums_cpp(
    K$p, K$j, K$x, K$Dim, symmetric = isTRUE(K$symmetric)
  )
  if (isTRUE(K$transposed) && !isTRUE(K$symmetric)) {
    list(
      rowSums = sums$colSums,
      colSums = sums$rowSums,
      sumSquares = sums$sumSquares
    )
  } else {
    sums
  }
}

#' Decode a float32 kernel into a standard sparse Matrix
#'
#' This diagnostic conversion materializes float64 values and should not be
#' used in a large-data computation.
#'
#' @param x A `CoProFloat32SparseMatrix`.
#' @return A `Matrix::dgCMatrix`, or a `Matrix::dsCMatrix` for a symmetric
#'   one-type kernel.
#' @export
asDoubleSparseMatrix <- function(x) {
  if (!.isFloat32SparseKernel(x)) {
    stop("x must be a CoProFloat32SparseMatrix.")
  }
  row_counts <- diff(x$p)
  decoded <- Matrix::sparseMatrix(
    i = rep.int(seq_len(x$Dim[1]), row_counts),
    j = x$j + 1L,
    x = float32_csr_values_cpp(x$x),
    dims = x$Dim,
    dimnames = if (
      isTRUE(x$transposed) && !isTRUE(x$symmetric)
    ) rev(x$Dimnames) else x$Dimnames,
    symmetric = isTRUE(x$symmetric)
  )
  if (isTRUE(x$transposed) && !isTRUE(x$symmetric)) t(decoded) else decoded
}

#' Materialize all encoded kernels in a CoPro object
#'
#' Converts every `CoProFloat32SparseMatrix` in `object@kernelMatrices` to a
#' standard double-precision sparse `Matrix`. This is a compatibility escape
#' hatch for third-party or older code that accesses the slot directly rather
#' than using [getKernelMatrix()]. It can temporarily require both encoded and
#' decoded copies, so it should be avoided in memory-critical analysis steps.
#'
#' @param object A `CoPro` object.
#' @param verbose Whether to report how many kernels were converted.
#' @return A copy of `object` whose kernel list contains standard sparse
#'   matrices.
#' @export
materializeFloat32Kernels <- function(object, verbose = TRUE) {
  if (!methods::is(object, "CoPro")) {
    stop("object must inherit from CoPro.")
  }
  encoded <- vapply(
    object@kernelMatrices, .isFloat32SparseKernel, logical(1)
  )
  if (!any(encoded)) return(object)
  object@kernelMatrices[encoded] <- lapply(
    object@kernelMatrices[encoded], asDoubleSparseMatrix
  )
  attr(object, "kernelNormalizerCache") <- NULL
  if (verbose) {
    message("Materialized ", sum(encoded), " float32 kernel(s).")
  }
  object
}

#' Compute block-streamed float32 sparse Gaussian kernels
#'
#' A large-data kernel path that never materializes a float64 sparse kernel.
#' Each cell-type block is enumerated once at the largest requested radius,
#' distances are retained temporarily as float32, and every sigma-specific
#' kernel is written directly as a row-compressed float32 byte buffer. The
#' block's neighbor data are released before the next block.
#'
#' The resulting kernels are consumed directly by [runSkrCCA()] through a
#' parallel float32 `X1' K X2` operator. Only the small PC-space result is
#' converted to float64. Unnormalized and globally normalized one-cell-type
#' kernels store only their strict upper triangle and are applied as symmetric
#' matrices without expansion. Row- or column-normalized self-kernels are
#' expanded because those operations make the represented matrix asymmetric.
#'
#' Both single- and multi-slide objects are supported, with any positive number
#' of cell types. Multi-slide construction streams one slide/cell-type block at
#' a time and schedules the largest block first. Older consumers can request a
#' temporary standard sparse matrix through [getKernelMatrix()] or
#' [asDoubleSparseMatrix()].
#'
#' Global, row, and column normalization are applied directly to float32 values.
#' The ordinary centered Frobenius objective normalizer is also evaluated from
#' the encoded float32 kernel without materialization. Exact whitened Frobenius
#' normalization currently uses the compatibility path through temporary
#' double-precision sparse matrices when within-type operators are supplied.
#'
#' @param object A `CoProSingle` or `CoProMulti` object.
#' @param sigmaValues Positive Gaussian-kernel scale values.
#' @param lowerLimit Kernel support threshold.
#' @param upperQuantile Quantile at which large kernel values are clipped.
#' @param normalizeKernel Whether to divide all entries by the median positive
#'   row sum.
#' @param minAveCellNeighor Minimum average represented neighbors.
#' @param rowNormalizeKernel,colNormalizeKernel Whether to normalize each
#'   nonempty row or column to sum to one. They cannot both be `TRUE`.
#' @param distType `"Euclidean2D"` or `"Euclidean3D"`.
#' @param xDistScale,yDistScale,zDistScale Per-axis coordinate scales.
#' @param normalizeDistance Whether to normalize the low distance percentile.
#' @param normalizeTarget Target low distance percentile after normalization.
#' @param truncateLowDist Whether to floor very small distances.
#' @param overwrite Whether to replace existing kernel matrices.
#' @param verbose Whether to report progress.
#' @return The object with encoded float32 kernels in `@kernelMatrices`.
#' @family spatial-pipeline
#' @seealso [computeSparseKernel()], [runSkrCCA()], [asDoubleSparseMatrix()]
#' @export
setGeneric(
  "computeSparseKernelFloat32",
  function(
      object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
      normalizeKernel = FALSE, minAveCellNeighor = 2,
      rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
      distType = c("Euclidean2D", "Euclidean3D"),
      xDistScale = 1, yDistScale = 1, zDistScale = 1,
      normalizeDistance = TRUE, normalizeTarget = 0.01,
      truncateLowDist = TRUE, overwrite = TRUE,
      verbose = TRUE) {
    standardGeneric("computeSparseKernelFloat32")
  }
)

#' Build a list of float32 construction blocks
#' @noRd
.float32KernelBlocks <- function(object, cts, is_multi) {
  slides <- if (is_multi) getSlideList(object) else NA_character_
  if (is_multi && length(slides) == 0L) {
    stop("No slides found in multi-slide object.")
  }
  blocks <- list()
  for (slide_key in slides) {
    slide <- if (is.na(slide_key)) NULL else slide_key
    if (length(cts) == 1L) {
      n <- .countSlideCellType(object, slide, cts)
      if (n <= if (is_multi) 5L else 1L) next
      blocks[[length(blocks) + 1L]] <- list(
        slide = slide, cellType1 = cts, cellType2 = cts,
        n1 = n, n2 = n, symmetric = TRUE
      )
    } else {
      pairs <- utils::combn(cts, 2L)
      for (pair_index in seq_len(ncol(pairs))) {
        cell_type_1 <- pairs[1L, pair_index]
        cell_type_2 <- pairs[2L, pair_index]
        n1 <- .countSlideCellType(object, slide, cell_type_1)
        n2 <- .countSlideCellType(object, slide, cell_type_2)
        if (is_multi && (n1 <= 5L || n2 <= 5L)) next
        if (n1 == 0L || n2 == 0L) next
        blocks[[length(blocks) + 1L]] <- list(
          slide = slide,
          cellType1 = cell_type_1, cellType2 = cell_type_2,
          n1 = n1, n2 = n2, symmetric = FALSE
        )
      }
    }
  }
  if (length(blocks) == 0L) {
    stop("No slide/cell-type blocks have enough cells to compute kernels.")
  }
  block_sizes <- vapply(blocks, function(block) {
    represented <- as.numeric(block$n1) * as.numeric(block$n2)
    if (block$symmetric) represented / 2 else represented
  }, numeric(1))
  blocks[order(block_sizes, decreasing = TRUE)]
}

#' Block-streamed float32 construction shared by single and multi-slide data
#' @noRd
.computeSparseKernelFloat32Core <- function(
    object, sigmaValues, lowerLimit, upperQuantile,
    normalizeKernel, minAveCellNeighor,
    rowNormalizeKernel, colNormalizeKernel,
    distType, xDistScale, yDistScale, zDistScale,
    normalizeDistance, normalizeTarget, truncateLowDist, overwrite,
    verbose, is_multi) {
  cts <- .checkInputSparseKernel(
    object, sigmaValues, lowerLimit, upperQuantile,
    minAveCellNeighor, rowNormalizeKernel, colNormalizeKernel, distType
  )
  if (!is.logical(overwrite) || length(overwrite) != 1L ||
      is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  blocks <- .float32KernelBlocks(object, cts, is_multi)
  normalization <- if (rowNormalizeKernel) {
    2L
  } else if (colNormalizeKernel) {
    3L
  } else if (normalizeKernel) {
    1L
  } else {
    0L
  }
  get_coordinates <- function(block, which_type) {
    cell_type <- block[[which_type]]
    .getCoordinateMatrix(
      object, cell_type, distType,
      xDistScale, yDistScale, zDistScale,
      slideID = block$slide
    )
  }

  # Pass 1 retains only one scalar per block. The global minimum matches the
  # distance normalization used by the ordinary sparse single/multi paths.
  need_percentile <- truncateLowDist || normalizeDistance
  percentiles <- rep(NA_real_, length(blocks))
  if (need_percentile) {
    for (index in seq_along(blocks)) {
      block <- blocks[[index]]
      coordinate_1 <- get_coordinates(block, "cellType1")
      coordinate_2 <- if (block$symmetric) {
        NULL
      } else {
        get_coordinates(block, "cellType2")
      }
      probability <- if (block$symmetric) {
        1e-4
      } else {
        .pairPercentileProb(nrow(coordinate_1), nrow(coordinate_2))
      }
      percentiles[index] <- .lowPercentileBlock(
        coordinate_1, coordinate_2, probability
      )$percentile
    }
  }
  scaling_factor <- if (normalizeDistance) {
    normalizeTarget / min(percentiles, na.rm = TRUE)
  } else {
    1
  }
  if (!is.finite(scaling_factor) || scaling_factor <= 0) {
    stop("Unable to determine a finite positive distance scale.")
  }

  kernel_matrices <- if (
    overwrite || length(object@kernelMatrices) == 0L
  ) list() else object@kernelMatrices
  invalid_sigmas <- stats::setNames(
    logical(length(sigmaValues)), as.character(sigmaValues)
  )
  valid_multi_sigmas <- stats::setNames(
    logical(length(sigmaValues)), as.character(sigmaValues)
  )
  written_names <- character()

  # Pass 2 builds all sigma values for a block, then releases its float32
  # neighbor temporary before the next persistent block accumulates.
  for (index in seq_along(blocks)) {
    block <- blocks[[index]]
    slide_suffix <- if (is.null(block$slide)) {
      ""
    } else {
      paste0(" on ", block$slide)
    }
    if (verbose) {
      message(sprintf(
        "Float32 kernel block %s -> %s%s (%d x %d)",
        block$cellType1, block$cellType2, slide_suffix,
        block$n1, block$n2
      ))
    }
    coordinate_1 <- get_coordinates(block, "cellType1")
    coordinate_2 <- if (block$symmetric) {
      coordinate_1
    } else {
      get_coordinates(block, "cellType2")
    }
    built <- float32_csr_gaussian_kernels_cpp(
      coordinate_1, coordinate_2, sigmaValues,
      percentiles[index], scaling_factor,
      lowerLimit, upperQuantile, truncateLowDist,
      symmetric = block$symmetric,
      normalization = normalization
    )

    for (sigma_index in seq_along(sigmaValues)) {
      sigma <- sigmaValues[sigma_index]
      stored_nnz <- built$nonzeros[sigma_index]
      represented_nnz <- if (block$symmetric) {
        2 * stored_nnz
      } else {
        stored_nnz
      }
      valid <- represented_nnz >=
        minAveCellNeighor * min(block$n1, block$n2)
      if (!valid) {
        if (!is_multi) invalid_sigmas[as.character(sigma)] <- TRUE
        warning(sprintf(
          paste0(
            "Float32 kernel %s -> %s%s at sigma=%g is too sparse ",
            "(%g represented entries)."
          ),
          block$cellType1, block$cellType2, slide_suffix,
          sigma, represented_nnz
        ))
        if (is_multi) next
      } else if (is_multi) {
        valid_multi_sigmas[as.character(sigma)] <- TRUE
      }

      kernel <- .newFloat32SparseKernel(
        built$kernels[[sigma_index]],
        dimnames = list(rownames(coordinate_1), rownames(coordinate_2))
      )
      name <- .createKernelMatrixName(
        sigma, block$cellType1, block$cellType2, slide = block$slide
      )
      kernel_matrices[[name]] <- kernel
      written_names <- c(written_names, name)
    }
    rm(built, coordinate_1, coordinate_2)
    invisible(gc(FALSE))
  }

  if (is_multi) {
    invalid_sigmas <- !valid_multi_sigmas
  }
  if (any(invalid_sigmas)) {
    invalid_values <- as.numeric(names(invalid_sigmas)[invalid_sigmas])
    if (!is_multi && length(sigmaValues) == 1L) {
      stop("The requested sigma produced a float32 kernel that was too sparse.")
    }
    remove <- vapply(names(kernel_matrices), function(name) {
      name %in% written_names &&
        .parseKernelMatrixName(name)$sigma %in% invalid_values
    }, logical(1))
    kernel_matrices <- kernel_matrices[!remove]
    sigmaValues <- sigmaValues[!invalid_sigmas]
  }
  if (length(sigmaValues) == 0L) {
    stop("No requested sigma produced a valid float32 sparse kernel.")
  }

  object@kernelMatrices <- kernel_matrices
  object@sigmaValues <- sigmaValues
  object@distanceScaleFactor <- scaling_factor
  object@distances <- list()
  object
}

#' @rdname computeSparseKernelFloat32
#' @export
setMethod(
  "computeSparseKernelFloat32", "CoProSingle",
  function(
      object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
      normalizeKernel = FALSE, minAveCellNeighor = 2,
      rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
      distType = c("Euclidean2D", "Euclidean3D"),
      xDistScale = 1, yDistScale = 1, zDistScale = 1,
      normalizeDistance = TRUE, normalizeTarget = 0.01,
      truncateLowDist = TRUE, overwrite = TRUE,
      verbose = TRUE) {
    .computeSparseKernelFloat32Core(
      object, sigmaValues, lowerLimit, upperQuantile,
      normalizeKernel, minAveCellNeighor,
      rowNormalizeKernel, colNormalizeKernel,
      match.arg(distType), xDistScale, yDistScale, zDistScale,
      normalizeDistance, normalizeTarget, truncateLowDist, overwrite,
      verbose, is_multi = FALSE
    )
  }
)

#' @rdname computeSparseKernelFloat32
#' @export
setMethod(
  "computeSparseKernelFloat32", "CoProMulti",
  function(
      object, sigmaValues, lowerLimit = 1e-7, upperQuantile = 0.85,
      normalizeKernel = FALSE, minAveCellNeighor = 2,
      rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE,
      distType = c("Euclidean2D", "Euclidean3D"),
      xDistScale = 1, yDistScale = 1, zDistScale = 1,
      normalizeDistance = TRUE, normalizeTarget = 0.01,
      truncateLowDist = TRUE, overwrite = TRUE,
      verbose = TRUE) {
    .computeSparseKernelFloat32Core(
      object, sigmaValues, lowerLimit, upperQuantile,
      normalizeKernel, minAveCellNeighor,
      rowNormalizeKernel, colNormalizeKernel,
      match.arg(distType), xDistScale, yDistScale, zDistScale,
      normalizeDistance, normalizeTarget, truncateLowDist, overwrite,
      verbose, is_multi = TRUE
    )
  }
)
