# =============================================================================
# Fixed-radius near-neighbor search via grid bucketing (dependency-free)
# -----------------------------------------------------------------------------
# These helpers enumerate, without ever forming a dense n x n distance matrix,
# every pair of cells whose Euclidean distance is <= r. They are the engine
# behind the sparse Gaussian-kernel construction in 12b_sparse_kernel.R.
#
# The Gaussian kernel exp(-0.5 (d/sigma)^2) drops below `lowerLimit` when
#   d > r = sigma * sqrt(-2 * log(lowerLimit)),
# so enumerating only pairs with d <= r reproduces the dense-then-threshold
# result exactly (it is lossless), while cost scales with the number of
# near pairs rather than n^2.
#
# Spatial transcriptomics coordinates are 2-D or 3-D, so a uniform grid with
# cell side = r is an exact and simple accelerator: two points within distance
# r differ by at most one grid cell per axis, so each query cell only needs to
# inspect its own cell and the surrounding 1-ring (3^d cells).
# =============================================================================

#' Encode integer cell indices (0-based, per axis) into a single integer key
#'
#' @param cells integer matrix (n x d) of per-axis cell indices in [0, G_axis-1]
#' @param G integer vector length d of per-axis grid dimensions
#' @return integer vector of length n
#' @noRd
.encodeCellKeys <- function(cells, G) {
  d <- ncol(cells)
  key <- cells[, 1]
  if (d >= 2) {
    mult <- 1
    for (ax in 2:d) {
      mult <- mult * G[ax - 1]
      key <- key + cells[, ax] * mult
    }
  }
  key
}

#' Assign points to grid cells of side `r` using a shared origin/grid
#'
#' @param coords numeric matrix (n x d)
#' @param origin numeric vector length d (min corner, shared across A and B)
#' @param r positive cell side length
#' @param G integer vector length d (per-axis grid dimensions, shared)
#' @return integer matrix (n x d) of per-axis cell indices, clamped to range
#' @noRd
.assignCells <- function(coords, origin, r, G) {
  d <- ncol(coords)
  cells <- matrix(0L, nrow = nrow(coords), ncol = d)
  for (ax in seq_len(d)) {
    idx <- as.integer(floor((coords[, ax] - origin[ax]) / r))
    # clamp into [0, G_axis - 1] to stay within the encoded grid
    idx[idx < 0L] <- 0L
    idx[idx > (G[ax] - 1L)] <- G[ax] - 1L
    cells[, ax] <- idx
  }
  cells
}

#' Fixed-radius neighbor search (grid-bucketing), 2-D or 3-D.
#'
#' Returns every pair (a in A, b in B) with Euclidean distance <= r, as parallel
#' index/distance vectors. Each qualifying pair is emitted exactly once.
#'
#' @param A numeric matrix (nA x d) of query coordinates.
#' @param B numeric matrix (nB x d) of index coordinates, or NULL for the
#'   within-set case (B = A). In the within-set case, self pairs (a == a) are
#'   excluded and both directions (a,b) and (b,a) are emitted, so the assembled
#'   matrix is symmetric with a zero diagonal.
#' @param r positive radius (same units as the coordinates).
#' @return list(i = integer row index into A, j = integer index into B (or A),
#'   d = numeric Euclidean distance). May be empty.
#' @noRd
.frnnGrid <- function(A, B = NULL, r) {
  within <- is.null(B)
  if (within) B <- A

  nA <- nrow(A)
  nB <- nrow(B)
  d <- ncol(A)
  empty <- list(i = integer(0), j = integer(0), d = numeric(0))
  if (nA == 0L || nB == 0L) return(empty)
  if (!is.finite(r) || r <= 0) stop("radius r must be a positive finite number")

  # Shared origin and grid extent across both point sets.
  origin <- pmin(apply(A, 2, min), apply(B, 2, min))
  upper  <- pmax(apply(A, 2, max), apply(B, 2, max))
  span <- upper - origin
  # number of cells per axis (at least 1); +1 guards the upper edge
  G <- as.integer(pmax(1L, floor(span / r) + 1L))

  cellA <- .assignCells(A, origin, r, G)
  cellB <- .assignCells(B, origin, r, G)

  keyB <- .encodeCellKeys(cellB, G)
  ukeys <- sort(unique(keyB))
  bucketList <- split(seq_len(nB), factor(keyB, levels = ukeys))

  # 1-ring neighbor offsets in {-1,0,1}^d
  offset_grid <- as.matrix(expand.grid(rep(list(c(-1L, 0L, 1L)), d)))

  i_acc <- vector("list", nrow(offset_grid))
  j_acc <- vector("list", nrow(offset_grid))
  d_acc <- vector("list", nrow(offset_grid))

  r2 <- r * r
  for (oi in seq_len(nrow(offset_grid))) {
    o <- offset_grid[oi, ]
    # candidate cell for each A point under this offset
    cand <- cellA + matrix(o, nrow = nA, ncol = d, byrow = TRUE)
    # keep only candidates that fall inside the grid on every axis
    inb <- rep(TRUE, nA)
    for (ax in seq_len(d)) {
      inb <- inb & cand[, ax] >= 0L & cand[, ax] <= (G[ax] - 1L)
    }
    if (!any(inb)) next
    aq <- which(inb)
    qk <- .encodeCellKeys(cand[aq, , drop = FALSE], G)
    m <- match(qk, ukeys)
    hit <- !is.na(m)
    if (!any(hit)) next
    aq <- aq[hit]
    mb <- m[hit]
    bcounts <- lengths(bucketList)[mb]
    if (sum(bcounts) == 0) next
    a_rep <- rep(aq, bcounts)
    b_idx <- unlist(bucketList[mb], use.names = FALSE)
    # squared Euclidean distance for the candidate pairs
    diff <- A[a_rep, , drop = FALSE] - B[b_idx, , drop = FALSE]
    dd2 <- rowSums(diff * diff)
    keep <- dd2 <= r2
    if (within) keep <- keep & (a_rep != b_idx)
    if (!any(keep)) next
    i_acc[[oi]] <- a_rep[keep]
    j_acc[[oi]] <- b_idx[keep]
    d_acc[[oi]] <- sqrt(dd2[keep])
  }

  list(
    i = unlist(i_acc, use.names = FALSE),
    j = unlist(j_acc, use.names = FALSE),
    d = unlist(d_acc, use.names = FALSE)
  )
}

#' Heuristic seed radius expected to capture roughly `target` near pairs
#'
#' Used only to seed the radius-doubling loop in the percentile pre-pass; the
#' loop corrects any under-estimate, so this only needs to be in the right
#' ballpark.
#' @noRd
.seedRadius <- function(A, B, target, within) {
  d <- ncol(A)
  nA <- nrow(A)
  Bc <- if (within || is.null(B)) A else B
  nB <- nrow(Bc)
  origin <- pmin(apply(A, 2, min), apply(Bc, 2, min))
  upper  <- pmax(apply(A, 2, max), apply(Bc, 2, max))
  span <- upper - origin
  span[span <= 0] <- .Machine$double.eps
  vol <- prod(span)
  # as.numeric() avoids integer overflow: the pair count exceeds
  # .Machine$integer.max (~2.1e9) once a block has more than ~46k cells.
  pairs_total <- if (within) as.numeric(nA) * (as.numeric(nA) - 1) else
    as.numeric(nA) * as.numeric(nB)
  if (pairs_total <= 0) return(max(span))
  frac <- min(1, max(target, 1) / pairs_total)
  # invert ball-volume fraction: c_d * r^d / vol = frac  (c_2 = pi, c_3 = 4pi/3)
  c_d <- if (d == 2) pi else if (d == 3) 4 * pi / 3 else 2^d
  r <- (frac * vol / c_d)^(1 / d)
  if (!is.finite(r) || r <= 0) r <- max(span) / max(2, nA^(1 / d))
  r
}

#' Exact R type-7 quantile of the lowest values from a retained subset
#'
#' The dense pipeline computes `quantile(finite_distances, p)` over ALL N finite
#' non-zero distances. For tiny `p` the relevant order statistics are the very
#' smallest distances, which a radius search retains in full. Given the sorted
#' retained distances and the TRUE total count N, this reproduces R's default
#' (type 7) quantile exactly.
#'
#' @param sorted_small numeric vector, ascending, of the retained distances
#'   (must contain at least the floor(h)+1 smallest, where h is computed below).
#' @param N integer true total count of finite non-zero distances.
#' @param p numeric probability in [0, 1].
#' @return numeric quantile value.
#' @noRd
.exactLowQuantile <- function(sorted_small, N, p) {
  if (N < 1) stop("N must be >= 1 for quantile computation")
  if (N == 1) return(sorted_small[1])
  h <- (N - 1) * p + 1
  lo <- floor(h)
  frac <- h - lo
  if (lo > length(sorted_small)) {
    stop("Not enough retained distances to reproduce the dense percentile.")
  }
  if (frac == 0) {
    return(sorted_small[lo])
  }
  if (lo + 1L > length(sorted_small)) {
    stop("Not enough retained distances to reproduce the dense percentile.")
  }
  sorted_small[lo] + frac * (sorted_small[lo + 1L] - sorted_small[lo])
}

#' Compute the low distance percentile for one cell-type block, exactly,
#' without forming the dense distance matrix.
#'
#' Mirrors `.processDistanceMatrix()` in 11_compute_distance.R: zero distances
#' are replaced by the smallest non-zero distance before the quantile is taken,
#' and the quantile probability matches the dense path (`min(1e-3, 2/max(dim))`
#' for cross-type pairs, a fixed value such as 1e-4 for within-type).
#'
#' @param A,B coordinate matrices (B = NULL for within-type).
#' @param p quantile probability used by the dense path for this block.
#' @return list(percentile = numeric, min_nonzero = numeric or NA).
#' @noRd
.lowPercentileBlock <- function(A, B, p) {
  within <- is.null(B)
  nA <- nrow(A)
  nB <- if (within) nA else nrow(B)
  # as.numeric() avoids integer overflow of the true pair count N (used only as
  # a count for the type-7 quantile rank), which would otherwise become NA past
  # ~46k cells and silently break the percentile reconstruction.
  N <- if (within) as.numeric(nA) * (as.numeric(nA) - 1) else
    as.numeric(nA) * as.numeric(nB)
  if (N < 1) stop("Block has too few cells to compute a distance percentile.")

  h <- (N - 1) * p + 1
  need <- floor(h) + 1L

  Bc <- if (within) A else B
  origin <- pmin(apply(A, 2, min), apply(Bc, 2, min))
  upper  <- pmax(apply(A, 2, max), apply(Bc, 2, max))
  max_span <- sqrt(sum((upper - origin)^2))

  r <- .seedRadius(A, B, target = 2 * need, within = within)
  d <- NULL
  repeat {
    tri <- .frnnGrid(A, B, r)
    d <- tri$d
    if (length(d) >= need) break
    if (r >= max_span) {
      # radius now spans the whole block; we have every pair there is
      break
    }
    r <- min(r * 2, max_span)
  }

  min_nonzero <- NA_real_
  if (length(d) == 0) {
    stop("No non-zero distances found while computing distance percentile.")
  }
  if (any(d == 0)) {
    nz <- d[d > 0]
    if (length(nz) == 0) {
      stop("All retained distances are zero; cannot compute percentile.")
    }
    min_nonzero <- min(nz)
    d[d == 0] <- min_nonzero
  }

  ds <- sort(d)
  list(
    percentile = .exactLowQuantile(ds, N, p),
    min_nonzero = min_nonzero
  )
}
