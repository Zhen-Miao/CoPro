# Direct correctness tests for the fixed-radius near-neighbor (FRNN) engine that
# underlies the sparse kernel. We compare .frnnGrid() against an O(n^2) brute
# force enumeration of all pairs within radius r. This is the ground-truth check
# for the engine's exactness claim and, crucially, the only coverage of the 3D
# code path (the package advertises distType = "Euclidean3D").

# Brute-force every pair (a in A, b in B) with Euclidean distance <= r. For the
# within-set case (B = NULL) self pairs are excluded and both directions kept,
# matching the .frnnGrid contract.
.bruteFrnn <- function(A, B = NULL, r) {
  within <- is.null(B)
  if (within) B <- A
  oi <- integer(0); oj <- integer(0); od <- numeric(0)
  for (a in seq_len(nrow(A))) {
    dv <- sqrt(colSums((t(B) - A[a, ])^2))
    keep <- which(dv <= r)
    if (within) keep <- keep[keep != a]
    oi <- c(oi, rep(a, length(keep)))
    oj <- c(oj, keep)
    od <- c(od, dv[keep])
  }
  ord <- order(oi, oj)
  list(i = oi[ord], j = oj[ord], d = od[ord])
}

.canonFrnn <- function(g) {
  ord <- order(g$i, g$j)
  list(i = g$i[ord], j = g$j[ord], d = g$d[ord])
}

# Assert .frnnGrid reproduces the brute-force pair set (indices) and distances.
expect_frnn_matches <- function(A, B, r) {
  g <- .canonFrnn(CoPro:::.frnnGrid(A, B, r))
  e <- .bruteFrnn(A, B, r)
  expect_identical(g$i, e$i)
  expect_identical(g$j, e$j)
  expect_equal(g$d, e$d, tolerance = 1e-12)
}

test_that(".frnnGrid matches brute force in 2D (cross, within, clustered)", {
  set.seed(3)
  A <- cbind(runif(80), runif(80))
  B <- cbind(runif(60), runif(60))

  expect_frnn_matches(A, B, 0.25)      # cross-type
  expect_frnn_matches(A, NULL, 0.25)   # within-type (symmetric, no diagonal)

  # two tight clusters: most mass is intra-cluster, exercises uneven bucket loads
  clustered <- rbind(cbind(rnorm(40, 0, 0.05), rnorm(40, 0, 0.05)),
                     cbind(rnorm(40, 1, 0.05), rnorm(40, 1, 0.05)))
  expect_frnn_matches(clustered, NULL, 0.2)
})

test_that(".frnnGrid matches brute force in 3D (cross, within)", {
  set.seed(5)
  A <- cbind(runif(80), runif(80), runif(80))
  B <- cbind(runif(60), runif(60), runif(60))

  expect_frnn_matches(A, B, 0.45)      # cross-type, 3D
  expect_frnn_matches(A, NULL, 0.45)   # within-type, 3D
})

test_that(".frnnGrid handles coincident points (d = 0 pairs kept, self excluded)", {
  set.seed(3)
  A <- cbind(runif(80), runif(80))
  # force exact duplicates: cells 5 and 6 land on cell 1
  co <- A
  co[5, ] <- co[6, ] <- A[1, ]
  g <- .canonFrnn(CoPro:::.frnnGrid(co, NULL, 0.25))
  e <- .bruteFrnn(co, NULL, 0.25)
  expect_identical(g$i, e$i)
  expect_identical(g$j, e$j)
  expect_equal(g$d, e$d, tolerance = 1e-12)
  # the coincident, non-self pairs are present with distance 0 ...
  expect_true(any(g$d == 0))
  # ... but no point is ever paired with itself
  expect_false(any(g$i == g$j))
})

test_that(".lowPercentileBlock matches dense quantile with many coincident cells", {
  # Many coincident cells produce a large block of zero distances. The dense path
  # replaces them with the smallest non-zero distance before quantiling, so the
  # sparse reproducer must expand its radius past the zeros to recover that value
  # rather than stopping short on an all-zero retained set.
  set.seed(42)
  A <- rbind(matrix(rep(c(0.3, 0.7), each = 150), ncol = 2),   # 150 stacked cells
             cbind(runif(50), runif(50)))                      # + 50 spread out
  p <- 1e-4

  D <- as.matrix(dist(A)); diag(D) <- NA
  full <- as.numeric(D[!is.na(D)])
  mn <- min(full[full > 0]); full[full == 0] <- mn
  dense_q <- as.numeric(stats::quantile(full, p, type = 7, names = FALSE))

  res <- CoPro:::.lowPercentileBlock(A, NULL, p)
  expect_true(is.finite(res$percentile))
  expect_equal(res$percentile, dense_q, tolerance = 1e-12)
  expect_equal(res$min_nonzero, mn, tolerance = 1e-12)
})

test_that(".frnnGrid returns an empty result when no pair is within radius", {
  A <- cbind(c(0, 10), c(0, 10))
  g <- CoPro:::.frnnGrid(A, NULL, 0.5)
  expect_length(g$i, 0L)
  expect_length(g$j, 0L)
  expect_length(g$d, 0L)
})

test_that("compiled and R fixed-radius engines return identical triplets", {
  set.seed(808)
  old_options <- options(CoPro.useRcppFRNN = FALSE)
  on.exit(options(old_options), add = TRUE)

  for (dimensions in 2:3) {
    A <- matrix(runif(240 * dimensions), ncol = dimensions)
    B <- matrix(runif(180 * dimensions), ncol = dimensions)
    radius <- if (dimensions == 2) 0.18 else 0.3

    for (within in c(FALSE, TRUE)) {
      B_current <- if (within) NULL else B
      r_result <- .canonFrnn(CoPro:::.frnnGrid(A, B_current, radius))
      options(CoPro.useRcppFRNN = TRUE)
      cpp_result <- .canonFrnn(CoPro:::.frnnGrid(A, B_current, radius))
      options(CoPro.useRcppFRNN = FALSE)

      expect_identical(cpp_result$i, r_result$i)
      expect_identical(cpp_result$j, r_result$j)
      expect_equal(cpp_result$d, r_result$d, tolerance = 1e-12)
    }
  }
})
