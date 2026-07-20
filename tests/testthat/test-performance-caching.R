# Regression tests for cached permutation setup and scoring. These compare the
# optimized paths against deliberately direct versions of the historical
# calculations so future performance work cannot silently change results.

.legacy_resample_spatial <- function(location_data, num_bins_x, num_bins_y,
                                     match_quantile) {
  original_order <- paste(location_data$x, location_data$y, sep = "_")
  rownames(location_data) <- original_order
  if (!all(c("x_bin", "y_bin") %in% colnames(location_data))) {
    location_data$x_bin <- cut(location_data$x, num_bins_x, labels = FALSE)
    location_data$y_bin <- cut(location_data$y, num_bins_y, labels = FALSE)
  }
  location_data$x_bin[is.na(location_data$x_bin)] <- 1
  location_data$y_bin[is.na(location_data$y_bin)] <- 1
  location_data$bin_id <- paste(location_data$x_bin, location_data$y_bin,
                                sep = "_")
  unique_bins <- unique(location_data$bin_id)
  mapping <- stats::setNames(sample(unique_bins), unique_bins)
  bin_coords <- unique(location_data[, c("bin_id", "x_bin", "y_bin")])
  out <- vector("list", length(unique_bins))

  for (i in seq_along(unique_bins)) {
    orig_bin <- unique_bins[i]
    target_bin <- mapping[[orig_bin]]
    orig_idx <- location_data$bin_id == orig_bin
    target_idx <- location_data$bin_id == target_bin
    orig_points <- location_data[orig_idx, , drop = FALSE]
    candidate_points <- location_data[target_idx, , drop = FALSE]
    n_points <- sum(orig_idx)

    if (nrow(candidate_points) < n_points) {
      neighbors <- CoPro:::.get_neighbor_bins(
        bin_coords, target_bin, num_bins_x, num_bins_y
      )
      for (neighbor in setdiff(neighbors, target_bin)) {
        neighbor_idx <- location_data$bin_id == neighbor
        if (any(neighbor_idx)) {
          candidate_points <- rbind(
            candidate_points,
            location_data[neighbor_idx, , drop = FALSE]
          )
        }
        if (nrow(candidate_points) >= n_points) break
      }
    }

    sampled <- if (match_quantile) {
      CoPro:::.match_by_quantile_position(
        orig_points, candidate_points, n_points
      )
    } else {
      sample(nrow(candidate_points), n_points,
             replace = nrow(candidate_points) < n_points)
    }
    orig_points$cell_ID <- candidate_points$cell_ID[sampled]
    orig_points$bin_id <- candidate_points$bin_id[sampled]
    out[[i]] <- orig_points
  }

  out <- do.call(rbind, out)
  rownames(out) <- paste(out$x, out$y, sep = "_")
  out[original_order, , drop = FALSE]
}


test_that("prepared spatial resampling preserves historical draws", {
  set.seed(710)
  location_data <- data.frame(
    x = runif(600), y = runif(600),
    cell_ID = paste0("cell_", seq_len(600))
  )

  for (match_quantile in c(FALSE, TRUE)) {
    set.seed(711)
    expected <- .legacy_resample_spatial(
      location_data, 9, 8, match_quantile
    )
    set.seed(711)
    actual <- resample_spatial(location_data, 9, 8, match_quantile)

    expect_identical(actual$cell_ID, expected$cell_ID)
    expect_identical(actual$bin_id, expected$bin_id)
    expect_identical(actual[, c("x", "y")], expected[, c("x", "y")])
  }
})


.legacy_permutation_ncorr <- function(object) {
  cts <- object@cellTypesOfInterest
  PCmats <- CoPro:::.getAllPCMats(object@pcaGlobal, object@scalePCs)
  pairs <- combn(cts, 2)
  sigma <- object@sigmaValueChoice
  nCC <- object@nCC
  out <- vector("list", object@nPermu)

  for (tt in seq_len(object@nPermu)) {
    values <- numeric(ncol(pairs) * nCC)
    for (pp in seq_len(ncol(pairs))) {
      ct1 <- pairs[1, pp]
      ct2 <- pairs[2, pp]
      K <- getKernelMatrix(object, sigma, ct1, ct2, verbose = FALSE)
      Rx <- tryCatch(getKernelMatrix(object, sigma, ct1, ct1, verbose = FALSE),
                     error = function(e) NULL)
      Ry <- tryCatch(getKernelMatrix(object, sigma, ct2, ct2, verbose = FALSE),
                     error = function(e) NULL)
      kernel_norm <- CoPro:::.whitenedFrobNorm(K, Rx, Ry)
      A <- PCmats[[ct1]][object@cellPermu[[ct1]][, tt], , drop = FALSE]
      B <- PCmats[[ct2]][object@cellPermu[[ct2]][, tt], , drop = FALSE]

      for (cc in seq_len(nCC)) {
        w1 <- object@skrCCAPermuOut[[tt]][[ct1]][, cc, drop = FALSE]
        w2 <- object@skrCCAPermuOut[[tt]][[ct2]][, cc, drop = FALSE]
        a <- A %*% w1
        b <- B %*% w2
        values[pp + (cc - 1L) * ncol(pairs)] <-
          as.numeric(crossprod(a, K %*% b)) /
          (sqrt(sum(a^2)) * sqrt(sum(b^2)) * kernel_norm)
      }
    }
    out[[tt]] <- values
  }
  names(out) <- paste0("permu_", seq_len(object@nPermu))
  out
}


test_that("batched permutation scoring matches direct pair/component loops", {
  skip_on_cran()
  obj <- create_test_copro_single(n_cells = 180, n_genes = 30,
                                  n_cell_types = 3, seed = 720)
  obj <- subsetData(obj, cellTypesOfInterest =
                      c("CellTypeA", "CellTypeB", "CellTypeC"))
  obj <- computePCA(obj, nPCA = 8, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = 0.1, method = "dense",
                             dropDistances = FALSE, verbose = FALSE)
  obj <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100, tol = 1e-6)
  ))
  obj <- suppressMessages(computeNormalizedCorrelation(obj))
  set.seed(721)
  obj <- suppressMessages(suppressWarnings(
    runSkrCCAPermu(obj, nPermu = 3, permu_method = "global",
                   maxIter = 100, tol = 1e-6, verbose = FALSE)
  ))

  expected <- .legacy_permutation_ncorr(obj)
  invisible(utils::capture.output(scored <- computeNormalizedCorrelationPermu(obj)))
  actual <- lapply(scored@normalizedCorrelationPermu,
                   function(x) x$normalizedCorrelation)
  expect_equal(actual, expected, tolerance = 1e-10)
})
