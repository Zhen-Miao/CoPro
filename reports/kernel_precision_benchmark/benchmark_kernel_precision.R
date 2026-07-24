#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Rcpp)
  devtools::load_all(".", quiet = TRUE)
})

set.seed(20260723)

out_dir <- "reports/kernel_precision_benchmark"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
Rcpp::sourceCpp(file.path(out_dir, "encoded_sparse_xky.cpp"))

data_path <- "data-raw/vignette_data/copro_colon_d3.rds"
cell_types <- c("Epithelial", "Fibroblast", "Immune")
sigma_values <- c(0.005, 0.01, 0.02, 0.05, 0.1)
n_pca <- 40L
n_cc <- 4L
lower_limit <- 1e-7
upper_quantile <- 0.85
normalize_target <- 0.01
speed_sigmas <- c(0.01, 0.05, 0.1)
speed_reps <- 3L

elapsed <- function(expr) {
  gc()
  timing <- system.time(value <- force(expr))
  list(value = value, seconds = unname(timing[["elapsed"]]))
}

median_elapsed <- function(fun, reps = speed_reps) {
  times <- numeric(reps)
  for (rep_index in seq_len(reps)) {
    gc()
    times[rep_index] <- unname(system.time(invisible(fun()))[["elapsed"]])
  }
  c(
    median_sec = median(times),
    min_sec = min(times),
    max_sec = max(times)
  )
}

kernel_name <- function(sigma, cell_type_1, cell_type_2) {
  paste("kernel", paste0("sigma", sigma), cell_type_1, cell_type_2, sep = "|")
}

encode_bins <- function(values, n_bins = 100L) {
  stopifnot(n_bins <= 256L, length(values) > 0L)
  boundaries <- as.numeric(stats::quantile(
    values,
    probs = seq_len(n_bins - 1L) / n_bins,
    names = FALSE,
    type = 7
  ))
  codebook <- as.numeric(stats::quantile(
    values,
    probs = (seq_len(n_bins) - 0.5) / n_bins,
    names = FALSE,
    type = 7
  ))
  codes <- findInterval(values, boundaries)
  list(values = as.raw(codes), codebook = codebook)
}

decode_bins <- function(encoded) {
  encoded$codebook[as.integer(encoded$values) + 1L]
}

encode_sparse <- function(kernel, encoding) {
  value_encoding <- switch(
    encoding,
    fp32 = list(values = cpp_pack_float32(kernel@x)),
    fp16 = list(values = cpp_pack_half(kernel@x)),
    bins100 = encode_bins(kernel@x, 100L),
    stop("Unknown encoding: ", encoding)
  )
  c(
    list(
      p = kernel@p,
      i = kernel@i,
      dims = kernel@Dim,
      encoding = encoding
    ),
    value_encoding
  )
}

decode_sparse_values <- function(encoded) {
  switch(
    encoded$encoding,
    fp32 = cpp_unpack_float32(encoded$values),
    fp16 = cpp_unpack_half(encoded$values),
    bins100 = decode_bins(encoded),
    stop("Unknown encoding: ", encoded$encoding)
  )
}

decode_sparse <- function(encoded) {
  Matrix::sparseMatrix(
    i = encoded$i + 1L,
    p = encoded$p,
    x = decode_sparse_values(encoded),
    dims = encoded$dims
  )
}

encoded_xky <- function(encoded, x_left, x_right,
                        float_accumulation = FALSE) {
  switch(
    encoded$encoding,
    fp32 = cpp_xky_float32(
      encoded$p, encoded$i, encoded$values, encoded$dims,
      x_left, x_right, float_accumulation
    ),
    fp16 = cpp_xky_half(
      encoded$p, encoded$i, encoded$values, encoded$dims,
      x_left, x_right, float_accumulation
    ),
    bins100 = cpp_xky_bins(
      encoded$p, encoded$i, encoded$values, encoded$codebook, encoded$dims,
      x_left, x_right, float_accumulation
    ),
    stop("Unknown encoding: ", encoded$encoding)
  )
}

encode_kernel_list <- function(kernel_list, encoding) {
  stats::setNames(
    lapply(kernel_list, encode_sparse, encoding = encoding),
    names(kernel_list)
  )
}

y_from_kernel_list <- function(kernel_list, x_list) {
  result <- stats::setNames(vector("list", length(cell_types)), cell_types)
  for (cell_type in cell_types) {
    result[[cell_type]] <- stats::setNames(
      vector("list", length(cell_types)), cell_types
    )
  }
  pairs <- utils::combn(cell_types, 2L)
  for (pair_index in seq_len(ncol(pairs))) {
    cell_type_1 <- pairs[1L, pair_index]
    cell_type_2 <- pairs[2L, pair_index]
    kernel <- kernel_list[[kernel_name(
      attr(kernel_list, "sigma"), cell_type_1, cell_type_2
    )]]
    y <- crossprod(
      x_list[[cell_type_1]],
      kernel %*% x_list[[cell_type_2]]
    )
    result[[cell_type_1]][[cell_type_2]] <- as.matrix(y)
    result[[cell_type_2]][[cell_type_1]] <- t(y)
  }
  result
}

y_from_encoded_list <- function(encoded_list, x_list, sigma,
                                float_accumulation = FALSE) {
  result <- stats::setNames(vector("list", length(cell_types)), cell_types)
  for (cell_type in cell_types) {
    result[[cell_type]] <- stats::setNames(
      vector("list", length(cell_types)), cell_types
    )
  }
  pairs <- utils::combn(cell_types, 2L)
  for (pair_index in seq_len(ncol(pairs))) {
    cell_type_1 <- pairs[1L, pair_index]
    cell_type_2 <- pairs[2L, pair_index]
    encoded <- encoded_list[[kernel_name(
      sigma, cell_type_1, cell_type_2
    )]]
    y <- encoded_xky(
      encoded,
      x_list[[cell_type_1]],
      x_list[[cell_type_2]],
      float_accumulation
    )
    result[[cell_type_1]][[cell_type_2]] <- y
    result[[cell_type_2]][[cell_type_1]] <- t(y)
  }
  result
}

solve_from_y <- function(y_original, x_list, seed = 20260723L) {
  set.seed(seed)
  feature_counts <- stats::setNames(
    vapply(x_list[cell_types], ncol, integer(1)), cell_types
  )
  weights <- CoPro:::initialize_weights_svd(x_list, cell_types)
  invisible(capture.output({
    weights <- CoPro:::bilinear_w_from_Y_resi(
      w_list_new = weights,
      Y_resi = y_original,
      n_features = feature_counts,
      max_iter = 500L,
      tol = 1e-5,
      step_size = 1,
      sdev2_list = NULL
    )
  }))
  for (cell_type in cell_types) {
    if (!is.matrix(weights[[cell_type]])) {
      weights[[cell_type]] <- matrix(weights[[cell_type]], ncol = 1L)
    }
  }

  y_residual <- y_original
  if (n_cc > 1L) {
    for (component in seq_len(n_cc - 1L)) {
      y_residual <- CoPro:::apply_deflation(
        y_residual,
        weights,
        component,
        cell_types,
        sdev2_list = NULL,
        deflation = "projection"
      )
      next_weights <- CoPro:::initialize_next_component(
        y_residual, cell_types
      )
      invisible(capture.output({
        next_weights <- CoPro:::bilinear_w_from_Y_resi(
          w_list_new = next_weights,
          Y_resi = y_residual,
          n_features = feature_counts,
          max_iter = 500L,
          tol = 1e-5,
          step_size = 1,
          sdev2_list = NULL
        )
      }))
      for (cell_type in cell_types) {
        weights[[cell_type]] <- cbind(
          weights[[cell_type]], next_weights[[cell_type]]
        )
      }
    }
  }
  weights
}

cell_scores_from_weights <- function(weights, x_list) {
  stats::setNames(
    lapply(cell_types, function(cell_type) {
      x_list[[cell_type]] %*% weights[[cell_type]]
    }),
    cell_types
  )
}

y_metric_rows <- function(variant, sigma, baseline, candidate) {
  pairs <- utils::combn(cell_types, 2L)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(pair_index) {
    cell_type_1 <- pairs[1L, pair_index]
    cell_type_2 <- pairs[2L, pair_index]
    reference <- baseline[[cell_type_1]][[cell_type_2]]
    alternative <- candidate[[cell_type_1]][[cell_type_2]]
    delta <- alternative - reference
    data.frame(
      variant = variant,
      sigma = sigma,
      pair = paste(cell_type_1, cell_type_2, sep = "|"),
      relative_frobenius_error =
        sqrt(sum(delta^2)) / sqrt(sum(reference^2)),
      max_absolute_error = max(abs(delta)),
      stringsAsFactors = FALSE
    )
  }))
}

score_metric_rows <- function(variant, sigma, baseline_weights,
                              candidate_weights, baseline_scores,
                              candidate_scores) {
  rows <- list()
  row_index <- 0L
  for (cell_type in cell_types) {
    for (component in seq_len(n_cc)) {
      reference <- baseline_scores[[cell_type]][, component]
      alternative <- candidate_scores[[cell_type]][, component]
      sign_use <- if (sum(reference * alternative) < 0) -1 else 1
      aligned <- alternative * sign_use
      delta <- aligned - reference
      weight_reference <- baseline_weights[[cell_type]][, component]
      weight_alternative <- candidate_weights[[cell_type]][, component] * sign_use
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        variant = variant,
        sigma = sigma,
        cell_type = cell_type,
        component = component,
        cell_score_correlation = stats::cor(reference, aligned),
        cell_score_nrmse = sqrt(mean(delta^2)) / stats::sd(reference),
        cell_score_max_absolute_error = max(abs(delta)),
        weight_cosine = sum(weight_reference * weight_alternative) /
          sqrt(sum(weight_reference^2) * sum(weight_alternative^2)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

kernel_metric_rows <- function(variant, baseline_kernels,
                               candidate_provider) {
  rows <- vector("list", length(baseline_kernels))
  for (index in seq_along(baseline_kernels)) {
    name <- names(baseline_kernels)[index]
    reference <- baseline_kernels[[name]]
    alternative <- candidate_provider(name)
    delta <- alternative - reference
    rows[[index]] <- data.frame(
      variant = variant,
      matrix = name,
      nnz_baseline = length(reference@x),
      nnz_candidate = length(alternative@x),
      relative_frobenius_error =
        sqrt(sum(delta@x^2)) / sqrt(sum(reference@x^2)),
      max_absolute_error =
        if (length(delta@x) == 0L) 0 else max(abs(delta@x)),
      stringsAsFactors = FALSE
    )
    rm(alternative, delta)
  }
  do.call(rbind, rows)
}

make_y_by_sigma_matrix <- function(kernel_list, x_list) {
  stats::setNames(lapply(sigma_values, function(sigma) {
    sigma_kernels <- kernel_list
    attr(sigma_kernels, "sigma") <- sigma
    y_from_kernel_list(sigma_kernels, x_list)
  }), as.character(sigma_values))
}

make_y_by_sigma_encoded <- function(encoded_list, x_list,
                                    float_accumulation) {
  stats::setNames(lapply(sigma_values, function(sigma) {
    y_from_encoded_list(
      encoded_list, x_list, sigma, float_accumulation
    )
  }), as.character(sigma_values))
}

make_kernel_list_from_distances <- function(blocks, distance_values) {
  kernels <- list()
  for (sigma in sigma_values) {
    for (block_index in seq_along(blocks)) {
      block <- blocks[[block_index]]
      triplets <- block$triplets
      triplets$dscaled <- distance_values[[block_index]]
      raw_kernel <- CoPro:::.sparseKernelFromTriplets(
        triplets, sigma, lower_limit, symmetric = FALSE
      )
      kernels[[kernel_name(sigma, block$cell_type_1, block$cell_type_2)]] <-
        CoPro:::.processSparseKernelMatrix(
          raw_kernel,
          lower_limit,
          upper_quantile,
          normalizeKernel = FALSE,
          rowNormalizeKernel = FALSE,
          colNormalizeKernel = FALSE
        )
    }
  }
  kernels
}

storage_rows <- list()
xky_speed_rows <- list()
kernel_build_rows <- list()
kernel_metric_all <- list()
y_metric_all <- list()
score_metric_all <- list()

cat("Loading colon D3 data\n")
data <- readRDS(data_path)
object <- newCoProSingle(
  normalizedData = data$normalizedData,
  locationData = data$locationData,
  metaData = data$metaData,
  cellTypes = data$cellTypes
)
object <- subsetData(object, cellTypesOfInterest = cell_types)

cat("Computing 40-PC input matrices\n")
pca_timing <- elapsed(
  computePCA(object, nPCA = n_pca, center = TRUE, scale. = TRUE)
)
object <- pca_timing$value

cat("Computing baseline fixed-radius sparse kernels\n")
kernel_timing <- elapsed(
  computeKernelMatrix(
    object,
    sigmaValues = sigma_values,
    method = "sparse",
    lowerLimit = lower_limit,
    upperQuantile = upper_quantile,
    normalizeKernel = FALSE,
    dropDistances = TRUE,
    verbose = FALSE
  )
)
object <- kernel_timing$value
baseline_kernels <- object@kernelMatrices
x_list <- CoPro:::.prepareDataMatrices(
  object = object,
  is_multi = FALSE,
  scalePCs = TRUE,
  cts = cell_types
)$PCmats

cat("Reconstructing the fixed-radius temporary distance triplets\n")
coordinate_list <- stats::setNames(lapply(cell_types, function(cell_type) {
  CoPro:::.getCoordinateMatrix(
    object, cell_type, "Euclidean2D", 1, 1, 1
  )
}), cell_types)
pairs <- utils::combn(cell_types, 2L)
percentiles <- vapply(seq_len(ncol(pairs)), function(pair_index) {
  cell_type_1 <- pairs[1L, pair_index]
  cell_type_2 <- pairs[2L, pair_index]
  coordinate_1 <- coordinate_list[[cell_type_1]]
  coordinate_2 <- coordinate_list[[cell_type_2]]
  probability <- CoPro:::.pairPercentileProb(
    nrow(coordinate_1), nrow(coordinate_2)
  )
  CoPro:::.lowPercentileBlock(
    coordinate_1, coordinate_2, probability
  )$percentile
}, numeric(1))
distance_scale <- normalize_target / min(percentiles)

blocks <- lapply(seq_len(ncol(pairs)), function(pair_index) {
  cell_type_1 <- pairs[1L, pair_index]
  cell_type_2 <- pairs[2L, pair_index]
  list(
    cell_type_1 = cell_type_1,
    cell_type_2 = cell_type_2,
    triplets = CoPro:::.buildBlockTriplets(
      coordinate_list[[cell_type_1]],
      coordinate_list[[cell_type_2]],
      percentiles[pair_index],
      distance_scale,
      max(sigma_values),
      lower_limit,
      truncateLowDist = TRUE
    )
  )
})
distance_double <- lapply(blocks, function(block) block$triplets$dscaled)

baseline_rebuilt_timing <- elapsed(
  make_kernel_list_from_distances(blocks, distance_double)
)
baseline_rebuilt <- baseline_rebuilt_timing$value
rebuild_error <- kernel_metric_rows(
  "double_distance_rebuild_check",
  baseline_kernels,
  function(name) baseline_rebuilt[[name]]
)
if (max(rebuild_error$max_absolute_error) > 1e-12 ||
    any(rebuild_error$nnz_baseline != rebuild_error$nnz_candidate)) {
  stop("Rebuilt double-distance kernels do not match the pipeline baseline")
}
rm(baseline_rebuilt)
gc()

storage_rows[[length(storage_rows) + 1L]] <- data.frame(
  object = "kernel_matrices",
  encoding = "double_Matrix_current",
  bytes = as.numeric(object.size(baseline_kernels)),
  persistent_in_current_pipeline = TRUE,
  stringsAsFactors = FALSE
)

distance_triplet_double <- lapply(blocks, function(block) {
  list(
    i = block$triplets$i,
    j = block$triplets$j,
    values = block$triplets$dscaled
  )
})
storage_rows[[length(storage_rows) + 1L]] <- data.frame(
  object = "temporary_distance_triplets",
  encoding = "double_current_temporary",
  bytes = as.numeric(object.size(distance_triplet_double)),
  persistent_in_current_pipeline = FALSE,
  stringsAsFactors = FALSE
)

cat("Computing baseline Y operators and final cell scores\n")
baseline_y <- make_y_by_sigma_matrix(baseline_kernels, x_list)
baseline_weights <- stats::setNames(lapply(sigma_values, function(sigma) {
  solve_from_y(baseline_y[[as.character(sigma)]], x_list)
}), as.character(sigma_values))
baseline_scores <- stats::setNames(lapply(sigma_values, function(sigma) {
  cell_scores_from_weights(
    baseline_weights[[as.character(sigma)]], x_list
  )
}), as.character(sigma_values))

for (sigma in speed_sigmas) {
  name <- kernel_name(sigma, "Epithelial", "Fibroblast")
  kernel <- baseline_kernels[[name]]
  timing <- median_elapsed(function() {
    crossprod(x_list$Epithelial, kernel %*% x_list$Fibroblast)
  })
  xky_speed_rows[[length(xky_speed_rows) + 1L]] <- data.frame(
    sigma = sigma,
    nnz = length(kernel@x),
    method = "Matrix_double_f64",
    median_sec = timing[["median_sec"]],
    min_sec = timing[["min_sec"]],
    max_sec = timing[["max_sec"]],
    stringsAsFactors = FALSE
  )
}

record_variant <- function(variant, candidate_y) {
  for (sigma in sigma_values) {
    key <- as.character(sigma)
    candidate_weights <- solve_from_y(candidate_y[[key]], x_list)
    candidate_scores <- cell_scores_from_weights(candidate_weights, x_list)
    y_metric_all[[length(y_metric_all) + 1L]] <<- y_metric_rows(
      variant, sigma, baseline_y[[key]], candidate_y[[key]]
    )
    score_metric_all[[length(score_metric_all) + 1L]] <<- score_metric_rows(
      variant,
      sigma,
      baseline_weights[[key]],
      candidate_weights,
      baseline_scores[[key]],
      candidate_scores
    )
  }
}

benchmark_kernel_encoding <- function(encoding) {
  cat("Encoding kernel values as", encoding, "\n")
  encode_timing <- elapsed(encode_kernel_list(baseline_kernels, encoding))
  encoded <- encode_timing$value
  storage_rows[[length(storage_rows) + 1L]] <<- data.frame(
    object = "kernel_matrices",
    encoding = encoding,
    bytes = as.numeric(object.size(encoded)),
    persistent_in_current_pipeline = FALSE,
    stringsAsFactors = FALSE
  )
  kernel_metric_all[[length(kernel_metric_all) + 1L]] <<-
    kernel_metric_rows(
      paste0("kernel_", encoding),
      baseline_kernels,
      function(name) decode_sparse(encoded[[name]])
    )

  y_double <- make_y_by_sigma_encoded(
    encoded, x_list, float_accumulation = FALSE
  )
  record_variant(
    paste0("kernel_", encoding, "_f64_accumulation"), y_double
  )
  rm(y_double)

  y_float <- make_y_by_sigma_encoded(
    encoded, x_list, float_accumulation = TRUE
  )
  record_variant(
    paste0("kernel_", encoding, "_f32_accumulation"), y_float
  )
  rm(y_float)

  for (sigma in speed_sigmas) {
    name <- kernel_name(sigma, "Epithelial", "Fibroblast")
    encoded_kernel <- encoded[[name]]
    reference_kernel <- baseline_kernels[[name]]

    timing_direct_double <- median_elapsed(function() {
      encoded_xky(
        encoded_kernel,
        x_list$Epithelial,
        x_list$Fibroblast,
        float_accumulation = FALSE
      )
    })
    timing_direct_float <- median_elapsed(function() {
      encoded_xky(
        encoded_kernel,
        x_list$Epithelial,
        x_list$Fibroblast,
        float_accumulation = TRUE
      )
    })
    timing_decode_matrix <- median_elapsed(function() {
      decoded_kernel <- reference_kernel
      decoded_kernel@x <- decode_sparse_values(encoded_kernel)
      crossprod(
        x_list$Epithelial,
        decoded_kernel %*% x_list$Fibroblast
      )
    })

    speed_specs <- list(
      direct_f64 = timing_direct_double,
      direct_f32 = timing_direct_float,
      decode_then_Matrix_f64 = timing_decode_matrix
    )
    for (method_name in names(speed_specs)) {
      timing <- speed_specs[[method_name]]
      xky_speed_rows[[length(xky_speed_rows) + 1L]] <<- data.frame(
        sigma = sigma,
        nnz = length(reference_kernel@x),
        method = paste(encoding, method_name, sep = "_"),
        median_sec = timing[["median_sec"]],
        min_sec = timing[["min_sec"]],
        max_sec = timing[["max_sec"]],
        stringsAsFactors = FALSE
      )
    }
  }

  data.frame(
    encoding = encoding,
    encode_all_kernel_matrices_sec = encode_timing$seconds,
    stringsAsFactors = FALSE
  )
}

encoding_time_rows <- do.call(
  rbind,
  lapply(c("fp32", "fp16", "bins100"), benchmark_kernel_encoding)
)
gc()

benchmark_distance_encoding <- function(encoding) {
  cat("Encoding temporary distances as", encoding, "\n")
  encode_one <- switch(
    encoding,
    fp32 = function(values) list(values = cpp_pack_float32(values)),
    fp16 = function(values) list(values = cpp_pack_half(values)),
    bins100 = function(values) encode_bins(values, 100L),
    stop("Unknown distance encoding")
  )
  decode_one <- switch(
    encoding,
    fp32 = function(value) cpp_unpack_float32(value$values),
    fp16 = function(value) cpp_unpack_half(value$values),
    bins100 = decode_bins,
    stop("Unknown distance encoding")
  )
  encoded_distance <- lapply(distance_double, encode_one)
  encoded_triplets <- lapply(seq_along(blocks), function(block_index) {
    list(
      i = blocks[[block_index]]$triplets$i,
      j = blocks[[block_index]]$triplets$j,
      values = encoded_distance[[block_index]]$values,
      codebook = encoded_distance[[block_index]]$codebook
    )
  })
  storage_rows[[length(storage_rows) + 1L]] <<- data.frame(
    object = "temporary_distance_triplets",
    encoding = encoding,
    bytes = as.numeric(object.size(encoded_triplets)),
    persistent_in_current_pipeline = FALSE,
    stringsAsFactors = FALSE
  )

  kernel_build_timing <- elapsed({
    decoded_distance <- lapply(encoded_distance, decode_one)
    make_kernel_list_from_distances(blocks, decoded_distance)
  })
  candidate_kernels <- kernel_build_timing$value
  kernel_build_rows[[length(kernel_build_rows) + 1L]] <<- data.frame(
    distance_encoding = encoding,
    decode_and_build_all_kernels_sec = kernel_build_timing$seconds,
    stringsAsFactors = FALSE
  )
  kernel_metric_all[[length(kernel_metric_all) + 1L]] <<-
    kernel_metric_rows(
      paste0("distance_", encoding),
      baseline_kernels,
      function(name) candidate_kernels[[name]]
    )
  candidate_y <- make_y_by_sigma_matrix(candidate_kernels, x_list)
  record_variant(
    paste0("distance_", encoding, "_kernel_f64"), candidate_y
  )

  if (encoding == "fp16") {
    cat("Evaluating combined fp16 distance and fp16 kernel storage\n")
    encoded_candidate_kernels <- encode_kernel_list(
      candidate_kernels, "fp16"
    )
    combined_y <- make_y_by_sigma_encoded(
      encoded_candidate_kernels,
      x_list,
      float_accumulation = FALSE
    )
    record_variant(
      "distance_fp16_plus_kernel_fp16_f64_accumulation",
      combined_y
    )
    rm(encoded_candidate_kernels, combined_y)
  }
  rm(candidate_kernels, candidate_y, encoded_distance, encoded_triplets)
  gc()
  invisible(NULL)
}

kernel_build_rows[[length(kernel_build_rows) + 1L]] <- data.frame(
  distance_encoding = "double_triplets",
  decode_and_build_all_kernels_sec = baseline_rebuilt_timing$seconds,
  stringsAsFactors = FALSE
)
invisible(lapply(
  c("fp32", "fp16", "bins100"), benchmark_distance_encoding
))

storage <- do.call(rbind, storage_rows)
xky_speed <- do.call(rbind, xky_speed_rows)
kernel_build_speed <- do.call(rbind, kernel_build_rows)
kernel_accuracy <- do.call(rbind, kernel_metric_all)
y_accuracy <- do.call(rbind, y_metric_all)
cell_score_accuracy <- do.call(rbind, score_metric_all)

score_summary <- do.call(rbind, lapply(
  split(cell_score_accuracy, cell_score_accuracy$variant),
  function(values) {
    data.frame(
      variant = values$variant[1L],
      min_cell_score_correlation =
        min(values$cell_score_correlation),
      max_cell_score_nrmse = max(values$cell_score_nrmse),
      max_cell_score_absolute_error =
        max(values$cell_score_max_absolute_error),
      min_weight_cosine = min(values$weight_cosine),
      stringsAsFactors = FALSE
    )
  }
))
y_summary <- do.call(rbind, lapply(
  split(y_accuracy, y_accuracy$variant),
  function(values) {
    data.frame(
      variant = values$variant[1L],
      max_y_relative_frobenius_error =
        max(values$relative_frobenius_error),
      max_y_absolute_error = max(values$max_absolute_error),
      stringsAsFactors = FALSE
    )
  }
))
accuracy_summary <- merge(y_summary, score_summary, by = "variant")

write.csv(storage, file.path(out_dir, "storage_memory.csv"), row.names = FALSE)
write.csv(xky_speed, file.path(out_dir, "xky_speed.csv"), row.names = FALSE)
write.csv(
  kernel_build_speed,
  file.path(out_dir, "distance_kernel_build_speed.csv"),
  row.names = FALSE
)
write.csv(
  encoding_time_rows,
  file.path(out_dir, "kernel_encoding_speed.csv"),
  row.names = FALSE
)
write.csv(
  kernel_accuracy,
  file.path(out_dir, "kernel_accuracy.csv"),
  row.names = FALSE
)
write.csv(y_accuracy, file.path(out_dir, "y_accuracy.csv"), row.names = FALSE)
write.csv(
  cell_score_accuracy,
  file.path(out_dir, "cell_score_accuracy.csv"),
  row.names = FALSE
)
write.csv(
  accuracy_summary,
  file.path(out_dir, "accuracy_summary.csv"),
  row.names = FALSE
)

run_metadata <- data.frame(
  item = c(
    "date", "R", "Matrix", "platform", "data", "cells", "cell_types",
    "n_pca", "n_cc", "sigmas", "lower_limit", "upper_quantile",
    "pca_elapsed_sec", "fixed_radius_kernel_elapsed_sec",
    "fixed_radius_distance_matrices_persisted",
    "sparse_kernel_class", "sparse_kernel_value_R_type"
  ),
  value = c(
    as.character(Sys.Date()),
    R.version.string,
    as.character(packageVersion("Matrix")),
    R.version$platform,
    data_path,
    nrow(object@metaDataSub),
    paste(cell_types, collapse = ","),
    n_pca,
    n_cc,
    paste(sigma_values, collapse = ","),
    lower_limit,
    upper_quantile,
    pca_timing$seconds,
    kernel_timing$seconds,
    length(object@distances),
    class(baseline_kernels[[1L]])[1L],
    typeof(baseline_kernels[[1L]]@x)
  ),
  stringsAsFactors = FALSE
)
write.csv(run_metadata, file.path(out_dir, "run_metadata.csv"), row.names = FALSE)

cat("\nAccuracy summary\n")
print(accuracy_summary)
cat("\nStorage\n")
print(transform(storage, mebibytes = bytes / 1024^2))
cat("\nX' K X speed\n")
print(xky_speed)
cat("\nDistance decode and kernel build speed\n")
print(kernel_build_speed)
cat("\nFinished. Results are in", out_dir, "\n")
