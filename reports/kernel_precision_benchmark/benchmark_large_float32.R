#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CoPro)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: benchmark_large_float32.R ",
    "<prepare|float64|float32> <input.rds> [output-prefix] [threads]"
  )
}

mode <- args[[1L]]
input_path <- args[[2L]]
output_prefix <- if (length(args) >= 3L) args[[3L]] else
  sub("\\.rds$", "", input_path)
n_threads <- if (length(args) >= 4L) as.integer(args[[4L]]) else 8L

sigmas <- c(0.005, 0.01, 0.02, 0.05, 0.1)
lower_limit <- 1e-7
upper_quantile <- 0.85
cell_types <- c("Epithelial", "Fibroblast", "Immune")

kernel_name <- function(sigma, cell_type_1, cell_type_2) {
  paste(
    "kernel", paste0("sigma", sigma),
    cell_type_1, cell_type_2, sep = "|"
  )
}

median_time <- function(fun, reps = 3L) {
  times <- numeric(reps)
  for (index in seq_len(reps)) {
    gc()
    times[index] <- unname(system.time(invisible(fun()))[["elapsed"]])
  }
  c(median = median(times), minimum = min(times), maximum = max(times))
}

if (mode == "prepare") {
  target_cells <- if (length(args) >= 3L) as.integer(args[[3L]]) else 200000L
  source_path <- "data-raw/vignette_data/copro_colon_d3.rds"
  data <- readRDS(source_path)
  object <- newCoProSingle(
    normalizedData = data$normalizedData,
    locationData = data$locationData,
    metaData = data$metaData,
    cellTypes = data$cellTypes
  )
  object <- subsetData(object, cellTypesOfInterest = cell_types)
  object <- computePCA(
    object, nPCA = 40L, center = TRUE, scale. = TRUE
  )
  base_coordinates <- stats::setNames(
    lapply(cell_types, function(cell_type) {
      CoPro:::.getCoordinateMatrix(
        object, cell_type, "Euclidean2D", 1, 1, 1
      )
    }),
    cell_types
  )
  base_x <- CoPro:::.prepareDataMatrices(
    object, is_multi = FALSE, scalePCs = TRUE, cts = cell_types
  )$PCmats

  pairs <- utils::combn(cell_types, 2L)
  base_percentiles <- vapply(seq_len(ncol(pairs)), function(index) {
    cell_type_1 <- pairs[1L, index]
    cell_type_2 <- pairs[2L, index]
    coordinate_1 <- base_coordinates[[cell_type_1]]
    coordinate_2 <- base_coordinates[[cell_type_2]]
    CoPro:::.lowPercentileBlock(
      coordinate_1,
      coordinate_2,
      CoPro:::.pairPercentileProb(
        nrow(coordinate_1), nrow(coordinate_2)
      )
    )$percentile
  }, numeric(1))
  scaling_factor <- 0.01 / min(base_percentiles)
  max_radius <- CoPro:::.kernelSupportMultiplier(lower_limit) *
    max(sigmas) / scaling_factor

  all_base_coordinates <- do.call(rbind, base_coordinates)
  coordinate_min <- apply(all_base_coordinates, 2L, min)
  coordinate_span <- apply(all_base_coordinates, 2L, function(values) {
    diff(range(values))
  })
  stride <- coordinate_span + max_radius + 1

  base_counts <- vapply(base_coordinates, nrow, integer(1))
  target_counts <- floor(
    target_cells * base_counts / sum(base_counts)
  )
  target_counts[length(target_counts)] <- target_cells -
    sum(target_counts[-length(target_counts)])

  tiled_coordinates <- stats::setNames(vector("list", length(cell_types)),
                                       cell_types)
  tiled_x <- stats::setNames(vector("list", length(cell_types)), cell_types)
  for (cell_type in cell_types) {
    n_base <- base_counts[[cell_type]]
    n_target <- target_counts[[cell_type]]
    linear_index <- seq_len(n_target) - 1L
    base_index <- linear_index %% n_base + 1L
    tile_index <- linear_index %/% n_base
    tile_column <- tile_index %% 5L
    tile_row <- tile_index %/% 5L

    coordinate <- base_coordinates[[cell_type]][base_index, , drop = FALSE]
    coordinate[, 1L] <- coordinate[, 1L] - coordinate_min[1L] +
      tile_column * stride[1L]
    coordinate[, 2L] <- coordinate[, 2L] - coordinate_min[2L] +
      tile_row * stride[2L]
    rownames(coordinate) <- NULL
    tiled_coordinates[[cell_type]] <- coordinate

    x <- base_x[[cell_type]][base_index, , drop = FALSE]
    rownames(x) <- NULL
    tiled_x[[cell_type]] <- unclass(x)
  }

  blocks <- lapply(seq_len(ncol(pairs)), function(index) {
    list(
      cellType1 = pairs[1L, index],
      cellType2 = pairs[2L, index],
      percentile = base_percentiles[index]
    )
  })
  block_sizes <- vapply(blocks, function(block) {
    as.numeric(nrow(tiled_coordinates[[block$cellType1]])) *
      as.numeric(nrow(tiled_coordinates[[block$cellType2]]))
  }, numeric(1))
  blocks <- blocks[order(block_sizes, decreasing = TRUE)]

  saveRDS(
    list(
      target_cells = target_cells,
      coordinates = tiled_coordinates,
      x = tiled_x,
      blocks = blocks,
      scaling_factor = scaling_factor,
      sigmas = sigmas,
      lower_limit = lower_limit,
      upper_quantile = upper_quantile,
      source = source_path,
      construction = paste(
        "Colon D3 tiled with gaps larger than the largest support radius;"
      )
    ),
    input_path,
    compress = FALSE
  )
  cat(
    "Prepared", target_cells, "cells at", input_path,
    "size", format(file.info(input_path)$size, big.mark = ","), "bytes\n"
  )
  quit(save = "no", status = 0L)
}

input <- readRDS(input_path)
coordinates <- input$coordinates
x_matrices <- input$x
blocks <- input$blocks
scaling_factor <- input$scaling_factor

if (mode == "float64") {
  gc()
  build_timing <- system.time({
    block_triplets <- lapply(blocks, function(block) {
      CoPro:::.buildBlockTriplets(
        coordinates[[block$cellType1]],
        coordinates[[block$cellType2]],
        block$percentile,
        scaling_factor,
        max(sigmas),
        lower_limit,
        truncateLowDist = TRUE
      )
    })
    kernels <- list()
    for (sigma in sigmas) {
      for (block_index in seq_along(blocks)) {
        block <- blocks[[block_index]]
        raw_kernel <- CoPro:::.sparseKernelFromTriplets(
          block_triplets[[block_index]],
          sigma, lower_limit, symmetric = FALSE
        )
        kernels[[kernel_name(
          sigma, block$cellType1, block$cellType2
        )]] <- CoPro:::.processSparseKernelMatrix(
          raw_kernel, lower_limit, upper_quantile,
          normalizeKernel = FALSE,
          rowNormalizeKernel = FALSE,
          colNormalizeKernel = FALSE
        )
      }
    }
  })
  temporary_bytes <- as.numeric(object.size(block_triplets))
  rm(block_triplets)
  gc()

  selected_name <- kernel_name(0.1, "Epithelial", "Fibroblast")
  selected_kernel <- kernels[[selected_name]]
  xky_timing <- median_time(function() {
    crossprod(
      x_matrices$Epithelial,
      selected_kernel %*% x_matrices$Fibroblast
    )
  })
  y <- crossprod(
    x_matrices$Epithelial,
    selected_kernel %*% x_matrices$Fibroblast
  )
  result <- data.frame(
    mode = mode,
    cells = input$target_cells,
    threads = 1L,
    build_elapsed_sec = unname(build_timing[["elapsed"]]),
    kernel_bytes = as.numeric(object.size(kernels)),
    temporary_neighbor_bytes = temporary_bytes,
    selected_nnz = length(selected_kernel@x),
    xky_median_sec = xky_timing[["median"]],
    xky_min_sec = xky_timing[["minimum"]],
    xky_max_sec = xky_timing[["maximum"]]
  )
} else if (mode == "float32") {
  gc()
  peak_temporary_bytes <- 0
  build_timing <- system.time({
    kernels <- list()
    for (block in blocks) {
      built <- CoPro:::float32_csr_gaussian_kernels_cpp(
        coordinates[[block$cellType1]],
        coordinates[[block$cellType2]],
        sigmas,
        block$percentile,
        scaling_factor,
        lower_limit,
        upper_quantile,
        TRUE
      )
      peak_temporary_bytes <- max(
        peak_temporary_bytes, built$temporary_bytes
      )
      for (sigma_index in seq_along(sigmas)) {
        name <- kernel_name(
          sigmas[sigma_index],
          block$cellType1,
          block$cellType2
        )
        kernels[[name]] <- CoPro:::.newFloat32SparseKernel(
          built$kernels[[sigma_index]], NULL
        )
      }
      rm(built)
      gc(FALSE)
    }
  })

  selected_name <- kernel_name(0.1, "Epithelial", "Fibroblast")
  selected_kernel <- kernels[[selected_name]]
  thread_values <- unique(c(1L, 2L, 4L, n_threads))
  timing_rows <- lapply(thread_values, function(threads) {
    timing <- median_time(function() {
      CoPro:::.kernelXKY(
        x_matrices$Epithelial,
        selected_kernel,
        x_matrices$Fibroblast,
        n_threads = threads
      )
    })
    data.frame(
      threads = threads,
      median = timing[["median"]],
      minimum = timing[["minimum"]],
      maximum = timing[["maximum"]]
    )
  })
  timing_table <- do.call(rbind, timing_rows)
  chosen <- timing_table[timing_table$threads == n_threads, , drop = FALSE]
  y <- CoPro:::.kernelXKY(
    x_matrices$Epithelial,
    selected_kernel,
    x_matrices$Fibroblast,
    n_threads = n_threads
  )
  result <- data.frame(
    mode = mode,
    cells = input$target_cells,
    threads = n_threads,
    build_elapsed_sec = unname(build_timing[["elapsed"]]),
    kernel_bytes = as.numeric(object.size(kernels)),
    temporary_neighbor_bytes = peak_temporary_bytes,
    selected_nnz = CoPro:::.float32KernelNnz(selected_kernel),
    xky_median_sec = chosen$median,
    xky_min_sec = chosen$minimum,
    xky_max_sec = chosen$maximum
  )
  write.csv(
    timing_table,
    paste0(output_prefix, "_float32_thread_scaling.csv"),
    row.names = FALSE
  )
} else {
  stop("Unknown mode: ", mode)
}

write.csv(result, paste0(output_prefix, "_", mode, ".csv"), row.names = FALSE)
saveRDS(y, paste0(output_prefix, "_", mode, "_y.rds"), compress = FALSE)
print(result)
