#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CoPro)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: benchmark_large_multislide_float32.R ",
    "<prepare|float64|float32> <input.rds> [output-prefix] ",
    "[target-cells] [slides]"
  )
}

mode <- args[[1L]]
input_path <- args[[2L]]
output_prefix <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  sub("\\.rds$", "", input_path)
}
target_cells <- if (length(args) >= 4L) as.integer(args[[4L]]) else 200000L
n_slides <- if (length(args) >= 5L) as.integer(args[[5L]]) else 8L
sigmas <- c(0.005, 0.01, 0.02, 0.05, 0.1)
cell_types <- c("Epithelial", "Fibroblast", "Immune")
lower_limit <- 1e-7

if (mode == "prepare") {
  source_path <- "data-raw/vignette_data/copro_colon_d3.rds"
  source <- readRDS(source_path)
  keep <- source$cellTypes %in% cell_types
  source_coordinates <- as.matrix(source$locationData[keep, c("x", "y")])
  source_types <- source$cellTypes[keep]
  base_coordinates <- stats::setNames(lapply(cell_types, function(cell_type) {
    source_coordinates[source_types == cell_type, , drop = FALSE]
  }), cell_types)

  pairs <- utils::combn(cell_types, 2L)
  base_percentiles <- vapply(seq_len(ncol(pairs)), function(index) {
    coordinate_1 <- base_coordinates[[pairs[1L, index]]]
    coordinate_2 <- base_coordinates[[pairs[2L, index]]]
    CoPro:::.lowPercentileBlock(
      coordinate_1, coordinate_2,
      CoPro:::.pairPercentileProb(
        nrow(coordinate_1), nrow(coordinate_2)
      )
    )$percentile
  }, numeric(1))
  scaling_factor <- 0.01 / min(base_percentiles)
  support_radius <- CoPro:::.kernelSupportMultiplier(lower_limit) *
    max(sigmas) / scaling_factor
  all_coordinates <- do.call(rbind, base_coordinates)
  coordinate_min <- apply(all_coordinates, 2L, min)
  coordinate_span <- apply(all_coordinates, 2L, function(values) {
    diff(range(values))
  })
  stride <- coordinate_span + support_radius + 1

  cells_per_slide <- rep(target_cells %/% n_slides, n_slides)
  cells_per_slide[seq_len(target_cells %% n_slides)] <-
    cells_per_slide[seq_len(target_cells %% n_slides)] + 1L
  proportions <- vapply(base_coordinates, nrow, integer(1))
  proportions <- proportions / sum(proportions)

  coordinate_chunks <- list()
  type_chunks <- list()
  slide_chunks <- list()
  id_chunks <- list()
  chunk_index <- 1L
  for (slide_index in seq_len(n_slides)) {
    slide <- paste0("Slide", slide_index)
    type_counts <- floor(cells_per_slide[slide_index] * proportions)
    type_counts[length(type_counts)] <- cells_per_slide[slide_index] -
      sum(type_counts[-length(type_counts)])
    for (cell_type in cell_types) {
      base <- base_coordinates[[cell_type]]
      n_target <- type_counts[[cell_type]]
      linear_index <- seq_len(n_target) - 1L
      base_index <- linear_index %% nrow(base) + 1L
      tile_index <- linear_index %/% nrow(base)
      tile_column <- tile_index %% 5L
      tile_row <- tile_index %/% 5L
      coordinate <- base[base_index, , drop = FALSE]
      coordinate[, 1L] <- coordinate[, 1L] - coordinate_min[1L] +
        tile_column * stride[1L]
      coordinate[, 2L] <- coordinate[, 2L] - coordinate_min[2L] +
        tile_row * stride[2L]
      ids <- paste0(
        slide, "_", cell_type, "_", seq_len(n_target)
      )
      rownames(coordinate) <- ids
      coordinate_chunks[[chunk_index]] <- coordinate
      type_chunks[[chunk_index]] <- rep(cell_type, n_target)
      slide_chunks[[chunk_index]] <- rep(slide, n_target)
      id_chunks[[chunk_index]] <- ids
      chunk_index <- chunk_index + 1L
    }
  }

  coordinates <- do.call(rbind, coordinate_chunks)
  cell_type_vector <- unlist(type_chunks, use.names = FALSE)
  slide_vector <- unlist(slide_chunks, use.names = FALSE)
  cell_ids <- unlist(id_chunks, use.names = FALSE)
  normalized_data <- matrix(
    0, nrow = target_cells, ncol = 1L,
    dimnames = list(cell_ids, "placeholder")
  )
  metadata <- data.frame(
    cell_id = cell_ids,
    slideID = slide_vector,
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )
  saveRDS(
    list(
      normalizedData = normalized_data,
      locationData = as.data.frame(coordinates),
      metaData = metadata,
      cellTypes = cell_type_vector,
      slideID = slide_vector,
      target_cells = target_cells,
      n_slides = n_slides,
      source = source_path
    ),
    input_path,
    compress = FALSE
  )
  cat(
    "Prepared", target_cells, "cells across", n_slides,
    "slides at", input_path, "\n"
  )
  quit(save = "no", status = 0L)
}

input <- readRDS(input_path)
object <- newCoProMulti(
  normalizedData = input$normalizedData,
  locationData = input$locationData,
  metaData = input$metaData,
  cellTypes = input$cellTypes,
  slideID = input$slideID
)
object <- subsetData(object, cellTypesOfInterest = cell_types)
rm(input)
invisible(gc())

build_time <- system.time({
  if (mode == "float64") {
    object <- computeKernelMatrix(
      object, sigmaValues = sigmas, method = "sparse",
      distType = "Euclidean2D", lowerLimit = lower_limit,
      verbose = FALSE
    )
  } else if (mode == "float32") {
    object <- computeSparseKernelFloat32(
      object, sigmaValues = sigmas,
      distType = "Euclidean2D", lowerLimit = lower_limit,
      verbose = FALSE
    )
  } else {
    stop("Unknown mode: ", mode)
  }
})

kernel_signatures <- lapply(
  names(object@kernelMatrices),
  function(kernel_name) {
    kernel <- object@kernelMatrices[[kernel_name]]
    if (inherits(kernel, "CoProFloat32SparseMatrix")) {
      sums <- CoPro:::.float32KernelSums(kernel)
      nnz <- CoPro:::.float32KernelNnz(kernel)
    } else {
      values <- kernel@x
      nnz <- length(values)
      sums <- list(
        rowSums = Matrix::rowSums(kernel),
        sumSquares = sum(values^2)
      )
      if (inherits(kernel, "symmetricMatrix")) {
        nnz <- 2 * nnz
        sums$sumSquares <- 2 * sums$sumSquares
      }
    }
    data.frame(
      kernel = kernel_name,
      represented_nnz = nnz,
      value_sum = sum(sums$rowSums),
      value_sum_squares = sums$sumSquares,
      stringsAsFactors = FALSE
    )
  }
)
kernel_signatures <- do.call(rbind, kernel_signatures)
result <- data.frame(
  mode = mode,
  cells = length(object@cellTypesSub),
  slides = length(getSlideList(object)),
  kernels = length(object@kernelMatrices),
  build_elapsed_sec = unname(build_time[["elapsed"]]),
  kernel_bytes = as.numeric(object.size(object@kernelMatrices)),
  represented_nnz = sum(kernel_signatures$represented_nnz)
)
write.csv(
  result, paste0(output_prefix, "_", mode, ".csv"),
  row.names = FALSE
)
write.csv(
  kernel_signatures,
  paste0(output_prefix, "_", mode, "_signatures.csv"),
  row.names = FALSE
)
print(result)
