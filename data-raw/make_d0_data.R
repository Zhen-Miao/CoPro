#!/usr/bin/env Rscript
# Prepare D0 multi-slide vignette data
# Run from repo root: Rscript data-raw/make_d0_data.R

library(Matrix)

data_dir <- "/Users/zhenmiao/Library/CloudStorage/Dropbox/DIALOGUE_plus project/Data"
out_dir <- "data-raw/vignette_data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("=== Preparing Colon D0 (multi-slide) ===")

d0_ra <- readRDS(file.path(data_dir, "ra_filtered_healthy_colon_all_slides.rds"))
d0_meta <- readRDS(file.path(data_dir, "meta_filtered_healthy_colon_all_slides.rds"))
d0_meta <- as.data.frame(d0_meta)

d0_bin <- d0_ra
d0_bin@x <- as.numeric(d0_bin@x > 0)
nz_genes <- Matrix::colSums(d0_bin)
d0_ra <- d0_ra[, nz_genes >= 0.008 * nrow(d0_bin)]

nz_cells <- Matrix::rowSums(d0_bin)
d0_ra <- d0_ra[nz_cells >= 20, ]

rownames(d0_meta) <- d0_meta$Cell_ID
d0_meta <- d0_meta[rownames(d0_ra), ]
rm(d0_bin); gc()

quant_98 <- quantile(d0_ra@x, probs = 0.98)
d0_ra@x[d0_ra@x > quant_98] <- quant_98
message("  98th percentile cap: ", quant_98)

selected_slides <- c("062921_D0_m3a_1_slice_1",
                      "082421_D0_m6_1_slice_1",
                      "082421_D0_m7_1_slice_1")
slide_counts <- table(d0_meta$Slice_ID)
selected_slides <- intersect(selected_slides, names(slide_counts))
message("  Selected slides: ", paste(selected_slides, collapse = ", "))

keep <- d0_meta$Slice_ID %in% selected_slides
d0_ra <- d0_ra[keep, ]
d0_meta <- d0_meta[keep, ]

ct_keep <- d0_meta$Tier1 %in% c("Epithelial", "Fibroblast", "Immune")
d0_ra <- d0_ra[ct_keep, ]
d0_meta <- d0_meta[ct_keep, ]

d0_ra <- as.matrix(d0_ra)

cell_ids_d0 <- paste0("cell_", seq_len(nrow(d0_ra)))
rownames(d0_ra) <- cell_ids_d0
rownames(d0_meta) <- cell_ids_d0

message("  Total cells: ", nrow(d0_ra))
for (sl in selected_slides) {
  n_sl <- sum(d0_meta$Slice_ID == sl)
  message("    ", sl, ": ", n_sl, " cells")
}
message("  Genes: ", ncol(d0_ra))

d0_location <- data.frame(
  x = d0_meta$x / 10,
  y = d0_meta$y / 10,
  row.names = cell_ids_d0
)

d0_data <- list(
  normalizedData = d0_ra,
  locationData = d0_location,
  metaData = d0_meta,
  cellTypes = d0_meta$Tier1,
  slideID = d0_meta$Slice_ID,
  selectedSlides = selected_slides,
  description = "Colon Day 0 healthy organoid (3 slides from different regions). Log-normalized seqFISH data. Gene filter >= 0.8%, cell filter >= 20 genes, 98th percentile cap. Coordinates /10.",
  source = "Miao et al. CoPro manuscript"
)

out_path <- file.path(out_dir, "copro_colon_d0_multi.rds")
saveRDS(d0_data, out_path, compress = "xz")
message("  Saved: ", round(file.size(out_path) / 1e6, 1), " MB")
message("Done.")
