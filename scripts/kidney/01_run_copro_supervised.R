# =============================================================
# Kidney: Supervised CoPro along the nephron axis
# =============================================================
#
# Uses known nephron segment ordering (proximal-to-distal) to
# guide the tubular cell type's axis. CoPro then finds the
# vascular gene program that co-varies with this axis.
#
# Inputs:
#   - Ctrl1_object_with_xy.rds, Ctrl2_object_with_xy.rds,
#     Ctrl3_object_with_xy.rds (Seurat objects with spatial coords)
#
# Outputs:
#   - CoPro objects per control replicate (RDS)
#   - Normalized correlation, in situ, and gene weight plots
#
# Dependencies: CoPro, SeuratObject, ggplot2, viridis

library(CoPro)
library(SeuratObject)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/kidney_seqFISH"  # UPDATE THIS
OUT_DIR  <- "output/kidney"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_pca <- 15  # fewer PCs for targeted spatial panels (~1300 genes)
sigma_choice <- c(0.04, 0.08, 0.1, 0.15)
sigma_opt <- 0.1

# Cell type grouping
tubular_types <- c("PTS1", "PTS2", "PTS3", "LOH-TL-C", "LOH-TL-JM",
                   "TAL_1", "TAL_2", "TAL_3", "DCT-CNT")
vasc_types <- c("Vasc_1", "Vasc_2", "Vasc_3")

# Nephron segment ordering for supervised mode
segment_order <- c(PTS1 = 1, PTS2 = 2, PTS3 = 3,
                   "LOH-TL-C" = 4, "LOH-TL-JM" = 4,
                   TAL_1 = 5, TAL_2 = 5, TAL_3 = 5,
                   "DCT-CNT" = 6)

ctrl_names <- c("Ctrl1", "Ctrl2", "Ctrl3")

# --- Process each control replicate ---
for (ctrl in ctrl_names) {
  message("=== Processing ", ctrl, " ===")

  # Load Seurat object
  seu <- readRDS(file.path(DATA_DIR,
    paste0(ctrl, "_object_with_xy.rds")))
  ra <- t(as.matrix(seu[["RNA"]]$data))  # cells x genes
  meta <- seu@meta.data

  # Assign grouped cell types
  meta$grouped_ct <- NA
  meta$grouped_ct[meta$celltype %in% tubular_types] <- "Tubular"
  meta$grouped_ct[meta$celltype %in% vasc_types] <- "Vascular"
  keep <- !is.na(meta$grouped_ct)
  ra   <- ra[keep, ]
  meta <- meta[keep, ]

  # Spatial coordinates (um -> mm)
  location_data <- data.frame(
    x = meta$x_um / 1000,
    y = meta$y_um / 1000,
    row.names = rownames(ra)
  )

  # Create CoPro object
  obj <- newCoProSingle(
    normalizedData = ra,
    locationData = location_data,
    metaData = meta,
    cellTypes = meta$grouped_ct
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("Tubular", "Vascular"))

  # Pipeline
  obj <- computePCA(obj, nPCA = n_pca, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D")
  obj <- computeKernelMatrix(obj, sigmaValues = sigma_choice)
  obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500, nCC = 2)
  obj <- computeNormalizedCorrelation(obj)
  obj <- computeGeneAndCellScores(obj)
  obj <- computeRegressionGeneScores(obj)

  # Save
  saveRDS(obj, file.path(OUT_DIR,
    paste0("kidney_", ctrl, "_CoPro_object.rds")))

  # In situ plot
  cs <- getCellScoresInSitu(obj, sigmaValueChoice = sigma_opt)
  pdf(file.path(OUT_DIR,
    paste0("insitu_CC1_", ctrl, ".pdf")), width = 8, height = 5)
  print(
    ggplot(cs) +
      geom_point(aes(x = x, y = y, color = cellScores), size = 0.8) +
      scale_color_viridis_c() +
      facet_wrap(~ cellType) +
      coord_fixed() +
      ggtitle(paste(ctrl, "- CC1 cell scores")) +
      theme_classic()
  )
  dev.off()

  message("  Tubular cells: ", sum(meta$grouped_ct == "Tubular"))
  message("  Vascular cells: ", sum(meta$grouped_ct == "Vascular"))
}
