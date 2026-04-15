# =============================================================
# Colon Day 3: Cross-cell-type co-progression analysis
# =============================================================
#
# Demonstrates CoPro with three cell types on a single colon D3
# organoid slide. Produces normalized correlation plots, in situ
# cell score maps, and gene weight bar charts.
#
# Inputs:
#   - ra_filtered_DSS3_colon_all_slides.rds  (expression matrix)
#   - meta_filtered_DSS3_colon_all_slides.rds (metadata)
#
# Outputs:
#   - CoPro object (RDS)
#   - Normalized correlation PDF
#   - In situ cell score PDFs (CC1, CC2)
#   - Gene weight bar plots
#
# Dependencies: CoPro, ggplot2, viridis

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/data"  # UPDATE THIS
OUT_DIR  <- "output/colon_D3"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

n_pca <- 40
cell_types <- c("Epithelial", "Fibroblast", "Immune")
sigma_choice <- c(0.002, 0.005, 0.01, 0.05, 0.1, 0.2)
region_id <- "092421_D3_m1_1_slice_1"

# --- Load data ---
meta <- readRDS(file.path(DATA_DIR, "meta_filtered_DSS3_colon_all_slides.rds"))
ra   <- readRDS(file.path(DATA_DIR, "ra_filtered_DSS3_colon_all_slides.rds"))

# --- Subset to region ---
keep <- meta$Slice_ID == region_id
meta <- meta[keep, ]
ra   <- ra[keep, ]

# --- Preprocess ---
# Filter genes: keep if >= 0.8% of cells have expression > 0
gene_frac <- colMeans(ra > 0)
ra <- ra[, gene_frac >= 0.008]

# Filter cells: keep if >= 20 genes expressed
cell_frac <- rowSums(ra > 0)
keep_cells <- cell_frac >= 20
ra   <- ra[keep_cells, ]
meta <- meta[keep_cells, ]

# Cap expression at 99th percentile
quant_99 <- apply(ra, 2, quantile, 0.99)
for (g in seq_len(ncol(ra))) {
  ra[, g] <- pmin(ra[, g], quant_99[g])
}

# Prepare spatial coordinates
location_data <- data.frame(
  x = meta$center_x / 10,
  y = meta$center_y / 10,
  row.names = rownames(ra)
)

# --- Create CoPro object ---
obj <- newCoProSingle(
  normalizedData = ra,
  locationData = location_data,
  metaData = meta,
  cellTypes = meta$Celltype
)
obj <- subsetData(obj, cellTypesOfInterest = cell_types)

# --- Run pipeline ---
obj <- computePCA(obj, nPCA = n_pca, center = TRUE, scale. = TRUE)
obj <- computeDistance(obj, distType = "Euclidean2D")
obj <- computeKernelMatrix(obj, sigmaValues = sigma_choice)
obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500, nCC = 4)
obj <- computeNormalizedCorrelation(obj)
obj <- computeGeneAndCellScores(obj)
obj <- computeRegressionGeneScores(obj)

# Save object
saveRDS(obj, file.path(OUT_DIR,
  paste0("colon_D3_", region_id, "_CoPro_object.rds")))

# --- Plot normalized correlation ---
ncorr <- getNormCorr(obj)
pdf(file.path(OUT_DIR, "norm_corr_D3.pdf"), width = 15, height = 10)
print(
  ggplot(ncorr, aes(x = sigmaValues, y = normalizedCorrelation)) +
    geom_point() + geom_line() +
    facet_wrap(~ ct12 + CC_index, scales = "free_y") +
    xlab("Sigma") + ylab("Normalized Correlation") +
    ggtitle("Colon D3: Normalized correlation") +
    theme_minimal()
)
dev.off()

# --- In situ cell scores ---
sigma_opt <- 0.01
for (cc_ind in 1:2) {
  cs <- getCellScoresInSitu(obj, sigmaValueChoice = sigma_opt,
                             ccIndex = cc_ind)
  for (ct in cell_types) {
    pdf(file.path(OUT_DIR,
      paste0("insitu_CC", cc_ind, "_", ct, ".pdf")),
      width = 6, height = 5)
    print(
      ggplot(cs[cs$cellType == ct, ]) +
        geom_point(aes(x = x, y = y, color = cellScores), size = 0.8) +
        scale_color_viridis_c() +
        coord_fixed() +
        ggtitle(paste("CC", cc_ind, "-", ct)) +
        theme_classic()
    )
    dev.off()
  }
}

# --- Gene weight bar plots (regression, top 20) ---
for (ct in cell_types) {
  key <- paste0("geneScores|sigma", sigma_opt, "|", ct)
  gs <- obj@geneScoresRegression[[key]]
  if (is.null(gs)) next

  gs_cc1 <- gs[, 1]
  top <- head(sort(abs(gs_cc1), decreasing = TRUE), 20)
  df <- data.frame(
    gene = factor(names(top), levels = rev(names(top))),
    weight = gs_cc1[names(top)]
  )
  df$dir <- ifelse(df$weight > 0, "pos", "neg")

  pdf(file.path(OUT_DIR,
    paste0("gene_weight_CC1_", ct, "_regression.pdf")),
    width = 3.5, height = 3.5)
  print(
    ggplot(df, aes(x = gene, y = weight, fill = dir)) +
      geom_col() + coord_flip() +
      ggtitle(paste(ct, "- CC1 gene weights")) +
      theme_classic() + theme(legend.position = "none")
  )
  dev.off()
}
