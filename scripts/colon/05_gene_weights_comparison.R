# =============================================================
# Gene weights comparison: PCA vs regression (Colon D0 example)
# =============================================================
#
# Compares PCA back-projection and regression-based gene weights
# on the colon D0 dataset. Generates bar plots and correlation
# scatter plots between the two methods.
#
# Inputs:
#   - Pre-computed CoPro object (colon D0)
#
# Outputs:
#   - Gene weight bar plots (PCA and regression, top 20 per cell type)
#   - Scatter plots comparing PCA vs regression weights
#   - Supplementary table CSV
#
# Dependencies: CoPro, ggplot2

library(CoPro)
library(ggplot2)

# --- Configuration ---
DATA_DIR <- "path/to/data"  # UPDATE THIS
OUT_DIR  <- "output/gene_weights_comparison"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cell_types <- c("Epithelial", "Fibroblast", "Immune")
sigma_opt <- 0.01

# --- Load pre-computed CoPro object ---
obj <- readRDS(file.path(DATA_DIR,
  "colon_healthy_all_all_slides_CoPro_object_region_3.rds"))

# --- Recompute gene scores with corrected formula ---
obj <- computeGeneAndCellScores(obj)            # PCA (corrected 1/sdev)
obj <- computeRegressionGeneScores(obj)         # Regression

# --- Compare methods per cell type ---
supp_table <- data.frame()

for (ct in cell_types) {
  key <- paste0("geneScores|sigma", sigma_opt, "|", ct)
  gs_pca <- obj@geneScores[[key]]
  gs_reg <- obj@geneScoresRegression[[key]]
  if (is.null(gs_pca) || is.null(gs_reg)) next

  # CC1 weights
  pca_cc1 <- gs_pca[, 1]
  reg_cc1 <- gs_reg[, 1]

  # Correlation
  common_genes <- intersect(names(pca_cc1), names(reg_cc1))
  r <- cor(pca_cc1[common_genes], reg_cc1[common_genes])
  message(ct, ": PCA-regression CC1 correlation = ", round(r, 3))

  # Scatter plot
  df <- data.frame(PCA = pca_cc1[common_genes],
                   Regression = reg_cc1[common_genes])
  pdf(file.path(OUT_DIR, paste0("PCA_vs_regression_", ct, ".pdf")),
      width = 4, height = 4)
  print(
    ggplot(df, aes(x = PCA, y = Regression)) +
      geom_point(size = 0.5, alpha = 0.3) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      ggtitle(paste(ct, "- r =", round(r, 3))) +
      theme_classic()
  )
  dev.off()

  # Bar plots (regression, top 20)
  top <- head(sort(abs(reg_cc1), decreasing = TRUE), 20)
  df_bar <- data.frame(
    gene = factor(names(top), levels = rev(names(top))),
    weight = reg_cc1[names(top)]
  )
  df_bar$dir <- ifelse(df_bar$weight > 0, "pos", "neg")

  pdf(file.path(OUT_DIR,
    paste0("top20_", ct, "_regression.pdf")),
    width = 3.5, height = 3.5)
  print(
    ggplot(df_bar, aes(x = gene, y = weight, fill = dir)) +
      geom_col() + coord_flip() +
      ggtitle(paste(ct, "CC1 - regression")) +
      theme_classic() + theme(legend.position = "none")
  )
  dev.off()

  # Supplementary table
  ct_table <- data.frame(
    Gene = common_genes,
    PCA_CC1 = pca_cc1[common_genes],
    Regression_CC1 = reg_cc1[common_genes]
  )
  names(ct_table)[2:3] <- paste0(ct, "_", names(ct_table)[2:3])
  if (nrow(supp_table) == 0) {
    supp_table <- ct_table
  } else {
    supp_table <- merge(supp_table, ct_table, by = "Gene", all = TRUE)
  }
}

write.csv(supp_table, file.path(OUT_DIR,
  "Supplementary_Table_Gene_Weights.csv"), row.names = FALSE)
