# Generate a small synthetic dataset for examples and testing
#
# This script creates a tiny CoPro-ready dataset with a spatial pattern
# embedded in the expression data. It is bundled in inst/extdata/ so
# that man page examples can run without downloading external data.
#
# Run this script from the package root:
#   Rscript data-raw/make_toy_data.R

set.seed(42)

n_cells <- 200
n_genes <- 80

# Two cell types: "Epithelial" and "Fibroblast"
cell_types <- rep(c("Epithelial", "Fibroblast"), each = n_cells / 2)
cell_ids <- paste0("cell_", seq_len(n_cells))
gene_names <- paste0("Gene", seq_len(n_genes))

# Spatial coordinates on a unit square
x <- runif(n_cells, 0, 1)
y <- runif(n_cells, 0, 1)

# Create a smooth spatial gradient (left-to-right)
spatial_score <- scale(x + 0.3 * y + rnorm(n_cells, sd = 0.1))[, 1]

# Generate expression data with a spatial pattern in the first 15 genes
normalizedData <- matrix(
  rnorm(n_cells * n_genes, mean = 2, sd = 0.5),
  nrow = n_cells, ncol = n_genes
)
for (g in 1:15) {
  sign_g <- ifelse(g %% 2 == 1, 1, -1)
  normalizedData[, g] <- normalizedData[, g] + spatial_score * sign_g * 0.8
}
normalizedData[normalizedData < 0] <- 0
rownames(normalizedData) <- cell_ids
colnames(normalizedData) <- gene_names

locationData <- data.frame(x = x, y = y, row.names = cell_ids)

metaData <- data.frame(
  cell_id = cell_ids,
  cell_type = cell_types,
  row.names = cell_ids
)

toy_data <- list(
  normalizedData = normalizedData,
  locationData = locationData,
  metaData = metaData,
  cellTypes = cell_types
)

saveRDS(toy_data, file = "inst/extdata/toy_copro_data.rds",
        compress = "xz")

message("Saved inst/extdata/toy_copro_data.rds (",
        round(file.size("inst/extdata/toy_copro_data.rds") / 1024, 1),
        " KB)")
