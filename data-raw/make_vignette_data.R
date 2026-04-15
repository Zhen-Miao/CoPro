# Prepare subsampled datasets for CoPro vignettes
#
# This script loads real datasets, subsamples them to manageable sizes
# (~5-30 MB), and saves them for upload to GitHub Releases via piggyback.
#
# After running this script, upload the files with:
#   library(piggyback)
#   pb_new_release("Zhen-Miao/CoPro", tag = "data-v1")
#   pb_upload("data-raw/vignette_data/copro_colon_d3.rds",
#             repo = "Zhen-Miao/CoPro", tag = "data-v1")
#   pb_upload("data-raw/vignette_data/copro_colon_d9.rds",
#             repo = "Zhen-Miao/CoPro", tag = "data-v1")
#   pb_upload("data-raw/vignette_data/copro_kidney.rds",
#             repo = "Zhen-Miao/CoPro", tag = "data-v1")
#   pb_upload("data-raw/vignette_data/copro_organoid.rds",
#             repo = "Zhen-Miao/CoPro", tag = "data-v1")
#   pb_upload("data-raw/vignette_data/copro_brain_merfish.rds",
#             repo = "Zhen-Miao/CoPro", tag = "data-v1")
#
# Prerequisites: Access to the original data files on Dropbox.

library(SeuratObject)
library(Matrix)

dir.create("data-raw/vignette_data", showWarnings = FALSE)
DATA_ROOT <- "/Users/zhenmiao/Library/CloudStorage/Dropbox/DIALOGUE_plus project/Data"

# ============================================================
# 1. Colon Day 3 -- cross-cell-type co-progression
#    Preprocessing matches: CoPro colon D3 all slides_by_region 2.R
# ============================================================
message("=== Preparing Colon D3 ===")

d3_ra <- readRDS(file.path(DATA_ROOT,
  "ra_filtered_DSS3_colon_all_slides.rds"))
d3_meta <- readRDS(file.path(DATA_ROOT,
  "meta_filtered_DSS3_colon_all_slides.rds"))
d3_meta <- as.data.frame(d3_meta)

# Gene filter: >= 0.8% cells with expression > 0 (matching original)
d3_bin <- d3_ra
d3_bin@x <- as.numeric(d3_bin@x > 0)
nz_genes <- Matrix::colSums(d3_bin)
d3_ra <- d3_ra[, nz_genes >= 0.008 * nrow(d3_bin)]

# Cell filter: >= 20 genes with expression > 0
nz_cells <- Matrix::rowSums(d3_bin)
d3_ra <- d3_ra[nz_cells >= 20, ]

# Match meta to data
rownames(d3_meta) <- d3_meta$Cell_ID
d3_meta <- d3_meta[rownames(d3_ra), ]
rm(d3_bin); gc()

# Cap at 99th percentile (matching original D3 script)
quant_99 <- quantile(d3_ra@x, probs = 0.99)
d3_ra@x[d3_ra@x > quant_99] <- quant_99
message("  99th percentile cap: ", quant_99)

# Subset to one slide and cell types
slide_id <- "092421_D3_m1_1_slice_1"
keep <- d3_meta$Slice_ID == slide_id
d3_ra <- d3_ra[keep, ]
d3_meta <- d3_meta[keep, ]

ct_keep <- d3_meta$Tier1 %in% c("Epithelial", "Fibroblast", "Immune")
d3_ra <- d3_ra[ct_keep, ]
d3_meta <- d3_meta[ct_keep, ]

# Convert to dense
d3_ra <- as.matrix(d3_ra)

# Create consistent cell IDs
cell_ids_d3 <- paste0("cell_", seq_len(nrow(d3_ra)))
rownames(d3_ra) <- cell_ids_d3
rownames(d3_meta) <- cell_ids_d3

message("  D3 slide: ", slide_id)
message("  Cells: ", nrow(d3_ra),
        " (Epi: ", sum(d3_meta$Tier1 == "Epithelial"),
        ", Fibro: ", sum(d3_meta$Tier1 == "Fibroblast"),
        ", Immune: ", sum(d3_meta$Tier1 == "Immune"), ")")
message("  Genes: ", ncol(d3_ra))

d3_location <- data.frame(
  x = d3_meta$x / 10,
  y = d3_meta$y / 10,
  row.names = cell_ids_d3
)

d3_data <- list(
  normalizedData = d3_ra,
  locationData = d3_location,
  metaData = d3_meta,
  cellTypes = d3_meta$Tier1,
  description = "Colon Day 3 organoid (slide 092421_D3_m1_1_slice_1). Log-normalized seqFISH data. Gene filter >= 0.8%, cell filter >= 20 genes, 99th percentile cap. Coordinates /10.",
  source = "Miao et al. CoPro manuscript"
)

saveRDS(d3_data, "data-raw/vignette_data/copro_colon_d3.rds", compress = "xz")
message("  Saved: ", round(file.size("data-raw/vignette_data/copro_colon_d3.rds") / 1e6, 1), " MB")

rm(d3_ra, d3_meta)
gc()

# ============================================================
# 2. Colon Day 9 -- multi-slide transfer
#    Preprocessing matches: CoPro colon D9 all slides_by_region 2.R
# ============================================================
message("\n=== Preparing Colon D9 ===")

d9_ra <- readRDS(file.path(DATA_ROOT,
  "ra_filtered_DSS9_colon_all_slides.rds"))
d9_meta <- readRDS(file.path(DATA_ROOT,
  "meta_filtered_DSS9_colon_all_slides.rds"))
d9_meta <- as.data.frame(d9_meta)

# Gene filter: >= 0.8% cells with expression > 0 (matching original)
d9_bin <- d9_ra
d9_bin@x <- as.numeric(d9_bin@x > 0)
nz_genes <- Matrix::colSums(d9_bin)
d9_ra <- d9_ra[, nz_genes >= 0.008 * nrow(d9_bin)]

# Cell filter: >= 20 genes with expression > 0
nz_cells <- Matrix::rowSums(d9_bin)
d9_ra <- d9_ra[nz_cells >= 20, ]

# Match meta to data
rownames(d9_meta) <- d9_meta$Cell_ID
d9_meta <- d9_meta[rownames(d9_ra), ]
rm(d9_bin); gc()

# Cap at 1.0 (matching original D9 script: daty@x[daty@x > 1] <- 1)
d9_ra@x[d9_ra@x > 1] <- 1
message("  Capped at 1.0")

# Select 3 slides with reasonable cell counts (>500 cells each)
slide_counts <- table(d9_meta$Slice_ID)
good_slides <- names(slide_counts[slide_counts > 500])
selected_slides <- head(good_slides, 3)
message("  Selected slides: ", paste(selected_slides, collapse = ", "))

keep <- d9_meta$Slice_ID %in% selected_slides
d9_ra <- d9_ra[keep, ]
d9_meta <- d9_meta[keep, ]

# Subset to cell types
ct_keep <- d9_meta$Tier1 %in% c("Epithelial", "Fibroblast", "Immune")
d9_ra <- d9_ra[ct_keep, ]
d9_meta <- d9_meta[ct_keep, ]

# Convert to dense
d9_ra <- as.matrix(d9_ra)

# Consistent cell IDs
cell_ids_d9 <- paste0("cell_", seq_len(nrow(d9_ra)))
rownames(d9_ra) <- cell_ids_d9
rownames(d9_meta) <- cell_ids_d9

message("  Total cells: ", nrow(d9_ra))
message("  Genes: ", ncol(d9_ra))

d9_location <- data.frame(
  x = d9_meta$x / 10,
  y = d9_meta$y / 10,
  row.names = cell_ids_d9
)

d9_data <- list(
  normalizedData = d9_ra,
  locationData = d9_location,
  metaData = d9_meta,
  cellTypes = d9_meta$Tier1,
  slideID = d9_meta$Slice_ID,
  selectedSlides = selected_slides,
  description = "Colon Day 9 organoid (3 slides). Log-normalized seqFISH data. Gene filter >= 0.8%, cell filter >= 20 genes, capped at 1.0. Coordinates /10.",
  source = "Miao et al. CoPro manuscript"
)

saveRDS(d9_data, "data-raw/vignette_data/copro_colon_d9.rds", compress = "xz")
message("  Saved: ", round(file.size("data-raw/vignette_data/copro_colon_d9.rds") / 1e6, 1), " MB")

rm(d9_ra, d9_meta, d9_ra_sub, d9_meta_sub)
gc()

# ============================================================
# 3. Kidney -- supervised/guided gradient
# ============================================================
message("\n=== Preparing Kidney ===")

kidney_obj <- readRDS(file.path(DATA_ROOT,
  "kidney_seqFISH/Ctrl2_object_with_xy.rds"))

# Extract data from Seurat object (assay is "Spatial", not "RNA")
kidney_ra <- as.matrix(SeuratObject::LayerData(kidney_obj[["Spatial"]], layer = "data"))
kidney_ra <- t(kidney_ra)  # genes x cells -> cells x genes
kidney_meta <- kidney_obj@meta.data

# Cell type column is "Celltype" (capital C), stored as factor
kidney_meta$celltype <- as.character(kidney_meta$Celltype)

# Define cell type groups (matching script 14/30)
tubular_types <- c("PTS1", "PTS2", "PTS3", "LOH-TL-C", "LOH-TL-JM",
                   "TAL_1", "TAL_2", "TAL_3", "DCT-CNT")
vasc_types <- c("Vasc_1", "Vasc_2", "Vasc_3")

# Assign grouped cell types
kidney_meta$grouped_celltype <- NA
kidney_meta$grouped_celltype[kidney_meta$celltype %in% tubular_types] <- "Tubular"
kidney_meta$grouped_celltype[kidney_meta$celltype %in% vasc_types] <- "Vascular"

# Keep only tubular and vascular cells
keep_kidney <- !is.na(kidney_meta$grouped_celltype)
kidney_ra <- kidney_ra[keep_kidney, ]
kidney_meta <- kidney_meta[keep_kidney, ]

message("  Cells: ", nrow(kidney_ra),
        " (Tubular: ", sum(kidney_meta$grouped_celltype == "Tubular"),
        ", Vascular: ", sum(kidney_meta$grouped_celltype == "Vascular"), ")")
message("  Genes: ", ncol(kidney_ra))

# Spatial coordinates (um -> mm)
kidney_location <- data.frame(
  x = kidney_meta$x_um / 1000,
  y = kidney_meta$y_um / 1000,
  row.names = rownames(kidney_ra)
)

# Segment ordering for supervised mode
segment_order <- c(PTS1 = 1, PTS2 = 2, PTS3 = 3,
                   "LOH-TL-C" = 4, "LOH-TL-JM" = 4,
                   TAL_1 = 5, TAL_2 = 5, TAL_3 = 5,
                   "DCT-CNT" = 6)

## ---- scRNA-seq data for transfer demonstration ----
## Include full-transcriptome expression matrices (sparse, filtered to
## >= 1% expressed genes), UMAP coordinates, and cell type annotations.
## Transfer uses only the shared genes, but the full transcriptome is
## needed for post-transfer regression to discover axis-associated genes
## beyond the spatial panel.
sc_dir <- "/Users/zhenmiao/Dropbox/kidney_sex_project/sex_RNA data/kidney_cello/"
seqfish_dir <- file.path(DATA_ROOT, "kidney_seqFISH")

library(Biobase)
library(Matrix)
eset_sc <- readRDS(file.path(sc_dir, "eset.rds"))
sc_meta_full <- eset_sc@phenoData@data
expr_full <- assayData(eset_sc)[["norm_exprs"]]  # genes x cells (sparse)

# Identify cell populations
vasc_cells <- rownames(sc_meta_full)[sc_meta_full$Lineages == "Vasculature"]
neph_cells <- rownames(sc_meta_full)[sc_meta_full$Lineages == "Nephron epithelium"]

# Full-transcriptome sparse matrices (genes x cells), filtered >= 1% expressed
sp_vasc <- expr_full[, vasc_cells]
pct_v <- Matrix::rowMeans(sp_vasc > 0)
sp_vasc <- sp_vasc[pct_v >= 0.01, ]
message("  Vascular scRNA-seq (full, filtered): ", nrow(sp_vasc), " genes x ",
        ncol(sp_vasc), " cells")

sp_neph <- expr_full[, neph_cells]
pct_n <- Matrix::rowMeans(sp_neph > 0)
sp_neph <- sp_neph[pct_n >= 0.01, ]
message("  Nephron scRNA-seq (full, filtered): ", nrow(sp_neph), " genes x ",
        ncol(sp_neph), " cells")

# Shared genes between seqFISH and scRNA-seq (for reference)
seqfish_genes <- colnames(kidney_ra)
shared_genes <- intersect(seqfish_genes, rownames(expr_full))
message("  Shared genes (seqFISH ∩ scRNA-seq): ", length(shared_genes))

# UMAP coordinates (pre-computed)
neph_umap <- readRDS(file.path(seqfish_dir, "sc_neph_umap_df.rds"))
vasc_umap <- readRDS(file.path(seqfish_dir, "sc_vasc_umap_df.rds"))

# Ontology_ID -> aggregated nephron label
ont_to_spatial <- c(
  "1" = NA, "2" = NA,
  "3" = "PTS1", "4" = "PTS1",
  "5" = "PTS2", "6" = "PTS2",
  "7" = "PTS3", "8" = "PTS3",
  "9A" = "LOH-TL-C", "9B" = "LOH-TL-JM",
  "10" = "LOH-TL-JM", "11" = "LOH-TL-JM", "12" = "LOH-TL-JM", "13" = "LOH-TL-JM",
  "14" = "TAL", "15" = "TAL", "16" = "TAL", "17" = "TAL",
  "18" = "DCT-CNT", "19" = "DCT-CNT"
)

# Add cell type annotations to UMAP dataframes
neph_mask <- sc_meta_full$Lineages == "Nephron epithelium"
ont_id <- sc_meta_full$Ontology_ID[neph_mask]
names(ont_id) <- rownames(sc_meta_full)[neph_mask]
neph_umap$ont_id <- ont_id[rownames(neph_umap)]
neph_umap$agg_label <- ont_to_spatial[as.character(neph_umap$ont_id)]

vasc_mask <- sc_meta_full$Lineages == "Vasculature"
vasc_cluster <- as.character(sc_meta_full$res.1[vasc_mask])
names(vasc_cluster) <- rownames(sc_meta_full)[vasc_mask]
vasc_umap$cluster <- vasc_cluster[rownames(vasc_umap)]

kidney_data <- list(
  normalizedData = kidney_ra,
  locationData = kidney_location,
  metaData = kidney_meta,
  cellTypes = kidney_meta$grouped_celltype,
  segmentOrder = segment_order,
  tubularTypes = tubular_types,
  vascTypes = vasc_types,
  ## scRNA-seq data for transfer workflow
  ## Sparse matrices (genes x cells) — full transcriptome, filtered >= 1% expressed.
  ## Transfer uses only shared genes; full transcriptome used for post-transfer regression.
  scRNA_vasc_expr = sp_vasc,
  scRNA_neph_expr = sp_neph,
  scRNA_neph_umap = neph_umap,
  scRNA_vasc_umap = vasc_umap,
  scRNA_ont_to_spatial = ont_to_spatial,
  description = "Kidney seqFISH control sample (Ctrl2). Tubular and Vascular cells. Includes segment ordering for supervised analysis and scRNA-seq full-transcriptome sparse expression matrices (genes x cells, filtered >= 1% expressed) with UMAP coordinates for transfer demonstration.",
  source = "Miao et al. CoPro manuscript; seqFISH data; scRNA-seq from Barry et al."
)

saveRDS(kidney_data, "data-raw/vignette_data/copro_kidney.rds", compress = "xz")
message("  Saved: ", round(file.size("data-raw/vignette_data/copro_kidney.rds") / 1e6, 1), " MB")

rm(kidney_obj, kidney_ra, kidney_meta)
gc()

# ============================================================
# 4. Organoid (72hr, single cell type)
# ============================================================
message("\n=== Preparing Organoid ===")

organoid_loc <- "/Users/zhenmiao/Library/CloudStorage/Dropbox/DIALOGUE_plus project/Raj_lab_data/72hr/roi1/"
if (dir.exists(organoid_loc)) {
  organoid_ra <- read.csv(file.path(organoid_loc,
    "output/cell_by_gene/cell_by_gene.csv"))
  organoid_ra <- organoid_ra[, 2:ncol(organoid_ra)]
  organoid_meta <- read.csv(file.path(organoid_loc,
    "output/attributes/cell_attributes.csv"))

  organoid_meta$label <- paste0("cell_", organoid_meta$label)
  rownames(organoid_ra) <- organoid_meta$label
  rownames(organoid_meta) <- organoid_meta$label

  # Cap at 95th percentile per gene (matching outlier_corrected script)
  filter_out <- function(x) {
    upper <- quantile(x, 0.95)
    x[x > upper] <- upper
    return(x)
  }
  organoid_ra <- apply(organoid_ra, 2, filter_out)

  # Filter cells >= 150 total counts
  organoid_meta <- organoid_meta[rowSums(organoid_ra) >= 150, ]
  organoid_ra <- organoid_ra[rowSums(organoid_ra) >= 150, ]
  message("  Cells after filtering (>= 150 counts): ", nrow(organoid_ra))

  # DESeq2 size-factor normalization + log1p
  counts_matrix <- t(organoid_ra) + 1
  mode(counts_matrix) <- "integer"
  condition <- factor(rep("condition", ncol(counts_matrix)))
  col_data <- data.frame(row.names = colnames(counts_matrix),
                         condition = condition)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_matrix, colData = col_data, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  normalized_counts <- t(organoid_ra) / DESeq2::sizeFactors(dds)
  organoid_ra <- t(log1p(normalized_counts))

  organoid_location <- data.frame(
    x = organoid_meta$center_x / 5000,
    y = organoid_meta$center_y / 5000,
    row.names = rownames(organoid_ra)
  )

  organoid_data <- list(
    normalizedData = organoid_ra,
    locationData = organoid_location,
    metaData = organoid_meta,
    cellTypes = rep("Epithelial", nrow(organoid_ra)),
    description = "72hr intestinal organoid culture (ROI-1). Single cell type (Epithelial). seqFISH data. 95th percentile cap per gene, cells >= 150 total counts, DESeq2 size-factor normalized + log1p. Coordinates /5000.",
    source = "Heyman et al. bioRxiv 2025.11.14.688372; Raj Lab organoid data"
  )

  saveRDS(organoid_data, "data-raw/vignette_data/copro_organoid.rds",
          compress = "xz")
  message("  Saved: ", round(file.size("data-raw/vignette_data/copro_organoid.rds") / 1e6, 1), " MB")
} else {
  message("  SKIPPED: organoid data not found at ", organoid_loc)
}

# ============================================================
# 5. Brain MERFISH (D1/D2 neurons)
# ============================================================
message("\n=== Preparing Brain MERFISH ===")

brain_data_path <- "/Users/zhenmiao/Dropbox/Zhuang-ABCA-1.054_1/Zhuang_ABCA_1.054_subset_data.rds"
brain_meta_path <- "/Users/zhenmiao/Dropbox/Zhuang-ABCA-1.054_1/Zhuang_ABCA_1.054_subset_metadata.rds"

if (file.exists(brain_data_path)) {
  brain_ra <- readRDS(brain_data_path)
  brain_meta <- readRDS(brain_meta_path)

  brain_types <- c("061 STR D1 Gaba", "062 STR D2 Gaba")
  brain_keep <- brain_meta$subclass %in% brain_types
  brain_ra <- brain_ra[brain_keep, ]
  brain_meta <- brain_meta[brain_keep, ]
  rownames(brain_meta) <- brain_meta$cell_label
  brain_meta$x <- as.numeric(brain_meta$x)
  brain_meta$y <- as.numeric(brain_meta$y)

  # Convert to dense if sparse
  if (inherits(brain_ra, "sparseMatrix")) {
    brain_ra <- as.matrix(brain_ra)
  }

  brain_location <- data.frame(
    x = brain_meta$x,
    y = brain_meta$y,
    row.names = rownames(brain_ra)
  )

  brain_out <- list(
    normalizedData = brain_ra,
    locationData = brain_location,
    metaData = brain_meta,
    cellTypes = brain_meta$subclass,
    description = "Brain MERFISH data. D1 and D2 GABAergic neurons from the striatum.",
    source = "Zhang et al. Nature 624, 343-354 (2023)"
  )

  saveRDS(brain_out, "data-raw/vignette_data/copro_brain_merfish.rds",
          compress = "xz")
  message("  Saved: ", round(file.size("data-raw/vignette_data/copro_brain_merfish.rds") / 1e6, 1), " MB")
} else {
  message("  SKIPPED: brain data not found at ", brain_data_path)
}

message("\n=== Done! ===")
message("Files in data-raw/vignette_data/:")
for (f in list.files("data-raw/vignette_data", full.names = TRUE)) {
  message("  ", basename(f), " (", round(file.size(f) / 1e6, 1), " MB)")
}
message("\nNext steps:")
message("  1. Review the files")
message("  2. Create a GitHub Release: pb_new_release('Zhen-Miao/CoPro', tag = 'data-v1')")
message("  3. Upload: pb_upload('data-raw/vignette_data/<file>.rds', repo = 'Zhen-Miao/CoPro', tag = 'data-v1')")
