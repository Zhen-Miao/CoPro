# CoPro object of spatial transcriptomics data

CoPro object of spatial transcriptomics data

## Slots

- `normalizedData`:

  A `matrix` object to store normalized data.

- `normalizedDataSub`:

  A `matrix` object to store the subset of the normalized data. The
  subset only contain relevant cell types of interest specified by the
  user.

- `locationData`:

  A `data.frame` object to store the location. It should either contain
  two columns named by "x" and "y", or three columns named by "x", "y",
  and "z". No other names allowed

- `locationDataSub`:

  A `data.frame` object to store the subset of the location data.

- `metaData`:

  A `data.frame` object to store metadata for each cell.

- `metaDataSub`:

  A `data.frame` object to store the subset of the meta data.

- `cellTypes`:

  A `vector` object with elements being character. It should match the
  number of cells in the data matrix and each represents a cell type
  label of a cell.

- `cellTypesSub`:

  A `vector` object with elements being character. It stores the subset
  of the cell type labels.

- `cellTypesOfInterest`:

  A `vector` object with elements being character. Specifies the cell
  types of interest. It will be used to subset the dataset

- `pcaResults`:

  A `list` object storing PCA results after integration. Recommended
  structure: `list(slideID = list(cellType = pc_matrix))`.

- `pcaGlobal`:

  A `list` object storing PCA results for each cell type.

- `distances`:

  A `list` object to store the pairwise distances between any two cell
  types of interest.

- `geneList`:

  A `vector` object with elements being character. To store the gene
  names.

- `kernelMatrices`:

  A `list` object. To store the kernel matrix generated from the
  distance matrices.

- `sigmaValues`:

  A `vector` object with elements being numeric. To store a set of sigma
  values used for generating the kernel matrix.

- `nPCA`:

  A single numeric value. Number of PCs to retain for downstream
  analyses.

- `nCC`:

  A single numeric value. Number of canonical components to retain for
  downstream analyses.

- `scalePCs`:

  A `logical` value. Whether to scale each PC before computing skrCCA

- `skrCCAOut`:

  A `list` object. Output from the skrCCA.

- `skrCCAPermuOut`:

  A `list` object. Output from the skrCCA after permutation. This helps
  establish the null distribution

- `cellPermu`:

  A `list` object that stores the cell permutation labels

- `nPermu`:

  A `numeric` value specifying the number of permutations conducted.

- `cellScores`:

  A `matrix` object. Cell scores for each cell type.

- `geneScores`:

  A `matrix` object. Gene scores for each cell type.

- `geneScoresRegression`:

  A `list` object. Regression-based gene scores. For each gene, the
  regression coefficient of gene expression on the cell score is used as
  the gene weight. This avoids collinearity issues present in the PCA
  back-projection approach stored in `geneScores`.

- `geneScoresTest`:

  A `list` object. Tested gene scores

- `normalizedCorrelation`:

  A `list` object. Normalized correlation values for each sigma value.

- `bidirCorrelation`:

  A `list` object. Bidirectional correlation values for each sigma
  value.

- `normalizedCorrelationPermu`:

  A `list` object. Normalized correlation values for each sigma value
  after permutation

- `bidirCorrelationPermu`:

  A `list` object. Bidirectional correlation values for each sigma value
  after permutation

- `sigmaValueChoice`:

  A `numeric` value. The optimal sigma squared based on the median
  normalized correlation value.
