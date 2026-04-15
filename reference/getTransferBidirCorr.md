# Compute Bidirectional Correlation from Transferred Cell Scores

Given a target object and a list of transferred cell scores (per cell
type), compute the bidirectional correlation for each pair of cell types
using the provided kernel matrices at a selected sigma value.

## Usage

``` r
getTransferBidirCorr(
  tar_obj,
  transfer_cell_scores,
  sigma_choice,
  calculationMode = NULL,
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
  K_row_sum_cutoff = 0.005,
  K_col_sum_cutoff = 0.005,
  sigma_choice_tar = NULL,
  verbose = TRUE
)
```

## Arguments

- tar_obj:

  A `CoProSingle` or `CoProMulti` object containing kernel matrices and
  metadata needed for alignment.

- transfer_cell_scores:

  A named list of matrices, with one entry per cell type (names must be
  the cell type names). Each matrix should be cells-by-CCs, where rows
  are cell IDs and columns are `CC_1`, `CC_2`, ..., as returned by
  `getTransferCellScores(agg_cell_type = FALSE)`.

- sigma_choice:

  Numeric scalar specifying the sigma value of the kernel to use.

- calculationMode:

  For `CoProMulti` objects only, either "perSlide" or "aggregate".
  Ignored for `CoProSingle`. Default "perSlide" if `tar_obj` is
  multi-slide.

- normalize_K:

  Character; method for normalizing the kernel matrix, one of
  "row_or_col", "sinkhorn_knopp", or "none". Default "row_or_col".

- filter_kernel:

  Logical; whether to filter the kernel matrix. Default TRUE.

- K_row_sum_cutoff:

  Numeric; cutoff for row sums when normalizing kernel matrix. Default
  5e-3.

- K_col_sum_cutoff:

  Numeric; cutoff for column sums when normalizing kernel matrix.
  Default 5e-3.

- sigma_choice_tar:

  Numeric; sigma value for target object kernel matrices. If NULL
  (default), uses sigma_choice. Not recommended for general use.

- verbose:

  Logical; whether to print progress messages.

## Value

A list with one element named `paste0("sigma_", sigma_choice)`, whose
value is a data.frame of results. For single-slide objects, the
data.frame has columns `sigmaValue`, `cellType1`, `cellType2`,
`CC_index`, `bidirCorrelation`. For multi-slide objects in `perSlide`
mode, the data.frame additionally includes `slideID`. For `aggregate`
mode, the correlation column is named `aggregateCorrelation`.

## Details

This function mirrors the correlation calculation used in
[`computeBidirCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeBidirCorrelation.md)
but operates on precomputed cell scores (e.g., obtained from
`getTransferCellScores(agg_cell_type = FALSE)`).

The bidirectional correlation is computed as the mean of two
correlations: cor(t(K) %*% A_w1, B_w2) and cor(A_w1, K %*% B_w2), where
A_w1 and B_w2 are the cell score vectors (for the same CC index) for
cell types A and B respectively, and K is the kernel matrix between the
two cell types at the chosen sigma.

## Examples

``` r
# Assuming `tar_obj` is prepared and `trans_scores` was computed with
# getTransferCellScores(..., agg_cell_type = FALSE)
# res <- getTransferBidirCorr(tar_obj, trans_scores, sigma_choice = 2.0)
```
