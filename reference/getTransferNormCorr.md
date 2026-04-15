# Compute Normalized Correlation from Transferred Cell Scores

Given a target object and a list of transferred cell scores (per cell
type), compute the normalized correlation for each pair of cell types
using the provided kernel matrices at a selected sigma value.

## Usage

``` r
getTransferNormCorr(
  tar_obj,
  transfer_cell_scores,
  sigma_choice,
  tol = 1e-04,
  calculationMode = NULL,
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

- tol:

  Numeric tolerance passed to the truncated SVD for spectral norm
  estimation (default 1e-4).

- calculationMode:

  For `CoProMulti` objects only, either "perSlide" or "aggregate".
  Ignored for `CoProSingle`. Default "perSlide" if `tar_obj` is
  multi-slide.

- sigma_choice_tar:

  Numeric; sigma value for target object kernel matrices. If NULL
  (default), uses sigma_choice. Not recommended for general use.

- verbose:

  Logical; whether to print progress messages.

## Value

A list with one element named `paste0("sigma_", sigma_choice)`, whose
value is a data.frame of results. For single-slide objects, the
data.frame has columns `sigmaValue`, `cellType1`, `cellType2`,
`CC_index`, `normalizedCorrelation`. For multi-slide objects in
`perSlide` mode, the data.frame additionally includes `slideID`. For
`aggregate` mode, the correlation column is named
`aggregateCorrelation`.

## Details

This function mirrors the correlation calculation used in
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)
but operates on precomputed cell scores (e.g., obtained from
`getTransferCellScores(agg_cell_type = FALSE)`).

The normalized correlation is computed as: numerator = t(A_w1) %*% K %*%
B_w2 denominator = sqrt(sum(A_w1^2)) \* sqrt(sum(B_w2^2)) \*
\|\|K\|\|\_2 where A_w1 and B_w2 are the cell score vectors (for the
same CC index) for cell types A and B respectively, K is the kernel
matrix between the two cell types at the chosen sigma, and \|\|K\|\|\_2
is the spectral norm of K.

## Examples

``` r
# Assuming `tar_obj` is prepared and `trans_scores` was computed with
# getTransferCellScores(..., agg_cell_type = FALSE)
# res <- getTransferNormCorr(tar_obj, trans_scores, sigma_choice = 2.0)
```
