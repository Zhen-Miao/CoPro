# Transfer cell scores between matrices using gene weights

Given reference (`mat_A`) and target (`mat_B`) cell-by-gene matrices
with the same set and order of genes (columns), this function optionally
quantile normalizes the target to match the reference per feature,
standardizes the target using the reference mean and standard deviation,
filters small-magnitude gene weights, and computes cell scores via
matrix multiplication \\B\_{standardized} \\times W\\.

## Usage

``` r
transfer_scores(
  mat_A,
  mat_B,
  gs_ct,
  use_quantile_normalization = TRUE,
  gs_weight_threshold = 0,
  verbose = TRUE
)
```

## Arguments

- mat_A:

  Numeric matrix of shape cells-by-genes (reference). Column order must
  correspond to genes in `gs_ct`.

- mat_B:

  Numeric matrix of shape cells-by-genes (target). Must have the same
  genes (columns) and order as `mat_A`.

- gs_ct:

  Numeric matrix of gene weights of shape genes-by-K, where K is the
  number of signatures or cell types. Column names, if present, are used
  in verbose messages.

- use_quantile_normalization:

  Logical; whether to quantile-normalize `mat_B` to `mat_A` per feature
  before standardization (default `TRUE`).

- gs_weight_threshold:

  Numeric; absolute threshold used to zero out small gene weights in
  `gs_ct` (default `0`).

- verbose:

  Logical; whether to print progress messages (default `TRUE`).

## Value

A numeric matrix of cell scores with shape (nrow(`mat_B`) x
ncol(`gs_ct`)).

## Details

- If `use_quantile_normalization = TRUE`, `mat_B` is quantile-normalized
  to the distribution of `mat_A` for each gene.

- Columns of the normalized `mat_B` are then centered and scaled using
  the column means and standard deviations computed from `mat_A` (with a
  small safeguard for near-zero standard deviations).

- Gene weights in `gs_ct` are thresholded by absolute value per column;
  values with `abs(weight) < gs_weight_threshold` are set to 0.

- Final cell scores are obtained as the matrix product of the
  standardized `mat_B` with the filtered weight matrix `gs_ct`.

This function assumes that `mat_A`, `mat_B`, and `gs_ct` are aligned on
the same set and ordering of genes (columns of `mat_A`/`mat_B`, rows of
`gs_ct`). Callers are responsible for subsetting/reordering if needed.

## See also

[getTransferCellScores](https://zhen-miao.github.io/CoPro/reference/getTransferCellScores.md)
