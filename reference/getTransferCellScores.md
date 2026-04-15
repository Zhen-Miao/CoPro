# Get cell score by transferring gene weights from another slide By default, quantile normalization is used to ensure distribution match

Get cell score by transferring gene weights from another slide By
default, quantile normalization is used to ensure distribution match

## Usage

``` r
getTransferCellScores(
  ref_obj,
  tar_obj,
  sigma_choice,
  use_quantile_normalization = TRUE,
  agg_cell_type = FALSE,
  gs_weight_threshold = 0,
  sigma_choice_tar = NULL,
  gene_score_type = c("PCA", "regression"),
  verbose = TRUE
)
```

## Arguments

- ref_obj:

  Reference object (where the gene weights will be obtained)

- tar_obj:

  Target object (where the cell scores will be obtained)

- sigma_choice:

  Sigma value to be used

- use_quantile_normalization:

  Logical; apply quantile normalization of target to reference
  distribution (default TRUE)

- agg_cell_type:

  Logical; if TRUE, returns a single matrix aggregated across cell types
  (default FALSE)

- gs_weight_threshold:

  Numeric; absolute threshold used to zero out small gene weights in the
  gene score matrix prior to transfer (default `0`).

- sigma_choice_tar:

  Numeric; sigma value for target object. If NULL (default), uses
  sigma_choice. Not recommended for general use.

- gene_score_type:

  Character; which gene score slot to use for transfer. `"PCA"`
  (default) uses the PCA back-projection weights in `@geneScores`.
  `"regression"` uses the regression-based weights in
  `@geneScoresRegression`, which avoids collinearity issues and produces
  more robust transfers. The regression slot must have been populated by
  calling
  [`computeRegressionGeneScores`](https://zhen-miao.github.io/CoPro/reference/computeRegressionGeneScores.md)
  on the reference object.

- verbose:

  verbose

## Value

cell scores as a matrix
