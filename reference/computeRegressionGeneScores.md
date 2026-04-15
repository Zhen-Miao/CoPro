# Compute regression-based gene scores

For each gene, regresses its expression on the cell score and uses the
regression coefficient (beta) as the gene weight. This evaluates each
gene independently and avoids collinearity issues present in the PCA
back-projection approach used by
[`computeGeneAndCellScores`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md).

## Usage

``` r
computeRegressionGeneScores(object, sigma = NULL, verbose = TRUE)

# S4 method for class 'CoPro'
computeRegressionGeneScores(object, sigma = NULL, verbose = TRUE)

# S4 method for class 'CoProMulti'
computeRegressionGeneScores(object, sigma = NULL, verbose = TRUE)
```

## Arguments

- object:

  A CoPro object with cell scores already computed via
  [`computeGeneAndCellScores`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md).

- sigma:

  Optional numeric vector of sigma values to process. If `NULL`
  (default), all sigma values in `object@sigmaValues` are used.

- verbose:

  Logical; print progress messages. Default `TRUE`.

## Value

The input object with the `@geneScoresRegression` slot populated with
regression-based gene weights (beta coefficients). The format mirrors
`@geneScores`: a flat list keyed by
`"geneScores|sigma<value>|<cellType>"`, each entry a genes x nCC matrix.

## Details

Results are stored in the `@geneScoresRegression` slot, leaving the
original PCA-based `@geneScores` slot untouched.

For each cell type, sigma value, and canonical component (CC):

1.  Retrieves the cell score vector from `@cellScores`.

2.  Subsets `@normalizedDataSub` to cells of that type.

3.  Computes `beta_g = cov(gene_g, cellScore) / var(cellScore)` for
    every gene *g* (equivalent to simple linear regression).

4.  Stores the beta vector in `@geneScoresRegression`.

## See also

[`computeGeneAndCellScores`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)
for the PCA-based gene scores,
[`testGeneGLM`](https://zhen-miao.github.io/CoPro/reference/testGeneGLM.md)
for statistical testing with covariates.
