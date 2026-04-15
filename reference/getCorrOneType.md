# Retrieve the Correlation within one cell types

Retrieve the Correlation within one cell types

## Usage

``` r
getCorrOneType(object, cellTypeA, ccIndex = 1, sigmaValueChoice = NULL)
```

## Arguments

- object:

  A `CoPro` object

- cellTypeA:

  Cell type label for the cell type of interest

- ccIndex:

  Canonical vector index, default = 1

- sigmaValueChoice:

  A particular sigma squared value for the correlation

## Value

A data.frame with two columns, AK and B, where AK represents the cell
score of cell type A times the kernel matrix, and B represents the cell
score of cell type B.
