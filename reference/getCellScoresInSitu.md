# Get cell score and location information as a data.frame

Get cell score and location information as a data.frame

## Usage

``` r
getCellScoresInSitu(object, sigmaValueChoice, ccIndex = 1)

# S4 method for class 'CoProSingle'
getCellScoresInSitu(object, sigmaValueChoice, ccIndex = 1)

# S4 method for class 'CoProMulti'
getCellScoresInSitu(object, sigmaValueChoice, ccIndex = 1)
```

## Arguments

- object:

  A `CoPro` object

- sigmaValueChoice:

  A value to specify the sigma squared to use for selecting the
  particular cell score information

- ccIndex:

  Canonical vector index, default = 1

## Value

A data.frame object with cell scores and their locations

## See also

[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)

Other accessors:
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md),
[`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)
