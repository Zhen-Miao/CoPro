# Get Distance Matrix

A unified accessor function for distance matrices that uses flat
structure for efficient access.

## Usage

``` r
getDistMat(
  object,
  cellType1,
  cellType2,
  slide = NULL,
  returnTranspose = FALSE,
  verbose = TRUE
)

# S4 method for class 'CoProSingle'
getDistMat(
  object,
  cellType1,
  cellType2,
  slide = NULL,
  returnTranspose = FALSE,
  verbose = TRUE
)

# S4 method for class 'CoProMulti'
getDistMat(
  object,
  cellType1,
  cellType2,
  slide = NULL,
  returnTranspose = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A CoPro object (CoProSingle or CoProMulti)

- cellType1:

  First cell type name

- cellType2:

  Second cell type name

- slide:

  Slide ID (required for CoProMulti objects, ignored for CoProSingle)

- returnTranspose:

  If TRUE, forces return of transpose when accessing symmetric matrices

- verbose:

  Whether to print detailed error messages

## Value

Distance matrix as a numeric matrix

## See also

[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)

Other accessors:
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md),
[`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)
