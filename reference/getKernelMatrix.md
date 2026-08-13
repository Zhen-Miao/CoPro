# Get Kernel Matrix

A unified accessor function for kernel matrices that handles the complex
nested structure and provides symmetric access when needed.

## Usage

``` r
getKernelMatrix(
  object,
  sigma,
  cellType1,
  cellType2,
  slide = NULL,
  returnTranspose = FALSE,
  verbose = TRUE,
  materialize = TRUE
)

# S4 method for class 'CoProSingle'
getKernelMatrix(
  object,
  sigma,
  cellType1,
  cellType2,
  slide = NULL,
  returnTranspose = FALSE,
  verbose = TRUE,
  materialize = TRUE
)

# S4 method for class 'CoProMulti'
getKernelMatrix(
  object,
  sigma,
  cellType1,
  cellType2,
  slide = NULL,
  returnTranspose = FALSE,
  verbose = TRUE,
  materialize = TRUE
)
```

## Arguments

- object:

  A CoPro object (CoProSingle or CoProMulti)

- sigma:

  Sigma value for kernel selection

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

- materialize:

  For float32 kernels, whether to return a temporary standard
  double-precision sparse `Matrix`. The default is `TRUE` for
  compatibility with plotting, transfer, and user code that performs
  ordinary matrix algebra. Set to `FALSE` to retain the encoded kernel.

## Value

A standard kernel matrix, or an encoded float32 kernel when
`materialize = FALSE`.

## See also

[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md)

Other accessors:
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md),
[`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md),
[`getDistanceGeometry()`](https://zhen-miao.github.io/CoPro/reference/getDistanceGeometry.md)
