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
  verbose = TRUE
)

# S4 method for class 'CoProSingle'
getKernelMatrix(
  object,
  sigma,
  cellType1,
  cellType2,
  slide = NULL,
  returnTranspose = FALSE,
  verbose = TRUE
)

# S4 method for class 'CoProMulti'
getKernelMatrix(
  object,
  sigma,
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

## Value

Kernel matrix as a numeric matrix
