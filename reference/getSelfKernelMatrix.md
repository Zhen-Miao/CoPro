# Get Self-Kernel Matrix

Convenience function to retrieve self-kernel matrices computed by
computeSelfKernel.

## Usage

``` r
getSelfKernelMatrix(
  object,
  sigma,
  cellType,
  slide = NULL,
  verbose = TRUE,
  materialize = TRUE
)
```

## Arguments

- object:

  A CoPro object

- sigma:

  Sigma value

- cellType:

  Cell type name

- slide:

  Slide ID (for CoProMulti objects)

- verbose:

  Whether to print error messages

- materialize:

  For encoded float32 kernels, whether to return a temporary standard
  double-precision sparse matrix.

## Value

Self-kernel matrix for the specified parameters
