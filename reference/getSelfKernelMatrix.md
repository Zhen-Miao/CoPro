# Get Self-Kernel Matrix

Convenience function to retrieve self-kernel matrices computed by
computeSelfKernel.

## Usage

``` r
getSelfKernelMatrix(object, sigma, cellType, slide = NULL, verbose = TRUE)
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

## Value

Self-kernel matrix for the specified parameters
