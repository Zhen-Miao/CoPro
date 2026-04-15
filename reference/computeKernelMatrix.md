# Compute Kernel Matrix for CoPro

This method calculates the kernel matrices for pairs of cell types based
on their distances and a range of sigma values. The formula of
calculating kernel matrix is: \$\$K(x, y) =
\exp\left(-\frac{\\x-y\\^2}{2 \sigma^2}\right)\$\$ The matrices are
adjusted by clipping the upper quantile of the values to reduce the
effect of outliers. The results are stored within the object.

## Usage

``` r
computeKernelMatrix(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  verbose = TRUE
)

# S4 method for class 'CoProSingle'
computeKernelMatrix(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  verbose = TRUE
)

# S4 method for class 'CoProMulti'
computeKernelMatrix(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A `CoPro` object.

- sigmaValues:

  A vector of sigma values used for kernel calculation.

- lowerLimit:

  The lower limit for the kernel function, default is 1e-7.

- upperQuantile:

  The quantile used for clipping the kernel values, default is 0.85.

- normalizeKernel:

  Whether to normalize the kernel matrix? Default = FALSE. Note that
  normalization will not affect any downstream analyses, it is for
  numerical stability and easier interpretation only.

- minAveCellNeighor:

  What is the minimum average number of cell in the neighbor? This step
  is to help set up the expected sparsity of the kernel matrix. If a
  kernel sigma value is too small, this result in too few neighbors for
  most cells, resulting in an overly-sparse matrix that makes the
  parameter estimation hard. Thus, the sigma values that results in an
  overly-sparse matrix will be removed for later analysis.

- rowNormalizeKernel:

  Whether the kernel matrix will be row-wise normalized? Note that row
  or column wise normalization will result in an asymmetric result in
  skrCCA inference.

- colNormalizeKernel:

  Whether the kernel matrix will be column-wise normalized? Note that
  row or column wise normalization will result in an asymmetric result
  in skrCCA inference.

- verbose:

  Whether to output the progress and related information

## Value

The `CoPro` object with computed kernel matrices added. The kernel
matrices are organized into a three-layer nested list object. The first
layer is indexed by the sigma value, and the second and the third layers
are cell types

## Note

To-do: Shall we include row or column normalization of the kernel?
