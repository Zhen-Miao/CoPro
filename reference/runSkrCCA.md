# runSkrCCA

runSkrCCA

## Usage

``` r
runSkrCCA(
  object,
  scalePCs = TRUE,
  nCC = 2,
  tol = 1e-05,
  transferred_weight_1 = NULL,
  maxIter = 200,
  sigmaChoice = NULL,
  n_cores = 1,
  step_size = 1
)

# S4 method for class 'CoPro'
runSkrCCA(
  object,
  scalePCs = TRUE,
  nCC = 2,
  tol = 1e-05,
  transferred_weight_1 = NULL,
  maxIter = 200,
  sigmaChoice = NULL,
  n_cores = 1,
  step_size = 1
)

# S4 method for class 'CoProMulti'
runSkrCCA(
  object,
  scalePCs = TRUE,
  nCC = 2,
  tol = 1e-05,
  transferred_weight_1 = NULL,
  maxIter = 200,
  sigmaChoice = NULL,
  n_cores = 1,
  step_size = 1
)
```

## Arguments

- object:

  A CoPro object

- scalePCs:

  Whether to scale each PCs to a uniform variance before running the
  program

- nCC:

  Number of canonical vectors to compute, default = 2

- tol:

  Tolerance for termination, default = 1e-5

- transferred_weight_1:

  If we use cross-slide weight transfer function, the transferred weight
  on each PC. Otherwise, the value should be set to NULL.

- maxIter:

  Maximum iterations

- sigmaChoice:

  Specific sigma value to use (CoProMulti only, ignored for CoPro)

- n_cores:

  Number of cores for parallel processing (CoProMulti only, ignored for
  CoPro)

- step_size:

  Step size for damped power iteration. Default 1 (standard power
  iteration). Values in (0,1) blend old and new weights for smoother
  convergence, which can help with many cells or many CCs.

## Value

CoPro object with skrCCA results computed

## See also

[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)

## Examples

``` r
# \donttest{
toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
obj <- newCoProSingle(
  normalizedData = toy$normalizedData,
  locationData   = toy$locationData,
  metaData       = toy$metaData,
  cellTypes      = toy$cellTypes
)
obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
obj <- computePCA(obj, nPCA = 10)
#> Input is dense (matrixarray), performing irlba pca...
#> Input is dense (matrixarray), performing irlba pca...
obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
#> Running skrCCA [1/2] for sigma = 0.05 ...
#> [1] "Convergence reached at 4 iterations (Max diff = 7.744e-06 )"
#> [1] "Convergence reached at 0 iterations (Max diff = 1.080e-14 )"
#> Running skrCCA [2/2] for sigma = 0.1 ...
#> [1] "Convergence reached at 4 iterations (Max diff = 4.003e-07 )"
#> [1] "Convergence reached at 0 iterations (Max diff = 1.943e-15 )"
#> skrCCA finished 2 sigma value(s) in 0.0 s.
#> Optimization succeeded for 2 sigma value(s): sigma_0.05, sigma_0.1
# }
```
