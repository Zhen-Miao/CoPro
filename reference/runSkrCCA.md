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
