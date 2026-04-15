# Compute Normalized Correlation from Ground Truth Scores

Computes the normalized correlation using user-provided cell scores
(e.g., ground truth scores from simulation) and the pre-calculated
kernel matrix from a CoPro object.

## Usage

``` r
compute_ground_truth_ncorr(
  object,
  scores_ct1,
  scores_ct2,
  cellType1,
  cellType2,
  sigma = NULL,
  tol = 1e-04
)
```

## Arguments

- object:

  A CoPro object with kernel matrices computed

- scores_ct1:

  Numeric vector of scores for cell type 1. Must be in the same order as
  cells in the CoPro object.

- scores_ct2:

  Numeric vector of scores for cell type 2. Must be in the same order as
  cells in the CoPro object.

- cellType1:

  Name of cell type 1

- cellType2:

  Name of cell type 2

- sigma:

  Sigma value for kernel matrix (default: uses sigmaValueChoice)

- tol:

  Tolerance for SVD calculation (default: 1e-4)

## Value

Normalized correlation value

## Details

This function is useful for comparing:

- Ground truth normalized correlation (from simulated scores)

- CoPro's optimized normalized correlation (from CCA)

- Permutation distribution

Under null simulation (no cross-type coordination), the ground truth
normalized correlation should be close to 0, while CoPro's optimized
value will be higher because it searches for the maximum correlation.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running CoPro analysis on simulated data
# Get ground truth scores from metadata
meta <- br@metaDataSub
gt_scores_A <- meta$smoothed_score[meta$cell_type == "A"]
gt_scores_B <- meta$smoothed_score[meta$cell_type == "B"]

# Compute ground truth normalized correlation
gt_ncorr <- compute_ground_truth_ncorr(
  object = br,
  scores_ct1 = gt_scores_A,
  scores_ct2 = gt_scores_B,
  cellType1 = "A",
  cellType2 = "B"
)
} # }
```
