# Replicate-level inference for multi-slide CoPro analyses

Estimates an equal-replicate spatial coordination effect without
treating cells as biological replicates. For each replicate, canonical
weights and the bandwidth are learned from all other replicates and
evaluated only on the held-out replicate. The held-out effects are
combined with equal weight.

## Usage

``` r
runSlideLevelInference(
  object,
  cc_index = 1,
  sigma_values = NULL,
  replicate_id = NULL,
  alternative = "greater",
  n_resamples = 9999,
  conf_level = 0.95,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  A `CoProMulti` object after
  [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  and kernel construction. Exactly two cell types are currently
  supported.

- cc_index:

  Canonical axis to evaluate.

- sigma_values:

  Candidate sigma values. Defaults to all fitted kernel bandwidths.

- replicate_id:

  Optional slide-to-replicate mapping. Supply a vector named by
  `getSlideList(object)`; an unnamed vector is matched in slide order.
  Defaults to one independent replicate per slide.

- alternative:

  One of `"greater"`, `"less"`, or `"two.sided"`.

- n_resamples:

  Number of Monte-Carlo sign flips and bootstrap resamples.

- conf_level:

  Bootstrap confidence level.

- seed:

  Optional random seed.

- verbose:

  Print fold progress.

## Value

A list of class `CoProSlideInference` containing the equal-replicate
effect, confidence interval, replicate sign-flip p-value, Monte-Carlo
floor, and one held-out result per replicate.

## Details

`replicate_id` may group several slides from the same donor. Every fold
then leaves out all slides belonging to one donor. Sigma selection is
performed using training slides only, so the held-out effect is not
evaluated at a bandwidth selected from that replicate.

The p-value uses sign flips of the replicate-level held-out effects and
therefore assumes independent replicates and a symmetric null
distribution around zero. Up to 15 replicates all sign patterns are
enumerated exactly; larger analyses use `n_resamples` Monte-Carlo sign
flips with the Phipson–Smyth correction. The confidence interval is a
percentile bootstrap over replicates. Both summaries target the
equal-replicate mean.

## See also

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`selectSigmaByPermutation()`](https://zhen-miao.github.io/CoPro/reference/selectSigmaByPermutation.md)
