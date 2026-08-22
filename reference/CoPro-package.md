# CoPro: Detecting Co-Progression Between Cell Types in Spatial Transcriptomics Data

Detecting co-progression between cell types in spatial transcriptomics
data. The method works in both supervised and unsupervised manner, so it
can be used for explorative data analyses.

## Naming convention

The exported API deliberately uses two naming styles, split by layer:

- **`camelCase` – the object pipeline.** Anything whose first argument
  is a `CoProSingle` or `CoProMulti` object: constructors
  ([`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md),
  [`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)),
  the S4 generics that advance an analysis
  ([`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
  [`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
  [`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)),
  and the accessors that read results back out
  ([`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
  [`getNormCorr()`](https://zhen-miao.github.io/CoPro/reference/getNormalizedCorrelation.md),
  [`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)).
  These take an object and return an object, so they chain.

- **`snake_case` – the engine and utility layer.** Functions that work
  on plain matrices, lists, and data frames with no `CoPro` object
  involved: the CCA solvers
  ([`optimize_bilinear()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear.md),
  [`optimize_sumcor_pca()`](https://zhen-miao.github.io/CoPro/reference/optimize_sumcor_pca.md),
  [`optimize_genespace_avg_corr()`](https://zhen-miao.github.io/CoPro/reference/optimize_genespace_avg_corr.md)
  and their `_n` / `_multi_slides` variants), the spatial-null builders
  ([`resample_spatial()`](https://zhen-miao.github.io/CoPro/reference/resample_spatial.md),
  [`generate_toroidal_permutations()`](https://zhen-miao.github.io/CoPro/reference/generate_toroidal_permutations.md),
  [`diagnose_bin_distribution()`](https://zhen-miao.github.io/CoPro/reference/diagnose_bin_distribution.md)),
  and standalone helpers
  ([`quantile_normalize()`](https://zhen-miao.github.io/CoPro/reference/quantile_normalize.md),
  [`transfer_scores()`](https://zhen-miao.github.io/CoPro/reference/transfer_scores.md),
  [`copro_download_data()`](https://zhen-miao.github.io/CoPro/reference/copro_download_data.md)).
  Call these directly when you want the numerical core without the
  object wrapper.

A few exports sit outside that rule:

- [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
  and
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md)
  keep the `camelCase` stem of
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  with a `_Variant` suffix, so the three permutation entry points sort
  together.

- [`calculate_pvalue()`](https://zhen-miao.github.io/CoPro/reference/calculate_pvalue.md),
  [`calculate_pvalue_stepdown()`](https://zhen-miao.github.io/CoPro/reference/calculate_pvalue_stepdown.md),
  [`compute_ground_truth_ncorr()`](https://zhen-miao.github.io/CoPro/reference/compute_ground_truth_ncorr.md),
  and
  [`fit_score_reference()`](https://zhen-miao.github.io/CoPro/reference/fit_score_reference.md)
  take a `CoPro` object but are `snake_case`. They derive standalone
  inference or transfer results rather than advancing and returning the
  object pipeline.

## Package options

Behavior that used to be set through
[`options()`](https://rdrr.io/r/base/options.html) is now reachable as
function arguments; the options are still read to supply those
arguments' defaults, so existing scripts keep working.

- `CoPro.factorizePermutation` – default for the `factorize` argument of
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  and friends.

- `CoPro.compactPermutation` – default for their `compactPermutation`
  argument.

- `CoPro.float32Threads` – default for the `nThreads` argument of
  [`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md).

- `CoPro.useRcppFRNN` – default for the compiled fixed-radius neighbor
  engine used internally by the sparse kernel routes.

## See also

Useful links:

- <https://zhen-miao.github.io/CoPro/>

- <https://github.com/Zhen-Miao/CoPro>

- Report bugs at <https://github.com/Zhen-Miao/CoPro/issues>

## Author

**Maintainer**: Zhen Miao <zhenmiao@sas.upenn.edu>
([ORCID](https://orcid.org/0000-0002-3255-9517))
