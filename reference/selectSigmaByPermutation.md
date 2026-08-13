# Select the kernel bandwidth by studentized permutation (max-T)

Chooses the spatial bandwidth \\\sigma\\ by comparing the co-progression
statistic at each candidate bandwidth to *its own* permutation null,
rather than by taking the largest normalized correlation across
bandwidths.

## Usage

``` r
selectSigmaByPermutation(
  object,
  sigmaValues = NULL,
  ccIndex = 1L,
  nPermu = 199L,
  alternative = c("greater", "two.sided"),
  minSigma = "spacing",
  maxCells = 2000L,
  domain = NULL,
  alpha = 0.05,
  blockSize = 1024L,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  A `CoProSingle` object that has been through
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  and
  [`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
  so per-bandwidth cell scores exist.

- sigmaValues:

  Candidate bandwidths. Defaults to `object@sigmaValues`. Every value
  must have stored cell scores.

- ccIndex:

  Which canonical component to select on. Default 1.

- nPermu:

  Number of toroidal draws \\B\\. Default 199, giving a Monte-Carlo
  floor of 0.005 under the Phipson–Smyth correction. Must be at least 2,
  since a single draw has no spread to studentize by.

- alternative:

  `"greater"` (default, matching the sign the optimizer fixes) or
  `"two.sided"`.

- minSigma:

  Floor for the scanned grid. `"spacing"` (default) drops candidates
  below the median nearest-partner distance, measured on the full
  coordinates. A single positive number sets the floor explicitly;
  `NULL` scans the whole grid. See the note on small bandwidths below.

- maxCells:

  Cap on the number of cells sampled per cell type. The cost is
  quadratic in this, and the same sample is used for the observed
  statistic and for every draw, so the studentization stays internally
  consistent. `NULL` uses every cell. Default 2000.

- domain:

  Optional `list(lower = , upper = )` numeric vectors giving the wrap
  box per coordinate axis, on the object's distance scale. Defaults to
  the range of the sampled coordinates.

- alpha:

  Level used to report `plateau`. Default 0.05.

- blockSize:

  Rows per block when accumulating the bilinear statistic. Caps peak
  memory; does not change the result. Default 1024.

- seed:

  Optional integer RNG seed.

- verbose:

  Whether to report progress.

## Value

An object of class `CoProSigmaSelection`: a list with

- `sigmaValueChoice` — the selected bandwidth.

- `pValue` — the max-T permutation p-value at that bandwidth.

- `zMax` — the observed \\\max z\\.

- `plateau` — bandwidths with `pAdjusted <= alpha`.

- `perSigma` — a `data.frame` with `sigma`, `cellType1`, `cellType2`,
  `statistic` (\\T\\), `nullMean` (\\m\\), `nullSD` (\\s\\), `z`, and
  `pAdjusted`.

- `nullMax` — the length-`nPermu` null distribution of the maximum.

- `cells` — the sampled cell indices per cell type.

- `spacing` — the median nearest-partner distance per pair.

- `settings` — `nPermu`, `alternative`, `ccIndex`, `alpha`, `maxCells`,
  `null`, and `minSigma` (the floor actually applied).

## Details

### Why not just take the largest normalized correlation

[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)
divides the bilinear statistic \\T(\sigma) = a' K(\sigma) b\\ by a norm
of the kernel, and `object@sigmaValueChoice` is the bandwidth where that
ratio is largest. No available denominator makes the *null level* of
that ratio constant in \\\sigma\\: the un-whitened \\\\K_c\\\_F\\
ignores within-type spatial autocorrelation and drifts down with
\\\sigma\\, while whitened variants need a within-type correlation
operator that the data do not pin down (the principal components of one
cell type do not share a single correlation length). Comparing the ratio
across \\\sigma\\ therefore compares numbers whose null expectations
differ, and the argmax is biased toward whichever bandwidth happens to
have the highest floor.

### What this function does instead

It measures the floor rather than modelling it. For each \\\sigma\\ on
the grid it computes the observed statistic \\T(\sigma)\\ at the cell
scores CoPro already fitted at that \\\sigma\\, then estimates the mean
\\m(\sigma)\\ and the standard deviation \\s(\sigma)\\ of \\T(\sigma)\\
under a spatial null and selects

\$\$\hat\sigma = \arg\max\_\sigma\\ \frac{T(\sigma) -
m(\sigma)}{s(\sigma)}.\$\$

Because both the location and the scale are read off the null itself,
the studentized statistic \\z(\sigma)\\ has the same null level at every
bandwidth by construction, and it carries no tuning constant.

Centering is not cosmetic. Across two cell types \\m(\sigma)\\ is zero
in expectation — a wrap-around shift leaves every shifted cell equally
likely to land anywhere in the box, so each column of the kernel has the
same mean and a centered score vector annihilates it — and subtracting
the estimate costs only Monte-Carlo noise. *Within* one cell type it is
not zero: the \\w\_{ii} = 0\\ convention takes \\\sum_i a_i^2\\k(x_i,
x_i + \delta)\\ off every draw, and that term has a negative expectation
which grows with \\\sigma\\. On the package's toy object (\\B = 4000\\)
the null mean ran from \\-0.06\\s(\sigma)\\ at \\\sigma = 0.05\\ to
\\-0.36\\s(\sigma)\\ at \\\sigma = 0.2\\ — a \\\sigma\\-dependent floor
of exactly the kind this function exists to remove, and one that
dividing by \\s(\sigma)\\ alone leaves in place.

### The null, and why the same draws serve selection and inference

The null is a toroidal (rigid wrap-around) shift of the coordinates. A
rigid shift preserves each cell type's own spatial autocorrelation — the
structure that a plain label shuffle destroys, and whose destruction is
what makes shuffle-based nulls anti-conservative — while removing the
alignment between the score fields. It assumes spatial stationarity and
wrap-around, which is false at tissue edges; see
[`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
for the bin-wise alternative.

One draw is one offset *per cell type*, reused everywhere that type
appears. With three or more cell types a type sits in several pairs at
once, so drawing its offset afresh per pair would put the same cells in
two places inside one draw: that row of statistics would come from no
single configuration, and a max-T reference built from such rows is not
a null for the scan. The first cell type anchors the frame, which costs
no randomization — the offset *between* any two types is still uniform
over the box — and leaves one fewer cloud cut by the wrap seam. The
observed configuration is the zero offset, so it is itself an admissible
draw, which is what the \\+1\\ in the p-value below refers to. With a
single cell type that type is its own partner, so the shift is applied
to one copy while the other is held: a common offset would move both
sides together and leave the statistic where it started.

One pass of \\B\\ draws evaluates the null at *every* bandwidth on the
grid, so the draws are coupled across \\\sigma\\. That coupling is what
makes the accompanying p-value valid: the null distribution of
\\\max\_\sigma z(\sigma)\\ is available from the same draws, so
comparing the observed \\\max\_\sigma z(\sigma)\\ to it accounts for
having scanned the grid (single-step Westfall–Young max-T). The
selection is replicated inside the null, so the reported p-value is not
circular. Running
[`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
once per bandwidth would *not* give this: it redraws permutations on
every call and its default bin grid is itself \\\sigma\\-dependent, so
the per-bandwidth nulls are neither coupled nor comparable.

`perSigma$pAdjusted` applies the same max-T reference to each bandwidth
individually, giving the scales at which coordination is detectable
after adjusting for the scan. `plateau` collects those with
`pAdjusted <= alpha`. It answers "where is there detectable signal", not
"which bandwidths are indistinguishable from the best" — \\z(\sigma)\\
is typically flat-topped and \\\hat\sigma\\ is weakly identified, so
treat it as a representative scale within a band and check that
conclusions hold across the band.

### Scope, and an honest limitation

The canonical directions are held at their observed values inside each
draw; they are not re-optimized. That is what is wanted for *selection*
— it measures the \\\sigma\\-shape of the noise floor with the signal
held fixed — and it is what makes one \\O(B)\\ pass enough. It also
means the p-value inherits the mild anti-conservativeness of a
fixed-direction null, most visibly at small \\\sigma\\. For a
re-optimizing test at a chosen bandwidth use
[`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md),
or
[`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
for a re-optimizing max-over-sigma test of the raw normalized
correlation.

That small-\\\sigma\\ inflation is why `minSigma` defaults to
`"spacing"`. Below about one cell spacing a Gaussian kernel has almost
no mass off the diagonal, the statistic rests on a handful of pairs, and
the fixed-direction null understates the floor — so \\z\\ keeps climbing
as \\\sigma\\ shrinks and the argmax rails at whatever the smallest
candidate happens to be. Flooring the scan at the median nearest-partner
distance removes that regime. If the selected bandwidth still lands at
an end of the grid the function warns: the scan found the edge, not an
optimum.

The kernel used here is the plain Gaussian
\\\exp(-\tfrac{1}{2}(d/\sigma)^2)\\ on the same *Euclidean* distance
scale
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
used. Euclidean is a requirement, not a default: the kernel is rebuilt
from coordinates at every candidate bandwidth and under every shifted
configuration, and a `"Morphology-Aware"` metric is not a function of
the coordinates alone — its geodesic filter rests on a k-NN graph of the
unshifted tissue, which a shift invalidates. Such objects are refused
rather than silently rescored with a Euclidean kernel. Three shaping
steps that the kernel pipeline applies are deliberately *not*
replicated: the `upperQuantile` cap, the `lowerLimit` truncation, and
`computeDistance(truncateLowDist = TRUE)`'s floor on the bottom ~0.1% of
distances. The cap is a data-dependent quantile, so it would be
recomputed on every shifted configuration and the statistic would stop
meaning the same thing across draws; the other two act only where the
Gaussian is flat (at distances far below \\\sigma\\, or at kernel values
near zero) and so move \\T\\ negligibly. The selected bandwidth is
therefore a bandwidth for the same kernel family on the same distance
scale, and can be passed straight back to
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md).

Because the directions are fixed, the statistic does not depend on which
canonical criterion produced them: `sumcov` and `sumcor` weights are
scored the same way, and no criterion is re-run inside the null.

### Cell types

With one cell type this is the within-type (self) problem: `a` and `b`
are the same score vector, the kernel is the within-type kernel with a
zero diagonal (a cell is not its own neighbour, the \\w\_{ii} = 0\\
convention), and the null shifts one copy of the coordinates against the
other. With two or more, every pair is scanned and the max-T reference
is taken over the whole (bandwidth x pair) grid, so `pAdjusted` is
adjusted for both.

## References

Meinshausen N, Maathuis MH, Buhlmann P (2011). Asymptotic optimality of
the Westfall–Young permutation procedure for multiple testing under
dependence. *Annals of Statistics* 39(6):3369–3391.

Phipson B, Smyth GK (2010). Permutation p-values should never be zero.
*Statistical Applications in Genetics and Molecular Biology* 9(1):39.

## See also

[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md)
to build the candidate grid,
[`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
for a re-optimizing test at a fixed bandwidth,
[`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
for a re-optimizing max-over-sigma test.

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- computeKernelMatrix(obj, sigmaValues = detectSigmaRange(obj)$sigmaValues)
obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
obj <- computeNormalizedCorrelation(obj)
obj <- computeGeneAndCellScores(obj)

sel <- selectSigmaByPermutation(obj, nPermu = 199, seed = 1)
sel
sel$sigmaValueChoice
} # }
```
