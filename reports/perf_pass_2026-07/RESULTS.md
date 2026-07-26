# Performance and memory pass, July 2026

Measured by A/B-ing two checkouts of this repository — the pre-pass commit
`e5ddaa4` and the branch head — driven by the scripts in this folder. Session
machine: Darwin 25.1.0, Apple silicon, R 4.5.2, clang `-O2`.

Every timed arm is repeated and reported as the minimum, which is the standard
defence against a co-tenant process inflating one sample. Where a single-shot
figure is quoted it is called out as such.

Numeric equivalence is verified separately by
`reports/perf_pass_baseline/compare_baseline.R`, which re-runs 14 pipeline
scenarios and diffs every downstream result. See section 5 — there are two
distinct claims there, and merging them would misrepresent both.

---

## Headline: end-to-end pipeline

`benchmark_end_to_end.R` → `end_to_end_before.csv`, `end_to_end_after.csv`
(per-repetition detail in the `*_runs.csv` files)

60,000 cells, 2 cell types, 120 genes, `nPCA = 30`, 2 sigmas, 99 permutation
draws. The pipeline is stateful, so it is the *whole* pipeline that is repeated,
3 times per arm; seconds are the per-stage minimum and `peak MB` the per-stage
maximum of Vcells allocated above the resting set.

| stage | before | after | speedup | peak before | peak after |
|---|---|---|---|---|---|
| `computePCA` | 0.73 s | 0.63 s | 1.15x | 798 MB | **597 MB** |
| `computeSparseKernelFloat32` | 42.25 s | 31.61 s | **1.34x** | 1575 MB | 1576 MB |
| `runSkrCCA` | 3.50 s | 2.74 s | **1.28x** | 32 MB | 32 MB |
| `computeNormalizedCorrelation` | 6.44 s | 3.46 s | **1.86x** | 39 MB | 38 MB |
| `computeGeneAndCellScores` | 0.04 s | 0.04 s | 0.92x | 79 MB | 79 MB |
| `runSkrCCAPermu` (99 draws) | 4.00 s | 3.53 s | 1.13x | 787 MB | 788 MB |
| `computeNormalizedCorrelationPermu` | 6.05 s | 3.49 s | **1.73x** | 1272 MB | **825 MB** |
| `@cellPermu` stored size | — | — | — | 24 MB | **12 MB** |
| **total** | **63.0 s** | **45.5 s** | **1.38x** | | |

`computeGeneAndCellScores` at 0.92x is 4 ms against 4 ms on untouched code — it
is a null control, not a regression. `computePCA` and `runSkrCCAPermu` are
nearly-null controls that moved slightly because they *are* touched (per-column
SDs and permutation storage respectively).

The `@cellPermu` halving is at 99 draws and 30k cells per type. It scales with
`n * nPermu`: at 200,000 cells and 999 draws the held-fixed type alone goes
from 799 MB to a few dozen bytes.

---

## 1. Row-major packing of the float32 operands

`benchmark_float32_operators.R` → `operators_strided.csv`, `operators_packed.csv`

Kernels built from uniform coordinates on a 1000 x 1000 field, giving the
"hundreds of neighbours per cell" density that
`reports/kernel_precision_benchmark/RESULTS.md` reports for real data. Seconds
per call, minimum of 3 x 3 reps.

| case | nnz/row | nPC | threads | op | strided | packed | speedup |
|---|---|---|---|---|---|---|---|
| cross-type | 354 | 30 | 1 | `kernelXKY` | 16.512 | 4.658 | **3.54x** |
| cross-type | 354 | 30 | 4 | `kernelXKY` | 3.978 | 1.198 | **3.32x** |
| cross-type | 354 | 30 | 4 | `matMult` | 3.270 | 1.084 | **3.02x** |
| cross-type | 354 | 10 | 1 | `kernelXKY` | 2.280 | 1.493 | 1.53x |
| cross-type | 354 | 10 | 4 | `kernelXKY` | 0.503 | 0.373 | 1.35x |
| within-type | 374 | 30 | 1 | `kernelXKY` | 20.705 | 13.283 | 1.56x |
| within-type | 374 | 30 | 4 | `kernelXKY` | 10.485 | 6.314 | 1.66x |
| within-type | 374 | 10 | 4 | `kernelXKY` | 2.032 | 1.492 | 1.36x |
| within-type | 374 | 10 | 4 | `matMult` | 1.033 | 1.003 | 1.03x |

**3.0-3.5x on the cross-type operators at nPC = 30**, tapering to ~1.3-1.6x at
nPC = 10 or on the within-type kernel. The gain tracks nPC because that is how
many separate cache lines the old column-major gather touched per nonzero.

`matMult` is unchanged on the within-type kernel because symmetric input takes
the transposed branch, which splits over right-hand sides and already reads
contiguously; only the general branch was repacked.

### Variance, and a correction

`benchmark_symmetric_variance.R` re-measures the sweep's first cell, which
initially reported the packed build as *slower* (0.73x):

```
strided  p=10 threads=1 : min 6.886  median 8.575   all [6.89 8.04 8.54 8.61 8.78 9.83]
packed   p=10 threads=1 : min 4.870  median 5.455   all [4.87 5.32 5.38 5.53 6.40 6.83]
strided  p=30 threads=1 : min 28.532 median 44.254  all [28.53 36.51 42.69 45.82 46.99 54.81]
packed   p=30 threads=1 : min 16.608 median 16.685  all [16.61 16.67 16.68 16.69 16.71 16.98]
```

No regression — the packed build wins on both min and median; the apparent one
was cold-start noise in the first cell of the sweep.

The spread is itself a result. The strided build ranges 28.5-54.8 s on the same
input while the packed build holds 16.6-17.0 s. A loop touching nPC cache lines
per nonzero is memory-bandwidth-bound and hostage to whatever else is running;
the packed loop fetches roughly nPC times less and stays predictable. On a
shared HPC node that matters as much as the mean. It is also why every arm in
this document is repeated.

**Correction to the pre-implementation estimate.** A standalone microbenchmark
written before the change predicted 4.0-6.6x. It used 40 nonzeros per row; real
kernels have ~350-375, and at that density more of the working set stays
resident regardless of layout. The in-package numbers above supersede it.

---

## 2. Kernel construction

120k cells, symmetric, ~375 neighbours per cell. Bit-identical to the serial
build at every thread count and normalization mode.

| configuration | 1 thread | 4 threads | 8 threads |
|---|---|---|---|
| `normalization = 0`, 1 sigma | 11.31 s | 7.19 s (1.57x) | 6.09 s (**1.86x**) |
| `normalization = 0`, 3 sigmas | 28.43 s | 22.35 s (1.27x) | 20.51 s (1.39x) |
| `normalization = 2`, 1 sigma | 14.76 s | 10.88 s (1.36x) | 9.55 s (1.55x) |

Enumeration was ~3.7 s of the 12.3 s single-sigma build; parallelizing it is
worth proportionally less as sigmas are added, because the per-sigma passes are
serial and dominate.

### What parallel enumeration costs in memory

`benchmark_enumeration_memory.R` → `enumeration_memory.csv`

Enumerating in parallel means k private edge buffers that are then concatenated
into one contiguous array. Both are allocated at the moment of the first
insert, so the *bound* on simultaneously-allocated edge storage rises. R's
`gc()` cannot see any of this — it counts Vcells only — so peak RSS is sampled
from `ps` alongside the builder's own `temporary_bytes` figure, which reports
that bound.

120,000 cells, 264M candidate edges:

| threads | seconds | bound (`temporary_bytes`) | peak RSS |
|---|---|---|---|
| 1 | 61.0 s | 4594 MB | 4271 MB |
| 4 | 43.5 s | 6609 MB | 4277 MB |
| 8 | 37.6 s | 6609 MB | 4280 MB |

**The bound rises 44%; measured resident memory rises 0.1%.** Two reasons: a
reserved allocation is only made resident as it is written, and each private
buffer is released as it is consumed, so the copy that the bound assumes is
live all at once never is. With one thread there is nothing to concatenate and
the single buffer is moved into place, so the serial path allocates strictly
less than before this pass.

This deserved checking rather than assuming, because the pass rejected the
Gaussian-weight optimization below for a ~50% increase in the largest
temporary. That one would have been a *permanent* increase held across the
whole build; this is a transient during one concatenation, and it does not show
up in RSS. Eliminating it entirely would need enumeration and concatenation to
be pipelined — a real concurrency rewrite, not justified by a 0.1% measurement.

---

## 3. Whitened-Frobenius normalizer

`benchmark_whitened_frobenius.R` → `frobenius_before.csv`, `frobenius_after.csv`,
`frobenius_block_sweep.csv`

`(Rx %*% K) %*% Ry` fills in 11x, and only to produce a scalar. Streaming
`sum((Rx K) * (K Ry))` over column blocks is both faster and smaller — the
intermediates stay in cache instead of streaming a hundreds-of-MB product
through memory.

True A/B on the shipped `.whitenedFrobNorm()` in both checkouts, 150k cells,
nnz(K) 3.91M, minimum of 5 repetitions:

| arm | min | median | peak |
|---|---|---|---|
| before (`e5ddaa4`) | 6.16 s | 6.21 s | 2231 MB |
| after | **4.35 s** | **4.43 s** | **609 MB** |

**1.42x faster and 3.7x less peak.** The two arms agree to 1.1e-15 relative
(7108.1463071682192 against 7108.146307168211).

**This corrects an earlier version of this document**, which claimed 1.46x
against a *benchmark-local reimplementation* of the old path rather than the old
path itself. That stand-in omitted the low-rank cross terms of the
double-centering expansion, which the real function computes and which this
change did not speed up — so it credited the change with work it never did. The
in-package figure is 1.42x. The memory claim was, if anything, understated.

The block budget is a genuine parameter, swept on the new helper alone:

| block budget | blocks | min | peak |
|---|---|---|---|
| unblocked | 1 | 5.48 s | 2026 MB |
| `2e6` | 2 | 4.77 s | 1525 MB |
| `1e6` | 4 | 4.65 s | 970 MB |
| **`2e5`** (default) | **20** | **3.73 s** | **528 MB** |
| `5e4` | 79 | 4.16 s | 529 MB |
| `2e4` | 196 | 5.59 s | 529 MB |

A budget large enough to yield a single block is worse than not blocking (both
one-sided products live at once); very small blocks lose to per-block overhead.
The default sits at the knee.

**This also corrects the plan's estimate**, which predicted blocking would be
*slightly slower* and worth taking for memory alone. That was based on a
badly-sized trial block; at a sensible block size it wins on both.

Floating-point reassociation shifts the normalizer by ~3.7e-15 relative.
Weights, cell scores, gene scores and the selected sigma stay bit-identical.

---

## 4. Smaller results

- **Sparse-kernel quantile**, 30M stored values: 0.88 s / 1622 MB →
  0.43 s / 1022 MB (**2.1x**, 600 MB less). Bit-identical, verified over 366
  combinations of length, repetition count and probability.
- **Kernel signature probe**, 11.2M nonzeros: 0.339 s → 0.000167 s
  (**~2000x**). Moves one scan to construction (+11%, 3.19 s → 3.53 s) and
  removes six per (sigma, pair).
- **PC score storage**, 12k cells x 30 PCs, 4 slides: `@pcaResults` 3.89 MB →
  0.054 MB; total PC-score memory 7.77 MB → 3.93 MB.
- **Per-column SDs**, `apply(x, 2, sd)` → `matrixStats::colSds()`, minimum of 5:

  | matrix | apply | colSds | speedup |
  |---|---|---|---|
  | 30,000 x 120 | 0.030 s / 47 MB | 0.013 s / 0.5 MB | 2.31x |
  | 50,000 x 500 | 0.207 s / 332 MB | 0.088 s / 0.4 MB | 2.35x |
  | 5,000 x 2,000 | 0.095 s / 190 MB | 0.035 s / 0.1 MB | 2.71x |

  The memory column is the larger result and was not part of the original
  rationale: `apply()` transposes the whole matrix and allocates a list of
  column vectors, so it costs a full copy of the input per call. That is most of
  why `computePCA`'s peak drops 798 → 597 MB in the headline table.

---

## 5. Numeric equivalence: two claims, kept apart

`reports/perf_pass_baseline/` runs 14 scenarios — three kernel backends,
within-type, multi-slide at both `center_per_slide` and both `scalePCs`
settings, and **6 permutation combinations** (the `"global"` null at all three
`permu_which` values, plus `"bin"`, `"toroidal"` and `"pc"` at
`second_only`) — and diffs every downstream result.

Each scenario now computes cell scores, gene scores and regression gene weights.
Earlier it did not: `run_within_type()` and `run_permutation()` never called
`computeGeneAndCellScores()`, so their score slots were empty and compared
equal to each other while proving nothing — and `within_type_float32` is exactly
the scenario the Frobenius reassociation perturbs. `capture_baseline.R` now
refuses to record an empty slot, so that cannot recur silently.

**Claim one — the eight performance items** (`e5bc8c3` … `50e9095`), against a
pre-pass capture at `e5ddaa4` with the current harness:

> **13 of 14 scenarios bit-identical.** `within_type_float32` differs only in
> `normalizedCorrelation`, by ~3.7e-15 — the reassociation in section 3. Its
> cell scores, gene scores and regression gene weights are bit-identical, which
> is now actually tested rather than vacuously true.

**Claim two — adopting `matrixStats::colSds()`** is *not* bit-identical, and an
earlier commit message that said it was was wrong. `colSds()` and `sd()` use
different variance algorithms and disagree by 1 ulp on some columns (1.11e-16 on
Gene56 of 60 in the fixture). That is enough to flip the sign of a principal
component, and the CCA weight's coordinate on that PC flips with it.

> `@skrCCAOut` — the raw weight vectors, whose per-PC sign is a free convention
> — differs in 3 of 14 scenarios. **Every reported quantity agrees to 5.4e-11
> relative with zero sign changes**, and the selected sigma is identical in all
> 14.

Worst case over all 14 scenarios, for the branch head against the pre-pass
capture (so these include the section-3 reassociation as well):

| quantity | max relative difference | sign flips |
|---|---|---|
| `cellScores` | 1.0e-11 | 0 |
| `geneScores` | 5.4e-11 | 0 |
| `geneScoresRegression` | 1.5e-11 | 0 |
| `normalizedCorrelation` | 5.2e-15 | 0 |
| null correlations / p-values | 3.7e-14 | 0 |
| `sigmaValueChoice` | identical | — |
| `@skrCCAOut` (not reported to users) | 2.0 — 4 coordinates of 64, in 3 scenarios | 4 |

`compare_baseline.R` now reports sign flips in reported quantities as a separate
line and fails on any, at any tolerance. A tolerance alone cannot express this
property: a lone flip shows up as a relative difference of 2, which no sensible
tolerance passes, but a *cancelling pair* of flips leaves the tolerance clean
while an intermediate looks catastrophically different. Both are worth knowing
about, so both are reported.

### What this does and does not say about sign determinism

1.1.2 made skrCCA independent of the RNG state: the same input must give the
same answer, run to run, session to session, sequential or PSOCK. That guarantee
is intact — two consecutive captures at branch HEAD are bit-identical across all
14 scenarios.

It never was, and could not be, a guarantee that a PC sign survives a 1-ulp
change in the *input data*. A different BLAS or R build moves it under the old
code too. `test-pca-workflow.R` therefore compares the two implementations
after per-component sign alignment for scores and gene weights, and directly for
the sign-invariant quantities. That is a deliberate choice, not an oversight:
holding the pipeline to bit-identity under a 1-ulp input perturbation would be
asserting something that was never true and that nothing downstream needs.
`test-cca-determinism.R` continues to hold the guarantee that does matter.

The gene involved was also not marginal — `Gene56`, SD rank 34 of 60, 99.4%
nonzero — so "it only moved a negligible weight" is not the explanation. The
explanation is that a PC sign and the CCA coordinate on it are the same
convention seen twice.

`matrixStats` costs 13 ms to load, declares no `Imports` of its own, and was
already a CoPro dependency before this pass.

---

## Measured and rejected

**One Gaussian-weight pass per edge in the kernel builder** — the default
`normalization = 0` path already makes only two passes, and collapsing to one
requires retaining the weights (4 bytes/edge) alongside the edge array
(8 bytes/edge), a permanent 50% increase in the largest temporary. The
four-pass case is `normalization != 0` only, which is not the default. Wrong
trade for a pass whose purpose is bounding memory.

**Pipelining enumeration against concatenation** — see section 2. The bound it
would remove does not appear in measured RSS.
