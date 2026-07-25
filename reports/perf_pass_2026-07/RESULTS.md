# Performance and memory pass, July 2026

Measured by A/B-ing two checkouts of this repository — the pre-pass commit
`e5ddaa4` and the branch head — driven by the scripts in this folder. Session
machine: Darwin 25.1.0, Apple silicon, R 4.5.2, clang `-O2`.

Numeric equivalence is verified separately by
`reports/perf_pass_baseline/compare_baseline.R`, which re-runs 14 pipeline
scenarios and diffs every downstream result. Thirteen are bit-identical; the
fourteenth differs only by the 3.9e-15 normalizer reassociation documented in
section 3.

---

## Headline: end-to-end pipeline

`benchmark_end_to_end.R` → `end_to_end_before.csv`, `end_to_end_after.csv`

60,000 cells, 2 cell types, 120 genes, `nPCA = 30`, 2 sigmas, 99 permutation
draws. `peak MB` is Vcells allocated above the resting set during that stage.

| stage | before | after | speedup | peak before | peak after |
|---|---|---|---|---|---|
| `computePCA` | 0.95 s | 0.95 s | 1.00x | 273 MB | 271 MB |
| `computeSparseKernelFloat32` | 43.47 s | 33.20 s | **1.31x** | 1575 MB | 1576 MB |
| `runSkrCCA` | 4.16 s | 3.04 s | **1.37x** | 31 MB | 31 MB |
| `computeNormalizedCorrelation` | 6.63 s | 3.62 s | **1.83x** | 39 MB | 38 MB |
| `computeGeneAndCellScores` | 0.04 s | 0.04 s | 1.02x | 77 MB | 77 MB |
| `runSkrCCAPermu` (99 draws) | 3.98 s | 3.88 s | 1.03x | 786 MB | 785 MB |
| `computeNormalizedCorrelationPermu` | 6.43 s | 3.81 s | **1.69x** | 1271 MB | **824 MB** |
| `@cellPermu` stored size | — | — | — | 24 MB | **12 MB** |
| **total** | **65.7 s** | **48.5 s** | **1.35x** | | |

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
shared HPC node that matters as much as the mean.

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

---

## 3. Whitened-Frobenius normalizer

`benchmark_whitened_frobenius.R`

`(Rx %*% K) %*% Ry` fills in 11x, and only to produce a scalar. Streaming
`sum((Rx K) * (K Ry))` over column blocks is both faster and smaller — the
intermediates stay in cache instead of streaming a hundreds-of-MB product
through memory.

| n = 150k, nnz(K) 3.91M | seconds | peak |
|---|---|---|
| unblocked | 5.97 | 1996 MB |
| `block_nnz = 2e6` (2 blocks) | 5.11 | 1525 MB |
| `block_nnz = 1e6` (4 blocks) | 4.81 | 871 MB |
| **`block_nnz = 2e5` (20 blocks)** | **4.09** | **466 MB** |
| `block_nnz = 5e4` (79 blocks) | 4.27 | 467 MB |
| `block_nnz = 2e4` (196 blocks) | 5.64 | 467 MB |

1.46x faster, 4.3x less peak. A budget large enough to yield a single block is
worse than not blocking (both one-sided products live at once); very small
blocks lose to per-block overhead. The default sits at the knee.

This is the one change that moves a number: floating-point reassociation shifts
the normalizer by 3.9e-15 relative. Weights, cell scores, gene scores and the
selected sigma stay bit-identical.

**This corrects the plan's estimate**, which predicted blocking would be
*slightly slower* and worth taking for memory alone. That was based on a
badly-sized trial block; at a sensible block size it wins on both.

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

---

## `matrixStats::colSds` — rejected, then adopted after a closer look

Per-column SDs use `matrixStats::colSds()`, which is ~2.5x faster than
`apply(x, 2, sd)`.

This was initially backed out, on a mistaken reading. The two use different
variance algorithms and disagree by 1 ulp on some columns (1.11e-16 on 1 of 60
genes in the fixture); switching produced a relative difference of exactly 2 in
`@skrCCAOut`, which is the signature of a sign flip. I reported that as "flips
the sign of a PC, and with it the sign of the gene weights read off that PC"
and rejected the change on those grounds.

The second half of that was wrong, and the downstream numbers I had already
collected contradicted it. Decomposing the difference:

- Exactly one principal component flipped (CellTypeB's PC2). All others were
  identical to ~1e-13.
- The CCA weight's *coordinate on that PC* flipped with it. That is what
  produced the relative difference of 2 — one coordinate, not a negated vector.
- The two cancel. Every reported quantity was unchanged:

  | output | max difference | sign flips |
  |---|---|---|
  | `cellScores` | 5.2e-14 | none |
  | `geneScores` | 5.0e-15 | none |
  | `geneScoresRegression` | 7.6e-15 | none |
  | `normalizedCorrelation` | 4.7e-15 | — |

The gene involved was also not marginal — `Gene56`, SD rank 34 of 60, 99.4%
nonzero — so "it only moved a negligible weight" is not the explanation
either. The explanation is that a PC sign and the CCA coordinate on it are the
same convention seen twice.

Re-measured at the final branch state, the two implementations are
**bit-identical end-to-end** on this fixture: no PC sign flips, and zero
difference in cell scores, gene scores, regression gene scores and normalized
correlations. The full 14-scenario baseline is likewise bit-identical.

The flip observed earlier in the session is not reproducible at the final
state, which is itself the point: **the sign of a PC or CCA axis is
knife-edge** and a different BLAS, R build or 1-ulp input change can move it
under any implementation. That fragility is pre-existing and not something
`colSds` introduces. What matters is that everything downstream is invariant to
it. `test-pca-workflow.R` now runs the whole pipeline under both
implementations and asserts that invariance — sign-invariant outputs
(normalized correlation, selected sigma) compared directly, scores and gene
weights compared after sign alignment — rather than relying on bit-identity,
which would be the wrong thing to depend on.

`matrixStats` costs 13 ms to load, declares no `Imports` of its own, and was
already a CoPro dependency before this pass.

## Measured and rejected

**One Gaussian-weight pass per edge in the kernel builder** — the default
`normalization = 0` path already makes only two passes, and collapsing to one
requires retaining the weights (4 bytes/edge) alongside the edge array
(8 bytes/edge), a permanent 50% increase in the largest temporary. The
four-pass case is `normalization != 0` only, which is not the default. Wrong
trade for a pass whose purpose is bounding memory.
