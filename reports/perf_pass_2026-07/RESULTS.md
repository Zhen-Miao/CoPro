# Performance and memory pass, July 2026

Measurements taken on the session machine (Darwin 25.1.0, Apple silicon,
R 4.5.2, clang `-O2`) by A/B-ing two checkouts of this repository: the pre-pass
commit `e5ddaa4` and the branch head. Both were driven by the same scripts in
this folder.

Numeric equivalence is verified separately by
`reports/perf_pass_baseline/compare_baseline.R`, which re-runs 14 pipeline
scenarios and diffs every downstream result.

---

## 1. Row-major packing of the float32 operands

`benchmark_float32_operators.R` -> `operators_strided.csv`, `operators_packed.csv`

Kernels built from uniform coordinates on a 1000 x 1000 field, giving the
"hundreds of neighbours per cell" density that
`reports/kernel_precision_benchmark/RESULTS.md` reports for real data.
Seconds per call, minimum of 3 x 3 reps.

| case | nnz/row | nPC | threads | op | strided | packed | speedup |
|---|---|---|---|---|---|---|---|
| cross-type | 354 | 30 | 1 | `kernelXKY` | 16.512 | 4.658 | **3.54x** |
| cross-type | 354 | 30 | 4 | `kernelXKY` | 3.978 | 1.198 | **3.32x** |
| cross-type | 354 | 30 | 4 | `matMult`   | 3.270 | 1.084 | **3.02x** |
| cross-type | 354 | 10 | 1 | `kernelXKY` | 2.280 | 1.493 | 1.53x |
| cross-type | 354 | 10 | 4 | `kernelXKY` | 0.503 | 0.373 | 1.35x |
| cross-type | 354 | 10 | 4 | `matMult`   | 0.415 | 0.326 | 1.27x |
| within-type | 374 | 30 | 1 | `kernelXKY` | 20.705 | 13.283 | 1.56x |
| within-type | 374 | 30 | 4 | `kernelXKY` | 10.485 | 6.314 | 1.66x |
| within-type | 374 | 10 | 4 | `kernelXKY` | 2.032 | 1.492 | 1.36x |
| within-type | 374 | 10 | 4 | `matMult`   | 1.033 | 1.003 | 1.03x |

**Headline: 3.0-3.5x on the cross-type operators at nPC = 30, tapering to
~1.3-1.6x at nPC = 10 or on the within-type (symmetric) kernel.** The gain
scales with nPC because that is how many separate cache lines the old
column-major gather touched per nonzero.

`matMult` is unchanged (1.03x) on the within-type kernel because symmetric
input takes the transposed branch, which splits over right-hand sides and
already reads contiguously; only the general branch was repacked.

### A caveat on the first sweep

The sweep's first measured cell (within-type, nPC = 10, 1 thread) initially
reported the packed build as *slower* (0.73x). `benchmark_symmetric_variance.R`
re-measures that cell alone, with a warm-up and six reps:

```
strided  p=10 threads=1 : min 6.886  median 8.575   all [6.89 8.04 8.54 8.61 8.78 9.83]
packed   p=10 threads=1 : min 4.870  median 5.455   all [4.87 5.32 5.38 5.53 6.40 6.83]
strided  p=30 threads=1 : min 28.532 median 44.254  all [28.53 36.51 42.69 45.82 46.99 54.81]
packed   p=30 threads=1 : min 16.608 median 16.685  all [16.61 16.67 16.68 16.69 16.71 16.98]
```

There is no regression: the packed build wins on both min and median. The
apparent one was cold-start noise in the first cell of the sweep.

The spread is itself a result. The strided build ranges 28.5-54.8 s on the same
input while the packed build holds 16.6-17.0 s. A loop that touches nPC cache
lines per nonzero is memory-bandwidth-bound and therefore hostage to whatever
else is running; the packed loop fetches roughly nPC times less and stays
predictable. On a shared HPC node that stability is worth as much as the mean
speedup.

### Correction to the pre-implementation estimate

A standalone microbenchmark written before the change predicted 4.0-6.6x. That
harness used 40 nonzeros per row; real kernels have ~350-375, and at that
density more of the working set stays resident regardless of layout. The
in-package numbers above are the ones to quote.
