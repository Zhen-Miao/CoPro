# Sparse distance/kernel precision benchmark

## Implemented multi-slide and one-type extension

The float32 backend now supports `CoProSingle` and `CoProMulti` objects with
any positive number of cell types. A one-type analysis stores only the strict
upper triangle and applies it with symmetric CSR semantics. Multi-slide
construction schedules the largest slide/cell-type block first and releases
its neighbor temporary before continuing.

An additional production-API benchmark used 200,000 tiled colon cells split
across eight slides, three cell types, and five sigma values. This produced
120 kernels with 239,519,800 represented nonzeros. Both modes ran in fresh
processes using the public construction functions and an optimized installed
package.

| Metric | Current float64 sparse path | Streamed float32 path | Improvement |
|---|---:|---:|---:|
| Peak RSS | 8.86 GiB | 2.34 GiB | 73.6% lower |
| All 120 kernels | 2.68 GiB | 1.96 GiB | 27.0% lower |
| Build all kernels | 25.539 s | 12.207 s | 2.09x faster |

Every kernel had the same represented nonzero count in both modes. Across all
120 kernels, the maximum relative error in the sum of values was `1.74e-7`
and the maximum relative error in the sum of squared values was `4.11e-7`.

For compatibility, `getKernelMatrix()` and `getSelfKernelMatrix()` materialize
an encoded kernel as a temporary standard sparse matrix by default. Internal
CCA, normalized-correlation, permutation, and score paths explicitly retain
the encoded representation. `materializeFloat32Kernels()` converts the entire
kernel list for third-party code that reads the slot directly; it is
deliberately documented as a high-memory escape hatch.

Global, row, and column kernel normalization are applied directly in float32.
Row/column-normalized self-kernels are expanded to general CSR because the
result is asymmetric; unnormalized and globally normalized self-kernels keep
the single-triangle representation. Regression comparisons against the
float64 sparse path agree within `2e-6`.

The ordinary centered Frobenius objective norm is exact for the stored float32
values and does not decode the kernel: it uses row sums, column sums, and the
sum of squared entries. Native exact whitened Frobenius normalization is
deferred because it requires a substantially more complex sparse-sparse
contraction equivalent to `Rx %*% K %*% Ry`. When within-type whitening
operators are supplied, the compatibility path temporarily materializes the
encoded kernels as double sparse matrices and evaluates the existing exact
formula.

## Implemented 200,000-cell result

The recommended float32 design has now been implemented as an experimental
CoPro backend:

- `computeSparseKernelFloat32()` processes the largest cell-type block first
  and releases each block's neighbor data before moving to the next;
- distances are calculated in double and retained temporarily as float32;
- sigma-specific kernels are emitted directly as float32 bytes in CSR format,
  without constructing a float64 `Matrix`;
- the parallel CSR operator partitions rows between workers, forms
  `X1' K X2` using thread-local float32 `nPC x nPC` accumulators, and returns
  only the small result in float64;
- no `n_cells x nPC` `K %*% X2` temporary is allocated.

The 200,000-cell workload tiles the real colon D3 coordinates, cell-type
proportions, and 40-PC matrices into spatially separated tissue regions. Gaps
are larger than the maximum kernel support radius, so local neighbor density
and kernel values are preserved without artificial cross-tile neighbors. The
benchmark contains all three cell types, five sigma values, 15 kernels, and
137,831,652 nonzeros in the largest kernel.

Both modes ran in fresh R processes using an optimized `R CMD INSTALL` build.
Peak RSS was measured by `/usr/bin/time -l`.

| Metric | Current float64 sparse path | Streamed float32 path | Improvement |
|---|---:|---:|---:|
| Peak RSS | 8.37 GiB | 2.83 GiB | 66.2% lower |
| All 15 kernels | 2.75 GiB | 1.84 GiB | 33.1% lower |
| Neighbor temporary | 2.54 GiB | 1.09 GiB | 57.3% lower |
| Build all 15 kernels | 27.244 s | 9.745 s | 2.80x faster |
| Largest `X1' K X2` | 2.184 s | 0.478 s | 4.57x faster |

The 8-thread operator scales as follows:

| Threads | Median `X1' K X2` |
|---:|---:|
| 1 | 3.135 s |
| 2 | 1.576 s |
| 4 | 0.821 s |
| 8 | 0.478 s |

The single-thread float32 operator is slower than Matrix's optimized float64
multiply; the measured speed win comes from the race-free CSR row
parallelization as well as the smaller representation.

Numerical stability remained high:

- relative Frobenius error in the 40 x 40 `X1' K X2`: `1.235e-6`;
- worst sign-aligned NRMSE across the first four Epithelial/Fibroblast cell
  score components: `1.173e-6`;
- all eight cell-score correlations rounded to `1.000000000`.

This makes the float32 path worth pursuing for large single- and multi-slide
workloads. It supports one or more cell types, global/row/column normalization,
`runSkrCCA()`, `computeNormalizedCorrelation()`,
`computeGeneAndCellScores()`, transposed kernel access, and gene-space kernel
application. Transfer, plotting, bidirectional-correlation, and third-party
consumers that require a standard `Matrix` use the documented temporary
materialization path. Exact whitened Frobenius normalization likewise remains
on that compatibility path rather than using a native float32 sparse-sparse
operator.

## Decision

Keep distance calculation and the temporary fixed-radius triplets at their
current precision. Do not use 100 percentile bins for distances.

For persistent kernels:

1. Use float32 storage when final-score stability is the overriding
   requirement. Apply the sparse kernel in float32, cast the resulting small
   `nPC x nPC` `X1' K X2` matrix to float64, and keep all optimization,
   deflation, SVD/eigendecomposition, and cell-score construction in float64.
2. Float16 kernel storage is a reasonable opt-in, higher-compression mode. On
   this workload, float32 accumulation added no meaningful error beyond the
   float16 value rounding.
3. Do not decode the entire kernel back into an R `Matrix` before
   `X1' K X2`. That recreates an 8-byte value array, loses the working-memory
   benefit, and was slower.
4. Do not use 100 percentile bins for kernels when scores must be "very very
   close." It saved only another 9 percentage points of total sparse storage
   beyond float16 but increased worst score NRMSE by about 59-fold.

The current R `Matrix::dgCMatrix` kernel values are float64 (`typeof(K@x)` is
`"double"`), not float32. The fixed-radius sparse path persists zero distance
matrices; distances exist only in temporary triplets during kernel
construction.

## Peak memory and scaling

The first benchmark table below reported persistent object sizes. An isolated
peak-RSS measurement was therefore added: each fresh R process loaded all 15
kernels and the two PC matrices, then calculated the largest
`X1' K X2`.

| Representation | Maximum RSS | Reduction |
|---|---:|---:|
| Current float64 kernels and float64 `X1' K X2` | 353.89 MiB | — |
| float16 kernels and direct float32 `X1' K X2` | 269.80 MiB | 23.76% |

At this small size, approximately 180 MiB of R/runtime/package overhead is
unaffected by the representation. The live R vector heap immediately before
the multiply was 188.1 MiB versus 105.4 MiB, a difference matching the 82.7
MiB saved by the kernel list. The relative RSS benefit grows when kernels
become the dominant object.

The following projections scale the observed colon cell-type proportions and
neighbors per cell linearly. They exclude normalized expression data,
metadata, allocator overhead, and unrelated R objects, so they describe the
distance/kernel/PC portion of peak memory rather than total application RSS.

| Cells | Current `X1' K X2` working set | float32 kernel + float32 compute | float16 kernel + float32 compute |
|---:|---:|---:|---:|
| 200,000 | 2.88 GiB | 1.93 GiB | 1.47 GiB |
| 400,000 | 5.75 GiB | 3.86 GiB | 2.94 GiB |

Here the working set includes all 15 stored kernels, all 40-PC matrices, and
the largest `K %*% X2` temporary. At 400,000 cells:

- float16 kernel storage saves approximately 2.77 GiB;
- making `K %*% X2` float32 saves another 0.048 GiB;
- therefore kernel storage is the major gain, while float32 arithmetic is a
  useful but much smaller peak-memory gain.

This ratio follows from sparse-kernel memory being approximately
`nnz * (4-byte row index + value bytes)`, whereas the dense multiply temporary
is only `n_rows * nPC * element_bytes`. With hundreds of neighbors per cell
and 40 PCs, the sparse kernels dominate.

The projection assumes larger tissue area at approximately constant cell
density, so neighbors per cell stay stable and both memory and time are
`O(n)`. If cell count increases within the same physical area, neighbors per
cell also increase; `nnz` can approach `O(n^2)`. Compression then saves more
absolute memory but cannot prevent the calculation itself from eventually
becoming infeasible.

Peak memory will only fall if encoded kernels are created directly and remain
encoded during multiplication. Building ordinary double `dgCMatrix` kernels
and encoding them afterward temporarily holds both copies and can increase
peak memory. For large data, kernel construction should instead:

1. process one cell-type/slide block at a time, preferably largest first;
2. retain temporary distances as float32, which was numerically indistinguishable
   here, rather than float16;
3. generate float16/float32 encoded kernel values directly for every sigma;
4. release that block's neighbor triplets before moving to the next block.

With the current all-block triplet cache, the modeled distance/kernel build
working set at 400,000 cells is 10.65 GiB. Direct float16 kernels plus float32
distance triplets reduce it to 6.61 GiB before applying the additional
block-at-a-time scheduling improvement.

## Workload

- Apple arm64, R 4.5.2, Matrix 1.7.5
- Colon D3 vignette data: 11,666 cells
- Cell types: Epithelial, Fibroblast, Immune
- 40 PCs, 4 final components
- Sigma: 0.005, 0.01, 0.02, 0.05, 0.1
- 15 cross-cell-type sparse kernels
- Kernel support threshold: `1e-7`
- Current fixed-radius sparse kernel build: 3.673 seconds

All final-score metrics are the worst case over 5 sigma values, 3 cell types,
and 4 components (60 score vectors). Scores were sign-aligned because a
canonical component and its negative are equivalent. NRMSE is RMSE divided by
the baseline score standard deviation.

## Accuracy

| Variant | Max relative `X1' K X2` error | Min cell-score correlation | Max cell-score NRMSE | Max absolute score error |
|---|---:|---:|---:|---:|
| Kernel float32, float64 accumulation | 4.11e-8 | 1.000000000 | 2.94e-8 | 1.87e-7 |
| Kernel float32, float32 accumulation | 1.70e-6 | 1.000000000 | 2.58e-6 | 1.40e-5 |
| Kernel float16, float64 accumulation | 2.74e-4 | 0.999999950 | 3.15e-4 | 1.59e-3 |
| Kernel float16, float32 accumulation | 2.73e-4 | 0.999999951 | 3.14e-4 | 1.59e-3 |
| Kernel 100 percentile bins, float64 accumulation | 1.44e-2 | 0.999828739 | 1.85e-2 | 1.16e-1 |
| Distance float32, float64 kernel/compute | 1.22e-7 | 1.000000000 | 1.31e-7 | 6.49e-7 |
| Distance float16, float64 kernel/compute | 1.46e-3 | 0.999996440 | 2.66e-3 | 1.46e-2 |
| Distance float16 + kernel float16 | 1.63e-3 | 0.999995947 | 2.84e-3 | 1.56e-2 |
| Distance 100 percentile bins | 1.00e0 | 0.012531603 | 1.41e0 | 8.10e0 |

Float16 distance error is amplified by
`exp(-0.5 * (distance / sigma)^2)`. Equal-frequency distance bins are worse:
bins are determined from all neighbors needed by the largest sigma, so the
first 1% bin is far too coarse for the smallest sigma. At sigma 0.005, one
matrix changed from 8,833 to 15,867 nonzeros and its relative kernel error was
approximately 1.

## Memory

Measured `object.size()` includes compressed-column indices and R object
overhead, not only values.

| Object | Encoding | MiB | Reduction from float64 |
|---|---|---:|---:|
| All persistent kernels | Current float64 `dgCMatrix` | 165.38 | — |
| All persistent kernels | float32 bytes | 110.27 | 33.32% |
| All persistent kernels | float16 bytes | 82.72 | 49.98% |
| All persistent kernels | 100 one-byte bins + codebooks | 68.96 | 58.30% |
| Temporary neighbor triplets | Current float64 distances | 152.73 | — |
| Temporary neighbor triplets | float32 distances | 114.55 | 25.00% |
| Temporary neighbor triplets | float16 distances | 95.46 | 37.50% |
| Temporary neighbor triplets | 100 one-byte bins + codebooks | 85.92 | 43.75% |

Distance reductions are smaller because the two 32-bit index vectors remain
unchanged. More importantly, the current method does not retain these
distances after kernel construction, so distance encoding cannot reduce the
final object size.

## Speed

Median `X1' K X2` time for the largest cell-type pair:

| Sigma (nonzeros) | Current Matrix float64 | float16 direct, float32 accumulation | float16 decode then Matrix float64 |
|---|---:|---:|---:|
| 0.01 (185,798) | 0.004 s | 0.021 s | 0.005 s |
| 0.05 (2,732,143) | 0.057 s | 0.125 s | 0.066 s |
| 0.10 (8,090,652) | 0.131 s | 0.329 s | 0.163 s |

The original direct encoded operator was a correctness/measurement prototype,
not a parallel production sparse backend. It showed that encoding alone is
not a speed optimization: decoding and less optimized sparse loops outweighed
the lower memory traffic. The implemented CSR backend above supersedes that
prototype by eliminating the dense intermediate and parallelizing over
disjoint rows.

Rebuilding all 15 kernels from already-enumerated triplets took 1.176 seconds
with float64 distances, 1.148 seconds after float16 decode, and 1.215 seconds
after 100-bin decode. These differences are negligible relative to the full
3.673-second fixed-radius path and do not justify persistent distance storage.

The 200k result demonstrates the needed optimized path. With constant
neighbors per cell, both implementations remain `O(nnz * nPC)`, but the
float32 CSR implementation has a smaller constant, parallel row ownership,
and substantially lower allocation/copy pressure.

## Recommended compute boundary

For an optimized encoded implementation:

```text
uint16/float16 K values in compressed sparse storage
    -> decode each value to float32 inside sparse K %*% X2
    -> accumulate K %*% X2 in float32
    -> accumulate X1' %*% (...) in float32
    -> cast the small 40 x 40 Y matrix to float64
    -> optimize/deflate/SVD and form final cell scores in float64
```

The benchmark directly tested this boundary. For float16 kernels, float32
versus float64 accumulation changed worst score NRMSE from `3.146e-4` to
`3.139e-4`; storage rounding, not accumulation, dominated the error.

## Reproduction

Run from the repository root:

```sh
Rscript reports/kernel_precision_benchmark/benchmark_kernel_precision.R
```

Detailed results are in:

- `accuracy_summary.csv`
- `cell_score_accuracy.csv`
- `kernel_accuracy.csv`
- `y_accuracy.csv`
- `storage_memory.csv`
- `xky_speed.csv`
- `distance_kernel_build_speed.csv`
- `kernel_encoding_speed.csv`
- `run_metadata.csv`
- `peak_memory_observed.csv`
- `peak_memory_projection.csv`
- `large_200k_summary.csv`
- `large_200k_float32_thread_scaling.csv`
- `large_200k_accuracy.csv`
- `large_200k_cell_score_accuracy.csv`

`encoded_sparse_xky.cpp` contains the float16/float32 byte encoders and the
direct compressed-sparse `X1' K X2` prototype used by the benchmark.
`benchmark_large_float32.R` constructs and measures the 200k workload.
