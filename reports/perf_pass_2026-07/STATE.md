# Performance pass — state of play

Branch: `perf/memory-and-throughput-pass`, branched from `main` at `c7223bd`.

Everything below is committed. Nothing in this pass depends on network access:
`devtools::load_all()`, compilation, the test suite, `R CMD check` and the
baseline harness are all local.

Commit hashes here are the ones on the branch as pushed. An earlier revision of
this file quoted pre-rewrite hashes that no longer resolve; the branch was
rewritten with `git filter-branch` to drop a 256 MB benchmark artifact that had
been committed by accident.

## Status: all eight planned items landed, plus one adopted afterwards

| # | Item | Commit | Numerics |
|---|---|---|---|
| 1 | Row-major float32 operands | `e5bc8c3` | bit-identical |
| 3 | Compact held-fixed permutations | `45b047e` | bit-identical |
| 5 | Repeated-value type-7 quantile | `ddeb1d2` | bit-identical |
| 4 | Blocked whitened-Frobenius normalizer | `225fe1d` | normalizer ~3.7e-15; scores/weights/sigma bit-identical |
| 6 | Cached kernel content signature | `ed844db` | bit-identical |
| 7 | Call-site fixes (kernel decode, nonzero count, centering) | `1a52140` | bit-identical |
| 8 | Per-slide PC scores as views | `19f887a` | bit-identical |
| 2 | Parallel enumeration + flat grid buckets | `50e9095` | bit-identical |
| — | `matrixStats::colSds()` adoption | see below | **not** bit-identical; ≤5.4e-11 on reported quantities, no sign flips |

Plus `e5ddaa4` (baseline harness), `6da3072` (Rd link fix), `ef2a57a` /
`b9ce541` (write-ups) and `f5ca41b` (gitignore).

## The two numeric claims, kept separate

They are separate because they have genuinely different strengths, and merging
them would understate the first and overstate the second.

**Items 1–8**, measured against a pre-pass capture at `e5ddaa4` with the same
harness: **13 of 14 scenarios bit-identical**. The exception,
`within_type_float32`, differs only in `normalizedCorrelation` by ~3.7e-15 —
the documented reassociation from item 4. Its cell scores, gene scores and
regression gene weights are bit-identical.

**Adopting `colSds()`** additionally perturbs the column SDs by 1 ulp. That
flips the sign of a principal component and, with it, the CCA weight's
coordinate on that PC. The two cancel:

- `@skrCCAOut` — the raw weight vectors, a per-PC sign convention — differs in
  3 of 14 scenarios.
- Every reported quantity (cell scores, gene scores, regression gene weights,
  normalized correlations, null correlations, p-values) agrees to **5.4e-11
  relative with zero sign changes**, and the selected sigma is identical in all
  14 scenarios.
- Run-to-run reproducibility is untouched: two consecutive captures at branch
  HEAD are bit-identical across all 14 scenarios. The RNG-independence
  guarantee 1.1.2 established is about repeatability for fixed input, which
  this does not weaken. A PC sign was never invariant to a 1-ulp change in the
  input under any implementation.

An earlier commit message for this change claimed it was bit-identical
end-to-end. That was wrong, and `compare_baseline.R` now reports sign flips in
reported quantities explicitly so the claim cannot be made loosely again.

## Verification already run

- Full suite: **1457 passing, 0 failures, 0 errors**, 10 skips (all
  pre-existing: uninstalled Suggests — spatstat.geom, SingleCellExperiment,
  BPCells — and the PSOCK path that needs an installed CoPro).
- `R CMD check`: **`checking tests ... OK`**. Needs
  `_R_CHECK_FORCE_SUGGESTS_=false` locally because BPCells, RcppAnnoy,
  SingleCellExperiment, viridis and the spatstat packages are not installed
  here.
- `reports/perf_pass_baseline/compare_baseline.R 1e-9`: 11 of 14 scenarios
  clean; the 3 that differ do so only in `@skrCCAOut`, with the reported-
  quantity summary line confirming 0 sign flips.

## To resume or re-verify

```sh
# Numeric equivalence. The baseline was captured at e5ddaa4 with the current
# capture script; re-capture it there if you change what the script records.
Rscript reports/perf_pass_baseline/compare_baseline.R 1e-9

NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'
```

Snapshot tests need `NOT_CRAN=true`, otherwise `expect_snapshot_value()` skips.

The A/B benchmarks need a second checkout:

```sh
git worktree add /tmp/copro_base e5ddaa4
cp reports/perf_pass_2026-07/benchmark_*.R /tmp/copro_base/reports/perf_pass_2026-07/
cd /tmp/copro_base && Rscript reports/perf_pass_2026-07/benchmark_end_to_end.R before <outdir> 3
```

Copy the scripts across rather than using the ones committed at `e5ddaa4`: the
current versions repeat their timed arms, and the pre-pass copies do not.

## Deliberately not done

- Computing the Gaussian weight once per edge in the kernel builder — the
  default `normalization = 0` path already makes only two passes, and
  collapsing to one costs a permanent ~50% increase in the largest temporary.
  Rationale and measurements are in the `50e9095` commit message.
- Eliminating the transient second copy of the edge array when enumeration runs
  on more than one thread. It would need the enumeration and the concatenation
  to be pipelined, which is a real concurrency rewrite, and the measurement
  does not justify it: the bound rises 44% but peak RSS moves 0.1%, because
  reserved pages are only made resident as they are written and the private
  buffers are freed as they are consumed. See `enumeration_memory.csv`.
