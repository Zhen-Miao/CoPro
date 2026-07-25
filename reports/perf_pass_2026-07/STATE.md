# Performance pass — state of play

Branch: `perf/memory-and-throughput-pass`, branched from `main` at `c7223bd`.

Everything below is committed. Nothing in this pass depends on network access:
`devtools::load_all()`, compilation, the test suite, `R CMD check` and the
baseline harness are all local.

## Status: all eight planned items landed

| # | Item | Commit | Numerics |
|---|---|---|---|
| 1 | Row-major float32 operands | `e5bc8c3` | bit-identical |
| 3 | Compact held-fixed permutations | `664fb50` | bit-identical |
| 5 | Repeated-value type-7 quantile | `6466e5e` | bit-identical |
| 4 | Blocked whitened-Frobenius normalizer | `290a1ba` | normalizer 3.9e-15; scores/weights/sigma bit-identical |
| 6 | Cached kernel content signature | `dfa460d` | bit-identical |
| 7 | Call-site fixes (kernel decode, nonzero count, centering) | `5ffc001` | bit-identical |
| 8 | Per-slide PC scores as views | `d9db60a` | bit-identical |
| 2 | Parallel enumeration + flat grid buckets | `280d899` | bit-identical |

Plus `e5ddaa4` (baseline harness) and `c71d094` (Rd link fix).

## Verification already run

- Full suite: **1382 passing, 0 failures, 0 errors**, 10 skips (all
  pre-existing: uninstalled Suggests, and the PSOCK path that needs an
  installed CoPro).
- `R CMD check`: **`checking tests ... OK`**, 3 WARNINGs, all pre-existing or
  invocation artifacts — non-portable names under the untracked `poster/`
  directory, and two `inst/doc` warnings caused by building with
  `--no-build-vignettes`. Needs `_R_CHECK_FORCE_SUGGESTS_=false` locally
  because BPCells, RcppAnnoy, SingleCellExperiment, viridis and the spatstat
  packages are not installed here.
- `reports/perf_pass_baseline/compare_baseline.R`: 14 scenarios, 13
  bit-identical and 1 (`within_type_float32`) differing only by the documented
  3.9e-15 normalizer reassociation from item 4.

## To resume or re-verify

```sh
Rscript reports/perf_pass_baseline/compare_baseline.R        # bit-identical bar
Rscript reports/perf_pass_baseline/compare_baseline.R 1e-10  # tolerance bar
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'
```

Snapshot tests need `NOT_CRAN=true`, otherwise `expect_snapshot_value()` skips.

## Remaining work

Only the write-up of `RESULTS.md` and the "before" half of the end-to-end
A/B. The "after" numbers are in `end_to_end_after.csv`. To produce the
"before" half, check out the pre-pass commit in a worktree and run the same
script there:

```sh
git worktree add /tmp/copro_base e5ddaa4
cd /tmp/copro_base && Rscript reports/perf_pass_2026-07/benchmark_end_to_end.R before <outdir>
```

## Deliberately not done

- Computing the Gaussian weight once per edge in the kernel builder — the
  default `normalization = 0` path already makes only two passes, and
  collapsing to one costs a permanent ~50% increase in the largest temporary.
  Rationale and measurements are in the `280d899` commit message.
