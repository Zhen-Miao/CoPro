# Contributing to CoPro

Thanks for your interest in CoPro. This guide covers the minimum you
need to know to submit a well-formed pull request.

## Development setup

CoPro uses `renv` for pinned dependencies and expects R \>= 4.1.

``` r

install.packages("renv")
renv::restore()
devtools::load_all()
```

## Running tests

Tests use `testthat` 3rd edition and only depend on lightweight
synthetic data (no network, no piggyback downloads):

``` r

devtools::test()
```

The CI `R-CMD-check` workflow runs the same tests plus `R CMD check`
with vignette builds disabled (see below).

## R CMD check (CI-parity)

The GitHub Actions workflow invokes:

``` r

rcmdcheck::rcmdcheck(
  build_args = c("--no-build-vignettes"),
  args       = c("--no-vignettes", "--no-manual", "--ignore-vignettes")
)
```

Reproduce locally with the same flags before sending a PR.

## Vignettes

Vignette builds are intentionally **disabled in CI** because rendering
pulls datasets via `piggyback` from GitHub Releases, which is unreliable
on hosted runners. Please do not remove the `--no-build-vignettes` /
`--no-vignettes` flags from `.github/workflows/R-CMD-check.yaml`.

For end-to-end validation before a release, render vignettes locally
(piggyback will cache the data after the first run):

``` r

rmarkdown::render("vignettes/organoid_one_type.Rmd")
rmarkdown::render("vignettes/brain_merfish_two_type.Rmd")
rmarkdown::render("vignettes/colon_d3_cross_type.Rmd")
rmarkdown::render("vignettes/colon_d9_multi_slide.Rmd")
rmarkdown::render("vignettes/kidney_guided_gradient.Rmd")
```

Expect 5-15 minutes per vignette on a laptop-class machine.

## Regenerating documentation

Every change that touches roxygen blocks must be followed by:

``` r

devtools::document()
```

which regenerates `man/*.Rd` and `NAMESPACE`. Commit the regenerated
files together with the source change so the package builds cleanly from
a fresh clone.

## Linting

The lint workflow runs `lintr::lint_package()` with
`LINTR_ERROR_ON_LINT=true`. To reproduce locally:

``` r

lintr::lint_package()
```

A PR should not introduce new lints beyond the existing baseline.

## Commits and pull requests

- Branch from `main` with a short descriptive name.
- Keep each PR focused on one topic. Mixing unrelated refactors in a
  single PR slows review and makes reverts painful.
- Update `NEWS.md` under a “CoPro X.Y.Z (development)” heading for any
  user-visible change (new function, new argument, changed behavior, bug
  fix).
- If you add a new exported function, include at least one runnable
  `@examples` block that uses `inst/extdata/toy_copro_data.rds` (or an
  in-memory synthetic dataset) so `devtools::run_examples()` covers it.

## Reporting bugs

Please open an issue at <https://github.com/Zhen-Miao/CoPro/issues> with
a minimal reproducible example (ideally one that runs against the toy
dataset bundled in `inst/extdata/`).
