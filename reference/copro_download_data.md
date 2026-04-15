# Download example datasets for CoPro vignettes

Downloads pre-processed, subsampled datasets from GitHub Releases for
use in vignettes and tutorials. The datasets are attached to GitHub
releases via the `piggyback` package.

## Usage

``` r
copro_download_data(
  dataset = c("colon_d3", "colon_d9", "kidney", "organoid", "brain_merfish"),
  destdir = NULL,
  tag = "data-v1",
  overwrite = FALSE
)
```

## Arguments

- dataset:

  Character string specifying which dataset to download. One of
  `"colon_d3"`, `"colon_d9"`, `"kidney"`, `"organoid"`, or
  `"brain_merfish"`.

- destdir:

  Directory to save the downloaded file. Defaults to a package-specific
  cache directory via
  [`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html).

- tag:

  The GitHub release tag to download from. Defaults to `"data-v1"`.

- overwrite:

  Logical; if `TRUE`, re-download even if the file already exists
  locally. Default `FALSE`.

## Value

The file path to the downloaded RDS file (invisibly).

## Details

Available datasets:

- `colon_d3`:

  Colon Day 3 organoid data (Epithelial, Fibroblast, Immune).
  Demonstrates cross-cell-type co-progression with orthogonal CCA axes.

- `colon_d9`:

  Colon Day 9 organoid data (multiple slides). Demonstrates multi-slide
  analysis and score transfer.

- `kidney`:

  Kidney seqFISH data (tubular and vascular cells). Demonstrates
  supervised/guided spatial gradient detection.

- `organoid`:

  72hr organoid culture (single cell type). Demonstrates
  within-cell-type spatial pattern detection.

- `brain_merfish`:

  Brain MERFISH data (D1/D2 neurons). Demonstrates two-cell-type
  co-progression.

The data files are hosted as GitHub Release assets and are typically
5–30 MB each. They are subsampled from the full datasets to allow fast
vignette execution while preserving biological signal.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download and load the colon D3 dataset
path <- copro_download_data("colon_d3")
dat <- readRDS(path)
} # }
```
