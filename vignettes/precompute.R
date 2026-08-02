# Pre-compute expensive vignettes locally.
# Run from the repo root:
#   Rscript vignettes/precompute.R
#
# Knits .Rmd.orig files into .Rmd with output + figures embedded,
# so pkgdown renders them without re-executing the code.

library(knitr)

# Render against the SOURCE TREE, not whatever CoPro happens to be installed.
#
# Every .Rmd.orig calls library(CoPro). Without the load_all() below that
# resolves to the installed library, so the vignettes get rendered by an old
# build and their numbers silently stop describing this package. That is not
# hypothetical: the committed outputs were produced by CoPro 1.1.2 while the
# source tree was at 1.3.0, so they predated the whole 1.2.0 distance-unit
# rework (`normalizeMethod`, `detectSigmaRange`, the `normalizeDistance`
# default) and nothing in the workflow flagged it.
#
# load_all() attaches CoPro to the search path, so the library(CoPro) inside
# each vignette is a no-op and the dev version is what runs.
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("precompute.R needs devtools to load the source tree. ",
       "install.packages(\"devtools\")")
}
devtools::load_all(".", quiet = TRUE)

# Fail loudly rather than render a vignette with the wrong package. A silent
# version mismatch is exactly the failure this script existed to hide.
source_version <- as.character(read.dcf("DESCRIPTION", fields = "Version"))
loaded_version <- as.character(utils::packageVersion("CoPro"))
if (!identical(source_version, loaded_version)) {
  stop(sprintf(
    paste0("CoPro %s is loaded but DESCRIPTION says %s. The vignettes would be ",
           "rendered by the wrong build. Restart R and re-run from the repo root."),
    loaded_version, source_version
  ))
}
message("Rendering vignettes with CoPro ", loaded_version, " (source tree).")

orig_files <- list.files("vignettes", pattern = "\\.Rmd\\.orig$",
                         full.names = FALSE)

old_wd <- setwd("vignettes")
on.exit(setwd(old_wd))

for (orig in orig_files) {
  out <- sub("\\.orig$", "", orig)
  message("Knitting: ", orig, " -> ", out)
  knit(orig, output = out, envir = new.env(parent = globalenv()))
  message("Done: ", out)
}

# Note on figure paths: each .Rmd.orig sets `fig.path` in its own setup chunk,
# and it must point into `figures/` with a per-vignette prefix. Setting it here
# would not work -- the document's setup chunk runs later and wins -- so the
# convention is enforced there. Two reasons it matters:
#
#   * The knitr default `<base>_files/` is treated by `R CMD build` as vignette
#     build output and stripped from the tarball, so the shipped .Rmd referenced
#     images that were absent and every `R CMD check` rebuild failed with
#     "File <base>_files/plot-x-1.png not found in resource path" (pandoc 99).
#   * All six vignettes share `figures/`, and their chunk labels collide
#     (`plot-layout` exists in both brain_merfish_two_type and
#     colon_d3_cross_type), so the `<base>-` prefix stops them overwriting each
#     other's figures.
