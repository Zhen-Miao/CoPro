# Pre-compute expensive vignettes locally.
# Run from the repo root:
#   Rscript vignettes/precompute.R
#
# Knits .Rmd.orig files into .Rmd with output + figures embedded,
# so pkgdown renders them without re-executing the code.

library(knitr)

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
