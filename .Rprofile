# source("renv/activate.R")

# Automatically load development libraries when opening the CoPro project
if (interactive()) {
  # List of libraries to load automatically
  dev_libs <- c("devtools", "usethis", "testthat", "roxygen2")
  
  # Function to safely load libraries
  load_dev_libs <- function(libs) {
    for (lib in libs) {
      if (requireNamespace(lib, quietly = TRUE)) {
        suppressPackageStartupMessages(library(lib, character.only = TRUE))
        cat("✓ Loaded", lib, "\n")
      } else {
        cat("⚠ ", lib, "not available - install with install.packages('", lib, "')\n", sep = "")
      }
    }
  }
  
  # Load the libraries
  cat("Loading development libraries for CoPro project...\n")
  load_dev_libs(dev_libs)
  cat("Development environment ready! 🚀\n")
  
  # Optional: Set some useful options for development
  options(
    repos = c(CRAN = "https://cran.rstudio.com/"),
    browserNLdisabled = TRUE,
    deparse.max.lines = 2
  )
}
