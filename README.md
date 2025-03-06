
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CoPro

<!-- badges: start -->
<!-- badges: end -->

<img src="man/figures/copro-refined-logo-final.jpg" width="400" />

The goal of CoPro is to identify co-progression between cell types, in
either supervised or unsupervised manner. In the presence of tissue
structures, CoPro will also be able to identify tissue
structure-associated cellular programs.

## Installation

You can install the current version of CoPro from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Zhen-Miao/CoPro")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CoPro)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
