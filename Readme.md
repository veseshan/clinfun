
<!-- README.md is generated from README.Rmd. Please edit that file -->

# clinfun

<!-- badges: start -->

[![R-CMD-check](https://github.com/veseshan/clinfun/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/veseshan/clinfun/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/clinfun)](https://CRAN.R-project.org/package=clinfun)
<!-- badges: end -->

Utilities to make your clinical collaborations easier if not fun. It
contains functions for designing studies such as Simon 2-stage and group
sequential designs and for data analysis such as Jonckheere-Terpstra
test and estimating survival quantiles.

## Installation

You can install the released version of clinfun from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("clinfun")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("veseshan/clinfun")
```

## Example

Here is a basic example which shows how calculate CPS and Fisher’s exact
sample sizes with power estimates for given response rates of standard
and experiment groups.

``` r
library(clinfun)

fe.ssize(p1 = .8,
         p2 = .2,
         alpha=0.05,
         power=0.8)
#>              Group 1 Group 2 Exact Power
#> CPS               13      13   0.8688275
#> Fisher Exact      12      12   0.8115276
```
