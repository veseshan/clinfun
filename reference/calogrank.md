# Survival curves analysis of covariance

Logrank test to compare survival curves adjusting for covariates

## Usage

``` r
calogrank(ftime, fstatus, grp, cvt, strat=NULL)
```

## Arguments

- ftime:

  failure times

- fstatus:

  status indicator

- grp:

  group indicator

- cvt:

  continuous covariates used for adjusted analysis

- strat:

  stratification variable

## Details

calogrank is the covariate adjusted version of k-sample survdiff. The
function in its current form only does basic error checking.

## References

Heller G. and Venkatraman E.S. (2004) A nonparametric test to compare
survival distributions with covariate adjustment. *JRSS-B* 66, 719-733.

## Examples

``` r
if (FALSE)   library(survival)
  data(pbc)
#> Warning: data set ‘pbc’ not found
  pbc1 <- pbc
#> Error: object 'pbc' not found
  pbc1$trt[pbc1$trt == -9] <- NA
#> Error: object 'pbc1' not found
  pbc1$copper[pbc1$copper == -9] <- NA
#> Error: object 'pbc1' not found
  # only death (2) is considered; transplant(1) is censored
  calogrank(pbc1$time, pbc1$status==2, pbc1$trt, pbc1[,c("copper")])
#> Loading required namespace: survival
#> Error: object 'pbc1' not found
  calogrank(pbc1$time, pbc1$status==2, pbc1$trt,
                                  pbc1[,c("protime", "copper")]) # \dontrun{}
#> Error: object 'pbc1' not found
```
