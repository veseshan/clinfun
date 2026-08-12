# Survival time quantile as a function of covariate

Given a specified quantile level p, draws a curve between survival time
and specified covariate. The curve represents the p-th quantile of
survival as a function of covariate.

## Usage

``` r
coxphQuantile(phfit, xrange, p=0.5, whichx=1, otherx=NULL, ...)
```

## Arguments

- phfit:

  a `coxph` object.

- xrange:

  the range of covariate values for which the quantiles of survival
  times are computed.

- p:

  the probability level for the quantile (default is median).

- whichx:

  if there are more than one covariates in the Cox model, the one chosen
  for the quantile plot.

- otherx:

  the values for other covariates in the Cox model. If missing uses
  their average values.

- ...:

  additional parameters to be passed on to the lines command.

## Details

This function is used to draw quantile curves. It requires a plot of the
data (time & covariate of interest) to be present. See example.

It invisibly returns the observed failure times and the covariate values
at which the estimated survival probability is (exactly) p.

## References

Heller G. and Simonoff J.S. (1992) Prediction in censored survival data:
A comparison of the proportional hazards and linear regression models.
*Biometrics* 48, 101-115.

## Examples

``` r
if (FALSE)   library(survival)
data(pbc)
#> Warning: data set ‘pbc’ not found
pbcfit <- coxph(Surv(time, status==2) ~ trt + log(copper), pbc,
                      subset=(trt>0 & copper>0)) 
#> Error in coxph(Surv(time, status == 2) ~ trt + log(copper), pbc, subset = (trt >     0 & copper > 0)): could not find function "coxph"
plot(log(pbc$copper[pbc$trt>0 & pbc$copper>0]), pbc$time[pbc$trt>0 &
  pbc$copper>0], pch=c("o","x")[1+pbc$status[pbc$trt>0 & pbc$copper>0]], 
  xlab="log Copper", ylab="Survival time")
#> Error: object 'pbc' not found
coxphQuantile(pbcfit, c(2.5,6), whichx=2, otherx=1)
#> Error: object 'pbcfit' not found
coxphQuantile(pbcfit, c(2.5,6), p=0.75, whichx=2, otherx=2, col=2) # \dontrun{}
#> Error: object 'pbcfit' not found
```
