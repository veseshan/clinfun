# Concordance Probability Estimate for Cox model

Calculates the Concordance Probability Estimate for a Cox proportional
hazards model. Both the Gonen and Heller (2005) version for continuous
risk score and Heller and Mo (2016) for discrete risk score can be
calculated.

## Usage

``` r
coxphCPE(phfit, out.ties=FALSE)
```

## Arguments

- phfit:

  output from a proportional hazards fit.

- out.ties:

  binary flag to decide if pairs with tied risk scores should be used.

## Value

coxphCPE returns a vector with CPE, smooth.CPE and se.CPE which are the
estimate, the smoothed estimate and its standard error respectively.

## References

Gonen M and Heller G. (2005) Concordance probability and discriminatory
power in proportional hazards regression. *Biometrika* 92, 965-970.

Heller G and Mo Q. (2016). Estimating the concordance probability in a
survival analysis with a discrete number of risk groups. *Lifetime Data
Analysis*, 22,263-79.

## Examples

``` r
if (FALSE)  library(survival)
  data(pbc)
#> Warning: data set ‘pbc’ not found
  pbcfit <- coxph(Surv(time, status==2) ~ trt + log(copper), pbc,
    subset=(trt>0 & copper>0))
#> Error in coxph(Surv(time, status == 2) ~ trt + log(copper), pbc, subset = (trt >     0 & copper > 0)): could not find function "coxph"
  coxphCPE(pbcfit)
#> Error: object 'pbcfit' not found
 # \dontrun{}
```
