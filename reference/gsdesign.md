# Group Sequential Designs

Functions to calculate sample size for group sequential designs

## Usage

``` r
gsdesign.binomial(ifrac, pC, pE, r = 1, sig.level = 0.05, power = 0.8,
  delta.eb=0.5, delta.fb = NULL, alternative = c("two.sided",
  "one.sided"), pooled.variance = FALSE, CPS = TRUE, tol=0.00001, ...) 
gsdesign.normal(ifrac, delta, sd = 1, r = 1, sig.level = 0.05,
  power = 0.8, delta.eb = 0.5, delta.fb = NULL, alternative = 
  c("two.sided", "one.sided"), tol=0.00001, ...)
gsdesign.survival(ifrac, haz.ratio, r = 1, sig.level = 0.05, 
  power = 0.8, delta.eb = 0.5, delta.fb = NULL, alternative = 
  c("two.sided", "one.sided"), tol=0.00001, ...)
```

## Arguments

- ifrac:

  information fraction or the ratio of current sample size or number of
  events to the total sample size or number of events. This should be an
  increasing vector of numbers from 0 to 1 with the last one being 1. If
  just 1 is given a fixed sample design is derived.

- pC:

  prob of success of the standard therapy (for binomial data)

- pE:

  prob of success of the experimental therapy (for binomial data)

- delta:

  true difference in means (for normal data)

- sd:

  standard deviation (for normal data)

- haz.ratio:

  hazard ratio (for survival comparison)

- r:

  treatment allocation of r (default=1) experimental per 1 control.

- sig.level:

  significance level (type I error probability)

- power:

  power of test (1 minus type II error probability)

- delta.eb:

  power for efficacy boundary in the Pocock (=0) to O'Brien-Fleming
  (=0.5) family (default is 0.5)

- delta.fb:

  power for futility boundary in the Pocock (=0) to O'Brien-Fleming
  (=0.5) family (default is NULL i.e. no futility boundary is
  requested.)

- alternative:

  one- or two-sided test.

- pooled.variance:

  whether the test statistic is standardized by pooled
  (2\*pbar\*(1-pbar)) or unpooled variance (pC\*(1-pC) + pE\*(1-pE)).
  Default is unpooled variance.

- CPS:

  whether continuity correction is used for sample size calculation as
  in Casagrande, Pike & Smith. Default is to use it.

- tol:

  tolerance level for multivariate normal probability computation.

- ...:

  additional options passed on the pmvnorm function.

## Value

a list with ifrac, sig.level, power, alternative, delta.eb, delta.fb
and:

- efbdry:

  the critical value to use at the different looks. For two-sided
  alternative the absolute test statistic should exceed this.

- futbdry:

  the critical value to use at the different looks. For two-sided
  alternative the absolute test statistic should be below this.

- sample.size:

  the sample size per arm for binomial/normal data.

- num.events:

  the total number of failures which should be converted to number of
  subjects using censoring proportion.

## Details

The futility boundary is not returned when delta.fb is not specified
i.e. stopping for futility is not requested. The futility boundary is
non-binding. That is the significance level is not adjusted to account
for early stopping for futility. This makes the test a bit conservative
in that the true size is less than the nominal level.

If the alternative is two-sided by default the futility boundary will
also be two-sided i.e. continuation region is wedge shaped. However, if
the goal is to show the superiority of the experimental treatment then
futility boundary should be one sided. This can be achieved by deriving
the boundaries for one-sided alternative and significance level set at
half of the value used for two sided alternative. See the examples
section for a representative design for which the trial cannot be
stopped at the first look for futility.

The Casagrande-Pike-Smith type continuity correction is obtained using
the formula \$\$n\*\[1 + \sqrt{1+4/(\|pC-pE\|\*n)}\]^2\$\$ where n is
the uncorrected sample size.

## Examples

``` r
  gsdesign.normal(1:4/4, 0.25, sig.level=0.05, alt="t", delta.fb=0.5)
#> 
#>  Group sequential design for comparing normal data with delta = 0.25 , sd = 1 
#>    Treatment allocated at 1:1 (C:E) ratio 
#>    power family of boundary; 0 (Pocock) to 0.5 (O'Brien-Fleming) 
#> 
#>   sample sizes (by arm) = 292 292 
#>    information fraction = 0.25 0.50 0.75 1.00 
#>       efficacy boundary = 4.048 2.862 2.337 2.024 (power = 0.5) 
#>       futility boundary = -Inf 0.726 1.465 2.024 (power = 0.5) 
#>               sig.level = 0.05 
#>                   power = 0.8 
#>             alternative = two.sided 
#> 
  gsdesign.normal(1:4/4, 0.25, sig.level=0.025, alt="o", delta.fb=0.5)
#> 
#>  Group sequential design for comparing normal data with delta = 0.25 , sd = 1 
#>    Treatment allocated at 1:1 (C:E) ratio 
#>    power family of boundary; 0 (Pocock) to 0.5 (O'Brien-Fleming) 
#> 
#>   sample sizes (by arm) = 294.84 294.84 
#>    information fraction = 0.25 0.50 0.75 1.00 
#>       efficacy boundary = 4.048 2.862 2.337 2.024 (power = 0.5) 
#>       futility boundary = -0.506 0.716 1.461 2.024 (power = 0.5) 
#>               sig.level = 0.025 
#>                   power = 0.8 
#>             alternative = one.sided 
#> 
```
