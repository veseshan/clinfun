# Probability of selection under pick the winner rule

Calculates the probability of selecting the treatment with the higher
response rate under the pick the winner rule.

## Usage

``` r
pselect(n, p, min.diff=NULL, min.resp=NULL)
```

## Arguments

- n:

  sample size for each treatment arm. This is either a single integer or
  a vector of two integers for the special case of comparing two
  treatments with unequal sample sizes

- p:

  vector of response rates for the treatments.

- min.diff:

  this is the number of responses or the rate by which the best
  treatment should be superior to the others to be chosen. This must be
  a positive integer or a rate between 0 and 1. If missing it defaults
  to 1 for the equal sample size case but quits with a warning for the
  unequal sample size case.

- min.resp:

  the minimum number of responses in each treatment arm for it to be
  considered further. If missing defaults to 0.

## Value

the function returns a list with:

- prob.none.worthy:

  is the probability that no treatment has the minimum number of
  responses specified in min.resp. this element is present only if
  min.resp is greater than 0 for at least one arm.

- prob.inconclusive:

  this is the probability that the best treatment has the requisite
  min.resp responses but exceeds the second best by less than min.diff
  responses (rate) provided the second best also has at least min.resp
  responses.

- prob.selection:

  this is a matrix which for each treatment gives the response
  probability and the probability of selecting it i.e. the number of
  responses in the chosen arm is at least min.resp and either none of
  the remaining arms exceed the min.resp threshold or the chosen (best)
  arm is better than the second best by at least min.diff responses
  (rate).

## References

Simon R, Wittes RE, Ellenberg SS. (1985). Randomized phase II clinical
trials. *Cancer Treat Rep* 69, 1375-1381.

## Examples

``` r
  # selection when no diffrence i.e. type I error
  pselect(18, c(0.2, 0.2, 0.2))
#> $prob.inconclusive
#> [1] 0.2236879
#> 
#> $prob.selection
#>       sample.size resp.rate prob.selection
#> trt 1          18       0.2      0.2587707
#> trt 2          18       0.2      0.2587707
#> trt 3          18       0.2      0.2587707
#> 
  # selection probability
  pselect(18, c(0.2, 0.2, 0.4))
#> $prob.inconclusive
#> [1] 0.09676346
#> 
#> $prob.selection
#>       sample.size resp.rate prob.selection
#> trt 1          18       0.2     0.05136236
#> trt 2          18       0.2     0.05136236
#> trt 3          18       0.4     0.80051182
#> 
  pselect(26, c(0.2, 0.2, 0.4), min.diff=2, min.resp=3)
#> $prob.none.worthy
#> [1] 1.961946e-06
#> 
#> $prob.inconclusive
#> [1] 0.1706533
#> 
#> $prob.selection
#>       sample.size resp.rate prob.selection
#> trt 1          26       0.2     0.01536975
#> trt 2          26       0.2     0.01536975
#> trt 3          26       0.4     0.79860522
#> 
  # unequal sample size case
  pselect(c(27,54), c(0.5, 0.65), min.diff=0.05)
#> $prob.inconclusive
#> [1] 0.1405518
#> 
#> $prob.selection
#>       sample.size resp.rate prob.selection
#> trt 1          27      0.50     0.04570331
#> trt 2          54      0.65     0.81374489
#> 
  # unequal sample size case
  pselect(c(27,54), c(0.5, 0.65), min.diff=0.05, min.resp=c(14,27))
#> $prob.none.worthy
#> [1] 0.00405262
#> 
#> $prob.inconclusive
#> [1] 0.1340373
#> 
#> $prob.selection
#>       sample.size resp.rate prob.selection
#> trt 1          27      0.50     0.04625472
#> trt 2          54      0.65     0.81565539
#> 
```
