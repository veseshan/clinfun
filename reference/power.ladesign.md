# Power of k-sample rank test under Lehmann alternative

Functions to calculate the power of rank tests for animal studies.

## Usage

``` r
power.ladesign(gsize, odds.ratio, sig.level = 0.05, statistic =
     c("Kruskal-Wallis", "Jonckheere-Terpstra"), alternative =
     c("two.sided", "one.sided"), nrep=1e+6) 
  # S3 method for class 'ladesign'
print(x,...)
```

## Arguments

- gsize:

  sample size of the k (= length of vector) groups.

- odds.ratio:

  odds ratio parameters for the k-1 groups. The first group is
  considered the control.

- sig.level:

  the significance level of the test (default = 0.05)

- statistic:

  the test statistic for the k-group comparison. Is one of
  Kruskal-Wallis (default) or Jonckeere-Terpstra.

- alternative:

  one- or two-sided test. Valid only for the Jonckheere-Terpstra test.

- nrep:

  number of reps (default 1 million) for Monte Carlo.

- x:

  object of class ladesign returned by power.ladesign

- ...:

  arguments to be passed on left for S3 method consistency.

## Value

returns a list with objects group.size, odds.ratio, statistic, sig.level
and power. The "print" method formats the output.

## Details

Although the power for Jonckheere-Terpstra test is calculated for any
set of odds ratio, the test is meant for monotone alternative. Thus it
is preferable to specify odds ratios that are monotonically increasing
with all values larger than 1 or decreasing with all values smaller than
1.

## References

Heller G. (2006). Power calculations for preclinical studies using a
K-sample rank test and the Lehmann alternative hypothesis. *Statistics
in Medicine* 25, 2543-2553.

## Examples

``` r
  power.ladesign(c(9,7), 4, statistic="K")
#>          Number of groups = 2 
#>                group size = 9 7 
#> odds ratios w.r.t group 1 = 4 
#>            test statistic = Kruskal-Wallis 
#>               alternative = two.sided 
#>        significance level = 0.05 
#>                     power = 0.587 
  power.ladesign(c(9,7,9), c(2,4), statistic="J")
#>          Number of groups = 3 
#>                group size = 9 7 9 
#> odds ratios w.r.t group 1 = 2 4 
#>            test statistic = Jonckheere-Terpstra 
#>               alternative = two.sided 
#>        significance level = 0.05 
#>                     power = 0.66 
  power.ladesign(c(9,7,9), c(2,4), statistic="J", alt="o")
#>          Number of groups = 3 
#>                group size = 9 7 9 
#> odds ratios w.r.t group 1 = 2 4 
#>            test statistic = Jonckheere-Terpstra 
#>               alternative = one.sided 
#>        significance level = 0.05 
#>                     power = 0.788 
```
