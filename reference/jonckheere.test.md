# Exact/permutation version of Jonckheere-Terpstra test

Jonckheere-Terpstra test to test for ordered differences among classes

## Usage

``` r
jonckheere.test(x, g, alternative = c("two.sided", "increasing",
                 "decreasing"), nperm=NULL)
```

## Arguments

- x, g:

  data and group vector

- alternative:

  means are monotonic (two.sided), increasing, or decreasing

- nperm:

  number of permutations for the reference distribution. The default is
  null in which case the permutation p-value is not computed. Recommend
  that the user set nperm to be 1000 or higher if permutation p-value is
  desired.

## Details

jonckheere.test is the exact (permutation) version of the
Jonckheere-Terpstra test. It uses the statistic \$\$\sum\_{k\<l}
\sum\_{ij} I(X\_{ik} \< X\_{jl}) + 0.5 I(X\_{ik} = X\_{jl}),\$\$ where
\\i, j\\ are observations in groups \\k\\ and \\l\\ respectively. The
asymptotic version is equivalent to cor.test(x, g, method="k"). The
exact calculation requires that there be no ties and that the sample
size is less than 100. When data are tied and sample size is at most 100
permutation p-value is returned.

## References

Jonckheere, A. R. (1954). A distribution-free k-sample test again
ordered alternatives. *Biometrika* 41:133-145.

Terpstra, T. J. (1952). The asymptotic normality and consistency of
Kendall's test against trend, when ties are present in one ranking.
*Indagationes Mathematicae* 14:327-333.

## Examples

``` r
  set.seed(1234)
  g <- rep(1:5, rep(10,5))
  x <- rnorm(50)
  jonckheere.test(x+0.3*g, g)
#> 
#>  Jonckheere-Terpstra test
#> 
#> data:  
#> JT = 629, p-value = 0.02734
#> alternative hypothesis: two.sided
#> 
  x[1:2] <- mean(x[1:2]) # tied data
  jonckheere.test(x+0.3*g, g)
#> Warning: Sample size > 100 or data with ties 
#>  p-value based on normal approximation. Specify nperm for permutation p-value
#> 
#>  Jonckheere-Terpstra test
#> 
#> data:  
#> JT = 639, p-value = 0.01741
#> alternative hypothesis: two.sided
#> 
  jonckheere.test(x+0.3*g, g, nperm=5000)
#> 
#>  Jonckheere-Terpstra test
#> 
#> data:  
#> JT = 639, p-value = 0.0136
#> alternative hypothesis: two.sided
#> 
```
