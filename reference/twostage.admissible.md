# Admissible design options between Minimax and Optimal

Lists the admissible design options between

## Usage

``` r
twostage.admissible(x)
```

## Arguments

- x:

  output from `ph2simon` call

## Value

twostage.admissible returns design options that are admissible (Jung et
al, 2004). The output is a matrix with 8 columns: r1, n1, r, n, EN(p0),
PET(p0), qLo, qHi. The columns qLo and qHi give the range of probability
values for which the particular design is admissible.

## References

Jung SH, Lee T, Kim K, and George, SL. (2004). Admissible two-stage
designs for phase II cancer clinical trials. *Statistics in medicine*
23(4), 561-569.

## Examples

``` r
  oo = ph2simon(0.5, 0.7, 0.05, 0.1)
  twostage.admissible(oo)
#>            r1 n1  r  n   EN(p0)   PET(p0)   qLo   qHi
#> Minimax    14 27 32 53 36.11440 0.6494460 0.285 1.000
#> Admissible 12 23 34 57 34.51987 0.6611803 0.112 0.285
#> Optimal    13 24 36 61 34.01324 0.7293719 0.000 0.112
```
