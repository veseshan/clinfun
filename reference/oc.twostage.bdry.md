# Two-stage boundary operating characteristics

Calculates the operating characteristics of a two-stage boundary.

## Usage

``` r
oc.twostage.bdry(pu, pa, r1, n1, r, n)
```

## Arguments

- pu:

  unacceptable response rate

- pa:

  response rate that is desirable

- r1:

  first stage threshold to declare treatment undesirable

- n1:

  first stage sample size

- r:

  overall threshold to declare treatment undesirable

- n:

  total sample size

## Value

oc.twostage.bdry returns the type I and II error rates as well as the
probability of early temination and expected sample size under pu for a
specific boundary.
