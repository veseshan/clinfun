# Exact single stage Phase II design

Calculates the exact one stage Phase II design

## Usage

``` r
ph2single(pu,pa,ep1,ep2,nsoln=5)
```

## Arguments

- pu:

  unacceptable response rate

- pa:

  response rate that is desirable

- ep1:

  threshold for the probability of declaring drug desirable under p0

- ep2:

  threshold for the probability of rejecting the drug under p1

- nsoln:

  number of designs with given alpha and beta

## Value

ph2single returns a data frame with variables: n, r, and the Type I and
Type II errors. Treatment desirable if \>r responses seen.
