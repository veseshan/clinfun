# Permutation version of survdiff

Small sample survdiff using permutation reference distributions.

## Usage

``` r
permlogrank(formula, data, subset, na.action, rho=0, nperm=5000)
```

## Arguments

- nperm:

  number of permutations for the reference distribution

- formula, data, subset, na.action, rho:

  see survdiff for details

## Details

permlogrank is the permutation version of k-sample survdiff. see
survdiff in survival package for details.

## References

Heller G, Venkatraman ES. (1996). Resampling procedures to compare two
survival distributions in the presence of right censored data.
*Biometrics* 52:1204-1213.
