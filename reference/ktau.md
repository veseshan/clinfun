# Kendall's tau-b estimate

Calculates the Kendall's tau-b.

## Usage

``` r
ktau(x, y)
```

## Arguments

- x:

  first variable

- y:

  second variable

## Value

ktau returns Kendall's tau-b.

## Details

ktau computes the same quantity as cor(x, y, method="kendall"). It uses
a faster algorithm than pairwise comparisons used by cor.

## Examples

``` r
  set.seed(1234)
  x <- rnorm(10000); y <- x+rnorm(10000)
  cor(x, y, method="k")
#> [1] 0.4956963
  clinfun:::ktau(x,y)  
#> [1] 0.4956963
```
