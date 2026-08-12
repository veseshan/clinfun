# Nonparametric area under the ROC curve

Computes the nonparametric area under the ROC curve and its variance
based on U-statistic theory (DeLong et al, 1988).

## Usage

``` r
roc.area.test(markers, status)
  # S3 method for class 'roc.area.test'
print(x, ...)
```

## Arguments

- markers:

  The marker values for each subject. If there are more than one markers
  then this should be a matrix.

- status:

  binary disease status indicator

- x:

  object of class roc.area.test output from this function.

- ...:

  optional arguments to the print function.

## Value

a list with the following elements

- area:

  estimated area.

- var:

  estimated variance (matrix).

- stat:

  test statistic for equality of AUCs. Is not returned when only one
  diagnostic marker is present.

- p.value:

  the p-value for the test of equality (2-sided).

- df:

  the degrees of freedom of the chi-square.

The "print" method formats and returns the output.

## Details

It calculates the area and its variance. For more than one marker it
calculates the statistic to test for the equality of all AUCs. This
statistic has a standard normal reference distribution for two variables
and chi-square with number of variables minus 1.

## References

DeLong, E. R., D. M. DeLong, and D. L. Clarke-Pearson. 1988. Comparing
the areas under two or more correlated receiver operating characteristic
curves: A nonparametric approach. *Biometrics* 44:837-845.

## Examples

``` r
g <- rep(0:1, 50)
x <- rnorm(100) + g
y <- rnorm(100) + g
z <- rnorm(100) + g
roc.area.test(cbind(x,y), g)
#>   2 markers with AUC 0.7368 0.7368 
#>   test statistic = 0 
#>   p-value = 1 from standard normal reference 
roc.area.test(cbind(x,y,z), g)
#>   3 markers with AUC 0.7368 0.7368 0.7496 
#>   test statistic = 0.03830606 
#>   p-value = 0.9810292 from chi-square (df = 2) reference 
y1 <- y + 0.75*g
roc.area.test(cbind(x,y1), g)
#>   2 markers with AUC 0.7368 0.8628 
#>   test statistic = 2.212703 
#>   p-value = 0.02691816 from standard normal reference 
```
