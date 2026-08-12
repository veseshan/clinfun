# Permutation test to compare ROC curve

Computes the test statistic and permutation reference distribution for
comparing paired or unpaired ROC curves.

## Usage

``` r
roc.perm.test(marker, status, marker2=NULL, group=NULL,
              nperm=2500, mp=NULL)
# S3 method for class 'roc.perm.test'
print(x, ...)
# S3 method for class 'roc.perm.test'
plot(x, ...)
```

## Arguments

- marker:

  marker values for each subject.

- status:

  binary disease status indicator.

- marker2:

  second diagnostic marker for the same subjects (paired).

- group:

  indicator of which diagnostic test was used (unpaired).

- nperm:

  number of permutations for the reference distribution.

- mp:

  mixing proportion for the unpaired case when proportion of diseased
  subjects can differ.

- x:

  object of class roc.perm.test output from this function.

- ...:

  optional arguments to print and plot functions.

## Value

an object of class roc.perm.test with the following elements

- ostat:

  test statistic from the observed data.

- pstat:

  test statistic from permuted data.

- p.value:

  the p-value for the test of equality (2-sided).

The "print" method formats and returns the statistic and p-value. The
"plot" method plots the density from the permutation reference
distribution and marks the location of the observed statistic.

## Details

This function implements the permutation method described in the
Venkatraman and Begg (1996) paper for the paired case and the
Venkatraman (2000) paper for the unpaired case.

The function detects whether the data are paired or unpaired by testing
which of the options marker2 and group is specified. If both are missing
it will stop with an error message. At present exactly one should be
missing.

## References

Venkatraman, E.S. and Begg, C.B. (1996). A distribution-free procedure
for comparing receiver operating characteristic curves from a paired
experiment. *Biometrika* 83, 835-848.

Venkatraman, E.S. (2000) A permutation test to compare receiver
operating characteristic curves. *Biometrics* 56(4):1134-8.

## Examples

``` r
d <- rep(0:1, 50)
x <- rnorm(100) + 1.2*d
y <- rnorm(100) + 1.2*d
oo <- roc.perm.test(x, d, marker2=y)
plot(oo)

oo <- roc.perm.test(c(x,y), c(d,d), group=rep(1:2,each=100))
plot(oo)
```
