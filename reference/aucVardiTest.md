# Two-Sample Tests for Growth Curves under Dependent Right Censoring

Permutation test for comparing growth curves across tow groups under
dependent right censoring.

## Usage

``` r
aucVardiTest(meas, grp, tim=NULL, cgrps=NULL, nperm=5000)
```

## Arguments

- meas:

  Matrix of measurements where the rows are the subjects and columns the
  timepoints. At least one value should not be missing in each row. For
  example they can be tumor sizes measured over time.

- grp:

  Group indicator for each subject. There must be at least two different
  groups. This can represent each subject's treatment.

- tim:

  Times at which the measurements in `meas` are taken. If missing, the
  times are set to 1 through `ncol(meas)`.

- cgrps:

  The two groups that are being compared. If missing the first two
  groups will be compared.

- nperm:

  Number of permutations for the reference distribution.

## Value

returns a list with objects ostat, pstat and p.value which are the
observed test statistic for the two groups being compared, values of the
statistics when the group labels are permuted

## Details

The test statistic is defined as the sum of pairwise differences in the
partial areas under the growth curve. For each pair of subjects the
partial area is computed until the smaller of the maximum followup
times. For each subject, linear interpolation is is used to fill-in
missing values prior to the maximum followup time. The reference
distribution of obtained by permuting the group labels.

## References

Vardi Y., Ying Z. and Zhang C.H. (2001). Two-Sample Tests for Growth
Curves under Dependent Right Censoring. *Biometrika* 88, 949-960.

## Examples

``` r
  grp <- sample(1:3, 100, replace=TRUE)
  grp0 <- LETTERS[grp]
  maxfup <- sample(5:20, 100, replace=TRUE)
  meas <- matrix(NA, 100, 20)
  for(i in 1:100) {
    meas[i, 1:maxfup[i]] <- cumsum((3+0.04*grp[i]) + rnorm(maxfup[i]))
  }
  aucVardiTest(meas, grp)
#> statistic for 2 - 1  =  -2408.479 ; p-value =  0.5806839 
  aucVardiTest(meas, grp0, cgrps=c("C","B"))
#> statistic for B - C  =  -219.7759 ; p-value =  0.954809 
```
