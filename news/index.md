# Changelog

## clinfun 1.1.6 (08/11/2026)

- Issue in dwilcox call in jtpdf.f due to wilcox_free(); fortran code
  rewritten
- Rd files for some functions updated
- version on GitHub (12/06/2023) - moved twostage.admissible call from
  print method to ph2simon; added reference

## clinfun 1.1.5 (10/19/2023)

CRAN release: 2023-10-19

- clarified one-sided vs two-sided futility boundary for gsdesign
  functions

## clinfun 1.1.4 (10/12/2023)

CRAN release: 2023-10-12

- added check for n; since it is for stopping early at least two values
  should be given
- addressed the warning in the example for calogrank; status had both
  death and transplant

## clinfun 1.1.3 (07/11/2023)

CRAN release: 2023-07-13

- fixed the bug in ph2simon/twostage.admissible when minimax is also
  optimal
- ph2simon and toxbdry help have been made more clear

## clinfun 1.1.2 (06/19/2023)

CRAN release: 2023-06-19

- changed all dfloat in fortran to dble
- print method for ph2simon now lists admissible designs

#### New Functionality

- plot and lines method for roc.curve includes the option for
  precision-recall curve
- roc.curve returns PPV and NPV which are needed for precision-recall
  curve

## clinfun 1.1.1 (03/06/2023)

CRAN release: 2023-03-07

- changed class check to use inherits instead of comparing strings

#### New Functionality

- Added out.ties option to coxphCPE for discrete risk score case

## clinfun 1.1.0 (02/22/2022)

CRAN release: 2022-02-22

#### New Functionality

- Added added admissible two stage design function
- Added added futilbdry function for sequential futility stopping
- added 1:r sample allocation in gsdesign functions
- Fixed the warnings messages from deltaAUC.f

#### New Documentation

- Added a `NEWS.md` file to track changes to the package.
- deposited development version on GitHub
- created {pkgdown} website and github actions for building site

## clinfun 1.0.15 (04/13/2018)

CRAN release: 2018-04-13

- Added new function deltaAUC to test for difference in the area under
  the ROC curves from nested binary regression models.
- power.ladesign does exact null if number of combinations is small
- Added to help of fedesign, ph2simon and ph2single functions

## clinfun 1.0.14 (04/25/2017)

CRAN release: 2017-04-25

- added init.c to register native routines

## clinfun 1.0.13 (11/29/2016)

CRAN release: 2016-11-30

- permlogrank bug fix - strata option gave an error. Names of elements
  in survfit.object changed. (strata instead of ntimes.strata and n a
  vector).

## clinfun 1.0.12 (08/03/2016)

CRAN release: 2016-10-17

- added aucVardiTest to compare growth curves

## clinfun 1.0.11 (08/11/2015)

CRAN release: 2015-08-12

- check for misssing values in jonckheere.test

## clinfun 1.0.10 (05/18/2015)

CRAN release: 2015-05-18

- bug fix for pselect; nlen \> 3 should have been nlen \> 2

## clinfun 1.0.9 (02/24/2015)

CRAN release: 2015-02-24

- use requireNamespace to remove NOTES on functions using survival pkg

## clinfun 1.0.8 (02/19/2015)

- fix name swap in the result of mdrr in v1.0.7

## clinfun 1.0.7 (02/19/2015)

- new function mdrr added to calculate minimum detectable difference in
  response rates for given average response rate and class proportion

## clinfun 1.0.6 (06/10/2014)

CRAN release: 2014-07-22

- On rare occasions jonckheere.test gave a p-value bigger than 1.
  Sometimes 2*min(iPVAL, dPVAL, 1) can be larger than 1. Replace with
  2*min(iPVAL, dPVAL, 0.5) (Thanks to Drs. Shterev and Owzar of Duke).

## clinfun 1.0.5 (04/16/2013)

CRAN release: 2013-04-16

- Fixed the bound checks

## clinfun 1.0.4 (01/22/2013)

CRAN release: 2013-01-23

- Added the option to calculate continuity corrected sample size in the
  function gsdesign.binomial.

## clinfun 1.0.3 (10/15/2012)

CRAN release: 2012-10-15

- Fixed two-sided p-value \> 1 bug when statistic is exactly its mean.
- Fixed Rd file to address LaTeX warnings.

## clinfun 1.0.2 (09/25/2012)

CRAN release: 2012-09-26

- integer overflow in djonck for the Jonckheere-Terpstra statistic.
  replace with pdf calculation using Mark van de Wiel convolution.
- Create seperate help for functions permlogrank and jonckheere.test

## clinfun 1.0.1 (08/14/2012)

CRAN release: 2012-08-14

- Fix integer overflow because of n0*n1 == 0 in roc.curve and nn*nd == 0
  in roc.area.test
- In coxphQuantile eliminate times for which survival probability is 0
  or 1 from quantile computation.

## clinfun 1.0.0 (03/13/2012)

CRAN release: 2012-03-13

- ## clinfun number changed to 1.0.0 in preparation for R 2.15.

  CRAN release: 2012-03-13

- toxbdry now does the entire Pocock to O’Brien-Fleming range of
  boundaries. Added references for the method.

- Fixed linebreak in the help for coxphERR.

## clinfun 0.9.9 (01/18/2012)

CRAN release: 2012-01-18

- Added coxphERR to calculate Heller’s explained relative risk.
- Fixed NaN bug in toxbdry when priority=“alt” is used.

## clinfun 0.9.8 (09/13/2011)

CRAN release: 2011-09-14

- fixed bug in print.gsdesign for binomial case (p1,p2 instead of pC,pE)
- ROC functions now check that there are at least one each of status=0,1

## clinfun 0.9.7 (04/27/2011)

CRAN release: 2011-04-27

- fixed fortran code to address gfortran-4.6 warnings
- added ktau a faster implementation of cor(x, y, method=“k”). not in
  NAMESPACE. Should be called using clinfun:::ktau

## clinfun 0.9.7 (04/25/2011)

CRAN release: 2011-04-27

- bug fix: roc.area.test integer overflow for large nn\*nd.
- use sort function to speed up roc curve and area estimation.

## clinfun 0.9.6 (03/24/2011)

CRAN release: 2011-03-24

- bug fix: roc.area.test gave NaN as the statistic and p-value when the
  markers are identical. Changed it to 0.

## clinfun 0.9.5 (03/09/2011)

CRAN release: 2011-03-09

- bug fix: gsdesign funtions not returning the sample size / \# events.

## clinfun 0.9.4 (02/24/2011)

CRAN release: 2011-02-25

- twostage.inference for umvue, p-value and CI for 2 stage design.

## clinfun 0.9.3 (12/06/2010)

CRAN release: 2010-12-07

- ph2simon was testing whether dim is null for feasible solution.
  Replaced with nrow == 0 since it is now possible to have 0 rows.

## clinfun 0.9.2 (11/05/2010)

CRAN release: 2010-11-10

- Added functions to compute and plot the empirical ROC curve.

## clinfun 0.9.1 (11/03/2010)

CRAN release: 2010-11-05

- Added functions for the area and permutation tests to compare ROC.
- Checks that min.diff is greater than 0 in pselect.

## clinfun 0.9.0 (11/01/2010)

CRAN release: 2010-11-03

- Added a non-binding futility boundary to gsdesign

## clinfun 0.8.10 (04/16/2010)

CRAN release: 2010-04-16

- variable names for returned data.frame in coxphQuantile
- examples in coxphQuantile & coxphCPE use status==2

## clinfun 0.8.9 (04/09/2010)

CRAN release: 2010-04-10

- check R# clinfun so that coxphCPE works for any \# clinfun (see 0.8.9)

## clinfun 0.8.8 (04/07/2010)

CRAN release: 2010-04-08

- Change coxphCPE to reflect the fact that model.matrix.coxph doesn’t
  have an intercept term.

## clinfun 0.8.8 (02/25/2010)

CRAN release: 2010-04-08

- Added the function or2pcase

## clinfun 0.8.7 (11/20/2009)

CRAN release: 2009-11-22

- Fixed the 0/0 bug in the revised pselect

## clinfun 0.8.6 (11/17/2009)

CRAN release: 2009-11-19

- Changed Venkat’s affiliation to MSKCC.
- Fixed pselect to calculate the selection probability correctly when
  only one treatment exceeds the min.resp threshold.

## clinfun 0.8.5 (07/10/2009)

CRAN release: 2009-07-18

- Changed the Jonckheere-Terpstra statistic such that large value is
  indicative of increasing group locations and small for decreasing.
  Function warns that p-value is based on approximation for tied data.

## clinfun 0.8.4 (12/02/2008)

CRAN release: 2008-12-02

- Added functionality to pselect. It can do unequal sample size for the
  case of two treatments. min.diff can be specified as a rate instead of
  number of responses. Output element names changed to be more
  descriptive.

## clinfun 0.8.3 (11/18/2008 & 09/18/2008)

CRAN release: 2008-11-18

- Fixed the bug CPS.ssize. call inside used fixed alpha, power & r.
- Fixed the bug in the approximate one-sided p-value of Jonckheere test.
  wrong tail was used.

## clinfun 0.8.2 (06/20/2008)

CRAN release: 2008-07-01

- toxbdry allows the error threshold to prioritize when the sample size
  is too small to have both satisfied.

## clinfun 0.8.1 (06/17/2008)

CRAN release: 2008-06-17

- Fixed gsdesign to allow for fixed sample designs. Help file fixed.

## clinfun 0.8.0 (05/23/2008)

CRAN release: 2008-05-23

- New \# clinfun with one new function power.lehmann.design
