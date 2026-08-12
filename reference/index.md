# Package index

## Trial Design Functions

- [`ph2simon()`](ph2simon.md) [`print(`*`<ph2simon>`*`)`](ph2simon.md)
  [`plot(`*`<ph2simon>`*`)`](ph2simon.md) : Simon's 2-stage Phase II
  design
- [`ph2single()`](ph2single.md) : Exact single stage Phase II design
- [`twostage.admissible()`](twostage.admissible.md) : Admissible design
  options between Minimax and Optimal
- [`twostage.inference()`](twostage.inference.md) : Inference following
  a two-stage design for binary response
- [`toxbdry()`](toxbdry.md) [`futilbdry()`](toxbdry.md)
  [`bdrycross.prob()`](toxbdry.md)
  [`print(`*`<toxbdry>`*`)`](toxbdry.md)
  [`print(`*`<futilbdry>`*`)`](toxbdry.md) : Stopping rule for
  toxicity/futility monitoring
- [`gsdesign.binomial()`](gsdesign.md)
  [`gsdesign.normal()`](gsdesign.md)
  [`gsdesign.survival()`](gsdesign.md) : Group Sequential Designs
- [`oc.twostage.bdry()`](oc.twostage.bdry.md) : Two-stage boundary
  operating characteristics
- [`fe.ssize()`](fedesign.md) [`CPS.ssize()`](fedesign.md)
  [`fe.mdor()`](fedesign.md) [`mdrr()`](fedesign.md)
  [`fe.power()`](fedesign.md) [`or2pcase()`](fedesign.md) : Trial
  Designs Based On Fisher's Exact Test
- [`power.ladesign()`](power.ladesign.md)
  [`print(`*`<ladesign>`*`)`](power.ladesign.md) : Power of k-sample
  rank test under Lehmann alternative
- [`pselect()`](pselect.md) : Probability of selection under pick the
  winner rule

## Analysis Functions

### ROC Curves

- [`roc.area.test()`](roc.area.test.md)
  [`print(`*`<roc.area.test>`*`)`](roc.area.test.md) : Nonparametric
  area under the ROC curve
- [`roc.curve()`](roc.curve.md)
  [`print(`*`<roc.curve>`*`)`](roc.curve.md)
  [`plot(`*`<roc.curve>`*`)`](roc.curve.md)
  [`lines(`*`<roc.curve>`*`)`](roc.curve.md) : Empirical ROC curve
- [`roc.perm.test()`](roc.perm.test.md)
  [`print(`*`<roc.perm.test>`*`)`](roc.perm.test.md)
  [`plot(`*`<roc.perm.test>`*`)`](roc.perm.test.md) : Permutation test
  to compare ROC curve
- [`ROCanalysis`](ROCanalysis.md) [`rocanalysis`](ROCanalysis.md) :
  Functions to plot and compare ROC curves.
- [`deltaAUC()`](deltaAUC.md) : Comparing the AUC from ROC curves from
  nested binary regression

### Survival Analysis

- [`permlogrank()`](permlogrank.md) : Permutation version of survdiff
- [`calogrank()`](calogrank.md) : Survival curves analysis of covariance
- [`coxphCPE()`](coxphCPE.md) : Concordance Probability Estimate for Cox
  model
- [`coxphQuantile()`](coxphQuantile.md) : Survival time quantile as a
  function of covariate
- [`coxphERR()`](coxphERR.md) : Heller Explained Relative Risk
- [`jonckheere.test()`](jonckheere.test.md) : Exact/permutation version
  of Jonckheere-Terpstra test

### Miscellaneous

- [`aucVardiTest()`](aucVardiTest.md) : Two-Sample Tests for Growth
  Curves under Dependent Right Censoring
- [`ktau()`](ktau.md) : Kendall's tau-b estimate
