# Empirical ROC curve

Computes the empirical ROC curve for a diagnostic tool.

## Usage

``` r
roc.curve(marker, status, method=c("empirical"))
  # S3 method for class 'roc.curve'
print(x, ...)
  # S3 method for class 'roc.curve'
plot(x, PRC=FALSE, ...)
  # S3 method for class 'roc.curve'
lines(x, PRC=FALSE, ...)
```

## Arguments

- marker:

  the marker values for each subject.

- status:

  binary disease status indicator

- method:

  the method for estimating the ROC curve. Currently only the empirical
  curve is implemented.

- x:

  object of class roc.curve output from this function.

- PRC:

  flag to tell whether ROC or Precision-Recall curve plotted.

- ...:

  optional arguments to the print, plot and lines functions.

## Value

a list with the following elements

- marker:

  the diagnostic marker being studied.

- status:

  binary disease

- tpr:

  true positive rates for all thresholds.

- fpr:

  true positive rates for all thresholds.

- ppv:

  positive predictive values for all thresholds.

- npv:

  negative predictive values for all thresholds.

The "print" method returns the nonparametric AUC and its s.e.

The "plot" and "lines" methods can be used to draw a new plot and add to
an existing plot of ROC curve.

## Details

The computation is based on assuming that larger values of the marker is
indicative of the disease. So for a given threshold x0, TPR is P(marker
\>= x0\|status =1) and FPR is P(marker \>= x0\|status =0). This function
computes the empirical estimates of TPR and FPR.

## Examples

``` r
g <- rep(0:1, 50)
x <- rnorm(100) + g
y <- rnorm(100) + 1.5*g
o <- roc.curve(x, g)
plot(o)
lines(roc.curve(y, g), col=2)
```
