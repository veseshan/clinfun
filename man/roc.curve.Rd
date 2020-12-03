\name{roc.curve}
\title{Empirical ROC curve}
\alias{roc.curve}
\alias{print.roc.curve}
\alias{plot.roc.curve}
\alias{lines.roc.curve}
\description{
  Computes the empricial ROC curve for a diagnostic tool.
}
\usage{
  roc.curve(marker, status, method=c("empirical"))
  \method{print}{roc.curve}(x, \dots)
  \method{plot}{roc.curve}(x, \dots)
  \method{lines}{roc.curve}(x, \dots)
}
\arguments{
  \item{marker}{the marker values for each subject.}
  \item{status}{binary disease status indicator}
  \item{method}{the method for estimating the ROC curve.  Currently only
    the empirical curve is implemented.}
  \item{x}{object of class roc.area.test output from this function.}
  \item{...}{optional arguments to the print, plot and lines functions.}
}
\value{a list with the following elements
  \item{tpr}{true positive rates for all thresholds.}
  \item{fpr}{true positive rates for all thresholds.}
  \item{marker}{the diagnostic marker being studied.}
  \item{status}{binary disease }

  The "print" method returns the nonparametric AUC and its s.e.

  The "plot" and "lines" methods can be used to draw a new plot and add
  to an existing plot of ROC curve.
}
\details{
  The computation is based on assuming that larger values of the marker
  is indicative of the disease.  So for a given threshold x0, TPR is
  P(marker >= x0|status =1) and FPR is P(marker >= x0|status =0).  This
  function computes the empirical estimates of TPR and FPR.
}
\examples{
g <- rep(0:1, 50)
x <- rnorm(100) + g
y <- rnorm(100) + 1.5*g
o <- roc.curve(x, g)
plot(o)
lines(roc.curve(y, g), col=2)
}
