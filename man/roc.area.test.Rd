\name{roc.area.test}
\title{Nonparametric area under the ROC curve}
\alias{roc.area.test}
\alias{print.roc.area.test}
\keyword{htest}
\description{
  Computes the nonparametric area under the ROC curve and its variance
  based on U-statistic theory (DDCP). 
}
\usage{
  roc.area.test(markers, status)
  \method{print}{roc.area.test}(x, \dots)
}
\arguments{
  \item{markers}{The marker values for each subject.  If there are more
    than one markers then this should be a matrix.}
  \item{status}{binary disease status indicator}
  \item{x}{object of class roc.area.test output from this function.}
  \item{...}{optional arguments to the print function.}
}
\value{a list with the following elements
  \item{area}{estimated area.}
  \item{var}{estimated variance (matrix).}
  \item{stat}{test statistic for equality of AUCs. Is not returned when
    only one diagnostic marker is present.}
  \item{p.value}{the p-value for the test of equality (2-sided).}
  \item{df}{the degrees of freedom of the chi-square.}

  The "print" method formats and returns the output.
}
\details{
  It calculates the area and its variance.  For more than one marker it
  calculates the statistic to test for the equality of all AUCs.  This
  statistic has a standard normal reference distribution for two
  variables and chi-square with number of variables minus 1.
}
\examples{
g <- rep(0:1, 50)
x <- rnorm(100) + g
y <- rnorm(100) + g
z <- rnorm(100) + g
roc.area.test(cbind(x,y), g)
roc.area.test(cbind(x,y,z), g)
y1 <- y + 0.75*g
roc.area.test(cbind(x,y1), g)
}
\references{
  DeLong, E. R., D. M. DeLong, and D. L. Clarke-Pearson. 1988. Comparing
  the areas under two or more correlated receiver operating characteristic
  curves:  A nonparametric approach. \emph{Biometrics} 44:837-845.
}
