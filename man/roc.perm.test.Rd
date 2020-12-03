\name{roc.perm.test}
\title{Permutation test to compare ROC curve}
\alias{roc.perm.test}
\alias{print.roc.perm.test}
\alias{plot.roc.perm.test}
\keyword{htest}
\description{
  Computes the test statistic and permutation reference distribution for
  comparing paired or unpaired ROC curves.
}
\usage{
roc.perm.test(marker, status, marker2=NULL, group=NULL,
              nperm=2500, mp=NULL)
\method{print}{roc.perm.test}(x, \dots)
\method{plot}{roc.perm.test}(x, \dots)
}
\arguments{
  \item{marker}{marker values for each subject.}
  \item{status}{binary disease status indicator.}
  \item{marker2}{second diagnostic marker for the same subjects (paired).}
  \item{group}{indicator of which diagnostic test was used (unpaired).}
  \item{nperm}{number of permutations for the reference distribution.}
  \item{mp}{mixing proportion for the unpaired case when proportion of
    diseased subjects can differ.} 
  \item{x}{object of class roc.perm.test output from this function.}
  \item{...}{optional arguments to print and plot functions.}
}

\value{an object of class roc.perm.test with the following elements
  \item{ostat}{test statistic from the observed data.}
  \item{pstat}{test statistic from permuted data.}
  \item{p.value}{the p-value for the test of equality (2-sided).}

  The "print" method formats and returns the statistic and p-value.
  The "plot" method plots the density from the permutation reference
  distribution and marks the location of the observed statistic.
}

\details{
  This function implements the permutation method described in the
  Venkatraman and Begg (1996) paper for the paired case and the
  Venkatraman (2000) paper for the unpaired case.

  The function detects whether the data are paired or unpaired by
  testing which of the options marker2 and group is specified.  If both
  are missing it will stop with an error message.  At present exactly
  one should be missing.
}

\examples{
d <- rep(0:1, 50)
x <- rnorm(100) + 1.2*d
y <- rnorm(100) + 1.2*d
oo <- roc.perm.test(x, d, marker2=y)
plot(oo)
oo <- roc.perm.test(c(x,y), c(d,d), group=rep(1:2,each=100))
plot(oo)
}
\references{
  Venkatraman, E.S. and Begg, C.B. (1996). A distribution-free procedure
  for comparing receiver operating characteristic curves from a paired
  experiment. \emph{Biometrika} 83, 835-848.
  
  Venkatraman, E.S. (2000) A permutation test to compare receiver
  operating characteristic curves. \emph{Biometrics} 56(4):1134-8.
}
