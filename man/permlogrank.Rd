\name{permlogrank}
\title{Permutation version of survdiff}
\alias{permlogrank}
\description{
  Small sample survdiff using permutation reference distributions.
}
\usage{
  permlogrank(formula, data, subset, na.action, rho=0, nperm=5000)
}
\arguments{
  \item{nperm}{number of permutations for the reference distribution}
  \item{formula, data, subset, na.action, rho}{see survdiff for details}
}
\details{
  permlogrank is the permutation version of k-sample survdiff. see
  survdiff in survival package for details.
}
\references{
  Heller G, Venkatraman ES. (1996). Resampling procedures to compare two
  survival distributions in the presence of right censored data.
  \emph{Biometrics} 52:1204-1213.
}
\keyword{htest}
