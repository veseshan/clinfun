\name{calogrank}
\title{Survival curves analysis of covariance}
\alias{calogrank}
\keyword{htest}
\description{
  Logrank test to compare survival curves adjusting for covariates
}
\usage{
  calogrank(ftime, fstatus, grp, cvt, strat=NULL)
}
\arguments{
  \item{ftime}{failure times}
  \item{fstatus}{status indicator}
  \item{grp}{group indicator}
  \item{cvt}{continuous covariates used for adjusted analysis}
  \item{strat}{stratification variable}
}
\details{
  calogrank is the covariate adjusted version of k-sample survdiff.  The
  function in its current form only does basic error checking.
}
\examples{
\dontrun{  library(survival)
  data(pbc)
  pbc1 <- pbc
  pbc1$trt[pbc1$trt == -9] <- NA
  pbc1$copper[pbc1$copper == -9] <- NA
  calogrank(pbc1$time, pbc1$status, pbc1$trt, pbc1[,c("copper")])
  calogrank(pbc1$time, pbc1$status, pbc1$trt,
                                  pbc1[,c("protime", "copper")])}
}
\references{
  Heller G. and Venkatraman E.S. (2004) A nonparametric test to compare
  survival distributions with covariate adjustment. \emph{JRSS-B} 66,
  719-733.
}
\keyword{htest}
