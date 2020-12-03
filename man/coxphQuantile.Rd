\name{coxphQuantile}
\title{Survival time quantile as a function of covariate}
\alias{coxphQuantile}
\keyword{survival}
\description{
  Draws a quantile curve of survival distribution as a function of covariate.
}
\usage{
  coxphQuantile(phfit, xrange, p=0.5, whichx=1, otherx=NULL, ...)
}
\arguments{
  \item{phfit}{output from a proportional hazards fit.}
  \item{xrange}{the range of covariate values for which the quantiles of
    survival times are computed.}
  \item{p}{the probability level for the quantile (default is median).}
  \item{whichx}{if there are more than one covariates in the Cox model,
    the one chosen for the quantile plot.}
  \item{otherx}{the values for other covariates in the Cox model.  If
    missing uses their average values.}
  \item{...}{additional parameters to be passed on to the lines command.}
}
\details{
  This function is used to draw quantile curves.  It requires a plot of
  the data (time & covariate of interest) to be present.  See example.

  It invisibly returns the observed failure times and the covariate
  values at which the estimated survival probability is (exactly) p.
}
\examples{
\dontrun{  library(survival)
data(pbc)
pbcfit <- coxph(Surv(time, status==2) ~ trt + log(copper), pbc,
                      subset=(trt>0 & copper>0)) 
plot(log(pbc$copper[pbc$trt>0 & pbc$copper>0]), pbc$time[pbc$trt>0 &
  pbc$copper>0], pch=c("o","x")[1+pbc$status[pbc$trt>0 & pbc$copper>0]], 
  xlab="log Copper", ylab="Survival time")
coxphQuantile(pbcfit, c(2.5,6), whichx=2, otherx=1)
coxphQuantile(pbcfit, c(2.5,6), p=0.75, whichx=2, otherx=2, col=2)}
}
\references{
Heller G. and Simonoff J.S. (1992) Prediction in censored survival data:
A comparison of the proportional hazards and linear regression models.
\emph{Biometrics} 48, 101-115.
}
