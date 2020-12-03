\name{coxphCPE}
\title{Gonen \& Heller Concordance Probability Estimate}
\alias{coxphCPE}
\keyword{survival}
\description{
  Calculates the Concordance Probability Estimate for a Cox model.
}
\usage{
  coxphCPE(phfit)
}
\arguments{
  \item{phfit}{output from a proportional hazards fit.}
}
\value{
  coxphCPE returns a vector with CPE, smooth.CPE & se.CPE which are the
  estimate, the smoothed estimate and its standard error respectively.
}  
\examples{
\dontrun{ library(survival)
  data(pbc)
  pbcfit <- coxph(Surv(time, status==2) ~ trt + log(copper), pbc,
    subset=(trt>0 & copper>0))
  coxphCPE(pbcfit)
}}
\references{
  Gonen M and Heller G. (2005) Concordance probability and discriminatory
  power in proportional hazards regression. \emph{Biometrika} 92, 965-970.
}
