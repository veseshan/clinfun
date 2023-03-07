\name{coxphCPE}
\title{Concordance Probability Estimate for Cox model}
\alias{coxphCPE}
\keyword{survival}
\description{
  Calculates the Concordance Probability Estimate for a Cox proportional 
  hazards model. Both the Gonen and Heller (2005) version for continuous 
  risk score and Heller and Mo (2016) for discrete risk score can be
  calculated.
}
\usage{
  coxphCPE(phfit, out.ties=FALSE)
}
\arguments{
  \item{phfit}{output from a proportional hazards fit.}
  \item{out.ties}{binary flag to decide if pairs with tied risk scores
    should be used.}
}
\value{
  coxphCPE returns a vector with CPE, smooth.CPE and se.CPE which are the
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

  Heller G and Mo Q. (2016). Estimating the concordance probability in a
  survival analysis with a discrete number of risk groups.
  \emph{Lifetime Data Analysis}, 22,263-79.
}
