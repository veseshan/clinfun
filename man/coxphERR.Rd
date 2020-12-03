\name{coxphERR}
\title{Heller Explained Relative Risk}
\alias{coxphERR}
\keyword{survival}
\description{
  Calculates the contribution of a subset of covariates to the explained
  relative risk derived from the full Cox proportional hazards model.
}
\usage{
  coxphERR(phfit, ngamma=NULL)
}
\arguments{
  \item{phfit}{The output from a proportional hazards fit.}
  \item{ngamma}{A vector of indices corresponding to covariates of
    interest. If missing (default), the explained relative risk is
    computed for the full model.} 
}
\value{
  The function coxphERR returns the vector (ERR, se.ERR).
  The first component ERR represents the contribution of a subset of
  covariates to the explained relative risk estimate of the full
  model. If a set of covariates is not provided, then it computes the
  estimate of the full model. The second component se.ERR is the
  standard error of the estimate.
}
\details{
  The object phfit should be the result of a call to coxph with the
  option x=TRUE.
}
\examples{
\dontrun{
 library(survival)
 ovarianfit <- coxph(Surv(futime, fustat) ~ age + resid.ds + rx +
                                           ecog.ps, data=ovarian,x=T) 
 # Compute the explained relative risk (ERR) and
 # its standard error (se.ERR) for the full model.
 coxphERR(ovarianfit)       
 # Compute the contribution of age and ECOG performance status to
 # the explained relative risk. Age and ECOG performance status are 
 # the first and fourth covariates in the model.
 coxphERR(ovarianfit, c(1,4))
}}
\references{
  Heller G. (2012) A measure of explained risk in the proportional hazards
  model. \emph{Biostatistics}
}
