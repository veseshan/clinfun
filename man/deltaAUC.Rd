\name{deltaAUC}
\title{Comparing the AUC from ROC curves from nested binary regression}
\alias{deltaAUC}
\alias{ipmn}
\keyword{htest}
\description{
  Conducts the test 
}
\usage{
  deltaAUC(y, x, z)
}
\arguments{
  \item{y}{binary response variable}
  \item{x}{matrix of set of covariates that is the basis of the existing
    (reduced) model}
  \item{z}{matrix of set of covariates that are added to to get the new
    (full) model}
}
\value{
  It returns a list with the following elements
  \item{par.full}{the MRC estimate of parameters for the full model}
  \item{par.red}{the MRC estimate of parameters for the reduced model}
  \item{results}{matrix od results which gives the full reduced model
    AUCs along with the test statistic and p-value}
}
\details{
  The models are fit using maximum rank correlation (MRC) method which
  is an alternate approach to logistic regression. In MRC the area under
  the ROC curve (AUC) is maximized as opposed to the likelihood in
  logistic regression. Due to invariance of AUC to location and scale 
  shifts one of the parameters (anchor variable) is set to 1.

  The first variable (column) in x is used as the anchor variable.

  The IPMN data set used as an example in the paper below is included.
  The columns are high risk lesion (V1), recent weight loss (V2), main
  duct involvement (V4), presence of a solid component in imaging (V3),
  and lesion size (V5).
}
\examples{
  data(ipmn)
  deltaAUC(ipmn$V1, cbind(ipmn$V4, ipmn$V3, ipmn$V5), ipmn$V2)
}
\references{
  Heller G., Seshan V.E., Moskowitz C.S. and Gonen M. (2016) Inference
  for the difference in the area under the ROC curve derived from nested
  binary regression models. \emph{Biostatistics} 18, 260-274.
}
