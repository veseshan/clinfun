\name{twostage.inference}
\title{Inference following a two-stage design for binary response}
\alias{twostage.inference}
\keyword{design}
\description{
  Calculates the p-value, UMVUE and CI for the data from a study using a
  two stage design for response.
}
\usage{
twostage.inference(x, r1, n1, n, pu, alpha=0.05)
}
\arguments{
  \item{x}{number of responses observed at the end of the study}
  \item{r1}{first stage threshold to declare treatment undesirable}
  \item{n1}{first stage sample size}
  \item{n}{total sample size}
  \item{pu}{unacceptable response rate (null hypothesis)}
  \item{alpha}{the confidence level. For consistency with the design use
    the same value from the design. (default is 0.05)}
}
\value{
  twostage.inference returns the UMVUE (Jung & Kim, 2004), p-value and
  CI (Koyama & Chen, 2008). The CI has confidence level 1-2*alpha and
  the one-sided (1-alpha) interval consistent with the design is
  obtained by changing the upper confidence limit (UCL) to 1.
}
\references{
  Jung SH and Kim KM. (2004). On the estimation of the binomial
  probability in multistage clinical trials. \emph{Statistics in
    Medicine}  23, 881-896.

  Koyama T and Chen H. (2008). Proper inference from Simon's two-stage
  designs. \emph{Statistics in Medicine} 27, 3145-3154.
}
