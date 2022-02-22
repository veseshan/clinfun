\name{twostage.admissible}
\title{Admissible design options between Minimax and Optimal}
\alias{twostage.admissible}
\keyword{design}
\description{
  Lists the admissible design options between
}
\usage{
twostage.admissible(x)
}
\arguments{
  \item{x}{output from \code{ph2simon} call}
}
\value{
  twostage.admissible returns design options that are admissible (Jung
  et al, 2004). The output is a matrix with 8 columns: r1, n1, r, n,
  EN(p0), PET(p0), qLo, qHi. The columns qLo and qHi give the range of
  probability values for which the particular design is admissible. 
}
\references{
  Jung SH, Lee T, Kim K, and George, SL. (2004). Admissible two-stage
  designs for phase II cancer clinical trials. \emph{Statistics in 
  medicine} 23(4), 561-569.
}
\examples{
  oo = ph2simon(0.5, 0.7, 0.05, 0.1)
  twostage.admissible(oo)
}
