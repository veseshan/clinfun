\name{ph2simon}
\title{Simon's 2-stage Phase II design}
\alias{ph2simon}
\alias{print.ph2simon}
\alias{plot.ph2simon}
\keyword{design}
\description{
  Calculates Optimal and Minimax 2-stage Phase II designs given by
  Richard Simon
}
\usage{
ph2simon(pu, pa, ep1, ep2, nmax=100)
\method{print}{ph2simon}(x, \dots)
\method{plot}{ph2simon}(x, \dots)
}
\arguments{
  \item{pu}{unacceptable response rate}
  \item{pa}{response rate that is desirable}
  \item{ep1}{threshold for the probability of declaring drug desirable
    under p0}
  \item{ep2}{threshold for the probability of rejecting the drug under
    p1}
  \item{nmax}{maximum total sample size (default 100; can be at most
    500)}
  \item{x}{object returned by ph2simon}
  \item{...}{arguments to be passed onto plot and print commands called
    within}
}
\value{
  ph2simon returns a list with pu, pa, alpha, beta and nmax as above
  and:
  \item{out}{matrix of best 2 stage designs for each value of total
    sample size n.  the 6 columns are: r1, n1, r, n, EN(p0), PET(p0)}

  Trial is stopped early if <= r1 responses are seen in the first stage
  and treatment is considered desirable only when >r responses seen.

  The "print" method formats and returns the minimax and optimal
  designs.  The "plot" plots the expected sample size agains the maximum
  sample size as in Jung et al., 2001
}

\seealso{
  \code{\link{twostage.inference}}, \code{\link{oc.twostage.bdry}}
}

\examples{
  ph2simon(0.2, 0.4, 0.1, 0.1)
  ph2simon(0.2, 0.35, 0.05, 0.05)
  ph2simon(0.2, 0.35, 0.05, 0.05, nmax=150)
}

\references{
  Simon R. (1989).  Optimal Two-Stage Designs for Phase II Clinical
  Trials. \emph{Controlled Clinical Trials} 10, 1-10.

  Jung SH, Carey M and Kim KM. (2001). Graphical Search for Two-Stage
  Designs for Phase II Clinical Trials. \emph{Controlled Clinical
    Trials}  22, 367-372.
}

\keyword{design}
