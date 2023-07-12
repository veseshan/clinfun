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
  \item{pu}{unacceptable response rate; baseline response rate that needs to
    be exceeded for treatment to be deemed promising}
  \item{pa}{response rate that is desirable; should be larger than pu}
  \item{ep1}{threshold for the probability of declaring drug desirable
    under pu (target type 1 error rate); between 0 and 1}
  \item{ep2}{threshold for the probability of rejecting the drug under
    pa (target type 2 error rate); between 0 and 1}
  \item{nmax}{maximum total sample size (default 100; can be at most 1000)}
  \item{x}{object returned by ph2simon}
  \item{...}{arguments to be passed onto plot and print commands called
    within}
}
\value{
  ph2simon returns a list with pu, pa, alpha, beta and nmax as above
  and: \item{out}{matrix of best 2 stage designs for each value of total
    sample size n.  The 6 columns in the matrix are:
    \tabular{rl}{
      r1 \tab number of responses needed to exceeded in first stage \cr
      n1 \tab number of subjects treated in first stage \cr
      r \tab number of responses needed to exceeded at the end of trial \cr 
      n \tab total number of subjects to be treated in the trial \cr
      EN(pu) \tab expected number pf patients in the trial under pu \cr
      PET(pu) \tab probability of stopping after the first stage under pu 
    }}

  Trial is stopped early if <= r1 responses are seen in the first stage
  and treatment is considered desirable only when >r responses seen.
}
\section{Methods (by generic)}{
  \itemize{
    \item \code{print(ph2simon)}: formats and returns the minimax,
    optimal and any admissible designs.
    \item \code{plot(ph2simon)}: plots the expected sample size against
    the maximum sample size as in Jung et al., 2001
  }
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

\seealso{
  \code{\link{twostage.inference}}, \code{\link{oc.twostage.bdry}}
}

\keyword{design}
