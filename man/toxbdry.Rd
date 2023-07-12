\name{toxbdry}
\title{Stopping rule for toxicity/futility monitoring}
\alias{toxbdry}
\alias{futilbdry}
\alias{bdrycross.prob}
\alias{print.toxbdry}
\alias{print.futilbdry}
\keyword{design}
\description{
  Computes a stopping rule and its operating characteristics for
  toxicity monitoring based repeated significance testing.
}
\usage{
toxbdry(pLo, pHi, n, cP0=0.1, cP1=0.9, ngrid=6, niter=10, delta=0,
          priority=c("null","alt"))
futilbdry(rLo, rHi, n, size=0.1, power=0.9, ngrid=3, niter=10, delta=0.5)
bdrycross.prob(n, r, ptox)
\method{print}{toxbdry}(x, \dots)
\method{print}{futilbdry}(x, \dots)
}
\arguments{
  \item{pLo}{the toxicity rate that is acceptable.}
  \item{rLo}{baseline (null) response rate.}
  \item{pHi}{the toxicity rate that is too high and hence unacceptable.}
  \item{rHi}{desirable response rate. Stop when it is too unlikely.}
  \item{n}{vector of times (sample size) when toxicty/response is monitored.}
  \item{r}{vector of maximum acceptable toxicities (non-responders for
    futility) corresponding to n.}
  \item{ptox}{the toxicity rates for which the operating characteristics
    are calculated. For futility this is the non-response rate.}
  \item{cP0}{boundary crossing probability under pLo i.e. type I error
    or the probability of declaring a treatment with toxicity rate pLo
    unacceptable.}
  \item{cP1}{boundary crossing probability under pHi i.e. power or the
    probability of declaring a treatment with toxicity rate pHi
    unacceptable.}
  \item{size}{probability of calling drug effective if response rate is rLo.}
  \item{power}{probability of calling drug effective if response rate is rHi.}
  \item{ngrid}{the number of toxicity rates from pLo to pHi for which
    the operating characteristics are computed.}
  \item{niter}{the number of iterations run to obtain the boundary.}
  \item{delta}{power determining the shape of the boundary. Should be
    between 0 (default) and 0.5.}
  \item{priority}{the error threshold to prioritize when the max sample 
    size is too small to have both error thresholds satisfied. Default
    is the null i.e. error under pLo.}
  \item{x}{object returned by the function toxbdry.}
  \item{...}{additional arguments to print.}
}
\value{
  the function returns a list with:
  \item{looks}{when toxicty is monitored - same as input n.}
  \item{lo.bdry}{lower boundary is a vector of maximum acceptable number
    of toxicities corresponding the number of subjects in n.  The
    boundary crossing probability for this is slightly above cP0.}
  \item{hi.bdry}{upper boundary is a vector of maximum acceptable number
    of toxicities corresponding the number of subjects in n.  The
    boundary crossing probability for this is slightly below cP0.}
  \item{bdry.oc}{the operating characteristics i.e the toxicity
    rate, the probability of crossing, stopping (i.e. cross before the
    last observation) and the expected sample size for both the low (lo)
    and high (hi) boundaries.}
  \item{bdry.alpha}{the alpha levels for testing at each look for the
    two boundaries.}

  stopping for toxicity is done when the number of toxicities exceeds
  the boundary i.e. the boundary gives the maximum acceptable number.
}
\details{
  The shape parameter delta is used to determine the stopping boundary
  with 0 corresponding to the Pocock boundary where the same
  significance level is used for all looks and 0.5 corresponding to the
  O'Brien-Fleming boundary which has smaller probability of stopping at
  early looks.

  Default value of delta for toxicity monitoring is 0; value between 0.1
  and 0.2 is a reasonable choice to make it less likely to stop early.
  Default values of delta for futility stopping is 0.5.

  For toxicity monitoring two sets of probabilities - pstop and pcross -
  are given which correspond to probability of stopping early and
  probability of declaring the treatment too toxic with the full
  complement of study subjects accrued and treated.

  For futility monitoring instead two sets of probabilities - pstop and
  peffective - are given corresponding to the probability of stopping
  early for futility and probability of finishing the trial and
  declaring it a success. Note that peffective is the complement of
  pcross in toxicity monitoring.

  The futility boundary can have a -1 in earlier looks which means that
  even zero responses is not sufficient for stopping at that look.

  The exact calculations in this function are done along the lines of
  the method in Chapter 12 of Jennison and Turnbull (2000). Ivanova,
  Qaqish and Schell (2005) have an illustrative paper.
}
\examples{
  toxbdry(0.2, 0.35, c(20,40,60,75))
  toxbdry(0.2, 0.3, c(20,40,60,75), cP0=0.15, cP1=0.8)
  # continuous monitoring
  toxbdry(0.1, 0.3, 2:30)
  # prioritize cP1 error threshold
  toxbdry(0.1, 0.3, 2:25, priority="alt")
}
\references{
  Jennison C and Turnbull BW. (2000).  Group Sequential Methods with
  Applications to Clinical Trials.  \emph{Chapman and Hall/CRC}

  Ivanova A, Qaqish BF and Schell MJ. (2005). Continuous Toxicity
  Monitoring in Phase II Trials in Oncology. \emph{Biometrics} 61, 540-545.
}
