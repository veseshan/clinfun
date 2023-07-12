\name{fedesign}
\title{Trial Designs Based On Fisher's Exact Test}
\alias{fedesign}
\alias{fe.ssize}
\alias{CPS.ssize}
\alias{fe.mdor}
\alias{mdrr}
\alias{fe.power}
\alias{or2pcase}
\keyword{design}
\description{Calculates sample size, effect size and power based on
  Fisher's exact test}
\usage{
fe.ssize(p1, p2, alpha=0.05, power=0.8, r=1, npm=5, mmax=1000)
CPS.ssize(p1, p2, alpha=0.05, power=0.8, r=1)
fe.mdor(ncase, ncontrol, pcontrol, alpha=0.05, power=0.8)
mdrr(n, cprob, presp, alpha=0.05, power=0.8, niter=15)
fe.power(d, n1, n2, p1, alpha = 0.05)
or2pcase(pcontrol, OR)
}
\arguments{
  \item{p1}{response rate of standard treatment}
  \item{p2}{response rate of experimental treatment}
  \item{d}{difference = p2-p1}
  \item{pcontrol}{control group probability}
  \item{n1}{sample size for the standard treatment group}
  \item{n2}{sample size for the standard treatment group}
  \item{ncontrol}{control group sample size}
  \item{ncase}{case group sample size}
  \item{alpha}{size of the test (default 5\%)}
  \item{power}{power of the test (default 80\%)}
  \item{r}{treatments are randomized in 1:r ratio (default r=1)}
  \item{npm}{the sample size program searches for sample sizes in a
    range (+/- npm) to get the exact power}
  \item{mmax}{the maximum group size for which exact power is
    calculated}
  \item{n}{total number of subjects}
  \item{cprob}{proportion of patients who are marger positive}
  \item{presp}{probability of response in all subjects}
  \item{niter}{number of iterations in binary search}
  \item{OR}{odds-ratio}
}
\details{
CPS.ssize returns Casagrande, Pike, Smith sample size which is a very
close to the exact.  Use this for small differences p2-p1 (hence large
sample sizes) to get the result instantaneously.

Since Fisher's exact test orders the tables by their probability the
test is naturally two-sided.

fe.ssize return a 2x3 matrix with CPS and Fisher's exact sample sizes
with power.

fe.mdor return a 3x2 matrix with Schlesselman, CPS and Fisher's exact
minimum detectable odds ratios and the corresponding power.

fe.power returns a Kx2 matrix with probabilities (p2) and exact power.

mdrr computes the minimum detectable P(resp|marker+) and P(resp|marker-)
configurations when total sample size (n), P(response) (presp) and
proportion of subjects who are marker positive (cprob) are specified.

or2pcase give the probability of disease among the cases for a given
probability of disease in controls (pcontrol) and odds-ratio (OR).
}
\references{
  Casagrande, JT., Pike, MC. and Smith PG. (1978). An improved
  approximate formula for calculating sample sizes for comparing two
  binomial distributions. \emph{Biometrics} 34, 483-486.

  Fleiss, J. (1981) Statistical Methods for Rates and Proportions.

  Schlesselman, J. (1987) Re: Smallest Detectable Relative Risk With
  Multiple Controls Per Case.  \emph{Am. J. Epi.}
}
