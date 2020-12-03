\name{power.ladesign}
\title{Power of k-sample rank test under Lehmann alternative}
\alias{power.ladesign}
\alias{print.ladesign}
\description{
  Functions to calculate the power of rank tests for animal studies.
}
\usage{
  power.ladesign(gsize, odds.ratio, sig.level = 0.05, statistic =
     c("Kruskal-Wallis", "Jonckheere-Terpstra"), alternative =
     c("two.sided", "one.sided"), nrep=1e+6) 
  \method{print}{ladesign}(x,\dots)
}
\arguments{
  \item{gsize}{sample size of the k (= length of vector) groups.}
  \item{odds.ratio}{odds ratio parameters for the k-1 groups.  The first
    group is considered the control.}
  \item{sig.level}{the significance level of the test (default = 0.05)}
  \item{statistic}{the test statistic for the k-group comparison. Is one
    of Kruskal-Wallis (default) or Jonckeere-Terpstra.}
  \item{alternative}{one- or two-sided test.  Valid only for the
    Jonckheere-Terpstra test.}
  \item{nrep}{number of reps (default 1 million) for Monte Carlo.}
  \item{x}{object of class ladesign returned by power.ladesign}
  \item{...}{arguments to be passed on left for S3 method consistency.}
}
\value{
  returns a list with objects group.size, odds.ratio, statistic,
  sig.level and power.   The "print" method formats the output.
}
\examples{
  power.ladesign(c(9,7), 4, statistic="K")
  power.ladesign(c(9,7,9), c(2,4), statistic="J")
  power.ladesign(c(9,7,9), c(2,4), statistic="J", alt="o")
}
\details{
  Although the power for Jonckheere-Terpstra test is calculated for any
  set of odds ratio, the test is meant for monotone alternative.  Thus
  it is preferable to specify odds ratios that are monotonically
  increasing with all values larger than 1 or decreasing with all values
  smaller than 1.
}
\references{
  Heller G. (2006).  Power calculations for preclinical studies using a
  K-sample rank test and the Lehmann alternative hypothesis.
  \emph{Statistics in Medicine} 25, 2543-2553. 
}
\keyword{design}
