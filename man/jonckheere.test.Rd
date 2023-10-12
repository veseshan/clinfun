\name{jonckheere.test}
\title{Exact/permutation version of Jonckheere-Terpstra test}
\alias{jonckheere.test}
\description{
  Jonckheere-Terpstra test to test for ordered differences among classes
}
\usage{
  jonckheere.test(x, g, alternative = c("two.sided", "increasing",
                 "decreasing"), nperm=NULL)
}
\arguments{
  \item{x, g}{data and group vector}
  \item{alternative}{means are monotonic (two.sided), increasing, or
    decreasing}
  \item{nperm}{number of permutations for the reference distribution.
    The default is null in which case the permutation p-value is not
    computed. Recommend that the user set nperm to be 1000 or higher if
    permutation p-value is desired.}
}
\details{
  jonckheere.test is the exact (permutation) version of the
  Jonckheere-Terpstra test.  It uses the statistic
  \deqn{\sum_{k<l} \sum_{ij} I(X_{ik} < X_{jl}) + 0.5 I(X_{ik} =
    X_{jl}),} where \eqn{i, j} are observations in groups \eqn{k} and
  \eqn{l} respectively.  The asymptotic version is equivalent to
  cor.test(x, g, method="k"). The exact calculation requires that there
  be no ties and that the sample size is less than 100. When data are
  tied and sample size is at most 100 permutation p-value is returned.
}
\examples{
  set.seed(1234)
  g <- rep(1:5, rep(10,5))
  x <- rnorm(50)
  jonckheere.test(x+0.3*g, g)
  x[1:2] <- mean(x[1:2]) # tied data
  jonckheere.test(x+0.3*g, g)
  jonckheere.test(x+0.3*g, g, nperm=5000)
}
\references{
  Jonckheere, A. R. (1954). A distribution-free k-sample test again
  ordered alternatives. \emph{Biometrika} 41:133-145.

  Terpstra, T. J. (1952). The asymptotic normality and consistency of
  Kendall's test against trend, when ties are present in one ranking.
  \emph{Indagationes Mathematicae} 14:327-333.
}
\keyword{htest}
