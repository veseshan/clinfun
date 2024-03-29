
\name{gsdesign}
\title{Group Sequential Designs}
\alias{gsdesign}
\alias{gsdesign.binomial}
\alias{gsdesign.normal}
\alias{gsdesign.survival}
\description{
  Functions to calculate sample size for group sequential designs
}
\usage{
gsdesign.binomial(ifrac, pC, pE, r = 1, sig.level = 0.05, power = 0.8,
  delta.eb=0.5, delta.fb = NULL, alternative = c("two.sided",
  "one.sided"), pooled.variance = FALSE, CPS = TRUE, tol=0.00001, ...) 
gsdesign.normal(ifrac, delta, sd = 1, r = 1, sig.level = 0.05,
  power = 0.8, delta.eb = 0.5, delta.fb = NULL, alternative = 
  c("two.sided", "one.sided"), tol=0.00001, ...)
gsdesign.survival(ifrac, haz.ratio, r = 1, sig.level = 0.05, 
  power = 0.8, delta.eb = 0.5, delta.fb = NULL, alternative = 
  c("two.sided", "one.sided"), tol=0.00001, ...)
}
\arguments{
  \item{ifrac}{information fraction or the ratio of current sample size
    or number of events to the total sample size or number of events.
    This should be an increasing vector of numbers from 0 to 1 with the
    last one being 1.  If just 1 is given a fixed sample design is derived.}
  \item{pC}{prob of success of the standard therapy (for binomial data)}
  \item{pE}{prob of success of the experimental therapy (for binomial data)}
  \item{delta}{true difference in means (for normal data)}
  \item{sd}{standard deviation (for normal data)}
  \item{haz.ratio}{hazard ratio (for survival comparison)}
  \item{r}{treatment allocation of r (default=1) experimental per 1 control.}
  \item{sig.level}{significance level (type I error probability)}
  \item{power}{power of test (1 minus type II error probability)}
  \item{delta.eb}{power for efficacy boundary in the Pocock (=0) to
    O'Brien-Fleming (=0.5) family (default is 0.5) }
  \item{delta.fb}{power for futility boundary in the Pocock (=0) to
    O'Brien-Fleming (=0.5) family (default is NULL i.e. no futility
    boundary is requested.) }
  \item{alternative}{one- or two-sided test.}
  \item{pooled.variance}{whether the test statistic is standardized by
    pooled (2*pbar*(1-pbar)) or unpooled variance (pC*(1-pC) + pE*(1-pE)).
    Default is unpooled variance.}
  \item{CPS}{whether continuity correction is used for sample size
    calculation as in Casagrande, Pike & Smith. Default is to use it.}
  \item{tol}{tolerance level for multivariate normal probability
    computation.}
  \item{...}{additional options passed on the pmvnorm function.}
}
\value{
  a list with ifrac, sig.level, power, alternative, delta.eb, delta.fb and:
  \item{efbdry}{the critical value to use at the different looks. For
  two-sided alternative the absolute test statistic should exceed this.}
  \item{futbdry}{the critical value to use at the different looks. For
  two-sided alternative the absolute test statistic should be below this.}
  \item{sample.size}{the sample size per arm for binomial/normal data.}
  \item{num.events}{the total number of failures which should be
    converted to number of subjects using censoring proportion.}  
}
\details{
  The futility boundary is not returned when delta.fb is not specified
  i.e. stopping for futility is not requested.  The futility boundary is
  non-binding.  That is the significance level is not adjusted to account
  for early stopping for futility.  This makes the test a bit conservative 
  in that the true size is less than the nominal level.

  If the alternative is two-sided by default the futility boundary will
  also be two-sided i.e. continuation region is wedge shaped. However,
  if the goal is to show the superiority of the experimental treatment
  then futility boundary should be one sided. This can be achieved by
  deriving the boundaries for one-sided alternative and significance
  level set at half of the value used for two sided alternative. See the
  examples section for a representative design for which the trial
  cannot be stopped at the first look for futility.

  The Casagrande-Pike-Smith type continuity correction is obtained using
  the formula \deqn{n*[1 + \sqrt{1+4/(|pC-pE|*n)}]^2} where n is
  the uncorrected sample size.
}
\examples{
  gsdesign.normal(1:4/4, 0.25, sig.level=0.05, alt="t", delta.fb=0.5)
  gsdesign.normal(1:4/4, 0.25, sig.level=0.025, alt="o", delta.fb=0.5)
}
\keyword{design}
