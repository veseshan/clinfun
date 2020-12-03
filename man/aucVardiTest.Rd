\name{aucVardiTest}
\title{Two-Sample Tests for Growth Curves under Dependent Right Censoring}
\alias{aucVardiTest}
\description{
  Permutation test for comparing growth curves across tow groups under
  dependent right censoring.
}
\usage{
  aucVardiTest(meas, grp, tim=NULL, cgrps=NULL, nperm=5000)
}
\arguments{
  \item{meas}{Matrix of measurements where the rows are the subjects and
    columns the timepoints. At least one value should not be missing in
    each row. For example they can be tumor sizes measured over time.} 
  \item{grp}{Group indicator for each subject. There must be at least
    two different groups. This can represent each subject's treatment.}
  \item{tim}{Times at which the measurements in \code{meas} are taken.
    If missing, the times are set to 1 through \code{ncol(meas)}.} 
  \item{cgrps}{The two groups that are being compared. If missing the
    first two groups will be compared.}
  \item{nperm}{Number of permutations for the reference distribution.}
}
\value{
  returns a list with objects ostat, pstat and p.value which are the
  observed test statistic for the two groups being compared, values of
  the statistics when the group labels are permuted
}
\examples{
  grp <- sample(1:3, 100, replace=TRUE)
  grp0 <- LETTERS[grp]
  maxfup <- sample(5:20, 100, replace=TRUE)
  meas <- matrix(NA, 100, 20)
  for(i in 1:100) {
    meas[i, 1:maxfup[i]] <- cumsum((3+0.04*grp[i]) + rnorm(maxfup[i]))
  }
  aucVardiTest(meas, grp)
  aucVardiTest(meas, grp0, cgrps=c("C","B"))
}
\details{
  The test statistic is defined as the sum of pairwise differences in
  the partial areas under the growth curve. For each pair of subjects
  the partial area is computed until the smaller of the maximum followup
  times. For each subject, linear interpolation is is used to fill-in
  missing values prior to the maximum followup time. The reference
  distribution of obtained by permuting the group labels.
}
\references{
  Vardi Y., Ying Z. and Zhang C.H. (2001). Two-Sample Tests for Growth
  Curves under Dependent Right Censoring. \emph{Biometrika} 88, 949-960. 
}
\keyword{htest}
