\name{ph2single}
\title{Exact single stage Phase II design}
\alias{ph2single}
\keyword{design}
\description{
  Calculates the exact one stage Phase II design
}
\usage{
ph2single(pu,pa,ep1,ep2,nsoln=5)
}
\arguments{
  \item{pu}{unacceptable response rate}
  \item{pa}{response rate that is desirable}
  \item{ep1}{threshold for the probability of declaring drug desirable
    under p0}
  \item{ep2}{threshold for the probability of rejecting the drug under
    p1}
  \item{nsoln}{number of designs with given alpha and beta}
}
\value{
  ph2single returns a data frame with variables: n, r, and the Type I
  and Type II errors. Treatment desirable if >r respenses seen.
}
