\name{oc.twostage.bdry}
\title{Two-stage boundary operating characteristics}
\alias{oc.twostage.bdry}
\keyword{design}
\description{
  Calculates the operating characteristics of a two-stage boundary.
}
\usage{
oc.twostage.bdry(pu, pa, r1, n1, r, n)
}
\arguments{
  \item{pu}{unacceptable response rate}
  \item{pa}{response rate that is desirable}
  \item{r1}{first stage threshold to declare treatment undesirable}
  \item{n1}{first stage sample size}
  \item{r}{overall threshold to declare treatment undesirable}
  \item{n}{total sample size}
}
\value{
  oc.twostage.bdry returns the type I and II error rates as well as the
  probability of early temination and expected sample size under pu for
  a specific boundary.
}  
