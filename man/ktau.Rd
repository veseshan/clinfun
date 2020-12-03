\name{ktau}
\title{Kendall's tau-b estimate}
\alias{ktau}
\keyword{multivariate}
\description{
  Calculates the Kendall's tau-b.
}
\usage{
  ktau(x, y)
}
\arguments{
  \item{x}{first variable}
  \item{y}{second variable}
}
\value{
  ktau returns Kendall's tau-b.
}
\details{
  ktau computes the same quantity as cor(x, y, method="kendall").
  It uses a faster algorithm than pairwise comparisons used by cor.
}
\examples{
  set.seed(1234)
  x <- rnorm(10000); y <- x+rnorm(10000)
  cor(x, y, method="k")
  clinfun:::ktau(x,y)  
}
