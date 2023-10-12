\name{clinfun-internal}
\title{Internal clinfun functions}
\alias{compareROC.paired}
\alias{compareROC.unpaired}
\alias{exactNull}
\alias{getBoundary}
\alias{gsd.bdryconstant}
\alias{gsd.drift.efficacy}
\alias{gsd.drift.both}
\alias{gsd.drift}
\alias{mrcobj}
\alias{smmrcobj}
\alias{rocauc}
\alias{smrocauc}
\description{
 Internal functions used by the functions in the package.
}
\usage{
compareROC.paired(x,y,d,nperm=2500)
compareROC.unpaired(x, dx, y, dy, nperm=2500, mp=NULL)
exactNull(gsize, kw)
getBoundary(pLo, pHi, n, cP0=0.1, cP1=0.9, ngrid=6, niter=10, delta=0,
    priorityNull=TRUE) 
gsd.bdryconstant(ifrac, eprob = 0.05, delta = 0.5,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
gsd.drift.efficacy(ifrac, delta.eb, sig.level=0.05, pow = 0.8,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
gsd.drift.both(ifrac, delta.eb, delta.fb, sig.level = 0.05, pow=0.8,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
gsd.drift(ifrac, sig.level = 0.05, pow = 0.8, delta.eb = 0.5, delta.fb = NULL,
    alternative = c("two.sided", "one.sided"), tol = 0.00001, ...)
mrcobj(coefs, xmat, n, nn)
smmrcobj(coefs, xmat, n, nn, bw)
rocauc(marker, n, nn)
smrocauc(marker, n, nn)
}
\details{
  These functions are not to be called by the user.
}
\keyword{internal}
