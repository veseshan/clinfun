##  
##  coded by E. S. Venkatraman based on the R code by Gonen 02/08/2007
##  modified to add version with ties (Heller & Mo 2016, Lifetime Data Analysis)
##

coxphCPE <- function(phfit, out.ties=FALSE) {
  if (!inherits(phfit, "coxph")) stop("phfit shoud be coxph class object")
  n <- phfit$n
  betahat <- phfit$coefficients
  p <- length(phfit$coefficients)
  vbetahat <- phfit$var
  xbeta <- phfit$linear.predictors
  # sort xbeta to make check lesser/greater check redundant
  ii = order(xbeta, decreasing=TRUE)
  xbeta = xbeta[ii]
  # removed the R 2.9.1 check
  xMat <- as.matrix(model.matrix(phfit)[ii,])
  # bandwidth for kernel smoothing
  bw <- 0.5*sd(xbeta)*(n^(-1/3))
  # call Fortran code for CPE and variance terms
  if (out.ties) {
    zzz <- .Fortran("cpesubt",
                    as.integer(n),
                    as.integer(p),
                    as.double(xMat),
                    as.double(xbeta),
                    nPI=double(1),
                    CPE=double(1),
                    varDeriv=double(p),
                    uRowSum=double(n),
                    uSSQ=double(1),
                    PACKAGE="clinfun")
    CPE <- zzz$CPE
    CPEsmooth <- NA # no smooth counterpart in discrete case
    varTerm1 <- (sum(zzz$uRowSum^2) - zzz$uSSQ)/zzz$nPI^2
    varDeriv <- zzz$varDeriv/zzz$nPI
    varTerm2 <- t(varDeriv)%*%vbetahat%*%varDeriv
  } else {
    zzz <- .Fortran("cpesub",
                    as.integer(n),
                    as.integer(p),
                    as.double(xMat),
                    as.double(xbeta),
                    as.double(bw),
                    CPE=double(1),
                    CPEsmooth=double(1),
                    varDeriv=double(p),
                    uRowSum=double(n),
                    uSSQ=double(1),
                    PACKAGE="clinfun")
    CPE <- 2*zzz$CPE/(n*(n-1))
    CPEsmooth <- 2*zzz$CPEsmooth/(n*(n-1))
    varTerm1 <- 4*(sum((zzz$uRowSum+rep(0.5,n)-n*CPEsmooth)^2) - (zzz$uSSQ + n/4 - n*CPEsmooth - n*(n-2)*CPEsmooth^2))/(n*(n-1))^2
    varDeriv <- 2*zzz$varDeriv/(n*(n-1))
    varTerm2 <- t(varDeriv)%*%vbetahat%*%varDeriv
  }
  # estimated variance of CPE estimate
  varCPE <- varTerm1 + varTerm2
  # results
  out <- c(CPE, CPEsmooth, sqrt(varCPE))
  names(out) <- c("CPE", "smooth.CPE", "se.CPE")
  out
}
