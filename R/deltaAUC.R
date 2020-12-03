deltaAUC <- function(y, x, z) {
    n <- length(y)    # number of observations
    nn <- n - sum(y)  # number of subjects with y == 0
    res <- list()     # initialize output
    # logistic coefficients of full model
    coefs <- glm.fit(cbind(1, x, z), y, family=binomial())$coefficients[-1]
    names(coefs) <- NULL
    # anchor variable is the one with best AUC
    # xauc <- apply(x[order(y),], 2, rocauc, n, nn)
    # divide by anchor variable (first column of x)
    coefs.f <- coefs[-1]/coefs[1]
    # get the maximum rank correlation (MRC) fit
    xmat <- cbind(x,z)[order(y),]
    p <- ncol(xmat)
    mrcfull <- optim(coefs.f, fn=mrcobj, xmat=xmat, n=n, nn=nn, control=list(fnscale=-1))
    # optim message
    optim.msg <- list(full=mrcfull[c("convergence","message")])
    # compute the linear predictor
    xbfull <- c(xmat %*% c(1, mrcfull$par))
    # print(summary(xbfull))
    # calculate bandwidth for smoothed AUC, variance etc.
    seu <- sqrt(2*var(xbfull)*n/(n-1))
    bwcdf <- seu/n^(1/3)
    bwpdf <- seu/n^(1/5)
    # re-estimate coefficients using smoothed AUC
    smmrcfull <- optim(coefs.f, fn=smmrcobj, xmat=xmat, n=n, nn=nn, bw=bwcdf, control=list(fnscale=-1))
    # print(smmrcfull)
    # full model coefficients and AUC
    res$par.full <- c(1, smmrcfull$par)
    # print(res)
    # compute the linear predictor
    xbfull <- c(xmat %*% res$par.full)
    # empirical and smoothed AUC full model
    area.full <- rocauc(xbfull, n, nn)
    smarea.full <- smrocauc(xbfull/bwcdf, n, nn)
    # optim message
    optim.msg$smfull <- smmrcfull[c("convergence","message")]
    # logistic coefficients of reduced model (check if >1 covariate)
    xmat0 <- cbind(x)[order(y), , drop=FALSE]
    p0 <- ncol(xmat0)
    if (p0 > 1) {
        coefs <- glm.fit(cbind(1, x), y, family=binomial())$coefficients[-1]
        names(coefs) <- NULL
        # divide by anchor variable (first column of x)
        coefs.r <- coefs[-1]/coefs[1]
        # get the maximum rank correlation (MRC) fit
        if (p0 > 2) {
            mrcred <- optim(coefs.r, smmrcobj, xmat=xmat0, n=n, nn=nn, bw=bwcdf, control=list(fnscale=-1))
        } else {
            mrcred <- optim(coefs.r, smmrcobj, xmat=xmat0, n=n, nn=nn, bw=bwcdf, control=list(fnscale=-1), lower=-100, upper=100, method="Brent")
        }
        optim.msg$red <- mrcred[c("convergence","message")]
        # reduced model coefficients
        res$par.red <- c(1,mrcred$par)
        # compute the linear predictor
        xbred <- c(xmat0 %*% res$par.red)
        # print(summary(xbred))
    } else {
        # linear predictor is just the single variable in x
        xbred <- c(xmat0)
        # reduced model coefficients and AUC
        res$par.red <- 1
    }
    # empirical and smoothed AUC reduced model
    area.red <- rocauc(xbred, n, nn)
    smarea.red <- smrocauc(xbred/bwcdf, n, nn)
    # print(res)
    # D and W matrix for hypothesis test and variance under non-null
    p1 <- p - 1 # number of columns in full model minus 1
    zzz <- .Fortran("daucmats",
                    n=as.integer(n),
                    p=as.integer(p1),
                    n0=as.integer(nn),
                    x=as.double(xmat[,-1]/bwpdf),
                    xb=as.double(xbfull/bwpdf),
                    xb0=as.double(xbred/bwcdf),
                    delta=as.double(smarea.full-smarea.red),
                    dmat=double(p1^2),
                    wmat=double(p1^2),
                    valt=double(1))
    dmat <- matrix(zzz$dmat, p1, p1)
    dinv <- solve(dmat)
    vmat <- dinv %*% matrix(zzz$wmat, p1, p1) %*% dinv
    # p0 is variables in reduced model; so added variables start at p0+1
    # since first column is dropped in dmat, wmat etc. it's p0 instead
    # negative sign used for the matrix so that eigen values are correct
    vdinv <- -vmat[p0:p1, p0:p1] %*% solve(dinv[p0:p1, p0:p1])
    # test statistic an p-value using sum of chi-squares
    veigen <- eigen(vdinv)$values
    xstat <- rep(0, 99999)
    for (ve in veigen) xstat <- xstat + ve*rnorm(99999)^2
    tstat <- 2*n*c(area.full - area.red, smarea.full - smarea.red)
    pval <- c(1+sum(xstat >= tstat[1]), 1+sum(xstat >= tstat[2]))/1e5
    res$results <- cbind(c(area.full, smarea.full), c(area.red, smarea.red), tstat, pval)
    colnames(res$results) <- c("AUCfull", "AUCreduced", "statistic","p-value")
    rownames(res$results) <- c("Empirical","Smooth")
    res
}

mrcobj <- function(coefs, xmat, n, nn) {
    marker <- xmat %*% c(1, coefs)
    rocauc(marker, n, nn)
}

smmrcobj <- function(coefs, xmat, n, nn, bw) {
    marker <- c(xmat %*% c(1, coefs))/bw
    smrocauc(marker, n, nn)
}

rocauc <- function(marker, n, nn) {
    zzz <- .Fortran("rocauc",
                    as.integer(n),
                    as.integer(nn),
                    as.double(marker),
                    area=double(1))
    zzz$area
}

smrocauc <- function(marker, n, nn) {
    zzz <- .Fortran("smrocauc",
                    as.integer(n),
                    as.integer(nn),
                    as.double(marker),
                    area=double(1))
    zzz$area
}
