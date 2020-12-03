aucVardiTest <- function(meas, grp, tim=NULL, cgrps=NULL, nperm=5000) {
    # check that meas has at least one value for each subject
    if (any(apply(meas, 1, function(x) {sum(is.finite(x))}) == 0))
        stop("every row of meas should have at least one non missing value")
    # change grp into a factor variable
    grp0 <- as.factor(grp)
    # only two groups can be compared
    lgrp <- levels(grp0) # group labels
    ngrp <- length(lgrp) # number of groups
    if (ngrp < 2) stop("grp should have at least 2 levels")
    # if tim is not give assume they are equally spaced
    if (missing(tim)) tim <- 1:ncol(meas)
    # if missing compgrps the levels 1 & 2 are compared
    if (missing(cgrps)) cgrps <- lgrp[1:2]
    # exactly 2 groups are compared and should be in lgrp
    if (length(intersect(cgrps, lgrp)) != 2) stop("\n misspecified cgrps; specify two values from: ", paste(lgrp, collapse=", "))
    refgrp <- cgrps[1]
    # number of subjects with refgrp
    nref <- sum(grp == refgrp)
    # number of timepoints
    ntim <- length(tim)
    # add time 0
    tim0 <- c(0, tim)
    # restrict data to the 2 groups being compared
    meas0 <- meas[grp0 %in% cgrps, ]
    grp <- grp0[grp0 %in% cgrps]
    # interpolate missing values until the last observed time
    # and return partial auc contributions upto each timepoint
    pauc <- t(apply(cbind(0, meas0), 1, function(x, tim, dtim, ntim) {
        k <- max(which(is.finite(x)))
        if (length(is.na(x[1:k])) > 0) {
            ii <- which(is.finite(x[1:k]))
            jj <- (1:k)[-ii]
            # approx function does linear interpolation
            x[jj] <- approx(tim[ii], x[ii], tim[jj])$y
        }
        cumsum((x[1:ntim]+x[2:(ntim+1)])*dtim/2)
    }, tim0, diff(tim0), ntim))
    # max followup time for each subject
    maxfup <- apply(pauc, 1, function(x) {sum(is.finite(x))})
    # pairwise difference in partial auc
    n <- nrow(pauc)
    nn <- 1:n
    pwaucdiff <- rep(0, n*n)
    l <- 0
    for (i in 1:n) {
        for(j in 1:n) {
            l <- l+1
            k <- min(maxfup[i], maxfup[j])
            pwaucdiff[l] <- pauc[i, k] - pauc[j, k]
        }
    }
    # initialize index of relevant (non-ref, ref) pairs
    ij <- rep(0, nref*(n-nref))
    # locations of (non-ref, ref) pairs under original group allocation 
    jj <- which(grp == refgrp)
    l <- 0
    for (i in which(grp != refgrp)) {
        ij[l + 1:nref] <- (i-1)*n + jj
        l <- l + nref
    }
    # observed test statistic
    ostat <- sum(pwaucdiff[ij])
    # permute the grp variable and calculate statistic
    pgrp <- grp
    pstat <- rep(0, nperm)
    for (k in 1:nperm) {
        # locations of (non-ref, ref) pairs under permuted group allocation 
        jj <- sort(sample(nn, nref))
        l <- 0
        for (i in nn[-jj]) {
            ij[l + 1:nref] <- (i-1)*n + jj
            l <- l + nref
        }
        pstat[k] <- sum(pwaucdiff[ij])
    }
    p.value <- (sum(abs(pstat) >= abs(ostat))+1)/(nperm+1)
    cat("statistic for", cgrps[2], "-", cgrps[1], " = ", ostat,"; p-value = ", p.value,"\n")
    invisible(list(ostat=ostat, pstat=pstat, p.value=p.value))
}
