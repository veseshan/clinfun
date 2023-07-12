power.ladesign <- function(gsize, odds.ratio, sig.level = 0.05, statistic = c("Kruskal-Wallis", "Jonckheere-Terpstra"), alternative = c("two.sided", "one.sided"), nrep=1e+6) {
  ng <- length(gsize)
  nOR <- length(odds.ratio)
  if (nOR != ng-1) {
    errmesg <- paste(" There are", ng, "groups and", nOR,"odds ratios, there should be", ng-1)
    stop(message=errmesg)
  }
  statistic <- match.arg(statistic)
  kw <- statistic=="Kruskal-Wallis"
  alternative <- match.arg(alternative)
  if (kw & alternative == "one.sided") stop("Only two-sided alternative is allowed for Kruskal-Wallis test")
  ng <- length(gsize)
  # number of possible combinations
  ncomb <- exp(lgamma(sum(gsize)+1) - sum(lgamma(gsize+1)))
  if (ncomb > 1e5) {
      zz <- .Fortran("lehman",
                     as.integer(ng),
                     as.integer(gsize),
                     double(ng),
                     as.double(rep(1, ng)),
                     as.double(sum(gsize)),
                     double(ng),
                     as.logical(kw),
                     as.integer(nrep),
                     tstat=double(nrep))
  } else {
      zz <- exactNull(gsize, kw)
  }
  if (kw) {
      # round the statistic to 3 significant digits after decimal
      ldq <- quantile(round(zz$tstat, 3), 1-sig.level)
  } else {
      if (alternative=="one.sided") {
          ldq <- quantile(zz$tstat, c(sig.level, 1-sig.level))
      } else {
          ldq <- quantile(zz$tstat, c(sig.level/2, 1-sig.level/2))
      }
  }
  zz <- .Fortran("lehman",
                 as.integer(ng),
                 as.integer(gsize),
                 double(ng),
                 as.double(c(1,odds.ratio)),
                 as.double(sum(gsize*c(1,odds.ratio))),
                 double(ng),
                 as.logical(kw),                 
                 as.integer(nrep),
                 tstat=double(nrep))
  if (kw) {
      # round the statistic to 3 significant digits after decimal
      ldpow <- sum(round(zz$tstat, 3) >= ldq)/nrep
  } else {
      if (alternative=="one.sided") {
          ldpow <- max(sum(zz$tstat <= ldq[1]), sum(zz$tstat >= ldq[2]))/nrep
      } else {
          ldpow <- (sum(zz$tstat <= ldq[1]) + sum(zz$tstat >= ldq[2]))/nrep
      }
  }
  out <- list()
  out$group.size <- gsize
  out$odds.ratio <- odds.ratio
  out$statistic <- statistic
  out$alternative <- alternative
  out$sig.level <- sig.level
  out$power <- ldpow
  class(out) <- "ladesign"
  out
}

print.ladesign <- function(x, ...) {
  cat("         Number of groups =", length(x$group.size), "\n")
  cat("               group size =", x$group.size, "\n")
  cat("odds ratios w.r.t group 1 =", round(x$odds.ratio, 3), "\n")
  cat("           test statistic =", x$statistic, "\n")
  cat("              alternative =", x$alternative, "\n")
  cat("       significance level =", x$sig.level, "\n")
  cat("                    power =", round(x$power, 3), "\n")
  invisible(x)
}

# k groups of sizes n_1,...,n_k; Total is N
# number of combinations is factorial(N)/product(factorial(n_i))
# if it is small we can enumerate the whole set
exactNull <- function(gsize, kw) {
    ng <- length(gsize)
    # use function combn to recursively generate combinations
    allcombs <- list()
    N <- Ni <- sum(gsize)
    ncomb <- 1
    for (i in 1:(ng-1)) {
        allcombs[[i]] <- combn(Ni, gsize[i])
        ncomb <- ncomb*ncol(allcombs[[i]])
        Ni <- Ni - gsize[i]
    }
    # vector of test statistics
    tstat <- rep(0, ncomb)
    # ways of assigning sum(gsize) elements into ng groups of size gsize
    igrps <- rep(0, N)
    # get the rank sum for each assignment
    # for KW the observations are ranked globally
    # for JT group j is ranked among observations from groups j thru ng
    if (kw) {rsums <- rep(0, ng)} else {rsums <- rep(0, ng-1)}
    # loop through the columns of each element in allcombs in a nested way
    ids <- rep(1, ng-1)
    nmax <- unlist(lapply(allcombs, ncol))
    irow <- 0
    while (all(ids <= nmax)) {
        # increment irow
        irow <- irow + 1
        # the column index for each group is specified by ids[i]
        igrps <- rep(0, N)
        for (i in 1:(ng-1)) {
            if (kw) {
                ii <- which(igrps==0)[allcombs[[i]][, ids[i]]]
                igrps[ii] <- i
                rsums[i] <- sum(ii)
            } else {
                rsums[i] <- sum(allcombs[[i]][, ids[i]])
            }
        }
        # igrps[which(igrps==0)] <- ng
        # rank sum of the last group for KW
        if (kw) {rsums[ng] <- sum(which(igrps==0))}
        # compute KW or JT statistic for each of the row of grps
        if (kw) {
            tstat[irow] <- sum(rsums^2/gsize)
        } else {
            tstat[irow] <- sum(rsums)
        }
        # increment ids to the next possible combination
        # find the last index of ids below nmax
        jj <- which(ids < nmax)
        if (length(jj) > 0) {
            jj0 <- max(jj)
            # increment the value at jj0 and set the ones after it to 1
            ids[jj0] <- ids[jj0] + 1
            if (jj0 < ng-1) ids[(jj0+1):(ng-1)] <- 1
        } else {
            ids[1] <- nmax[1] + 1
        }
    }
    # convert the JT statistic from (Wilcoxon) rank sum to Mann-Whitney form
    if (!kw) tstat <- tstat - sum(gsize[-ng]*(gsize[-ng]+1)/2)
    # return the statistics
    list(tstat=tstat)
}
