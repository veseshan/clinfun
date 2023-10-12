roc.curve <- function(marker, status, method=c("empirical")) {
  if (any(!is.finite(marker))) stop("Marker values should be finite")
  if (any(!is.finite(status))) stop("All status should be finite")
  method <- match.arg(method)
  ii <- order(marker)
  nu <- length(unique(marker)) + 1
  n <- length(marker)
  n1 <- sum(status)
  n0 <- n - n1
  if (min(n0,n1) == 0) stop("Status vector should have least one each of 0 & 1")
  zzz <- .Fortran("roccurve",
                  as.integer(n),
                  as.integer(n0),
                  as.integer(n1),
                  as.double(marker[ii]),
                  as.integer(status[ii]),
                  as.integer(nu),
                  tpr=double(nu),
                  fpr=double(nu))
  out <- NULL
  out$marker <- marker
  out$status <- status
  out$tpr <- zzz$tpr
  out$fpr <- zzz$fpr
  # add PPV for precision(ppv)-recall(tpr) curve (and NPV for completeness)
  prev = n1/n
  out$ppv = prev*out$tpr/(prev*out$tpr + (1-prev)*out$fpr)
  out$npv = (1-prev)*(1-out$fpr)/(prev*(1-out$tpr) + (1-prev)*(1-out$fpr))
  class(out) <- "roc.curve"
  out
}

print.roc.curve <- function(x, ...) {
  out <- roc.area.test(x$marker, x$status)
  cat("  ROC curve with AUC =", out$area, "and s.e. =", sqrt(out$var), "\n")
}

plot.roc.curve <- function(x, PRC=FALSE, ...) {
  if (PRC) {
    plot(x$tpr, x$ppv, xlab="Recall (TPR)", ylab="Precision (PPV)", ylim=c(0,1), type="l", ...)
  } else {
    plot(x$fpr, x$tpr, xlab="False positive rate", ylab="True positive rate", type="l", ...)
  }
}

lines.roc.curve <- function(x, PRC=FALSE, ...) {
  if (PRC) {
    lines(x$tpr, x$ppv, ...)
  } else {
    lines(x$fpr, x$tpr, ...)
  }
}
