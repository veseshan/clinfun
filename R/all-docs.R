

# Helper Functions --------------------------------------------------------

# Since roxygen2 doesn't parse the `@usage` block, it doesn't know what the usage is, and hence it doesn't find any parameters to inherit, so we can't rely on
# `@inheritParams`, cf. https://github.com/r-lib/roxygen2/issues/836#issuecomment-513476356
#
# Instead we use `@eval` to minimize duplication, cf. https://roxygen2.r-lib.org/articles/rd.html#evaluating-arbitrary-code
shared_ph2_args <- function() {
  c('@param pu unacceptable response rate; baseline response rate that needs to be exceeded for treatment to be deemed promising',
    '@param pa response rate that is desirable; should be larger than pu',
    '@param ep1 threshold for the probability of declaring drug promising under pu (target type 1 error rate); between 0 and 1',
    '@param ep2 threshold for the probability of declaring the drug not promising under pa (target type 2 error rate); between 0 and 1')
  }


# ph2simon ---------------------------------------------------------------

#' @name ph2simon
#' @title Simon's two-stage Phase II design
#' @usage
#' ph2simon(pu, pa, ep1, ep2, nmax = 100)
#'
#' # S3 method for ph2simon
#' print(x, ...)
#'
#' # S3 method for ph2simon
#' plot(x, ...)
#'
#' @description
#' Calculates the sample size and decision rules for Optimal and Minimax
#' two-stage Phase II designs given by Richard Simon, as well as any admissible design options (Jung et al, 2004) under
#' the given set of parameters.
#'
#' The trial proceeds to the second stage only if a minimal number of
#' responses (> `r1`) is observed at the end of the first stage.
#' Trial is stopped early if <= `r1` responses are seen in the first stage.
#'
#' @aliases ph2simon print.ph2simon plot.ph2simon
#'
#' @eval shared_ph2_args()
#' @param nmax maximum total sample size (default is 100; can be at most 1000)
#' @param x object returned by ph2simon
#' @param ... arguments to be passed onto plot and print commands called
#' within
#'
#' @return
#' ph2simon returns a list with `pu`, `pa`, `alpha`, `beta` and `nmax` (as defined above) and:
#' \item{out}{matrix of best two-stage designs for each value of total sample size n.
#'  The 6 columns in the matrix are:
#'     \tabular{rl}{
#'       r1 \tab Number of responses to be exceeded in stage 1 in order to continue to stage 2 \cr
#'       n1 \tab Number of subjects treated in first stage \cr
#'       r \tab Number of total (stage 1 + stage 2) responses to be exceeded at the end of stage 2 for treatment to be deemed promising \cr
#'       n \tab Total number (stage 1 + stage 2) of subjects to be treated in the trial \cr
#'       EN(p0) \tab Expected number of trial patients accrued under pu, calculated as the weighted
#'        average of number of patients, weighted by the probability of early termination under pu \cr
#'       PET(p0) \tab Probability of stopping after the first stage under pu
#'        (i.e. probability of stopping at the end of stage 1 if the true response rate is pu)
#'     }}
#'
#' @details
#' Trial is stopped early if <= `r1` responses are seen in the first stage
#' and treatment is considered desirable only when > `r` responses seen at the end of the second stage.
#'
#' Function will return design characteristics for the Optimal design (minimizes expected sample size under the null hypothesis),
#' the minimax design (minimizes maximum sample size under the null), and any admissible designs
#' (minimizes weighted average of expected sample size and maximum sample size).
#'
#' The `print` method formats and returns characteristics for admissible designs (including minimax and optimal designs)
#' The `plot` method plots the expected sample size against the
#' maximum sample size as in Jung et al., 2001

#' @section Methods (by generic):
#' - `print(ph2simon)`: Formats and returns the minimax, optimal, and any admissible designs.
#' - `plot(ph2simon)`: Plots the expected sample size against the maximum sample size, as in Jung et al., 2001.
#'
#' @family phase2
#' @seealso
#' `twostage.inference`, `oc.twostage.bdry`
#'
#' @references
#' - Simon R. (1989). Optimal Two-Stage Designs for Phase II Clinical Trials. *Controlled Clinical Trials*, 10, 1-10.
#' - Jung SH, Carey M, and Kim KM. (2001). Graphical Search for Two-Stage Designs for Phase II Clinical Trials. *Controlled Clinical Trials*, 22, 367-372.
#' - Jung SH, Lee T, Kim K, and George, SL. (2004). Admissible Two-Stage Designs for Phase II Cancer Clinical Trials. *Statistics in medicine*, 23(4), 561-569.
#'
#' @keywords design
#' @examples
#'
#' ph2simon(0.2, 0.4, 0.1, 0.1)
#' ph2simon(0.2, 0.35, 0.05, 0.05)
#' result <- ph2simon(0.2, 0.35, 0.05, 0.05, nmax=150)
#' plot(result)
#'
NULL



# ph2single ---------------------------------------------------------------

#' @name ph2single
#' @title Exact single stage Phase II design
#'
#' @description
#' Calculates the sample size and decision rule for exact single stage
#' Phase II design
#'
#' @usage ph2single(pu,pa,ep1,ep2,nsoln=5)
#' @eval shared_ph2_args()
#' @param nsoln number of designs to be returned for given `ep1` and `ep2`
#' @return  ph2single returns a data frame with the following variables:
#'    - `n` - the total number of subjects treated (n)
#'    - `r` - the number of responses (r) needed to be exceeded to consider
#' treatment promising
#'    - `Type I error`
#'    - `Type II error`
#'
#' @details
#' The function assesses multiple designs and returns the `nsoln` possible
#' ones with the type-I and type-II errors closest to and below the target ones
#' given in arguments (`ep1`, `ep2`). The treatment is deemed promising if > `r`
#' responses are observed in the `n` patients.
#'
#' @keywords design
#' @family phase2
#' @examples
#' ph2single(0.2, 0.4, 0.1, 0.1)
#' ph2single(0.2, 0.35, 0.05, 0.05)
#' ph2single(0.2, 0.35, 0.05, 0.05, nsoln = 10)
#'
NULL



# calogrank ---------------------------------------------------------------

##' @name calogrank
##' @title Survival curves analysis of covariance
##'
##' @description
##' Log-rank test to compare survival curves adjusting for covariates
##'
##' @usage calogrank(ftime, fstatus, grp, cvt, strat=NULL)
##'
##' @details
##' calogrank is the covariate adjusted version of k-sample survdiff. The
##' function in its current form only does basic error checking.
##'
##' @param ftime failure times
##' @param fstatus status indicator
##' @param grp group indicator
##' @param cvt continuous covariates used for adjusted analysis
##' @param strat stratification variable
##' @references Heller G. and Venkatraman E.S. (2004) A nonparametric test to
##' compare survival distributions with covariate adjustment. \emph{JRSS-B} 66,
##' 719-733.
##' @keywords htest
##' @examples
##'
##' \dontrun{
##' library(survival)
##' data(pbc)
##' pbc1 <- pbc
##' pbc1$trt[pbc1$trt == -9] <- NA
##' pbc1$copper[pbc1$copper == -9] <- NA
##' # only death (2) is considered; transplant(1) is censored
##' calogrank(pbc1$time, pbc1$status==2, pbc1$trt, pbc1[,c("copper")])
##' calogrank(pbc1$time, pbc1$status==2, pbc1$trt,
##'                                   pbc1[,c("protime", "copper")])}
##'
NULL
