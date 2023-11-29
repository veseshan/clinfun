
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
#' two-stage Phase II designs given by Richard Simon.
#' The trial proceeds to the second stage only if a minimal number of
#' responses (> `r`) is observed at the end of the first stage.
#' Trial is stopped early if <= `r1` responses are seen in the first stage.
#'
#' @aliases ph2simon print.ph2simon plot.ph2simon
#'
#' @param pu unacceptable response rate; baseline response rate that needs to
#' be exceeded for treatment to be deemed promising
#' @param pa response rate that is desirable; should be larger than pu
#' @param ep1 threshold for the probability of declaring drug promising
#; under pu (target type 1 error rate); between 0 and 1
#' @param ep2 threshold for the probability of declaring the drug not
#' promising under pa (target type 2 error rate); between 0 and 1
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
#'       EN(pu) \tab Expected number of patients in the trial under pu \cr
#'       PET(pu) \tab Probability of stopping after the first stage under pu
#'     }}
#'
#' @details
#'
#' Trial is stopped early if <= `r1` responses are seen in the first stage
#' and treatment is considered desirable only when > `r` responses seen at the end of the second stage.
#'
#' The `print` method formats and returns the minimax and optimal
#' designs.  The `plot` method plots the expected sample size against the
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
#' Simon R. (1989). Optimal Two-Stage Designs for Phase II Clinical Trials. *Controlled Clinical Trials*, 10, 1-10.
#'
#' Jung SH, Carey M, and Kim KM. (2001). Graphical Search for Two-Stage Designs for Phase II Clinical Trials. *Controlled Clinical Trials*, 22, 367-372.
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
#' Calculates the sample size and decision rule exact single stage
#' Phase II design
#'
#' @usage ph2single(pu,pa,ep1,ep2,nsoln=5)
#' @inheritParams ph2simon
#' @param nsoln number of designs to be returned for given `ep1` and `ep2`
#' @return  ph2single returns a data frame with the following variables:
#'    - `n` - the total number of subjects treated (n)
#'    - `r` - the number of responses (r) needed to be exceeded to consider
#' treatment promising
#'    - `Type I error`
#'    - `Type II error`
#' @keywords design
#' @family phase2
#' @examples
#' ph2single(0.2, 0.4, 0.1, 0.1)
#' ph2single(0.2, 0.35, 0.05, 0.05)
#' ph2single(0.2, 0.35, 0.05, 0.05, nsoln = 10)
#'
NULL
