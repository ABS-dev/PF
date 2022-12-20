#' @title Mantel-Haenszel method, CI for common RR over strata or clusters with
#'   sparse data.
#' @description Estimates confidence intervals for the risk ratio or prevented
#'   fraction from clustered or stratified data, using a Mantel-Haenszel
#'   estimator for sparse data.
#' @details Based on the Mantel-Haenszel (1959) procedure for sparse data
#'   developed by Greenland and Robins (1985). The confidence limits are based
#'   on asymptotic normality of the log(risk ratio).  Agresti and Hartzel (2000)
#'   favor this procedure for small, sparse data sets, but they warn that it is
#'   less efficient than maximum likelihood for large data sets.
#' @param formula Formula of the form \code{cbind(y, n) ~ x + cluster(w)}, where
#'   \code{y} is the number positive, \code{n} is the group size, \code{x} is a
#'   factor with two levels of treatment, and \code{w} is a factor indicating
#'   the clusters.
#' @param data \code{data.frame} containing variables for formula
#' @param compare Text vector stating the factor levels: compare[1] is the
#'   vaccinate group to which compare[2] (control or reference) is compared.
#' @param Y Matrix of data, \eqn{K \times 4}{K x 4}. Each row is a stratum or
#'   cluster. The columns are \eqn{y1, n1, y2, n2}, where the y's are the number
#'   of positive in each group, and the n is the total in each group. Group 1
#'   corresponds to vaccinates and group 2 are controls or reference. If data
#'   entered by formula and dataframe, \code{Y} is generated automatically.
#' @param pf Estimate \emph{RR} or its complement \emph{PF}?
#' @param alpha Complement of the confidence level.
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#' @return An object of class \code{\link{rr1}}  with the following fields.
#'
#'   \item{estimate}{vector of point and interval estimates:  point estimate,
#'   lower confidence limit, upper confidence limit}
#'
#'   \item{estimator}{either \code{"PF"} or \code{"RR"}}
#'
#'   \item{y}{data.frame of restructured input}
#'
#'   \item{rnd}{how many digits to round the display}
#'
#'   \item{alpha}{complement of confidence level}
#'
#' @export
#' @note If either all y1's or all y2's are zero, a division by zero may occur,
#'   and a NaN returned for some values. \cr \cr Vignette \emph{Examples for
#'   Stratified Designs} forthcoming with more examples. \cr Call to this
#'   function may be one of two formats: (1) specify \code{data} and
#'   \code{formula} or (2) as a matrix \code{Y} \cr \cr \code{RRmh(formula,
#'   data, compare = c('b','a'), pf = TRUE, alpha = 0.05, rnd = 3)} \cr \cr
#'   \code{RRmh(Y, pf = TRUE, alpha = 0.05, rnd = 3)}
#' @references Mantel N, Haenszel W, 1959.  Statistical aspects of the analysis
#'   of data from retrospective studies of disease. \emph{Journal of the
#'   National Cancer Institute} 22:  719-748. \cr \cr Greenland S, Robins JM,
#'   1985.  Estimation of a common effect parameter from sparse follow-up data.
#'   \emph{Biometrics}  41:  55-68.  Errata, 45:  1323-1324. \cr \cr Agresti A,
#'   Hartzel J, 2000.  Strategies for comparing treatments on a binary response
#'   with multi-centre data.  \emph{Statistics in Medicine}  19:  1115-1139. \cr
#'   Lachin JM, 2000.  \emph{Biostatistical Methods:  The Assessment of Relative
#'   Risks} (Wiley, New York), Sec. 4.3.1.
#' @author \link{PF-package}
#' @seealso \code{\link{rr1}}
#' @examples
#'
#' ## Table 1 from Gart (1985)
#' ##  as data frame
#'
#' # tx group "b" is control
#' RRmh(cbind(y,n) ~ tx + cluster(clus), Table6, compare = c('a', 'b'), pf = FALSE)
#'
#' # RR
#' # 95% interval estimates
#' #
#' #   RR   LL   UL
#' # 2.67 1.37 5.23
#' #
#'
#' ## or as matrix
#' RRmh(Y = table6, pf = FALSE)
#'
#' # RR
#' # 95% interval estimates
#' #
#' #   RR   LL   UL
#' # 2.67 1.37 5.23

##########################################################################################
#
# Mantel-Haenszel estimate of common risk ratio for K 2x2 tables.
# First draft 5 March 2010.  Updated 11 March 2010, 5 Nov 2010.
# Revised 28 Dec 2010 to bring input/output format consistent with David Siev's RRstr().
# Note that unlike earlier versions, this version does not check the data types and dimensions
# of the inputs.
#
##########################################################################################

RRmh <- function(formula = NULL, data = NULL, compare = c('vac', 'con'), Y, alpha = 0.05,
                 pf = TRUE, rnd = 3)
{
  # convert data to matrix
  if (!is.null(formula) & !is.null(data)) {
    Y <- .matricize(formula = formula, data = data, compare = compare)$Y
  }
  colnames(Y) <- c('y1', 'n1', 'y2', 'n2')
  rownames(Y) <- paste("Row", 1:nrow(Y), sep = "")
  # save data and empirical Rs
  Y <- cbind(Y, R.obs = (Y[, 1]/Y[, 2])/(Y[, 3]/Y[, 4]))

  ### Innards of the algorithm
  zcrit <- qnorm(1 - alpha/2)
  x1vec <- Y[, 1] # vac
  n1vec <- Y[, 2] # vac
  x2vec <- Y[, 3] # con or ref
  n2vec <- Y[, 4] # con or ref
  nk <- n1vec + n2vec

  # Point estimate (Greenland & Robins, eq. 4.  Lachin, eq. 4.17)
  numer <- sum(x1vec * n2vec / (nk))
  denom <- sum(x2vec * n1vec / (nk))
  rr.est <- numer / denom

  # variance of log RR (Greenland & Robins, eq. 13)
  numer <- ((n1vec * n2vec) * (x1vec + x2vec) - (x1vec * x2vec) * nk) / nk^2
  rk <- x1vec * n2vec / nk
  sk <- x2vec * n1vec / nk
  denom <- sum(rk) * sum(sk)
  var.log.rr <- sum(numer) / denom

  # Confidence limits
  rr.ci.lo <- exp(log(rr.est) - zcrit * sqrt(var.log.rr))
  rr.ci.hi <- exp(log(rr.est) + zcrit * sqrt(var.log.rr))

  int <- c(rr.est, rr.ci.lo, rr.ci.hi)

  if (!pf) {
    names(int) <- c("RR", "LL", "UL")
  } else {
    int <- 1 - int[c(1, 3, 2)]
    names(int) <- c("PF", "LL", "UL")
  }

  return(rr1$new(estimate = int, estimator = ifelse(pf, 'PF', 'RR'),
                 y = as.data.frame(Y), rnd = rnd, alpha = alpha))
}
