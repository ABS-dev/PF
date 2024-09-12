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
#' @param formula Formula of the form `cbind(y, n) ~ x + cluster(w)`, where `Y`
#'   is the number positive, `n` is the group size, `x` is a factor with two
#'   levels of treatment, and `w` is a factor indicating the clusters.
#' @param data `data.frame` containing variables for formula
#' @param compare Text vector stating the factor levels: `compare[1]` is the
#'   vaccinate group to which `compare[2]` (control or reference) is compared.
#' @param Y Matrix of data, \eqn{K \times 4}{K x 4}. Each row is a stratum or
#'   cluster. The columns are \eqn{y1, n1, y2, n2}, where the y's are the number
#'   of positive in each group, and the n is the total in each group. Group 1
#'   corresponds to vaccinates and group 2 are controls or reference. If data
#'   entered by formula and dataframe, `Y` is generated automatically.
#' @param pf Estimate *RR* or its complement *PF*?
#' @param alpha Complement of the confidence level.
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#' @returns An object of class [rr1]  with the following fields.
#' * `estimate`: vector of point and interval estimates:  point estimate, lower
#' confidence limit, upper confidence limit
#' * `estimator`: either `"PF"` or `"RR"`
#' * `y`: data.frame of restructured input
#' * `rnd`: how many digits to round the display
#' * `alpha`: complement of confidence level
#' @export
#' @note If either all y1's or all y2's are zero, a division by zero may occur,
#'   and a NaN returned for some values.
#'
#'   Vignette *Examples for Stratified Designs* forthcoming with more examples.
#'
#'   Call to this function may be one of two formats: (1) specify `data` and
#'   `formula` or (2) as a matrix `Y`
#'
#'   `RRmh(formula, data, compare = c("b", "a"), pf = TRUE, alpha = 0.05, rnd =
#'   3)`
#'
#'   `RRmh(Y, pf = TRUE, alpha = 0.05, rnd = 3)`
#' @references Mantel N, Haenszel W, 1959.  Statistical aspects of the analysis
#'   of data from retrospective studies of disease. *Journal of the National
#'   Cancer Institute* 22:719-748.
#'
#'   Greenland S, Robins JM, 1985.  Estimation of a common effect parameter from
#'   sparse follow-up data. *Biometrics*  41:  55-68.  Errata, 45: 1323-1324.
#'
#'   Agresti A, Hartzel J, 2000.  Strategies for comparing treatments on a
#'   binary response with multi-centre data.  *Statistics in Medicine*  19:
#'   1115-1139.
#'
#'   Lachin JM, 2000.  *Biostatistical Methods:  The Assessment of Relative
#'   Risks* (Wiley, New York), Sec. 4.3.1.
#' @author [PF-package]
#' @seealso [rr1]
#' @examples
#'
#' ## Table 1 from Gart (1985)
#' ##  as data frame
#'
#' # tx group "b" is control
#' RRmh(cbind(y, n) ~ tx + cluster(clus),
#'      Table6,
#'      compare = c("a", "b"), pf = FALSE)
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

################################################################################
#
# Mantel-Haenszel estimate of common risk ratio for K 2x2 tables.
# First draft 5 March 2010.  Updated 11 March 2010, 5 Nov 2010.
# Revised 28 Dec 2010 to bring input/output format consistent with David Siev's
# RRstr(). Note that unlike earlier versions, this version does not check the
# data types and dimensions of the inputs.
#
################################################################################

#' @importFrom stats qnorm
RRmh <- function(formula = NULL,
                 data = NULL,
                 compare = c("vac", "con"),
                 Y,
                 alpha = 0.05,
                 pf = TRUE,
                 rnd = 3) {
  # convert data to matrix
  if (!is.null(formula) && !is.null(data)) {
    Y <- .matricize(formula = formula, data = data, compare = compare)$Y
  }
  colnames(Y) <- c("y1", "n1", "y2", "n2")
  rownames(Y) <- paste("Row", seq_len(nrow(Y)), sep = "")

  ### Innards of the algorithm
  zcrit <- qnorm(1 - alpha / 2)
  x1vec <- Y[, 1] # vac
  n1vec <- Y[, 2] # vac
  x2vec <- Y[, 3] # con or ref
  n2vec <- Y[, 4] # con or ref
  R.obs <- (x1vec / n1vec) / (x2vec / n2vec)
  nk <- n1vec + n2vec
  rk <- x1vec * n2vec / nk
  sk <- x2vec * n1vec / nk
  sum_rk <- sum(rk)
  sum_sk <- sum(sk)

  # save data and empirical Rs
  Y <- cbind(Y, R.obs)

  # Point estimate (Greenland & Robins, eq. 4.  Lachin, eq. 4.17)
  rr.est <- sum_rk / sum_sk

  # variance of log RR (Greenland & Robins, eq. 13)
  numer <- sum((n1vec * rk + n2vec * sk - x1vec * x2vec) / nk)
  denom <- sum_rk * sum_sk
  var.log.rr <- numer / denom

  # Confidence limits
  radius <- exp(zcrit * sqrt(var.log.rr))
  int <- c(rr.est, rr.est / radius, rr.est * radius)

  if (pf) {
    int <- 1 - int[c(1, 3, 2)]
    names(int) <- c("PF", "LL", "UL")
  } else {
    names(int) <- c("RR", "LL", "UL")
  }

  return(rr1$new(estimate = int, estimator = ifelse(pf, "PF", "RR"),
                 y = as.data.frame(Y), rnd = rnd, alpha = alpha))
}
