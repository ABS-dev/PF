#' @title Wald confidence intervals for RR from matched pairs
#' @description Estimates confidence intervals for the risk ratio or prevented
#'   fraction from matched pairs.
#' @details Estimates confidence intervals for the risk ratio or prevented
#'   fraction from matched pairs. The response is the tetranomial vector `c(11,
#'   12, 21, 22)`, where the first index is the row and the the second index is
#'   the column when displayed as a 2x2 table. Wald type confidence intervals
#'   are found by applying the delta method to the multinomial variance. This
#'   method fails when there are no responders in one of the treatment groups.
#'   \cr \cr Alternative forms of data entry are illustrated by the output, say
#'   \cr \code{Y}, where \code{c(Y$xtable) = Y$freqvec = Y$multvec$Freq}. \cr
#'   \cr If RR = 0 (PF = 1), the function will return degenerate interval.
#' @name RRmpWald
#' @param formula Formula of the form \code{y ~ x + cluster(w)}, where y is the
#'   indicator for an individual's positive response, x is a factor with two
#'   levels of treatment, and w identifies the pairs.
#' @param data \code{data.frame} containing variables in formula
#' @param compare Text vector stating the factor levels: `compare[1]` is the
#'   vaccinate group to which `compare[2]` (control or reference) is compared.
#' @param affected Indicator for positive response
#' @param x Alternative data input. Instead of formula and data frame, data may
#'   be input as frequency vector. See example for how to order this vector.
#' @param alpha Complement of the confidence level
#' @param pf Estimate \emph{RR} or its complement \emph{PF}?
#' @param tdist Use t distribution?
#' @param df Degrees of freedom. When NULL, the function will default to
#'   \code{df = N - 2}, where N is the total number of pairs.
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#'
#' @return A \code{\link{rrmp}} object with the following fields:
#'
#'   \item{estimate}{vector of point and interval estimates - see details}
#'
#'   \item{estimator}{either \code{"PF"} or \code{"RR"}}
#'
#'   \item{compare}{text vector, same as input}
#'
#'   \item{alpha}{complement of confidence level}
#'
#'   \item{rnd}{how many digits to round the display}
#'
#'   \item{multvec}{data frame showing the multinomial representation of the
#'   data}
#'
#' @export
#' @author \link{PF-package}
#' @note Experimental functions for estimating profile likelihood intervals are
#'   in the CVBmisc package. \cr \cr Call to this function may be one of two
#'   formats: (1) specify \code{data} and \code{formula} or (2) as a vector
#'   \code{x} \cr \cr \code{RRmpWald(formula, data, compare = c('vac', 'con'),
#'   affected = 1, alpha = 0.05,} \cr \code{pf = TRUE, tdist = TRUE, df = NULL,
#'   rnd = 3)} \cr \cr \code{RRmpWald(x, compare = c('vac', 'con'), affected =
#'   1, alpha = 0, 05,} \cr \code{pf = TRUE, tdist = TRUE, df = NULL, rnd = 3)}
#' @examples
#' RRmpWald(pos ~ tx + cluster(cage), New, compare = c('vac', 'con'))
#'
#' # PF
#' # 95% interval estimates
#' #
#' #    PF    LL    UL
#' # 0.550 0.183 0.752
#'
#' require(magrittr)
#' thistable <- New %>%
#'   tidyr::spread(tx, pos) %>%
#'   dplyr::mutate(vac = factor(vac, levels = 1:0),
#'     con = factor(con, levels = 1:0)) %>%
#'   with(., table(vac, con))
#' thistable
#' #    con
#' # vac  1  0
#' #   1  7  2
#' #   0 13  4
#' as.vector(thistable)
#' # [1]  7 13  2  4
#'
#' RRmpWald(x = as.vector(thistable))
#'
#' # PF
#' # 95% interval estimates
#' #
#' #    PF    LL    UL
#' # 0.550 0.183 0.752
#' @importFrom dplyr "%>%"
#' @importFrom stats qnorm qt
RRmpWald <- function(formula = NULL, data = NULL, compare = c("vac", "con"),
                     affected = 1, x, alpha = 0.05, pf = TRUE, tdist = TRUE,
                     df = NULL, rnd = 3) {
  # CI for RR with matched pairs, based on asymptotic normality of log(RR) and
  # multinomial variance
  #
  # Data entry:
  #
  # formula of the form response ~ treatment + cluster(clustername) then it will
  # convert data to matrix and vector if entered as vector (x=) it must be
  # ordered by vac/con pairs: c(11, 01, 10, 00)
  multvec <- NULL
  if (!is.null(formula) && !is.null(data)) {
    Xx <- .twoby(formula = formula, data = data, compare = compare,
                 affected = affected)
    xtable <- Xx$xtable
    x <- Xx$freqvec
    multvec <- Xx$multvec
  } else if (is.matrix(x)) {
    # if (!all(dim(x) == 2)) {
    # 	stop("Table dimensions must be 2 x 2\n")
    # } else {
    # 	xtable <- x
    # 	x <- c(x)
    # }
    stop("RRmpWald: data input by matrix is deprecated.")
  } else if (is.vector(x)) {
    if (length(x) != 4) {
      stop("Vector length must be 4\n")
    }
    xtable <- matrix(x, 2, 2, byrow = FALSE)
    multvec <- xtable
    rownames(multvec) <- c("pos", "neg")
    colnames(multvec) <- c("pos", "neg")
    multvec <- multvec %>% as.table %>% as.data.frame
    colnames(multvec) <- c(compare, "Freq")

  }


  N <- sum(x)
  p <- x / N
  V <- (diag(p) - t(t(p)) %*% t(p)) / N
  p1 <- p[1] + p[2]
  p2 <- p[1] + p[3]
  R <- p2 / p1
  gradR <- c((p[2] - p[3]) / p1^2, -p2 / p1^2, 1 / p1, 0)
  logR <- log(p2) - log(p1)
  gradlogR <- c(1 / p2 - 1 / p1, -1 / p1, 1 / p2, 0)
  varR <- t(gradR) %*% V %*% t(t(gradR))
  varlogR <- t(gradlogR) %*% V %*% t(t(gradlogR))

  if (tdist) {
    if (is.null(df))	df <- N - 2
  }
  if (!is.null(df)) {
    q <- qt(c(0.5, alpha / 2, 1 - alpha / 2), df)
    what <- paste(100 * (1 - alpha), "% t intervals on ", df, " df\n", sep = "")
  } else {
    q <- qnorm(c(0.5, alpha / 2, 1 - alpha / 2))
    what <- paste(100 * (1 - alpha), "% gaussian interval\n", sep = "")
  }

  ci.dl <- exp(logR + q * sqrt(varlogR))

  # use the one based on log(R), unless R = 0
  if (R == 0) {
    ci.d <- R + q * sqrt(varR)
    ci.dl <- ci.d
  }

  int <- ci.dl
  if (!pf) {
    names(int) <- c("RR", "LL", "UL")
  } else {
    int <- 1 - int[c(1, 3, 2)]
    names(int) <- c("PF", "LL", "UL")
  }


  return(rrmp$new(estimate = int, estimator = ifelse(pf, "PF", "RR"),
                  compare = compare, alpha = alpha, rnd = rnd,
                  multvec =	multvec))

}
