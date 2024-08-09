#' @title IDR confidence interval.
#' @description Estimates confidence interval for the incidence density ratio or
#'   prevented fraction based on it.
#' @details The incidence density is the number of cases per subject-time; its
#'   distribution is assumed Poisson. `IDRsc` estimates a confidence interval
#'   for the incidence density ratio using Siev's formula based on the
#'   Poisson score statistic. \eqn{ IDR = \widehat{IDR}\left\{ 1 + \left(
#'   \frac{1}{{{y}_{1}}} + \frac{1}{{{y}_{2}}}
#'   \right)\frac{z_{\alpha / 2}^{2}}{2}\
#'   \ \pm \ \ \frac{z_{\alpha / 2}^{2}}{2{{y}_{1}}{{y}_{2}}}\sqrt{{{y}_{\bullet
#'   }}\left( {{y}_{\bullet }}z_{\alpha / 2}^{2} + 4{{y}_{1}}{{y}_{2}} \right)}
#'   \right\} }{(See
#'   \hyperref{https://www.aphis.usda.gov/animal_health/vet_biologics/
#'   publications/STATWI0007.pdf}{USDA CVB StatWI0007} for the formula.)}
#'
#'   The data may also be a matrix. In that case `y` would be entered as
#'
#'   `matrix(c(y1, n1 - y1, y2, n2 - y2), 2, 2, byrow = TRUE)`.
#' @param y Data vector c(y1, n1, y2, n2) where y are the positives, n are the
#'   total, and group 1 is compared to group 2 (control or reference).
#' @param formula Formula of the form cbind(y, n) ~ x, where y is the number
#'   positive, n is the group size, x is a factor with two levels of treatment.
#' @param data data.frame containing variables of formula.
#' @param vac_grp The name of the vaccinated group.
#' @param con_grp The name of the control group.
#' @param alpha Complement of the confidence level.
#' @param pf Estimate *IDR*, or its complement *PF*?
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#' @param compare `r badge("deprecated")`  Text vector stating the factor
#'   levels: `compare[1]` is the vaccinate group to which `compare[2]` (control
#'   or reference) is compared.
#' @returns A [rr1] object with the following elements.
#' * `estimate`: vector with point and interval estimate
#' * `estimator`: either *PF* or *IDR*
#' * `y`: data vector
#' * `rnd`: how many digits to round the display
#' * `alpha`: complement of confidence level
#' @references Siev D, 1994. Estimating vaccine efficacy in prospective studies.
#'   *Preventive Veterinary Medicine* 20:279-296, Appendix 1.
#'
#'   Graham PL, Mengersen K, Morton AP, 2003. Confidence limits for the ratio of
#'   two rates based on likelihood scores:non-iterative method *Statistics in
#'   Medicine* 22:2071-2083.
#'
#'   Siev D, 2004. Letter to the editor. *Statistics in Medicine* 23:693.
#'   (Typographical error in formula: replace the two final minus signs with
#'   subscript dots.)
#' @author [PF-package]
#' @seealso [IDRlsi]
#'
#' @examples
#' # All examples represent the same observation, with data entry by vector,
#' # matrix, and formula+data notation.
#'
#' y_vector <- c(26, 204, 10, 205)
#' IDRsc(y_vector, pf = FALSE)
#'
#' y_matrix <- matrix(c(26, 178, 10, 195), 2, 2, byrow = TRUE)
#' y_matrix
#'
#' IDRsc(y_matrix, pf = FALSE)
#'
#' require(dplyr)
#' data1 <- data.frame(group = rep(c("treated", "control"), each = 5),
#'             n = c(rep(41, 4), 40, rep(41, 5)),
#'             y = c(4, 5, 7, 6, 4, 1, 3, 3, 2, 1),
#'             cage = rep(paste("cage", 1:5), 2))
#' data2 <- data1 |>
#'   group_by(group) |>
#'   summarize(sum_y = sum(y),
#'   sum_n = sum(n))
#' IDRsc(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
#'     vac_grp = "treated", con_grp = "control", pf = FALSE)
#'
#' @importFrom stats qnorm
#' @importFrom lifecycle badge deprecate_warn is_present deprecated
#' @export
IDRsc <- function(y = NULL,
                  data = NULL,
                  formula = NULL,
                  vac_grp = "vac",
                  con_grp = "con",
                  alpha = 0.05,
                  pf = TRUE,
                  rnd = 3,
                  compare = deprecated()) {
  if (is_present(compare)) {
    deprecate_warn(
      "9.7.0",
      "IDRsc(compare)",
      details = "Please use the `vac_grp` and `con_grp` argumetns instead.")
    if (length(compare) != 2) {
      stop("`compare` must be a vector of length 2!")
    }
    vac_grp <- compare[1]
    con_grp <- compare[2]
  }

  ###########################################
  ## Error handling for input options
  ## - y can be matrix or vector (expects formula and data to be NULL)
  ## - if formula is specified, data is required (expects y is null)
  ###########################################
  .check_3input_cases_freq(data = data, formula = formula, y = y)

  ###########################################
  ## Data reshaping
  ## - y can be matrix or vector (expects formula and data to be NULL)
  ## - if formula is specified, data is required (expects y is null)
  ###########################################

  if (is.null(y)) {
    #extract from data+formula to vector c(y1, n1, y2, n2)
    y <- .extract_freqvec(formula, data, vac_grp, con_grp)

  } else if (is.matrix(y)) {
    y <- c(t(cbind(y[, 1], apply(y, 1, sum))))
  }
  y1 <- y[1.]
  s1 <- y[2.]
  y2 <- y[3.]
  s2 <- y[4.]
  y.dot <- y1 + y2
  z <- qnorm(1. - alpha / 2.)
  idr.hat <- (y1 * s2) / (y2 * s1)
  det <-
    (z * sqrt(y.dot * (y.dot * z^2. + 4. * y1 * y2))) / (2. * y1 * y2)
  det <- c(-det, det)
  ci <- idr.hat * (1. + (z^2. * (1. / y1 + 1. / y2)) / 2. + det)
  int <- c(idr.hat, ci)
  if (!pf) {
    names(int) <- c("IDR", "LL", "UL")
  } else {
    int <- 1 - int[c(1, 3, 2)]
    names(int) <- c("PF.IDR", "LL", "UL")
  }
  names(y) <- c("y1", "n1", "y2", "n2")

  return(rr1$new(
    estimate = int,
    estimator = ifelse(pf, "PF_IDR", "IDR"),
    y = as.data.frame(t(y)),
    rnd = rnd,
    alpha = alpha
  ))

}
