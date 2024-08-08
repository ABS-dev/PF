#' @title Wald confidence intervals for RR from matched pairs
#' @description Estimates confidence intervals for the risk ratio or prevented
#'   fraction from matched pairs.
#' @details Estimates confidence intervals for the risk ratio or prevented
#'   fraction from matched pairs. The response is the tetranomial vector `c(11,
#'   12, 21, 22)`, where the first index is the row and the the second index is
#'   the column when displayed as a 2x2 table. Wald type confidence intervals
#'   are found by applying the delta method to the multinomial variance. This
#'   method fails when there are no responders in one of the treatment groups.
#'
#'   Alternative forms of data entry are illustrated by the output, say `Y`,
#'   where `c(Y$xtable) = Y$freqvec = Y$multvec$Freq`.
#'
#'   If `RR = 0` (`PF = 1`), the function will return degenerate interval.
#' @name RRmpWald
#' @param formula Formula of the form `y ~ x + cluster(w)`, where `y` is the
#'   indicator for an individual's positive response, `x` is a factor with two
#'   levels of treatment, and `w` identifies the pairs.
#' @param data `data.frame` containing variables in formula
#' @param vac_grp The name of the vaccinated group.
#' @param con_grp The name of the control group.
#' @param affected Indicator for positive response
#' @param x Alternative data input. Instead of formula and data frame, data may
#'   be input as frequency vector. See example for how to order this vector.
#' @param alpha Complement of the confidence level
#' @param pf Estimate *RR* or its complement *PF*?
#' @param tdist Use t distribution?
#' @param df Degrees of freedom. When NULL, the function will default to `df = N
#'   - 2`, where N is the total number of pairs.
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#' @param compare `r badge("deprecated")`  Text vector stating the factor
#'   levels: `compare[1]` is the vaccinate group to which `compare[2]` (control
#'   or reference) is compared.
#' @returns A [rrmp] object with the following fields:
#' * `estimate`: vector of point and interval estimates - see details
#' * `estimator`: either `"PF"` or `"RR"`
#' * `compare`: text vector, same as input
#' * `alpha`: complement of confidence level
#' * `rnd`: how many digits to round the display
#' * `multvec`: data frame showing the multinomial representation of the data
#' @author [PF-package]
#' @note Experimental functions for estimating profile likelihood intervals are
#'   in the CVBmisc package.
#'
#'   Call to this function may be one of two formats: (1) specify `data` and
#'   `formula` or (2) as a vector `x`
#'
#'   `RRmpWald(formula, data, vac_grp = "vac", con_grp = "con", affected = 1,
#'   alpha = 0.05, pf = TRUE, tdist = TRUE, df = NULL, rnd = 3)`
#'
#'   `RRmpWald(x, vac_grp = "vac", con_grp = "con", affected = 1, alpha = 0, 05,
#'   pf = TRUE, tdist = TRUE, df = NULL, rnd = 3)`
#' @examples
#' RRmpWald(pos ~ tx + cluster(cage), New, vac_grp = "vac", con_grp = "con")
#'
#' thistable <- New |>
#'   tidyr::spread(tx, pos) |>
#'   tidyr::drop_na() |>
#'   dplyr::mutate(vac = factor(vac, levels = 1:0),
#'     con = factor(con, levels = 1:0)) |>
#'   with(table(vac, con))
#' thistable
#' as.vector(thistable)
#'
#' RRmpWald(x = as.vector(thistable))
#' @importFrom stats qnorm qt
#' @importFrom lifecycle badge deprecate_warn is_present deprecated
#' @export
RRmpWald <- function(formula = NULL,
                     data = NULL,
                     vac_grp = "vac",
                     con_grp = "con",
                     affected = 1,
                     x,
                     alpha = 0.05,
                     pf = TRUE,
                     tdist = TRUE,
                     df = NULL,
                     rnd = 3,
                     compare = deprecated()) {
  if (is_present(compare)) {
    deprecate_warn("9.7.0",
                   "RRmpWald(compare)",
                   "RRmpWald(vac_grp, con_grp)")
    if (length(compare) != 2) {
      stop("`compare` must be a vector of length 2!")
    }
    vac_grp <- compare[1]
    con_grp <- compare[2]
  }

  # CI for RR with matched pairs, based on asymptotic normality of log(RR) and
  # multinomial variance
  #
  # Data entry:
  #
  # formula of the form response ~ treatment + cluster(clustername) then it will
  # convert data to matrix and vector if entered as vector (x =) it must be
  # ordered by vac/con pairs: c(11, 01, 10, 00)
  multvec <- NULL
  if (!is.null(formula) && !is.null(data)) {
    x <- .twoby(formula = formula, data = data, vac_grp = vac_grp,
                con_grp = con_grp, affected = affected)
  } else if (is.matrix(x)) {
    stop("RRmpWald: data input by matrix is deprecated.")
  } else if (is.vector(x)) {
    if (length(x) != 4) {
      stop("Vector length must be 4\n")
    }
  }
  multvec <- matrix(x, 2, 2, byrow = FALSE)
  rownames(multvec) <- c("pos", "neg")
  colnames(multvec) <- c("pos", "neg")
  multvec <- multvec |> as.table() |> as.data.frame()
  colnames(multvec) <- c(vac_grp, con_grp, "Freq")

  N <- sum(x)
  p <- x / N
  V <- (diag(p) - t(t(p)) %*% t(p)) / N
  p1 <- p[1] + p[2]
  p2 <- p[1] + p[3]
  R <- p2 / p1
  gradR <- c((p[2] - p[3]) / p1^2, -p2 / p1^2, 1 / p1, 0)
  logR <- log(p2) - log(p1)
  gradlogR <- c(1 / p2 - 1 / p1, -1 / p1, 1 / p2, 0)
  varR <- as.numeric(t(gradR) %*% V %*% t(t(gradR)))
  varlogR <- as.numeric(t(gradlogR) %*% V %*% t(t(gradlogR)))

  if (tdist) {
    if (is.null(df))	df <- N - 2
  }
  if (!is.null(df)) {
    q <- qt(c(0.5, alpha / 2, 1 - alpha / 2), df)
  } else {
    q <- qnorm(c(0.5, alpha / 2, 1 - alpha / 2))
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
                  vac_grp = vac_grp, con_grp = con_grp,
                  alpha = alpha, rnd = rnd,
                  multvec =	multvec))
}

#' @importFrom stats model.frame terms
#' @importFrom data.table fifelse
#' @importFrom tidyr spread drop_na
.twoby <- function(formula, data, vac_grp, con_grp, affected) {
  tx <- pos <- NULL
  cluster <- function(x) {
    return(x)
  }
  data <- droplevels(data)
  Terms <- terms(formula, specials = "cluster", data = data)
  environment(Terms) <- environment()
  A <- model.frame(formula = Terms, data = data)
  names(A) <- c("pos", "tx", "cage")
  A <- tidyr::drop_na(A)
  A$pos <- fifelse(A$pos %in% affected, 1, 0)
  A$pos <- factor(A$pos, levels = 1:0)
  A$tx  <- fifelse(A$tx %in% vac_grp, "vac", "con")
  tbl <- A |>
    spread(tx, pos) |>
    drop_na() |>
    with(table(vac, con))
  as.vector(tbl)
}
