#' @title Rao-Scott weighting.
#' @description Rao-Scott weighting of clustered binomial observations.
#' @details Estimates the cluster design effect \eqn{{d}_{i}}{d_i} as the
#'   variance inflation due to clustering by the method of Rao and Scott.
#'   Observations are then weighted by the inverse of the \eqn{{d}_{i}}{d_i}.
#' @param fit A [glm] object.
#' @param subset.factor Factor for estimating phi by subset.  Will be converted
#'   to a factor if it is not a factor.
#' @param fit.only Return only the new fit? If FALSE, also returns the weights
#'   and phi estimates.
#' @returns A list with the following elements.
#' * `fit`: the new model fit, updated by the estimated weights
#' * `weights`: vector of weights
#' * `d`: vector of \eqn{{d}_{i}}{d_i} estimates
#'
#' @export
#' @references Rao JNK, Scott AJ, 1992. A simple method for the analysis of
#'   clustered binary data. *Biometrics* 48:577-585.
#' @author [PF-package]
#' @seealso [RRor], [rsb].
#' @examples
#' birdm.fit <- glm(cbind(y, n - y) ~ tx-1, binomial, birdm)
#' RRor(rsbWt(birdm.fit))
#' #
#' # 95% t intervals on 4 df
#' #
#' # PF
#' #     PF     LL     UL
#' #  0.479 -1.061  0.868
#' #
#' #       mu.hat    LL     UL
#' # txcon  0.768 0.968 0.2659
#' # txvac  0.400 0.848 0.0737
#' #
#' @importFrom stats update
rsbWt <- function(fit = NULL,
                  subset.factor = NULL,
                  fit.only = TRUE) {
  ## rsbWt uses rsb to refit model
  subset.factor <- .check_factor(subset.factor)
  fit <- update(fit, x = TRUE, y = TRUE)
  yovern <- fit$y
  n <- fit$prior.weights
  if (is.null(n))
    n <- rep(1, length(fit$y))
  y <- yovern * n
  if (is.null(subset.factor))
    subset.factor <- factor(rep("all", length(y)))
  rsbdw <- rsb(y, n, id = subset.factor)
  w <- rsbdw$w
  d <- rsbdw$d
  newfit <- update(fit, weights = w)
  if (fit.only)
    out <- newfit
  else
    out <- list(fit = newfit, weights = w, d = d)
  return(out)
}

#' @title Rao-Scott weights.
#' @details Estimates the cluster design effect \eqn{{d}_{i}}{d_i} as the
#'   variance inflation due to clustering by the method of Rao and Scott. `rsb`
#'   estimates the \eqn{{d}_{i}}{d_i} for use by `rsbWt` or other functions.
#' @description Rao-Scott weights.
#' @param y vector of number positive.
#' @param n vector of total number.
#' @param id vector of factor for estimating the weights by subset.
#' @param data data.frame containing variables of formula.
#' @param formula Formula of the form cbind(y, n) ~ id, where y is the number
#'   positive, n is the total number, id is a factor for estimating the weights
#'   by subset.
#' @returns A list with the following elements.
#' * `w`: vector of weights
#' * `d`: vector of \eqn{{d}_{i}}{d_i} estimates
#' @export
#' @references Rao JNK, Scott AJ, 1992. A simple method for the analysis of
#'   clustered binary data. *Biometrics* 48:577-585.
#' @author [PF-package]
#' @seealso [rsbWt].
#' @examples
#' # Weil's rat data (Table 1 of Rao and Scott)
#' rsb(rat$y, rat$n, id = rat$group)$d
#' #  control  treated
#' # 1.232495 3.952861
#' rsb(data = rat, formula = cbind(y, n) ~ group)$d
#' #  control  treated
#' # 1.232495 3.952861
#' @importFrom dplyr select ungroup
rsb <- function(y = NULL, n = NULL, formula = NULL, data = NULL, id = NULL) {
  # rsb returns d's and weights
  if (is.null(y) && is.null(n) && !(is.null(formula) && is.null(data))) {
    ## case when formula and data is input
    vars <- all.vars(formula)
    sumdata <- data |>
      ungroup() |>
      select(vars[1:3])
    colnames(sumdata) <- c("y", "n", "id")
    y <- sumdata$y
    n <- sumdata$n
    id <- sumdata$id
  } else if (!(is.null(y) && is.null(n)) && is.null(formula) && is.null(data)) {
    ## case when y and n is input
    if (is.null(id)) {
      id <- rep("all", length(y))
    }
  } else {
    stop("Cannot parse input. If \"y\" & \"n\" are provided, ",
         "\"formula\" & \"data\" must be NULL and vice-versa.")
  }

  # Rao-Scott design effect weights
  # Biometrics 48:577-585
  # based on raoscott.bin function

  y.i <- tapply(y, id, sum)
  n.i <- tapply(n, id, sum)
  m.i <- c(table(id))
  p.i <- y.i / n.i
  r.ij <- y - n * p.i[tapply(y, id)]
  v.i <- (tapply(r.ij^2, id, sum) * m.i) / ((m.i - 1) * n.i^2)
  d.i <- (n.i * v.i) / (p.i * (1 - p.i))
  w <- 1 / d.i[tapply(y, id)]
  return(list(weights = w, d = d.i))
}
