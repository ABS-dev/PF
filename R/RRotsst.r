#' @title RR exact CI, OTSST method.
#' @description Estimates confidence interval for the risk ratio or prevented
#'   fraction, exact method based on the score statistic (inverts one two-sided
#'   test).
#' @details Estimates confidence intervals based on the score statistic that are
#'   'exact' in the sense of accounting for discreteness. The score statistic is
#'   used to select tail area tables, and the binomial probability is estimated
#'   over the tail area by taking the maximum over the nuisance parameter.
#'   Algorithm is a simple step search.
#'
#'   The data may also be a matrix. In that case `Y` would be entered as
#'
#'   `matrix(c(y1, n1 - y1, y2, n2 - y2), 2, 2, byrow = TRUE)`.
#' @param y Data vector c(y1, n1, y2, n2) where y are the positives, n are the
#'   total, and group 1 is compared to group 2 (control or reference).
#' @param compare Text vector stating the factor levels: `compare[1]` is the
#'   vaccinate group to which `compare[2]` (control or reference) is compared.
#' @param data data.frame containing variables of the formula.
#' @param formula  Formula of the form cbind(y, n) ~ x, where y is the number
#'   positive, n is the group size, x is a factor with two levels of treatment.
#' @param alpha Complement of the confidence level.
#' @param pf Estimate *RR* or its complement *PF*?
#' @param trace.it Verbose tracking of the iterations?
#' @param iter.max Maximum number of iterations
#' @param converge Convergence criterion
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#' @param stepstart starting interval for step search
#' @param nuisance.points number of points over which to evaluate nuisance
#'   parameter
#' @param gamma parameter for Berger-Boos correction (restricts range of
#'   nuisance parameter evaluation)
#' @returns An object of class [rr1] with the following fields:
#' * `estimate`: vector with point and interval estimate
#' * `estimator`: either `"PF"` or `"RR"`
#' * `y`: `data.frame` with "y1", "n1", "y2", "n2" values.
#' * `rnd`: how many digits to round the display
#' * `alpha`: complement of confidence level
#' @export
#' @references Koopman PAR, 1984. Confidence intervals for the ratio of two
#'   binomial proportions. *Biometrics* 40:513-517.
#'
#'   Agresti A, Min Y, 2001.  On small-sample confidence intervals for
#'   parameters in discrete distribution. *Biometrics* 57: 963-971.
#'
#'   Berger RL, Boos DD, 1994. P values maximized over a confidence set for the
#'   nuisance parameter. *Journal of the American Statistical Association*
#'   89:214-220.
#' @author [PF-package]
#' @seealso [RRtosst], [rr1].
#'
#' @examples
#' # All examples represent the same observation, with data entry by multiple
#' # options.
#'
#' y_vector <- c(4, 24, 12, 28)
#' RRotsst(y_vector, rnd = 3)
#'
#' # PF
#' # 95% interval estimates
#'
#' #    PF     LL     UL
#' # 0.6111 0.0148 0.8519
#'
#' y_matrix <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
#' RRotsst(y_matrix, rnd = 3)
#'
#' # PF
#' # 95% interval estimates
#'
#' #    PF     LL     UL
#' # 0.6111 0.0148 0.8519
#'
#' require(dplyr)
#' data1 <- data.frame(group = rep(c("treated", "control"), each = 2),
#'   y = c(1, 3, 7, 5),
#'   n = c(12, 12, 14, 14),
#'   cage = rep(paste("cage", 1:2), 2))
#'
#' data2 <- data1 |>
#'   group_by(group) |>
#'   summarize(sum_y = sum(y),
#'     sum_n = sum(n))
#' RRotsst(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
#'    compare = c("treated", "control"))
#'
#' # PF
#' # 95% interval estimates
#' #
#' # PF     LL     UL
#' # 0.6111 0.0148 0.8519
#' @importFrom stats qbeta
RRotsst <- function(y = NULL,
                    data = NULL,
                    formula = NULL,
                    compare = c("vac", "con"),
                    alpha = 0.05,
                    pf = TRUE,
                    stepstart = .1,
                    iter.max = 36,
                    converge = 1e-6,
                    rnd = 3,
                    trace.it = FALSE,
                    nuisance.points = 120,
                    gamma = 1e-6) {
  # Estimates exact confidence interval by the OTSST method
  # Score statistic used to select tail area tables
  # Binomial probability estimated over the tail area
  #  by taking the maximum over the nuisance parameter

  # Written 9/17/07 by Siev
  #
  # Functions called by rrcix(): .rr.score.asymp - gets asymptotic interval
  # for starting value of upper bound found in this file below
  #
  # (if want to eliminate calling this function would have to search down from
  # r.max)
  #
  # binci - gets Clopper-Pearson intervals for Berger-Boos method
  #
  # included here now, but may be moved to another package

  ###########################################
  ## internal helper function
  ###########################################
  binci <- function(y,
                    n,
                    alpha = .05,
                    show.warnings = FALSE) {
    w <- 1 * show.warnings - 1
    options(warn = w)

    p <- y / n
    cpl <- ifelse(y > 0, qbeta(alpha / 2, y, n - y + 1), 0)
    cpu <- ifelse(y < n, qbeta(1 - alpha / 2, y + 1, n - y), 1)
    out <- cbind(y, n, p, cpl, cpu)
    dimnames(out) <-
      list(names(y), c("y", "n", "p.hat", "cp low", "cp high"))

    options(warn = 0)
    return(out)
  }
  ## END internal helpers
  ####

  ###########################################
  ## Error handling for input options
  ## - y can be matrix or vector (expects formula and data to be NULL)
  ## - if formula is specified, data is required (expects y is null)
  ###########################################
  .check_3input_cases_freq(data = data, formula = formula, y = y)

  ## END error handling
  ####

  ###########################################
  ## Data reshaping
  ## - y can be matrix or vector (expects formula and data to be NULL)
  ## - if formula is specified, data is required (expects y is null)
  ###########################################
  if (is.null(y)) {

    #extract from data+formula to vector c(y1, n1, y2, n2)
    #y1n1 from vaccinate group
    #y2n2 from control group
    y <- .extract_freqvec(formula, data, compare)

  } else if (is.matrix(y)) {
    # Data entry y = c(x2, n2, x1, n1) Vaccinates First (order same but
    # subscripts reversed)
    # data vector
    y <- c(t(cbind(y[, 1], apply(y, 1, sum))))
    # NOTE: the subscripts are reversed compared to the other functions
  }

  x2 <- y[1] ## vacc
  n2 <- y[2] ## vacc
  x1 <- y[3] ## control
  n1 <- y[4] ## control
  p1 <- x1 / n1
  p2 <- x2 / n2
  rho.mle <- p2 / p1

  # itemize all possible tables in omega (17.26)
  Y <- data.frame(y1 = rep(0:n1, (n2 + 1)),
                  y2 = rep(0:n2, rep(n1 + 1, n2 + 1)))
  observed <- (seq_len(nrow(Y)))[Y[, 1] == x1 & Y[, 2] == x2]
  Y$C <- choose(n1, Y$y1) * choose(n2, Y$y2)

  # score statistic - with pi.tilde by quadratic formula
  scst <- function(rho, y1, n1, y2, n2) {
    pih1 <- y1 / n1  # unrestricted MLE of current data
    pih2 <- y2 / n2
    if (y1 == 0 && y2 == 0)
      sc <- 0
    else if (y2 == n2) {
      sc <- 0
    } else {
      A <- rho * (n1 + n2)
      B <- -(rho * (y1 + n2) + y2 + n1)
      C <- y1 + y2
      pit1 <- (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
      pit2 <- rho * pit1
      sc <- (pih2 - rho * pih1) /
        sqrt(rho^2 * pit1 * (1 - pit1) / n1 + pit2 * (1 - pit2) / n2)
    }
    return(sc)
  }


  # get Clopper-Pearson intervals for Berger-Boos method
  cp <- binci(c(x1, x2), c(n1, n2), alpha = gamma)[, c("cp low", "cp high")]
  L1 <- cp[1, 1]
  U1 <- cp[1, 2]
  L2 <- cp[2, 1]
  U2 <- cp[2, 2]
  r.min <- L2 / U1
  r.max <- U2 / L1

  if (rho.mle == 0) {
    low <- 0
  } else {
    # search for lower endpoint
    iter <- 0
    step <- stepstart
    low <- max(0.0001, r.min) # start above 0 (for quadratic formula)
    repeat {
      iter <- iter + 1
      if (iter > iter.max)
        break
      if (iter > 1) {
        low <- low + step
      }
      scst_y <- rep(NA, nrow(Y))
      for (i in seq_along(scst_y))
        scst_y[i] <- scst(low, Y$y1[i], n1, Y$y2[i], n2)
      q_set <- Y[abs(scst_y) >= abs(scst_y[observed]), ]
      q_set$n1y1 <- n1 - q_set$y1
      q_set$n2y2 <- n2 - q_set$y2
      if (gamma > 0)
        # Berger-Boos method 17.164
        pn <- seq(max(L1, L2 / low), min(U1, U2 / low),
                  length = nuisance.points)
      else
        # simple method 17.138
        pn <- seq(0, min(1 / low, 1), length = nuisance.points)
      if (sum(pn > 1) > 0) {
        cat("\nIteration",
            iter,
            "nuisance parameter outside parameter space\n")
        next
      }
      fy <- rep(NA, nuisance.points)
      for (i in 1:nuisance.points) {
        pni <- pn[i]
        fy[i] <-
          sum(
            q_set$C * pni^q_set$y1 *
              (1 - pni)^q_set$n1y1 *
              (low * pni)^q_set$y2 * (1 - low * pni)^q_set$n2y2
          )
      }
      max.fy <- max(fy)
      if (trace.it)
        cat("\nIteration", iter, "rho.low", low, "tail", max.fy, "\n")
      if (abs(max.fy - (alpha - gamma / 2)) < converge)
        break
      if (max.fy > (alpha - gamma / 2)) {
        step <- step / 2
        low <- low - step * 2
      }
    } # end repeat
  } # end else

  # search for upper endpoint upward from just below asymptotic
  # rather than downward from r.max
  # get asymptotic interval for starting

  # koopman version (slightly narrower interval than mn)
  ci.asymp <- .rr.score.asymp(c(x2, n2, x1, n1))
  high <- ci.asymp[3] * .9

  iter <- 0
  step <- stepstart
  repeat {
    iter <- iter + 1
    if (iter > iter.max)
      break
    if (iter > 1) {
      high <- high + step
    }
    scst_y <- rep(NA, nrow(Y))
    for (i in seq_along(scst_y))
      scst_y[i] <- scst(high, Y$y1[i], n1, Y$y2[i], n2)
    p_set <- Y[abs(scst_y) >= abs(scst_y[observed]), ]
    p_set$n1y1 <- n1 - p_set$y1
    p_set$n2y2 <- n2 - p_set$y2
    if (gamma > 0)
      # Berger-Boos method 17.164
      pn <- seq(max(L1, L2 / high), min(U1, U2 / high),
                length = nuisance.points)
    else
      # Berger-Boos method 17.164
      pn <- seq(0, min(1 / high, 1), length = nuisance.points)
    if (sum(pn > 1) > 0) {
      cat("\nIteration",
          iter,
          "nuisance parameter outside parameter space\n")
      next
    }
    fy <- rep(NA, nuisance.points)
    for (i in 1:nuisance.points) {
      pni <- pn[i]
      fy[i] <-
        sum(p_set$C * pni^p_set$y1 * (1 - pni)^p_set$n1y1 *
              (high * pni)^p_set$y2 *  (1 - high * pni)^p_set$n2y2)
    }
    max.fy <- max(fy)
    if (trace.it)
      cat("\nIteration", iter, "rho.high", high, "tail", max.fy, "\n")
    if (abs(max.fy - (alpha - gamma / 2)) < converge)
      break
    if (max.fy < (alpha - gamma / 2)) {
      step <- step / 2
      high <- high - step * 2
    }
  } # end repeat

  int <- c(rho.hat = rho.mle,
           low = low,
           high = high)
  if (!pf) {
    names(int) <- c("RR", "LL", "UL")
  } else {
    int <- 1 - int[c(1, 3, 2)]
    names(int) <- c("PF", "LL", "UL")
  }
  y <- as.data.frame(t(y))
  names(y) <- c("y1", "n1", "y2", "n2")
  return(rr1$new(
    estimate = int,
    estimator = ifelse(pf, "PF", "RR"),
    y = y,
    rnd = rnd,
    alpha = alpha
  ))

}
