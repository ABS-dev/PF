#' @title RR score based asymptotic CI.
#' @name RRsc
#' @description Estimates confidence intervals for the risk ratio or prevented
#'   fraction based on the score statistic.
#' @details Estimates are returned for three estimators based on the score
#'   statistic. The score method was introduced by Koopman (1984). Gart and
#'   Nam's modification (1988) includes a skewness correction. The method of
#'   Miettinen and Nurminen (1985) is a version made slightly more conservative
#'   than Koopman's by including a factor of \eqn{(N - 1) / N}. The starting
#'   estimate for the DUD algorithm is obtained by the modified Katz method (log
#'   method with 0.5 added to each cell). Both forms of the Katz estimate may be
#'   retrieved from the returned object using `RRsc()$estimate`.
#'
#'   The data may also be a matrix. In that case `Y` would be entered as
#'
#'   `matrix(c(y1, n1-y1, y2, n2-y2), 2, 2, byrow = TRUE)`.
#' @param y Data vector c(y1, n1, y2, n2) where y are the positives, n are the
#'   total, and group 1 is compared to group 2 (control or reference group).
#' @param formula Formula of the form `cbind(y, n) ~ x`, where y is the number
#'   positive, n is the group size, x is a factor with two levels of treatment.
#' @param data data.frame containing variables of formula.
#' @param compare Text vector stating the factor levels: `compare[1]` is the
#'   vaccinate group to which `compare[2]` (control or reference) is compared.
#' @param alpha Complement of the confidence level.
#' @param pf Estimate *RR* or its complement *PF*?
#' @param trace.it Verbose tracking of the iterations?
#' @param iter.max Maximum number of iterations
#' @param converge Convergence criterion
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#' @returns A [rrsc] object with the following fields.
#' * `estimate`: matrix of point and interval estimates - see details
#' * `estimator`: either `"PF"` or `"RR"`
#' * `y`: data.frame with "y1", "n1", "y2", "n2" values.
#' * `rnd`: how many digits to round the display
#' * `alpha`: complement of confidence level
#' @export
#' @references Gart JJ, Nam J, 1988. Approximate interval estimation of the
#'   ratio of binomial parameters: a review and corrections for skewness.
#'   *Biometrics* 44:323-338.
#'
#'   Koopman PAR, 1984. Confidence intervals
#'   for the ratio of two binomial proportions. *Biometrics* 40:513-517.
#'
#'   Miettinen O, Nurminen M, 1985. Comparative analysis of two rates.
#'   *Statistics in Medicine* 4:213-226.
#'
#'   Ralston ML, Jennrich RI, 1978.
#'   DUD, A Derivative-Free Algorithm for Nonlinear Least Squares.
#'   *Technometrics* 20:7-14.
#' @author [PF-package]
#' @seealso [rrsc]
#'
#' @examples
#' # All examples represent the same observation, with data entry by using
#' # multiple notation options.
#'
#' y_vector <- c(4, 24, 12, 28)
#' RRsc(y_vector)
#'
#' # PF
#' # 95% interval estimates
#'
#' # 				PF     LL    UL
#' # MN method    0.611 0.0251 0.857
#' # score method 0.611 0.0328 0.855
#' # skew corr    0.611 0.0380 0.876
#'
#' y_matrix <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
#' #       [, 1] [, 2]
#' # [1, ]    4   20
#' # [2, ]   12   16
#'
#' RRsc(y_matrix)
#'
#' # PF
#' # 95% interval estimates
#'
#' # 				PF     LL    UL
#' # MN method    0.611 0.0251 0.857
#' # score method 0.611 0.0328 0.855
#' # skew corr    0.611 0.0380 0.876
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
#' RRsc(data = data2, formula = cbind(sum_y, sum_n) ~ group,
#'   compare = c("treated", "control"))
#'
#' # PF
#' # 95% interval estimates
#'
#' # 				PF     LL    UL
#' # MN method    0.611 0.0251 0.857
#' # score method 0.611 0.0328 0.855
#' # skew corr    0.611 0.0380 0.876
#' @importFrom stats qnorm
RRsc <- function(y = NULL,
                 data = NULL,
                 formula = NULL,
                 compare = c("vac", "con"),
                 alpha = 0.05,
                 pf = TRUE,
                 trace.it = FALSE,
                 iter.max = 18,
                 converge = 1e-6,
                 rnd = 3) {

  ###########################################
  ## Error handling for input options
  ## - y can be matrix or vector (expects formula and data to be NULL)
  ## - if formula is specified, data is required (expects y is null)
  ###########################################
  .check_3input_cases_freq(data = data, formula = formula, y = y)

  ## end error checking
  ###########################################

  ###########################################
  ##
  ## internal functions
  ###############################
  u_p <- function(p1, p2, n1, n2) {
    (1 - p1) / (n1 * p1) + (1 - p2) / (n2 * p2)
  }

  zsc.phi <- function(phi, x1, x2, n1, n2, u_p, root, za, MN = FALSE) {
    # for score interval in RRsc
    if (MN)
      mn <- sqrt((n1 + n2 - 1) / (n1 + n2))
    else
      mn <- 1
    p2 <- root(x1, x2, n1, n2, phi)
    p1 <- p2 * phi
    u <- u_p(p1, p2, n1, n2)
    u[u < 0] <- NA
    z <- ((x1 - n1 * p1) / (1. - p1)) * sqrt(u) * mn
    zd <- abs(z - za)
    zd[is.na(zd)] <- 2. * zd[!is.na(zd)]
    return(unique(z[zd == min(zd)]))
  }

  zsk.phi <- function(phi, x1, x2, n1, n2, u_p, root, za, MN = FALSE) {
    # for skewness-corrected interval in RRsc
    p2 <- root(x1, x2, n1, n2, phi)
    p1 <- p2 * phi
    q1 <- 1 - p1
    q2 <- 1 - p2
    u <- u_p(p1, p2, n1, n2)
    u[u < 0] <- NA
    z_ph <- ((x1 - n1 * p1) / (1 - p1)) * sqrt(u)
    g <- (q1 * (q1 - p1)) / (n1 * p1)^2 - (q2 * (q2 - p2)) / (n2 * p2) ^
      2
    g.ph <- g / u^1.5
    z_s <- z_ph - (g.ph * (za^2 - 1)) / 6
    zd <- abs(z_s - za)
    zd[is.na(zd)] <- 2 * zd[!is.na(zd)]
    return(unique(z_s[zd == min(zd)]))
  }

  root <- function(x1, x2, n1, n2, phi) {
    # in RRsc
    a <- phi * (n1 + n2)
    b <-  -(phi * (x2 + n1) + x1 + n2)
    cc <- x1 + x2
    det <- sqrt(b^2 - 4 * a * cc)
    r1 <- (-b + det) / (2 * a)
    r2 <- (-b - det) / (2 * a)
    if (r1 == 0)
      r1 <- 1e-006
    else if (r1 == 1)
      r1 <- 1 - 1e-006
    if (r2 == 0)
      r2 <- 1e-006
    else if (r2 == 1)
      r2 <- 1 - 1e-006
    if (r1 < 0 || r1 > 1)
      r1 <- r2
    if (r2 < 0 || r2 > 1)
      r2 <- r1
    return(c(r1, r2))
  }

  rr.opt <- function(z_phi, phi, za, trace.it, u_p, root, MN) {
    # optimizer function for RRsc
    # data from parent environment
    zz <- c(z_phi(phi[1], x1, x2, n1, n2, u_p, root, za, MN),
            z_phi(phi[2], x1, x2, n1, n2, u_p, root, za, MN))
    if (abs(za - zz[1]) > abs(za - zz[2]))
      phi <- rev(phi)
    phi_new <- phi[1]
    phi_old <- phi[2]
    z_old <- z_phi(phi_old, x1, x2, n1, n2, u_p, root, za, MN)
    if (trace.it)
      cat("\n\nR start", phi, "\n")
    iter <- 0
    repeat {
      iter <- iter + 1
      z_new <-
        z_phi(phi_new, x1, x2, n1, n2, u_p, root, za, MN)
      phi <-
        exp(log(phi_old) + log(phi_new / phi_old) *
              ((za - z_old) / (z_new - z_old)))
      phi_old <- phi_new
      z_old <- z_new
      phi_new <- phi
      if (trace.it)
        cat("iteration", iter, "  z", z_new, "phi", phi_new, "\n")
      if (abs(za - z_new) < converge)
        break
      if (iter == iter.max) {
        # no convergence
        cat("\nIteration limit reached without convergence\n")
        break
      }
    } # end repeat
    return(phi_new)
  }
  ###########################################
  ## Data reshaping
  ## - y can be matrix or vector (expects formula and data to be NULL)
  ## - if formula is specified, data is required (expects y is null)
  ###########################################

  if (is.null(y)) {
    # extract from data+formula to vector c(y1, n1, y2, n2)
    y <- .extract_freqvec(formula, data, compare)
  } else if (is.matrix(y)) {
    y <- c(t(cbind(y[, 1], apply(y, 1, sum))))
  }
  x1 <- y[1] ## vacc
  n1 <- y[2] ## vacc
  x2 <- y[3] ## control or ref
  n2 <- y[4] ## control or ref
  p1 <- x1 / n1
  p2 <- x2 / n2

  int <-
    matrix(NA, 6, 2, dimnames = list(
      c(
        "point",
        "log method",
        "0.5 method",
        "MN method",
        "score method",
        "skew corr"
      ),
      c("upper", "lower")
    ))
  al2 <- alpha / 2
  z_al2 <- qnorm(al2)
  z_ah2 <- qnorm(1 - al2)
  zv <- c(z_al2, z_ah2)
  int["point", ] <- rep((x1 / n1) / (x2 / n2), 2)
  p1 <- x1 / n1

  # log method
  p2 <- x2 / n2
  phi <- p1 / p2
  v <- sqrt(u_p(p1, p2, n1, n2))
  intv <- exp(v * zv + logb(phi))
  int["log method", ] <- intv

  # 0.5 log method
  p1 <- (x1 + 0.5) / (n1 + 0.5)
  p2 <- (x2 + 0.5) / (n2 + 0.5)
  phi <- p1 / p2
  v <- sqrt(u_p(p1, p2, (n1 + 0.5), (n2 + 0.5)))
  intv <- exp(v * zv + logb(phi))
  int["0.5 method", ] <- intv

  # which ends to estimate
  score.start <- rep(NA, 2)
  if (x1 > 0 && x2 > 0) {
    which <- 1:2
  } else if (x1 == 0) {
    score.start[1] <- 0
    which <- 2
  } else if (x2 == 0) {
    score.start[2] <- Inf
    which <- 1
  } else {
    message("Are you kidding?")
  }
  # MN method
  score <- score.start
  for (k in which) {
    if (trace.it)
      cat("\nMN", switch(k, "lower", "upper"))
    za <-  -zv[k]
    phi <- c(int["0.5 method", k], 0.9 * int["0.5 method", k])
    phi_new <-
      rr.opt(
        z_phi = zsc.phi,
        phi = phi,
        za = za,
        trace.it = trace.it,
        u_p = u_p,
        root = root,
        MN = TRUE
      )
    score[k] <- phi_new
  }
  int["MN method", ] <- score

  # Koopman score method
  score <- score.start
  for (k in which) {
    if (trace.it)
      cat("\nScore", switch(k, "lower", "upper"))
    za <-  -zv[k]
    phi <- c(int["0.5 method", k], 0.9 * int["0.5 method", k])
    phi_new <-
      rr.opt(
        z_phi = zsc.phi,
        phi = phi,
        za = za,
        trace.it = trace.it,
        u_p = u_p,
        root = root,
        MN = FALSE
      )
    score[k] <- phi_new
  }
  int["score method", ] <- score

  # skewness correction
  score <- score.start
  for (k in which) {
    if (trace.it)
      cat("\nSkew corr", switch(k, "lower", "upper"))
    za <-  -zv[k]
    phi <- c(int["0.5 method", k], 0.9 * int["0.5 method", k])
    phi_new <-
      rr.opt(
        z_phi = zsk.phi,
        phi = phi,
        za = za,
        trace.it = trace.it,
        u_p = u_p,
        root = root,
        MN = FALSE
      )
    score[k] <- phi_new
  }
  int["skew corr", ] <- score

  if (trace.it)
    cat("\n\n")

  int <- cbind(rep(int["point", 1], 5), int[-1, ])
  if (!pf) {
    dimnames(int)[[2]] <- c("RR", "LL", "UL")
  } else {
    int <- 1 - int[, c(1, 3, 2)]
    dimnames(int)[[2]] <- c("PF", "LL", "UL")
  }
  y <- data.frame(t(y))
  names(y) <- c("y1", "n1", "y2", "n2")
  return(rrsc$new(
    estimate = int,
    estimator = ifelse(pf, "PF", "RR"),
    y = y,
    rnd = rnd,
    alpha = alpha
  ))

}
