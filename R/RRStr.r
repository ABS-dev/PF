#' @title Gart-Nam method, CI for common RR over strata or clusters.
#' @description Estimates confidence intervals for the risk ratio or prevented
#'   fraction from clustered or stratified data.
#' @details Uses the DUD algorithm to estimate confidence intervals by the
#'   method of Gart.
#' @param formula Formula of the form `cbind(y, n) ~ x + cluster(w)`, where y is
#'   the number positive, n is the group size, x is a factor with two levels of
#'   treatment, and w is a factor indicating the clusters.
#' @param data data.frame containing variables of formula
#' @param compare Text vector stating the factor levels: `compare[1]` is the
#'   control or reference group to which `compare[2]` is compared
#' @param Y Matrix of data. Each row is a stratum or cluster. The columns are
#'   y2, n2, y1, n1. If data entered by formula and dataframe, Y is generated
#'   automatically.
#' @param pf Estimate *RR* or its complement *PF*?
#' @param alpha Size of the homogeneity test and complement of the confidence
#'   level.
#' @param trace.it verbose tracking of the iterations?
#' @param iter.max Maximum number of iterations
#' @param converge Convergence criterion
#' @param rnd Number of digits for rounding. Affects display only, not
#'   estimates.
#' @param multiplier internal control parameter for algorithm
#' @param divider internal control parameter for algorithm
#' @returns A [rrstr] object with the following fields:
#' * `estimate`: matrix of point and interval estimates - starting value, MLE,
#' and skewness corrected
#' * `hom`: list of homogeneity statistic, p-value, and degrees of freedom, or
#' error message if appropriate.
#' * `estimator`: either `"PF"` or `"RR"`
#' * `y`: `data.frame` of restructured input
#' * `compare`: groups compared
#' * `rnd`: how many digits to round the display
#' * `alpha`: size of test; complement of confidence level
#' @export
#' @note Vignette *Examples for Stratified Designs* forthcoming with more
#'   examples.
#' @references Gart JJ, 1985. Approximate tests and interval estimation of the
#'   common relative risk in the combination of \eqn{2 x 2} tables.
#'   *Biometrika* 72:673-677.
#'
#'   Gart JJ, Nam J, 1988. Approximate interval estimation of the ratio of
#'   binomial parameters: a review and corrections for skewness.
#'   *Biometrics* 44:323-338.
#'
#'   Ralston ML, Jennrich RI, 1978. DUD, A Derivative-Free Algorithm for
#'   Nonlinear Least Squares. *Technometrics* 20:7-14.
#' @author [PF-package]
#' @note Call to this function may be one of two formats: (1) specify `data` and
#'   `formula` or (2) as a matrix `Y`
#'
#'
#'   `RRstr(formula, data, compare = c("b", "a"), pf = TRUE, alpha = 0.05,
#'   trace.it = FALSE, iter.max = 24, converge = 1e-6, rnd = 3, multiplier =
#'   0.7, divider = 1.1)`
#'
#'   `RRstr(Y, compare = c("b", "a"), pf = TRUE, alpha = 0.05, trace.it = FALSE,
#'   iter.max = 24, converge = 1e-6, rnd = 3, multiplier = 0.7, divider = 1.1)`
#' @seealso [rrstr]
#' @examples
#' ## Table 1 from Gart (1985)
#' ##  as data frame
#' ## "b" is control group
#' RRstr(cbind(y, n) ~ tx + cluster(clus),
#'       Table6,
#'       compare = c("a", "b"), pf = FALSE)
#'
#' # Test of homogeneity across clusters
#'
#' # stat     0.954
#' # df       3
#' # p        0.812
#'
#' #  RR estimates
#'
#' # 		        RR   LL   UL
#' # starting   2.66 1.37 5.18
#' # mle        2.65 1.39 5.03
#' # skew corr  2.65 1.31 5.08
#'
#' ## or as matrix
#' RRstr(Y = table6, pf = FALSE)
#'
#' tst <- data.frame(y = c(0, 2, 0, 4, 0, 3, 0, 7),
#' n = rep(10, 8),
#' tx = rep(c("a", "b"), 4),
#' clus = rep(paste("Row", 1:4, sep = ""), each = 2))
#'
#' @importFrom stats pchisq qnorm
RRstr <- function(formula = NULL, data = NULL, compare = c("vac", "con"), Y,
                  alpha = 0.05, pf = TRUE, trace.it = FALSE, iter.max = 24,
                  converge = 1e-6, rnd = 3, multiplier = 0.7, divider = 1.1) {

  # define internal functions:
  #  u_p, zi_phi, zis_phi, root, rr.opt, matricize

  u_p <- function(p1, p2, n1, n2) {
    (1 - p1) / (n1 * p1) + (1 - p2) / (n2 * p2)
  }

  root <- function(y1, y2, n1, n2, phi) {
    # in RRstr
    a <- phi * (n1 + n2)
    b <-  -(phi * (y2 + n1) + y1 + n2)
    cc <- y1 + y2
    det <- sqrt(b^2 - 4 * a * cc)
    r1 <- (-b + det) / (2 * a)
    r2 <- (-b - det) / (2 * a)
    r1 <- round(r1, 8) ## round for comparison
    r2 <- round(r2, 8) ##
    if (r1 < 0 || r1 > 1) r1 <- r2
    if (r2 < 0 || r2 > 1) r2 <- r1
    return(c(r1, r2))
  }

  zi_phi <- function(phi, Y, u_p, root, za) {
    # for score interval in RRstr
    zi <- ui <- 0
    for (i in seq_len(nrow(Y))) {
      y1 <- Y[i, 1] # vac
      n1 <- Y[i, 2] # vac
      y2 <- Y[i, 3] # con or ref
      n2 <- Y[i, 4] # con or ref
      p2 <- root(y1, y2, n1, n2, phi)
      p1 <- p2 * phi
      u <- u_p(p1, p2, n1, n2)
      u[u < 0] <- NA
      zz <- ((y1 - n1 * p1) / (1 - p1))
      z <- zz * sqrt(u)
      zd <- abs(z - za)
      zd[is.na(zd)] <- 2 * zd[!is.na(zd)]
      zi <- zi + unique(zz[zd == min(zd)])
      ui <- ui + 1 / unique(u[zd == min(zd)])
    }
    return(zi / sqrt(ui))
  }

  zis_phi <- function(phi, Y, u_p, root, za) {
    # for skewness-corrected interval in RRstr
    zi <- ui <- gi <- 0
    for (i in seq_len(nrow(Y))) {
      y1 <- Y[i, 1]
      n1 <- Y[i, 2]
      y2 <- Y[i, 3]
      n2 <- Y[i, 4]
      p2 <- root(y1, y2, n1, n2, phi)
      p1 <- p2 * phi
      q1 <- 1 - p1
      q2 <- 1 - p2
      u <- u_p(p1, p2, n1, n2)
      u[u < 0] <- NA
      zz <- ((y1 - n1 * p1) / (1 - p1))
      z_ph <- zz * sqrt(u)
      g <- (q1 * (q1 - p1)) / (n1 * p1)^2 - (q2 * (q2 - p2)) / (n2 * p2)^2
      g.ph <- g / u^1.5
      z_s <- z_ph - (g.ph * (za^2 - 1)) / 6
      zd <- abs(z_s - za)
      zd[is.na(zd)] <- 2 * zd[!is.na(zd)]
      zi <- zi + unique(zz[zd == min(zd)])
      ui <- ui + 1 / unique(u[zd == min(zd)])
      gi <- gi + unique(g[zd == min(zd)]) / unique(u[zd == min(zd)])^3
    }
    return(zi / sqrt(ui) - (gi * (za^2 - 1)) / (6 * ui^1.5))
  }

  rr.opt <- function(z_phi, phi, za, divider, trace.it, u_p, root) {
    # optimizer function for RRstr
    zz <- c(z_phi(phi[1], Y, u_p, root, za), z_phi(phi[2], Y, u_p, root, za))
    if (abs(za - zz[1]) > abs(za - zz[2]))  phi <- rev(phi)
    phi_new <- phi[1]
    phi_old <- phi[2]
    z_old <- z_phi(phi_old, Y, u_p, root, za)
    if (trace.it) cat("\n\nR start", phi, "\n")
    iter <- 0
    repeat {
      iter <- iter + 1
      z_new <- z_phi(phi_new, Y, u_p, root, za)
      if (is.na(z_new)) f <- 1 ## adjust step
      while (is.na(z_new)) {
        phi_new <- phi_old + (phi_new - phi_old) / f
        z_new <- z_phi(phi_new, Y, u_p, root, za)
        f <- f * (1 / divider)
        if (trace.it) cat("Step adjustment =", 1 / f, "\n")
      }

      phi <- exp(log(phi_old) + log(phi_new / phi_old) *
                   ((za - z_old) / (z_new - z_old)))
      phi_old <- phi_new
      z_old <- z_new
      phi_new <- phi

      if (trace.it)
        cat("iteration", iter, "  z", z_new, "phi", phi_new, "\n")
      if (abs(za - z_new) < converge) break
      if (iter == iter.max) { # no convergence
        cat("\nIteration limit reached without convergence\n")
        break
      }
    } # end repeat
    return(phi_new)
  }

  # convert to matrix


  if (!is.null(formula) && !is.null(data)) {
    Y <- .matricize(formula = formula, data = data, compare = compare)$Y
  }
  rownames(Y) <- paste("Row", seq_len(nrow(Y)), sep = "")
  colnames(Y) <- c("y1", "n1", "y2", "n2")
  # save data and empirical Rs
  Y <- cbind(Y, R.obs = (Y[, 1] / Y[, 2]) / (Y[, 3] / Y[, 4]))

  zv <- qnorm(c(alpha / 2, 1 - alpha / 2))
  numer <- denom <- uu <- 0
  for (i in seq_len(nrow(Y))) {
    y1 <- Y[i, 1]
    n1 <- Y[i, 2]
    y2 <- Y[i, 3]
    n2 <- Y[i, 4]
    p1 <- y1 / n1
    p2 <- y2 / n2
    numer <- numer + (n2 * y1) / max((n1 + n2 - y1 - y2), 1) # added max
    denom <- denom + (n1 * y2) / max((n1 + n2 - y1 - y2), 1) #
    uu <- uu + ifelse(u_p(p1, p2, n1, n2) > 0, 1 / u_p(p1, p2, n1, n2), 0)
  } # end for i

  phi <- Phi <- numer / denom
  v <- sqrt(uu)
  int <- sort(exp(log(phi) + (log(phi) * zv) / v))
  int[int < 0] <- 0

  # if all zeros or ones estimate only one end of interval
  which <- 1:2
  score <- rep(0, 2)
  if (numer == 0) {
    int[1] <- score[1] <- 0
    int[2] <- 1
    which <- 2
  }
  if (denom == 0) {
    int[2] <- score[2] <- 1
    int[1] <- 0
    which <- 1
  }
  if (trace.it) {
    cat("\nstarting estimates:")
    print(int)
  }

  # get MLE of R
  # and test for heterogeneity
  if (Phi == 0 || Phi == 1) {
    Phi.ML <- Phi
    hom <- paste("MLE = ", Phi.ML, ", Homogeneity test not possible", sep = "")
  } else {
    za <-  0
    phi <- c(Phi, multiplier * Phi)
    if (trace.it) cat("\nMLE")
    phi_new <- rr.opt(z_phi = zi_phi, phi = phi, za = za, divider = divider,
                      trace.it = trace.it, u_p = u_p, root = root)
    Phi.ML <- phi_new
    # test for homogeneity of phi's
    test <- rep(0, nrow(Y))
    for (i in seq_len(nrow(Y))) {
      test[i] <- zi_phi(Phi.ML, t(as.matrix(Y[i, ])), u_p, root, za)^2
    }
    hom.test <- sum(test)
    df.test <- nrow(Y) - 1
    p.test <- 1 - pchisq(hom.test, df.test)
    hom <- list(stat = hom.test, df = df.test, p = p.test)
  }
  # Score intervals
  for (k in which) {
    if (trace.it) cat("\nScore", switch(k, "lower", "upper"))
    za <-  -zv[k]
    phi <- c(int[k], multiplier * int[k])
    phi_new <- rr.opt(z_phi = zi_phi, phi = phi, za = za, divider = divider,
                      trace.it = trace.it, u_p = u_p, root = root)
    score[k] <- phi_new
  }
  int <- rbind(int, score)

  # Skewness-corrected intervals
  for (k in which) {
    if (trace.it) cat("\nSkewness-corrected", switch(k, "lower", "upper"))
    za <-  -zv[k]
    phi <- c(int[1, k], multiplier * int[1, k])
    phi_new <- rr.opt(z_phi = zis_phi, phi = phi, za = za, divider = divider,
                      trace.it = trace.it, u_p = u_p, root = root)
    score[k] <- phi_new
  }   # end for k

  int <- rbind(int, score)
  int <- cbind(c(Phi, rep(Phi.ML, nrow(int) - 1)), int)
  if (trace.it) cat("\n\n")
  if (!pf) {
    dimnames(int) <- list(c("starting", "mle", "skew corr"),
                          c("RR", "LL", "UL"))
  } else {
    int <- 1 - int[, c(1, 3, 2)]
    dimnames(int) <- list(c("starting", "mle", "skew corr"),
                          c("PF", "LL", "UL"))
  }
  return(rrstr$new(estimate = int, hom = hom,
                   estimator = ifelse(pf, "PF", "RR"),
                   y = as.data.frame(Y), compare = compare,
                   rnd = rnd, alpha = alpha))
}
