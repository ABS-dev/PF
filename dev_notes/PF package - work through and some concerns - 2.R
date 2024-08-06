library(plyr)
library(data.table)
library(PF)

RRstr1 = function (formula = NULL, data = NULL, compare = c("vac", "con"), 
                   Y, alpha = 0.05, pf = TRUE, trace.it = FALSE, iter.max = 24, 
                   converge = 1e-06, rnd = 3, multiplier = 0.7, divider = 1.1) 
{
  u.p <- function(p1, p2, n1, n2) (1 - p1)/(n1 * p1) + (1 - p2) / (n2 * p2)
  
  zi.phi <- function(phi, Y, u.p, root, za) {
    zi <- ui <- 0
    for (i in 1:nrow(Y)) {
      y1 <- Y[i, 1]
      n1 <- Y[i, 2]
      y2 <- Y[i, 3]
      n2 <- Y[i, 4]
      p2 <- root(y1, y2, n1, n2, phi)
      p1 <- p2 * phi
      u <- u.p(p1, p2, n1, n2)
      u[u < 0] <- NA
      zz <- ((y1 - n1 * p1)/(1 - p1))
      z <- zz * sqrt(u)
      zd <- abs(z - za)
      zd[is.na(zd)] <- 2 * zd[!is.na(zd)]
      zi <- zi + unique(zz[zd == min(zd)])
      ui <- ui + 1/unique(u[zd == min(zd)])
    }
    return(zi/sqrt(ui))
  }
  
  zis.phi <- function(phi, Y, u.p, root, za) {
    zi <- ui <- gi <- 0
    for (i in 1:nrow(Y)) {
      y1 <- Y[i, 1]
      n1 <- Y[i, 2]
      y2 <- Y[i, 3]
      n2 <- Y[i, 4]
      p2 <- root(y1, y2, n1, n2, phi)
      p1 <- p2 * phi
      q1 <- 1 - p1
      q2 <- 1 - p2
      u <- u.p(p1, p2, n1, n2)
      u[u < 0] <- NA
      zz <- ((y1 - n1 * p1)/(1 - p1))
      z.ph <- zz * sqrt(u)
      g <- (q1 * (q1 - p1))/(n1 * p1)^2 - (q2 * (q2 - p2))/(n2 * p2)^2
      g.ph <- g/u^1.5
      z.s <- z.ph - (g.ph * (za^2 - 1))/6
      zd <- abs(z.s - za)
      zd[is.na(zd)] <- 2 * zd[!is.na(zd)]
      zi <- zi + unique(zz[zd == min(zd)])
      ui <- ui + 1/unique(u[zd == min(zd)])
      gi <- gi + unique(g[zd == min(zd)])/unique(u[zd == 
                                                     min(zd)])^3
    }
    return(zi/sqrt(ui) - (gi * (za^2 - 1))/(6 * ui^1.5))
  }
  
  root <- function(y1, y2, n1, n2, phi) {
    a <- phi * (n1 + n2)
    b <- -(phi * (y2 + n1) + y1 + n2)
    cc <- y1 + y2
    det <- sqrt(b^2 - 4 * a * cc)
    r1 <- (-b + det)/(2 * a)
    r2 <- (-b - det)/(2 * a)
    r1 <- round(r1, 8)
    r2 <- round(r2, 8)
    if (r1 < 0 || r1 > 1) 
      r1 <- r2
    if (r2 < 0 || r2 > 1) 
      r2 <- r1
    return(c(r1, r2))
  }
  
  rr.opt <- function(z.phi, phi, za, divider, trace.it, u.p, root) {
    zz <- c(z.phi(phi[1], Y, u.p, root, za), z.phi(phi[2], Y, u.p, root, za))
    if (abs(za - zz[1]) > abs(za - zz[2])) 
      phi <- rev(phi)
    phi.new <- phi[1]
    phi.old <- phi[2]
    z.old <- z.phi(phi.old, Y, u.p, root, za)
    if (trace.it) 
      cat("\n\nR start", phi, "\n")
    iter <- 0
    repeat {
      iter <- iter + 1
      z.new <- z.phi(phi.new, Y, u.p, root, za)
      if (is.na(z.new)) 
        f <- 1
      while (is.na(z.new)) {
        phi.new <- phi.old + (phi.new - phi.old)/f
        z.new <- z.phi(phi.new, Y, u.p, root, za)
        f <- f * (1/divider)
        if (trace.it) 
          cat("Step adjustment =", 1/f, "\n")
      }
      phi <- exp(log(phi.old) + log(phi.new/phi.old) * 
                   ((za - z.old)/(z.new - z.old)))
      phi.old <- phi.new
      z.old <- z.new
      phi.new <- phi
      if (trace.it) 
        cat("iteration", iter, "  z", z.new, "phi", phi.new, 
            "\n")
      if (abs(za - z.new) < converge) 
        break
      if (iter == iter.max) {
        cat("\nIteration limit reached without convergence\n")
        break
      }
    }
    return(phi.new)
  }
  
  this.call <- match.call() # This line seems to do nothing?
  
  if (!is.null(formula) & !is.null(data)) {
    cat("a\n")
    Y <- PF:::.matricize(formula = formula, data = data, compare = compare)$Y
  }
  rownames(Y) <- paste("Row", 1:nrow(Y), sep = "")
  colnames(Y) <- c("y1", "n1", "y2", "n2")
  Y <- cbind(Y, R.obs = (Y[, 1]/Y[, 2])/(Y[, 3]/Y[, 4]))
  zv <- qnorm(c(alpha/2, 1 - alpha/2))
  numer <- denom <- uu <- 0
  for (i in 1:nrow(Y)) {
    y1 <- Y[i, 1]
    n1 <- Y[i, 2]
    y2 <- Y[i, 3]
    n2 <- Y[i, 4]
    p1 <- y1/n1
    p2 <- y2/n2
    numer <- numer + (n2 * y1)/max((n1 + n2 - y1 - y2), 1)
    denom <- denom + (n1 * y2)/max((n1 + n2 - y1 - y2), 1)
    uu <- uu + ifelse(u.p(p1, p2, n1, n2) > 0, 1/u.p(p1, p2, n1, n2), 0)
  }
  phi <- Phi <- numer/denom
  v <- sqrt(uu)
  int <- sort(exp(log(phi) + (log(phi) * zv)/v))
  int[int < 0] <- 0
  which <- 1:2
  score <- rep(0, 2)
  cat("b\n")
  if (numer == 0) {
    cat("c\n")
    int[1] <- score[1] <- 0
    int[2] <- 1
    which <- 2
  }
  if (denom == 0) {
    cat("d\n")
    int[2] <- score[2] <- 1
    int[1] <- 0
    which <- 1
  }
  if (trace.it) {
    cat("e\n")
    cat("\nstarting estimates:")
    print(int)
  }
  if (Phi == 0 | Phi == 1) {
    cat("f\n")
    Phi.ML <- Phi
    hom <- paste("MLE = ", Phi.ML, ", Homogeneity test not possible", 
                 sep = "")
  } else {
    cat("g\n")
    za <- 0
    phi <- c(Phi, multiplier * Phi)
    if (trace.it) 
      cat("\nMLE")
    phi.new <- rr.opt(z.phi = zi.phi, phi = phi, za = za, 
                      divider = divider, trace.it = trace.it, u.p = u.p, 
                      root = root)
    Phi.ML <- phi.new
    test <- rep(0, nrow(Y))
    for (i in 1:nrow(Y)) {
      test[i] <- zi.phi(Phi.ML, t(as.matrix(Y[i, ])), u.p, 
                        root, za)^2
    }
    hom.test <- sum(test)
    df.test <- nrow(Y) - 1
    p.test <- 1 - pchisq(hom.test, df.test)
    hom <- list(stat = hom.test, df = df.test, p = p.test)
  }
  cat("h\n")
  for (k in which) {
    if (trace.it) 
      cat("\nScore", switch(k, "lower", "upper"))
    za <- -zv[k]
    phi <- c(int[k], multiplier * int[k])
    phi.new <- rr.opt(z.phi = zi.phi, phi = phi, za = za, 
                      divider = divider, trace.it = trace.it, u.p = u.p, 
                      root = root)
    score[k] <- phi.new
  }
  cat("i\n")
  int <- rbind(int, score)
  for (k in which) {
    if (trace.it) 
      cat("\nSkewness-corrected", switch(k, "lower", "upper"))
    za <- -zv[k]
    phi <- c(int[1, k], multiplier * int[1, k])
    phi.new <- rr.opt(z.phi = zis.phi, phi = phi, za = za, 
                      divider = divider, trace.it = trace.it, u.p = u.p, 
                      root = root)
    score[k] <- phi.new
  }
  cat("j\n")
  int <- rbind(int, score)
  int <- cbind(c(Phi, rep(Phi.ML, nrow(int) - 1)), int)
  if (trace.it) 
    cat("\n\n")
  cat("pf == ", pf,"\n")
  if (!pf) {
    dimnames(int) <- list(c("starting", "mle", "skew corr"), 
                          c("RR", "LL", "UL"))
  }
  else {
    int <- 1 - int[, c(1, 3, 2)]
    dimnames(int) <- list(c("starting", "mle", "skew corr"), 
                          c("PF", "LL", "UL"))
  }
  return(PF:::rrstr$new(estimate = int, 
                        hom = hom, 
                        estimator = ifelse(pf, "PF", "RR"), 
                        y = as.data.frame(Y), 
                        compare = compare, 
                        rnd = rnd, alpha = alpha))
}
cat("\014")

trial = data.table(animalid = 1:60,
                   group = rep(c("control", "vaccinate"), 30),
                   cage  = rep(c(1, 2, 3), each = 20),
                   sick = FALSE)

set.seed(6022)


trial[group == "control"   & cage == 1, sick := sample(c(TRUE, FALSE), 10, prob = c(.2, .8), replace = TRUE)]
trial[group == "control"   & cage == 2, sick := sample(c(TRUE, FALSE), 10, prob = c(.1, .9), replace = TRUE)]
trial[group == "control"   & cage == 3, sick := sample(c(TRUE, FALSE), 10, prob = c(.3, .7), replace = TRUE)]
trial[group == "vaccinate" & cage == 1, sick := sample(c(TRUE, FALSE), 10, prob = c(.8, .2), replace = TRUE)]
trial[group == "vaccinate" & cage == 2, sick := sample(c(TRUE, FALSE), 10, prob = c(.9, .1), replace = TRUE)]
trial[group == "vaccinate" & cage == 3, sick := sample(c(TRUE, FALSE), 10, prob = c(.7, .3), replace = TRUE)]


summary = trial[,  .N, by = .(group, cage, sick)]
wide = dcast(summary, group + cage ~ sick, value.var = "N", drop = FALSE)
setnames(wide, c("FALSE", "TRUE"), c("well", "sick"))
wide[is.na(well), well := 0]
wide[is.na(sick), sick := 0]
wide[, n := well + sick]
wide2 = wide
wide2[1, 1] = wide2[1, 1]
wide2[, group := ifelse(group == 'control', 'a', 'b')]

wide
wide2

RRstr1(cbind(sick, n) ~ group + cluster(cage), data.frame(wide),   compare = c("control", "vaccinate"))
RRstr1(cbind(sick, n) ~ group + cluster(cage), data.frame(wide2),  compare = c("a", "b"))

