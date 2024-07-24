#' @title Binomial dispersion parameter.
#' @description MME estimate of dispersion parameter phi.
#' @details Estimates binomial dispersion parameter \eqn{\phi} by the method of
#'   moments. Refits the model, weighting the observations by \eqn{1/\phi}. Uses
#'   `quasibinomial` family in `glm()`.
#' @param fit A [glm] object.
#' @param subset.factor Factor for estimating phi by subset.  Will be converted
#'   to a factor if it is not a factor.
#' @param fit.only Return only the new fit?  If FALSE, also returns the weights
#'   and phi estimates.
#' @param show.warns Show warnings
#' @returns A list with the following elements. `fit`: the new model fit,
#'   updated by the estimated weights `weights`: vector of weights `phi`: vector
#'   of phi estimates
#' @export
#' @references Wedderburn RWM, 1974. Quasi-likelihood functions, generalized
#'   linear models, and the Gauss-Newton method. *Biometrika* 61:439-447.
#' @author [PF-package]
#' @seealso [tauWt], [RRor].
#' @examples
#' birdm.fit <- glm(cbind(y, n - y) ~ tx-1, binomial, birdm)
#' RRor(phiWt(birdm.fit))
#' #
#' # 95% t intervals on 4 df
#' #
#' # PF
#' #     PF     LL     UL
#' #  0.479 -0.537  0.823
#' #
#' #       mu.hat   LL    UL
#' # txcon  0.768 0.95 0.367
#' # txvac  0.400 0.78 0.111
#' #
#' @importFrom stats glm update
phiWt <- function(fit,
                  subset.factor = NULL,
                  fit.only = TRUE,
                  show.warns = FALSE) {
  # Estimates weights = 1 / phi by MME where phi = dispersion parameter such
  # that var(y) = n * phi * mu * (1-mu) old family either binomial or poisson
  # newfamily is quasibinomial or quasipoisson
  subset.factor <- .check_factor(subset.factor)
  fit <- update(fit, x = TRUE, y = TRUE)
  y <- fit$y
  m <- fit$prior.weights
  oldfamily <- fit$family$family
  # works for binomial or poisson
  newfamily.name <- paste("quasi", oldfamily, sep = "")
  link <- fit$family$link
  newfamily <- get(newfamily.name)(link = link)
  if (is.null(subset.factor)) {
    subset.factor <- factor(rep("all", length(y)))
    w <- rep(1 / summary(update(fit, family = newfamily))$disp, length(y))
  } else {
    w <- rep(NA, length(y))
    for (lev in levels(subset.factor)) {
      xi <- rep(1, sum(subset.factor == lev))
      yi <- y[subset.factor == lev]
      mi <- m[subset.factor == lev]
      w[subset.factor == lev] <- 1 / summary(glm(yi ~ xi - 1,
                                                 family = newfamily,
                                                 weights = mi))$disp
    }
  }
  comment(w) <- paste(newfamily.name, "family,", link, "link, subsets:",
                      paste(levels(subset.factor), collapse = ", "))
  newfit <- update(fit, weights = w)
  phi <- 1 / tapply(w, subset.factor, unique)
  if (fit.only) out <- newfit
  else out <- list(fit = newfit, weights = w, phi = phi)
  return(out)
}
