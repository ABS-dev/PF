#' @title IDR likelihood support interval.
#' @description Estimates likelihood support interval for the incidence density ratio or prevented fraction based on it.
#' @details Estimates likelihood support interval for the incidence density ratio based on orthogonal factoring of reparameterized likelihood.
#' The incidence density is the number of cases per subject-time; its distribution is assumed Poisson.
#' \cr \cr Likelihood support intervals are usually formed based on the desired likelihood ratio, often 1/8 or 1/32. 
#' Under some conditions the log likelihood ratio
#' may follow the chi square distribution. If so, 
##  Latex version not good
## then \eqn{\alpha=1-{{F}_{{{\chi }^{2}}}}\left( 2\log (k),1 \right)}. 
#' then \eqn{\alpha=1-F(2log(k),1)}, where \eqn{F} is a chi-square CDF. 
#' \code{RRsc()} will make the conversion from \eqn{\alpha}
#' to \emph{k} if \code{use.alpha = TRUE}.
#' \cr \cr The data may also be a matrix. In that case \code{y} would be entered as \cr
#' \code{matrix(c(y1, n1 - y1, y2, n2 - y2), 2, 2, byrow = TRUE)}.
#' @param y Data vector c(y1, n1, y2, n2) where y are the positives, n are the total, 
#' and group 1 is compared to group 2 (control or reference).
#' @param compare  Text vector stating the factor levels: compare[1] is the 
#' vaccinate group to which compare[2] (control or reference) is compared.
#' @param data data.frame containing variables of the formula.
#' @param formula  Formula of the form cbind(y, n) ~ x, where y is the number 
#' positive, n is the group size, x is a factor with two levels of treatment.#' 
#' @param k Likelihood ratio criterion.
#' @param alpha Complement of the confidence level.
#' @param use.alpha Base choice of k on its relationship to alpha?
#' @param pf Estimate \emph{IDR} or its complement \emph{PF}? 
#' @param trace.it Verbose tracking of the iterations? 
#' @param iter.max Maximum number of iterations
#' @param converge Convergence criterion
#' @param rnd Number of digits for rounding. Affects display only, not estimates.
#' @param start describe here.
#' @return A \code{\link{rrsi}} object with the following elements.
#'  \item{estimate}{vector with point and interval estimate}
#'  \item{estimator}{either \emph{PF} or \emph{IDR}}
#'  \item{y}{data.frame with "y1", "n1", "y2", "n2" values.}
#'	\item{k}{Likelihood ratio criterion}
#'  \item{rnd}{how many digits to round the display}
#'  \item{alpha}{complement of confidence level}
#' @references Royall R. \emph{Statistical Evidence: A Likelihood Paradigm}. Chapman & Hall, Boca Raton, 1997. Section 7.2.
#' @author \link{PF-package}
#' @seealso \code{\link{IDRsc}}
#' @export
#' @examples
#' # Both examples represent the same observation, with data entry by vector
#' # and matrix notation.
#' 
#' y_vector <- c(26, 204, 10, 205)
#' IDRlsi(y_vector, pf = FALSE)
#' 
#' # 1/8 likelihood support interval for IDR 
#' 
#' # corresponds to 95.858% confidence
#' #   (under certain assumptions)
#' 
#' # IDR 
#' #  IDR   LL   UL 
#' # 2.61 1.26 5.88 
#' 
#' y_matrix <- matrix(c(26, 178, 10, 195), 2, 2, byrow = TRUE)
#' y_matrix
#' #      [,1] [,2]
#' # [1,]   26  178
#' # [2,]   10  195
#' 
#' IDRlsi(y_matrix, pf = FALSE)
#' 
#' # 1/8 likelihood support interval for IDR 
#' 
#' # corresponds to 95.858% confidence
#' #   (under certain assumptions)
#' 
#' # IDR 
#' #  IDR   LL   UL 
#' # 2.61 1.26 5.88 
#' 
#' data1 <- data.frame(group = rep(c("treated", "control"), each = 5),
#'              n = c(rep(41, 4), 40, rep(41, 5)),
#'              y = c(4, 5, 7, 6, 4, 1, 3, 3, 2, 1), 
#'              cage = rep(paste('cage', 1:5), 2))
#' IDRlsi(data = data1, formula = cbind(y, n) ~ group, 
#'                compare = c("treated", "control"), pf = FALSE)
#' 
#' # 1/8 likelihood support interval for IDR 
#' 
#' # corresponds to 95.858% confidence
#' #   (under certain assumptions)
#' 
#' # IDR 
#' #  IDR   LL   UL 
#' # 2.61 1.26 5.88 
#' 
#' require(dplyr)
#' data2 <- data1 %>%
#'   group_by(group) %>%
#'   summarize(sum_y = sum(y),
#'     sum_n = sum(n))
#'     
#' IDRlsi(data = data2, formula = cbind(sum_y, sum_n) ~ group, 
#'                compare = c("treated", "control"), pf = FALSE)
#' 
#' # 1/8 likelihood support interval for IDR 
#' 
#' # corresponds to 95.858% confidence
#' #   (under certain assumptions)
#' 
#' # IDR 
#' #  IDR   LL   UL 
#' # 2.61 1.26 5.88 
##---------------------------------------
## IDRlsi
##--------------------------------------- 

IDRlsi <-
  function(y = NULL,
    formula = NULL,
    data = NULL,
    alpha = 0.05,
    k = 8,
    use.alpha = FALSE,
    pf = TRUE,
    converge = 1e-8,
    rnd = 3,
    start = NULL,
    trace.it = FALSE,
    iter.max = 24,
    compare = c("con", "vac")) {
    
    ###########################################
    ## Error handling for input options
    ## - y can be matrix or vector (expects formula and data to be NULL)
    ## - if formula is specified, data is required (expects y is null)
    ###########################################
    .check_3input_cases_freq(data = data, formula = formula, y = y)

    
    ## END error handling
    ####
    
    ###########################################
    ## internal helper function
    ###########################################
    
    # support interval based on factoring orthogonal parameterization (Royall, p. 152)
    # Coded 1999
    L <- function(y, r) {
      y1 <- y[1]
      h1 <- y[2]
      y2 <- y[3]
      h2 <- y[4]
      h <- h1 / h2
      Lk <- (r * h) ^ y1 / (1 + r * h) ^ (y1 + y2)
      return(Lk)
    }
    
    ## END internal helpers
    ####

    
    ###########################################
    ## Data reshaping
    ## - y can be matrix or vector (expects formula and data to be NULL)
    ## - if formula is specified, data is required (expects y is null)
    ###########################################
    
    if (is.null(y)) {

      #extract from data+formula to vector c(y1, n1, y2, n2)
      y <- .extract_freqvec(formula, data, compare)
      
    } else if (is.matrix(y)) {
      y <- c(t(cbind(y[, 1], apply(y, 1, sum))))
    } 
 
    y1 <- y[1]
    h1 <- y[2]
    y2 <- y[3]
    h2 <- y[4]
    p1 <- y1 / h1
    p2 <- y2 / h2
    r.hat <- p1 / p2
    
    if (use.alpha) {
      k <- exp(qchisq(1 - alpha, 1) / 2)
    }	else {
      alpha <- 1 - pchisq(log(k) * 2, 1)
    }
    
    if (is.null(start)) {
      start <- r.hat * (cbind(c(.4, .5), c(1.5, 2)))
    }
    end <- L(y, r.hat) / k
    si <- rep(NA, 2)
    for (i in 1:2) {
      Q <- c(L(y, start[1, i]), L(y, start[2, i]))
      Q <- Q[rev(order(abs(Q - end)))]
      tt <- start[rev(order(abs(Q - end))), i]
      Q2 <- Q[2]
      Q1 <- Q[1]
      t2 <- tt[2]
      t1 <- tt[1]
      iter <- 0
      repeat {
        iter <- iter + 1
        if (trace.it) {
          cat('iter', iter, '\t')
        }
        t3 <- t2 + (t2 - t1) * (end - Q2) / (Q2 - Q1)
        if (trace.it) {
          cat('t321',
            t3,
            t2,
            t1,
            '321',
            L(y, t3),
            Q2,
            Q1,
            'converge',
            abs(t3 - t2) / t2,
            '\n')
        }
        if (iter > 1)
          if (abs((t3 - t2) / t2) < converge)
            break
        t1 <- t2
        t2 <- t3
        Q1 <- Q2
        Q2 <- L(y, t3)
        if (iter > iter.max)
          break
      }
      si[i] <- t3
    }
    int <- c(r.hat, si)
    if (!pf)
      names(int) <- c("IDR", "LL", "UL")
    else{
      int <- 1 - int[c(1, 3, 2)]
      names(int) <- c("PF.IDR", "LL", "UL")
    }
    
    y <- as.data.frame(t(y))
    names(y) <- c("y1", "n1", "y2", "n2")
    return(rrsi$new(
      estimate = int,
      estimator = ifelse(pf, 'PF_IDR', 'IDR'),
      y = y,
      rnd = rnd,
      k = k,
      alpha = alpha
    ))
    # out <- list(estimate = int, estimator = ifelse(pf, 'PF_IDR', 'IDR'), y = y, rnd = rnd, k=k, alpha = alpha)
    # class(out) <- 'rrsi'
    # return(out)
  }

