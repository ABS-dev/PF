#' @importFrom methods setRefClass new setMethod
setClassUnion("listORchar", c("list", "character"))

#' @title Data class pf
# @name pf-class
#' @description data class pf
#' @aliases pf
#' @rdname pf
#' @section Fields:
#' \itemize{
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' }
#' @seealso \code{\link{rr1}, \link{rrsi}, \link{rrsc}, \link{rrstr}}
#' @export
#' @author \link{PF-package}
pf <- setRefClass("pf",
                  fields = list(estimator = "character",
                                rnd = "numeric",
                                alpha = "numeric"))


#' @title Data class rr1
#' @description Data class rr1
#' @aliases rr1
#' @name rr1-class
#' @rdname rr1-class
#' @section Fields:
#' \itemize{
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#' \item[\code{y}]{data.frame with restructured input}
#' \item[\code{rnd}]{how many digits to round display}
#' \item[\code{alpha}]{complement of c.i.}
#' }
#' @seealso \code{\link{IDRsc}, \link{RRotsst}, \link{RRtosst}}
#' @exportClass rr1
#' @author \link{PF-package}
rr1 <- setRefClass("rr1",
                   contains = "pf",
                   fields = list(estimate = "numeric",
                                 y = "data.frame"))

#' @title Data class rror
#' @description data class rror
#' @aliases rror
#' @rdname rrorclass
#' @section Fields:
#'
#'   \itemize{
#'
#'   \item[\code{estimate}]{  vector with point and interval estimate}
#'
#'   \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#'
#'   \item[\code{y}]{  data vector}
#'
#'   \item[\code{rnd}]{  how many digits to round display}
#'
#'   \item[\code{alpha}]{  complement of c.i.}
#'
#'   \item[\code{norm}]{  logical indicating Gaussian or t interval}
#'
#'   \item[\code{degf}]{  degrees of freedom}
#'
#'   \item[\code{mu}]{  matrix with rows giving probability estimates for each
#'   of the groups}
#'
#'   }
#' @seealso \code{\link{RRor}}
#' @export
#' @author \link{PF-package}
rror <- setRefClass("rror",
                    contains = "rr1",
                    fields = list(norm = "logical",
                                  degf = "numeric",
                                  mu = "matrix"))

#' @title Data class rrsi
#' @description data class rrsi
#' @aliases rrsi
#' @rdname rrsi
#' @section Fields:
#'
#' \itemize{
#'
#' \item[\code{y}]{ data.frame with restructured input}
#'
#' \item[\code{k}]{  likelihood ratio criterion}
#'
#' \item[\code{rnd}]{  digits to round display}
#'
#' \item[\code{alpha}]{  complement of c.i.}
#'
#' \item[\code{estimate}]{  vector with point and interval estimate}
#'
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#'
#' }
#'
#' @seealso \code{\link{IDRlsi}, \link{RRlsi}}
#'
#' @export
#' @author \link{PF-package}
rrsi <- setRefClass("rrsi",
                    contains = "pf",
                    fields = list(y = "data.frame",
                                  k = "numeric",
                                  estimate = "numeric"))

#' @title Data class rrmp
#' @description data class rrmp
#' @aliases rrmp
#' @rdname rrmp
#' @section Fields:
#'
#'   \itemize{
#'
#'   \item[\code{estimate}]{  vector with point and interval estimate}
#'
#'   \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#'
#'   \item[\code{y}]{  data vector}
#'
#'   \item[\code{rnd}]{  how many digits to round display}
#'
#'   \item[\code{alpha}]{  complement of c.i.}
#'
#'   \item[\code{compare}]{  text vector, same as input}
#'
#'   \item[\code{multvec}]{  data.frame showing the multinomial representation
#'   of the data}
#'
#'   }
#' @seealso \code{\link{RRmpWald}}
#' @export
#' @author \link{PF-package}
rrmp <- setRefClass("rrmp",
                    contains = "rr1",
                    fields = list(compare = "character",
                                  multvec = "data.frame"))

#' @title Data class rrsc
#' @description data class rrsc
#' @rdname rrscclass
#' @section Fields:
#' \itemize{
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"RR"}}
#' \item[\code{y}]{ data.frame with restructured input}
#' }
#' @seealso \code{\link{RRsc}}
#' @export
#' @author \link{PF-package}
rrsc <- setRefClass("rrsc",
                    contains = "pf",
                    fields = list(estimate = "matrix",
                                  y = "data.frame"))

#' @title Data class rrstr
#' @description data class rrstr
#' @aliases rrstr
#' @rdname rrstrclass
#' @section Fields: \itemize{
#'
#'   \item[\code{estimate}]{  vector with point and interval estimate}
#'
#'   \item[\code{rnd}]{  how many digits to round display}
#'
#'   \item[\code{alpha}]{  complement of c.i.}
#'
#'   \item[\code{estimator}]{  either \code{"PF"} or \code{"RR"}}
#'
#'   \item[\code{hom}]{list of homogeneity statistic, p-value, and degrees of
#'   freedom. If \code{Phi == 0 | Phi == 1}, homogeneity test is not possible
#'   and error message displays}
#'
#'   \item[\code{y}]{data.frame of restructured input}
#'
#'   \item[\code{compare}]{  groups compared}
#'
#'   }
#' @seealso \code{\link{RRstr}}
#' @export
#' @author \link{PF-package}
rrstr <- setRefClass("rrstr",
                     contains = "pf",
                     fields = list(estimate = "matrix",
                                   hom = "listORchar",
                                   y = "data.frame",
                                   compare = "character"))
