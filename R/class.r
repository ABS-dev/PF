#' @importFrom methods setRefClass new setMethod
setClassUnion("listORchar", c("list", "character"))

#' @title Data class pf
# @name pf-class
#' @description data class pf
#' @aliases pf
#' @rdname pf
#' @section Fields:
#' * `estimator`:  either `"PF"` or `"IDR"`
#' * `rnd`:  how many digits to round display
#' * `alpha`:  complement of c.i.
#' @seealso [rr1], [rrsi], [rrsc], [rrstr]
#' @export
#' @author [PF-package]
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
#' * `estimate`:   vector with point and interval estimate
#' * `estimator`:   either `"PF"` or `"IDR"`
#' * `Y`: data.frame with restructured input
#' * `rnd`: how many digits to round display
#' * `alpha`: complement of c.i.
#' @seealso [IDRsc], [RRotsst], [RRtosst]
#' @exportClass rr1
#' @author [PF-package]
rr1 <- setRefClass("rr1",
                   contains = "pf",
                   fields = list(estimate = "numeric",
                                 y = "data.frame"))

#' @title Data class rror
#' @description data class rror
#' @aliases rror
#' @rdname rrorclass
#' @section Fields:
#' * `estimate`:  vector with point and interval estimate
#' * `estimator`:  either `"PF"` or `"IDR"`
#' * `Y`:  data vector
#' * `rnd`:  how many digits to round display
#' * `alpha`:  complement of c.i.
#' * `norm`:  logical indicating Gaussian or t interval
#' * `degf`:  degrees of freedom
#' * `mu`:  matrix with rows giving probability estimates for each of the groups
#' @seealso [RRor]
#' @export
#' @author [PF-package]
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
#' * `Y`: data.frame with restructured input
#' * `k`:  likelihood ratio criterion
#' * `rnd`:  digits to round display
#' * `alpha`:  complement of c.i.
#' * `estimate`:  vector with point and interval estimate
#' * `estimator`:  either `"PF"` or `"IDR"`
#' @seealso [IDRlsi], [RRlsi]
#'
#' @export
#' @author [PF-package]
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
#' * `estimate`:  vector with point and interval estimate
#' * `estimator`:  either `"PF"` or `"IDR"`
#' * `Y`:  data vector
#' * `rnd`:  how many digits to round display
#' * `alpha`:  complement of c.i.
#' * `compare`:  text vector, same as input
#' * `multvec`:  data.frame showing the multinomial representation of the data
#' @seealso [RRmpWald]
#' @export
#' @author [PF-package]
rrmp <- setRefClass("rrmp",
                    contains = "rr1",
                    fields = list(compare = "character",
                                  multvec = "data.frame"))

#' @title Data class rrsc
#' @description data class rrsc
#' @rdname rrscclass
#' @section Fields:
#' * `estimate`:  vector with point and interval estimate
#' * `rnd`:  how many digits to round display
#' * `alpha`:  complement of c.i.
#' * `estimator`:  either `"PF"` or `"RR"`
#' * `Y`: data.frame with restructured input
#' @seealso [rrsc]
#' @export
#' @author [PF-package]
rrsc <- setRefClass("rrsc",
                    contains = "pf",
                    fields = list(estimate = "matrix",
                                  y = "data.frame"))

#' @title Data class rrstr
#' @description data class rrstr
#' @aliases rrstr
#' @rdname rrstrclass
#' @section Fields:
#' * `estimate`:   vector with point and interval estimate
#' * `rnd`:   how many digits to round display
#' * `alpha`:   complement of c.i.
#' * `estimator`:   either `"PF"` or `"RR"`
#' * `hom`: list of homogeneity statistic, p-value, and degrees of freedom.
#' If `Phi == 0 | Phi == 1`, homogeneity test is not possible and error message
#' displays
#' * `Y`: data.frame of restructured input
#' * `compare`:   groups compared
#' @seealso [rrstr]
#' @export
#' @author [PF-package]
rrstr <- setRefClass("rrstr",
                     contains = "pf",
                     fields = list(estimate = "matrix",
                                   hom = "listORchar",
                                   y = "data.frame",
                                   compare = "character"))
