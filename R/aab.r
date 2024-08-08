##### Shared functions for use within package

### this function is used in RRstr and RRmh

#' @title Convert data.frame and formula to a matrix
#' @description This function converts a subset of columns in `data` as
#'   specified by `formula` into a matrix
#' @param formula the formula to use. of format cbind(y, n) ~ tx + cluster(clus)
#' @param vac_grp The name of the vaccinated group.
#' @param con_grp The name of the control group.
#' @returns list of: `A`: A data.frame containing only the variables of
#'   `formula` `Y`: a matrix where each compare element is the set of columns
#'   `(y, n)` and each `unique(clus)` is a row
#' @seealso [RRmh], [RRstr]
#' @importFrom plyr ddply
#' @importFrom stats model.frame terms
#' @importFrom lifecycle badge deprecate_warn is_present deprecated
#' @noRd
.matricize <- function(formula, data, vac_grp, con_grp) {
  # 1/18/2012 - added error checking for formula argument. mcv
  # goal: ensure that formula is of meaningful format. if formula is incorrect,
  #      the accessors to A will not work. A better solution would be to access
  # 	   via names.
  if (length(all.vars(formula)) != 4) {
    stop("matricize: formula argument must be of",
         " format cbind(y, n) ~ tx + cluster(clus)")
  } else {
    if (is.na(pmatch("cbind", strsplit(as.character(formula), "~")[[2]]))) {
      stop("matricize: left side of formula argument must be of",
           " format cbind(y, n)")
    }

    if (length(strsplit(strsplit(as.character(formula), "~")[[3]],
                        "+", fixed = TRUE)[[1]]) != 2) {
      stop("matricize: right side of formula argument must be of",
           " format tx + cluster(clus)")
    }
  }

  cluster <- function(x) {
    return(x)
  }
  environment(cluster) <- parent.env(environment())

  Terms <- terms(formula, data = data)
  environment(Terms) <- environment()
  A <- model.frame(formula = Terms, data = data)

  A <- data.frame(A[, 1], A[, 2:3]) # for easier subscripting
  A <- A[order(A[, 4], A[, 3]), ]

  counts <- ddply(A, names(A)[4], nrow)
  rmclus <- counts[counts$V1 != 2, 1]
  if (length(rmclus) > 0) {
    message(paste(".matricize: Cluster group(s):",
                  paste(rmclus, collapse = ", ", sep = ""),
                  " does not have both comparison treatment levels.",
                  " Removing from analysis.",
                  collapse = "", sep = ""))
    A <- droplevels(A[!A[, 4] %in% rmclus, ])
  }
  x <- as.factor(A[, 3])

  clus <- A[, 4]
  Y1 <- A[x == con_grp, 1:2]
  Y2 <- A[x == vac_grp, 1:2]
  Y <- as.matrix(cbind(Y2, Y1))
  dimnames(Y) <- list(levels(clus), c("y1", "n1", "y2", "n2"))
  return(list(A = A, Y = Y))
}


.check_factor <- function(x) {
  parameter_name <- deparse(substitute(x))
  if (!is.null(x) && !is.factor(x)) {
    x <- factor(x)
    cat("Converting paramter ", parameter_name,
            " to factor")
  }
  x
}
