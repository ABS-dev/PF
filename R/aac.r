## used by RRmpWald
#' @importFrom stats model.frame terms
.twoby <- function(formula, data, compare, affected) {
  cluster <- function(x) {
    return(x)
  }
  drop.levels <- function(x) {
    for (j in seq_len(ncol(x))) {
      if (is.factor(x[, j])) {
        x[, j] <- factor(as.character(x[, j]))
      }
    }
    return(x)
  }
  data <- drop.levels(data)
  Terms <- terms(formula, specials = "cluster", data = data)
  environment(Terms) <- environment()
  A <- model.frame(formula = Terms, data = data)
  dat <- A[, 1] # all.vars(Terms)[1], response variable y
  group <- A[, 2] # attr(Terms, "term.labels")[1], treatment variable x
  clusters <- A[, 3] # attr(Terms, "term.labels")[2], grouping variable w
  tbl <- table(clusters, dat, group)[, 2, 1:2]
  for (i in 1:2) {
    tbl[, i] <- ifelse(tbl[, i] == affected, "af", "un")
  }

  # compare[2] is control
  # compare[1] is vaccinate
  xtable <- table(tbl[, compare[1]], tbl[, compare[2]])
  # order table
  xtable <- xtable[c("af", "un"), c("af", "un")]
  names(dimnames(xtable)) <- compare
  dimnames(xtable) <- lapply(dimnames(xtable), function(x) {
    ifelse(x == "af", "pos", "neg")
  })
  freqvec <- c(xtable)

  names(freqvec) <- paste(rep(dimnames(xtable)[[1]], 2),
                          rep(dimnames(xtable)[[2]], c(2, 2)))
  multvec <- as.data.frame(xtable)
  return(list(xtable = matrix(xtable, 2, 2),
              freqvec = freqvec, multvec = multvec))
}
