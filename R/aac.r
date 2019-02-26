## used by RRmpWald
.twoby <- function(formula, data, compare, affected){
  cluster <- function(x) {return(x)}
  this.call <- match.call()
  drop.levels <- function(x){
    for (j in 1:ncol(x)) {
      if (is.factor(x[, j])) {
        x[, j] <- factor(as.character(x[, j]))
      }
    }
    return(x)
  }
  data <- drop.levels(data)
  Terms <- terms(formula, specials = 'cluster', data = data)
  environment(Terms) <- environment()
  A <- model.frame(formula = Terms, data = data)
  dat <- A[, 1]
  group <- A[, 2]
  clusters <- A[, 3]
  tbl <- table(clusters, dat, group)[, 2, 1:2]
  for (i in 1:2) {
    tbl[,i] <- ifelse(tbl[,i] == affected, 'af', 'un')
  }
  tbl <- data.frame(tbl)	
  # both levels must be present in both groups
  for (i in 1:2) {
    levels(tbl[,i]) <- c(levels(tbl[,i]), c('af', 'un')[!c('af', 'un') %in% 
        levels(tbl[,i])]) 
  }
  xtable <- table(tbl[,compare[2]], tbl[,compare[1]])
  # order table
  xtable <- xtable[c('af', 'un'), c('af', 'un')] 
  names(dimnames(xtable)) <- rev(compare)
  dimnames(xtable) <- lapply(dimnames(xtable), function(x){ifelse(x == 'af', 
    'pos', 'neg')})
  freqvec <- c(xtable)
  names(freqvec) <- paste(rep(dimnames(xtable)[[1]], 2), rep(dimnames(xtable)[[2]], 
    c(2, 2)))
  multvec <- as.data.frame(xtable)
  return(list(xtable = xtable, freqvec = freqvec, multvec = multvec))
}
