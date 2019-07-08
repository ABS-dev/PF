## ---- echo=FALSE, message=FALSE------------------------------------------
require(magrittr)
library(xtable)


## ---- warning=FALSE------------------------------------------------------
library(PF)


## ------------------------------------------------------------------------
input.mtx <- matrix(c(0, 3, 1, 3, 1, 3, 3, 3, 0, 2, 2, 2, 1, 3, 2, 3, 2, 3,
  3, 3, 1, 2, 3, 3, 0, 2, 3, 3, 0, 3, 2, 2), 8, 4, byrow = TRUE)


## ----echo=FALSE, eval=TRUE, results='asis'-------------------------------
x <- input.mtx
colnames(x) <- c('Positive vaccinates', 'Total vaccinates', 
  'Positive controls', 'Total controls')
xtable(x, digits = 0, label = "tab:8strata", 
  caption = "An example with 8 strata") %>%
  print(., comment = FALSE)


## ------------------------------------------------------------------------
rr.by.stratum <- (input.mtx[, 1] / input.mtx[, 2]) /
  (input.mtx[, 3] / input.mtx[, 4])
round(1 - rr.by.stratum, 2)


## ------------------------------------------------------------------------
control.prev <- input.mtx[, 3] / input.mtx[, 4]
round(control.prev, 2)


## ------------------------------------------------------------------------
RRmh(Y = input.mtx)


## ------------------------------------------------------------------------
a <- RRstr(Y = input.mtx)
a


## ------------------------------------------------------------------------
a$y


## ------------------------------------------------------------------------
round(1 - a$y[, 5], 2)


## ------------------------------------------------------------------------
input.mtx <- matrix(c(0, 3, 3, 3, 0, 1, 1, 1, 1, 2, 1, 2, 2, 3, 2, 2, 1, 1, 
  1, 1, 1, 3, 1, 2, 0, 3, 2, 3, 0, 4, 4, 5), 8, 4, byrow = TRUE)


## ---- echo=FALSE, eval=TRUE, results='asis'------------------------------
x <- input.mtx
colnames(x) <- 
	c('Positive vaccinates','Total vaccinates','Positive controls','Total controls')
xtable(x, digits = 0, label = "tab:inhomog", 
  caption = "An example with inhomogeneity") %>% 
  print(., comment = FALSE)


## ------------------------------------------------------------------------
rr.by.stratum <- (input.mtx[, 1] / input.mtx[, 2]) / 
  (input.mtx[, 3] / input.mtx[, 4])
round(1 - rr.by.stratum, 2)


## ------------------------------------------------------------------------
RRtosst(apply(input.mtx, 2, sum))


## ------------------------------------------------------------------------
RRmh(Y = input.mtx)
RRstr(Y = input.mtx)


## ------------------------------------------------------------------------
input.mtx2 <- matrix(c(9, 10, 4, 5, 1, 10, 5, 5), 2, 4, byrow = TRUE)


## ---- echo=FALSE, eval=TRUE, results='asis'------------------------------
x <- input.mtx2
colnames(x) <- c('Positive vaccinates', 'Total vaccinates', 
  'Positive controls', 'Total controls')
xtable(x, digits = 0, caption = "An example with severe inhomogeneity",
  label = "tab:trainwreck") %>%
  print(., comment = FALSE)



## ------------------------------------------------------------------------
rr.by.stratum <- (input.mtx2[, 1] / input.mtx2[, 2]) / 
  (input.mtx2[, 3] / input.mtx2[, 4])
round(1 - rr.by.stratum, 2)


## ------------------------------------------------------------------------
RRstr(Y = input.mtx2)


## ------------------------------------------------------------------------
RRtosst(apply(input.mtx2, 2, sum))


## ------------------------------------------------------------------------
RRmh(Y = input.mtx2)


## ------------------------------------------------------------------------
input.mtx <- matrix(c(2, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 
  3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 1, 2, 3, 3), 9, 4, 
  byrow = TRUE)


## ---- echo=FALSE, eval=TRUE, results='asis'------------------------------
x <- input.mtx
colnames(x) <- 
	c('Positive vaccinates', 'Total vaccinates', 'Positive controls', 'Total controls')
xtable(x, digits = 0, caption = "An example with 9 strata", 
  label = "tab:9strata") %>% 
  print(., comment = FALSE)


## ------------------------------------------------------------------------
rr.by.stratum <- (input.mtx[, 1] / input.mtx[, 2]) / 
  (input.mtx[, 3] / input.mtx[, 4])
round(1 - rr.by.stratum, 2)
summary(1 - rr.by.stratum)


## ------------------------------------------------------------------------
RRtosst(apply(input.mtx, 2, sum))


## ------------------------------------------------------------------------
RRmh(Y = input.mtx)


## ------------------------------------------------------------------------
RRstr(Y = input.mtx)


## ------------------------------------------------------------------------
input.mtx <- matrix(c(1, 3, 4, 4, 0, 6, 6, 6), 2, 4, byrow = TRUE)


## ---- echo=FALSE, eval=TRUE, results='asis'------------------------------
x <- input.mtx
colnames(x) <- 
	c('Positive vaccinates', 'Total vaccinates', 'Positive controls', 'Total controls')
xtable(x, digits = 0, caption = "First example with 2 strata", 
  label = "tab:2strata1") %>%
  print(., comment = FALSE)


## ------------------------------------------------------------------------
RRmh(Y = input.mtx)
RRstr(Y = input.mtx)


## ------------------------------------------------------------------------
RRtosst(apply(input.mtx,2,sum))


## ------------------------------------------------------------------------
input.mtx2 <- matrix(c(1, 8, 1, 7, 1, 8, 4, 7), 2, 4, byrow = TRUE)


## ---- echo=FALSE, eval=TRUE, results='asis'------------------------------
x <- input.mtx2
colnames(x) <- 
	c('Positive vaccinates','Total vaccinates','Positive controls','Total controls')
xtable(x, digits = 0, label = "tab:2strata2", 
  caption = "Second example with 2 strata") %>%
  print(., comment = FALSE)


## ------------------------------------------------------------------------
rr.by.stratum <- (input.mtx2[, 1] / input.mtx2[, 2]) / 
  (input.mtx2[, 3] / input.mtx2[, 4])
1 - round(rr.by.stratum, 2)


## ------------------------------------------------------------------------
RRmh(Y = input.mtx2)
RRstr(Y = input.mtx2)


## ------------------------------------------------------------------------
RRtosst(apply(input.mtx2, 2, sum))


## ------------------------------------------------------------------------
input.mtx <- matrix(c(0, 5, 0, 4, 0, 5, 1, 3, 0, 2, 4, 6), 3, 4, 
  byrow = TRUE)


## ---- echo=FALSE, eval=TRUE, results='asis'------------------------------
x <- input.mtx
colnames(x) <- 
	c('Positive vaccinates', 'Total vaccinates', 'Positive controls', 'Total controls')
xtable(x, digits = 0, label = "tab:nonInf", 
  caption = "An example with a non-informative stratum") %>%
  print(., comment = FALSE)


## ------------------------------------------------------------------------
RRmh(Y = input.mtx)


## ---- error =TRUE--------------------------------------------------------
input.mtx2 <- input.mtx[-1,]
RRstr(Y = input.mtx2)



## ------------------------------------------------------------------------
input.mtx3 <- input.mtx
input.mtx3[1, 3] <- 1


## ---- error =TRUE--------------------------------------------------------
RRstr(Y = input.mtx3)


## ------------------------------------------------------------------------
input.mtx <- matrix(c(29, 45, 14, 45, 37, 45, 24, 45), 2, 4, byrow = TRUE)


## ---- echo=FALSE, eval=TRUE, results='asis'------------------------------
x <- input.mtx
colnames(x) <- 
	c('Positive vaccinates', 'Total vaccinates', 'Positive controls', 'Total controls')
xtable(x, digits = 0, caption = "A large sample example", label = "tab:large") %>%
  print(., comment = FALSE)


## ------------------------------------------------------------------------
rr.by.stratum <- (input.mtx[, 1] / input.mtx[, 2]) / 
  (input.mtx[, 3] / input.mtx[, 4])
round(rr.by.stratum, 2)


## ------------------------------------------------------------------------
RRmh(Y = input.mtx, pf = FALSE)


## ------------------------------------------------------------------------
RRstr(Y = input.mtx, pf = FALSE)


## ------------------------------------------------------------------------
system.time(b <- RRtosst(apply(input.mtx, 2, sum), pf = FALSE))
b

