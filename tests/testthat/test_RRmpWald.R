context("RRmpWald")

test_that("examples work", {
  # warning will go away once issue #5 is resolved
  expect_warning(ex1 <- RRmpWald(pos ~ tx + cluster(cage), New, compare = c('con', 'vac')))
  thismultvec <- data.frame(vac = rep(c('pos', 'neg'), 2),
    con = rep(c('pos', 'neg'), each = 2), Freq = c(7, 13, 2, 4))
  thismultvec$vac  <- factor(thismultvec$vac, levels = c('pos', 'neg'))
  thismultvec$con  <- factor(thismultvec$con, levels = c('pos', 'neg'))
  
  expect_identical(ex1$estimate %>% signif(ex1$rnd) %>% unname, 
    c(0.550, 0.183, 0.752))
  expect_identical(names(ex1$estimate), c("PF", "LL", "UL"))
  expect_identical(ex1$estimator, "PF")
  expect_identical(ex1$compare, c("con", "vac"))
  expect_identical(ex1$alpha, 0.05)
  expect_equal(ex1$xtable, matrix(c(7, 2, 13, 4), 2, 2, byrow = TRUE))
  expect_equal(ex1$freqvec %>% unname, c(7, 13, 2, 4))
  expect_identical(ex1$freqvec %>% names, 
    c('pos pos', 'neg pos', 'pos neg', 'neg neg'))
  expect_equal(ex1$multvec, thismultvec )
  
  thistable <- New %>%
    tidyr::spread(tx, pos) %>%
    dplyr::mutate(vac = factor(vac, levels = 1:0),
      con = factor(con, levels = 1:0)) %>%
    with(., table(vac, con)) 
  # warning will go away once issue #5 is resolved
  expect_warning(ex2 <- RRmpWald(x = as.vector(thistable))) 
  expect_equal(ex2$estimate, ex1$estimate)
  expect_equal(ex2$estimator, ex1$estimator)  
  expect_equal(ex2$compare, ex1$compare)
  expect_equal(ex2$alpha, ex1$alpha)
  expect_equal(ex2$xtable, ex1$xtable) 
  # expect_equal(ex2$freqvec, ex1$freqvec) # need to resolve issue 16 first
})
  