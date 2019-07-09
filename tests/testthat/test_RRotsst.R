context("RRotsst")

test_that("examples work", {
  y1 <- c(4, 24, 12, 28)
  ex1 <- RRotsst(y1, rnd = 3)
  
  expect_s4_class(ex1, "rr1")
  expect_identical(ex1$estimator, "PF")
  expect_equal(ex1$estimate %>% signif(ex1$rnd) %>% unname, 
    c(0.6110, 0.0148, 0.8520))
  expect_equal(ex1$y %>% as.vector, y1)
  expect_equal(ex1$alpha, 0.05)
  
  y2 <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
  ex2 <- RRotsst(y2, rnd = 3)
  expect_equal(ex1, ex2)
})