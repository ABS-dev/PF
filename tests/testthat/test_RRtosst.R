context("RRtosst")

test_that("examples work", {
  y1 <- c(4, 24, 12, 28)
  ex1 <- RRtosst(y1)
  
  expect_s4_class(ex1, "rr1")
  expect_equal(ex1$estimate %>% signif(ex1$rnd) %>% unname, 
    c(0.611, 0.012, 0.902))
  expect_identical(ex1$estimator, "PF")
  expect_equal(ex1$y %>% as.vector, y1)
  expect_equal(ex1$alpha, 0.05)
  
  y2 <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
  ex2 <- RRtosst(y2)
  expect_equal(ex1, ex2)
})