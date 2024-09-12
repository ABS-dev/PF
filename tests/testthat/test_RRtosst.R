context("RRtosst")

test_that("examples work", {
  y1 <- c(4, 24, 12, 28)
  ex1 <- RRtosst(y1)
  
  expect_s4_class(ex1, "rr1")
  expect_equal(ex1$estimate %>% signif(ex1$rnd) %>% unname, 
    c(0.611, 0.012, 0.902))
  expect_identical(ex1$estimator, "PF")
  expect_equal(ex1$y %>% as.numeric, y1)
  expect_equal(ex1$alpha, 0.05)
  
  y2 <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
  ex2 <- RRtosst(y2)
  expect_equal(ex1, ex2)
  
  data1 <- data.frame(group = rep(c("treated", "control"), each = 2),
    y = c(1, 3, 7, 5),
    n = c(12, 12, 14, 14), 
    cage = rep(paste('cage', 1:2), 2))
  ex3 <- RRtosst(data = data1, formula = cbind(y, n) ~ group, 
    compare = c("treated", "control"))
  expect_equal(ex1, ex3)
  
  data2 <- data1 %>%
    group_by(group) %>%
    summarize(sum_y = sum(y),
      sum_n = sum(n))
  
  ex4 <- RRtosst(data = data2, formula =  cbind(sum_y, sum_n) ~ group, 
    compare = c("treated", "control"))
  expect_equal(ex1, ex4)
})