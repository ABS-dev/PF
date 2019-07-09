context("RRmh")

test_that("examples work", {
  ex1 <- RRmh(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = FALSE)
  thisy <- matrix(c(4, 2, 4, 1, 16, 16, 18, 15, 5, 3, 10, 3, 79, 87, 90, 82,
    3.95, 3.625, 2, 1.82222), nrow = 4, ncol = 5, byrow = FALSE )
  rownames(thisy) <- paste("Row", 1:4, sep = "")
  colnames(thisy) <- c("y1", "n1", "y2", "n2", "R.obs")
  
  expect_s4_class(ex1, "rr1")
  expect_identical(ex1$estimate %>% signif(ex1$rnd) %>% unname, 
    c(2.67, 1.37, 5.23))
  expect_identical(ex1$estimator, "RR")
  expect_equal(ex1$y, thisy, tolerance = 0.0002)
  
  ex2 <- RRmh(Y = table6, pf = FALSE)
  # expect_equal(ex2, ex1)  # resolve issue #15
})