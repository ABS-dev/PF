context("RRstr")

test_that("examples work", {
  ex1 <- RRstr(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = FALSE)
  thisestimate <- matrix(c(2.66090262285731, 1.36795706588105, 5.17589546114057,
    2.65200433797186, 1.39144872390597, 5.03137233509942, 2.65200433797186, 
    1.3107358581469, 5.07721615622715), nrow = 3, ncol = 3, byrow = TRUE)
  rownames(thisestimate) <- c("starting", "mle", "skew corr")
  colnames(thisestimate) <- c('RR', 'LL', 'UL')
  thisy <- matrix(c(4, 2, 4, 1, 16, 16, 18, 15, 5, 3, 10, 3, 79, 87, 90, 82,
    3.95, 3.625, 2, 1.82222222222222), nrow = 4, ncol = 5, byrow = FALSE)
  rownames(thisy) <- paste("Row", 1:4, sep = "")
  colnames(thisy) <- c("y1", "n1", "y2", "n2", "R.obs")
  
  expect_s4_class(ex1, "rrstr")
  expect_identical(ex1$estimator, "RR")
  expect_equal(ex1$estimate, thisestimate)
  expect_equal(ex1$hom, list(stat = 0.9540868, df = 3, p = 0.8123596))
  expect_equal(ex1$y, thisy)
  expect_identical(ex1$compare, c("b", "a"))
  
  ex2 <- RRstr(Y = table6, pf = FALSE)
  # expect_equal(ex1, ex2) # this will work when issue #17 is resolved
})