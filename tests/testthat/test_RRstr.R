context("RRstr")

test_that("examples work", {
  Table6$tx <- factor(Table6$tx, labels = c('vac', 'con'))
  ex1 <- RRstr(cbind(y,n) ~ tx + cluster(clus), Table6 ,
    pf = FALSE)
  thisestimate <- matrix(c(2.66090262285731, 1.36795706588105, 5.17589546114057,
    2.65200433797186, 1.39144872390597, 5.03137233509942, 2.65200433797186, 
    1.3107358581469, 5.07721615622715), nrow = 3, ncol = 3, byrow = TRUE)
  rownames(thisestimate) <- c("starting", "mle", "skew corr")
  colnames(thisestimate) <- c('RR', 'LL', 'UL')
  thisy <- data.frame(y1 = c(4, 2, 4, 1), 
    n1 = c(16, 16, 18, 15),
    y2 = c(5, 3, 10, 3),
    n2 = c(79, 87, 90, 82),
    R.obs = c(3.95, 3.625, 2, 1.82222222222222))
  rownames(thisy) <- paste("Row", 1:4, sep = "")
  
  
  expect_s4_class(ex1, "rrstr")
  expect_identical(ex1$estimator, "RR")
  expect_equal(ex1$estimate, thisestimate)
  expect_equal(ex1$hom, list(stat = 0.9540868, df = 3, p = 0.8123596))
  expect_equal(ex1$y, thisy)
  expect_identical(ex1$compare, c("vac", "con"))
  
  ex2 <- RRstr(Y = table6, pf = FALSE)
  expect_equal(ex1, ex2) 
})