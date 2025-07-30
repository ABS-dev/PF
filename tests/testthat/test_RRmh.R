test_that("examples work", {
  ex1 <- RRmh(cbind(y, n) ~ tx + cluster(clus),
              Table6,
              compare = c("a", "b"),
              pf = FALSE)
  thisy <- data.frame(y1 = c(4, 2, 4, 1),
                      n1 = c(16, 16, 18, 15),
                      y2 = c(5, 3, 10, 3),
                      n2 = c(79, 87, 90, 82),
                      R.obs = c(3.95, 3.625, 2, 1.82222))
  rownames(thisy) <- paste("Row", 1:4, sep = "")


  expect_s4_class(ex1, "rr1")
  expect_identical(ex1$estimate |> signif(ex1$rnd) |> unname(),
                   c(2.67, 1.37, 5.23))
  expect_identical(ex1$estimator, "RR")
  expect_equal(ex1$y, thisy, tolerance = 0.0002)

  ex2 <- RRmh(Y = table6, pf = FALSE)
  expect_equal(ex2, ex1)
})
