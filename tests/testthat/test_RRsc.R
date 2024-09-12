context("RRsc")

test_that("examples work", {
  y1 <- c(4, 24, 12, 28)
  ex1 <- RRsc(y1)
  thisestimate <- matrix(c(rep(0.6111111, 5),
                           -0.04824591, -0.06395343, 0.02505437, 0.03282894,
                           0.03796954, 0.855726059905179, 0.835168605625983,
                           0.856770623929744, 0.855455017510877,
                           0.876306176546387),
                         nrow = 5, ncol = 3, byrow = FALSE)
  rownames(thisestimate) <- c("log method", "0.5 method", "MN method",
                              "score method", "skew corr")
  colnames(thisestimate) <- c("PF", "LL", "UL")

  expect_s4_class(ex1, "rrsc")
  expect_identical(ex1$estimator, "PF")
  expect_equal(ex1$estimate, thisestimate)
  expect_equal(ex1$y %>% as.numeric, c(4, 24, 12, 28))

  y2 <- matrix(c(4, 20, 12, 16), 2, 2, byrow = TRUE)
  ex2 <- RRsc(y2)
  expect_equal(ex1, ex2)



  data1 <- data.frame(group = rep(c("treated", "control"), each = 2),
                      y = c(1, 3, 7, 5),
                      n = c(12, 12, 14, 14),
                      cage = rep(paste("cage", 1:2), 2))
  ex3 <- RRsc(data = data1, formula = cbind(y, n) ~ group,
              compare = c("treated", "control"))
  expect_equal(ex1, ex3)

  data2 <- data1 %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(sum_y = sum(y),
              sum_n = sum(n))

  ex4 <- RRsc(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
              compare = c("treated", "control"))
  expect_equal(ex1, ex4)
})
